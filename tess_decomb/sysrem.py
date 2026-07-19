#!/usr/bin/env python3
"""SysRem-style "de-combing" for TESS moving-target (asteroid) photometry.

PURPOSE
-------
Remove the shared per-sector instrumental systematics (momentum-dump comb at
328.8/n h, scattered-light ramps) from moving-target light curves BY
PROJECTION on an eigen-basis learned from ordinary field stars observed in
the same sector -- removal, instead of comb-proximity avoidance.

METHOD (PCA on a stellar ensemble + unweighted LSQ projection; deliberately
NOT full iterative SysRem with per-point weights -- see note below):

  1. build_eigenbasis(sector): download SAP-flux light curves for ~200
     ordinary SPOC 2-min targets from the same TESS sector (sector-wide, not
     camera/ccd-matched: the dominant comb and scattered-light components are
     shared across the sector), cached under sysrem_cache/.
  2. Bin every reference star (and, at fit time, the target) onto a common
     30-min time grid; build a (stars x time) matrix of median-normalized
     mag deviations; take its SVD; keep the top K right-singular (time-domain)
     vectors as "eigen-systematics".
  3. Fit an asteroid's UNBINNED light curve as offset + linear combination of
     the K eigen-systematics (interpolated to the asteroid's own timestamps),
     least squares; subtract the systematics part (not the offset).
  4. Smoke-validate: synthetic-signal injection/recovery on a quiet held-out
     reference star that never entered the SVD (`validate` subcommand).

NOTE ON THE SIMPLIFICATION. The unweighted, non-iterative projection with
fixed K = 5 is not a prototype shortcut: it is the exact estimator the
validation campaign calibrated. The 17,372-injection recovery rates, the
kill/survive asymmetry, the ROC-derived thresholds, the composition
bootstrap, and the 0.45x-baseline validity domain (validation/ in this
repository) all apply to THIS variant. Substituting full iterative SysRem
(or changing K) changes the estimator and voids those calibrations; any such
change requires re-validation.

Everything here builds from PUBLIC data: the eigen-basis is learned from SPOC
2-min field-star light curves downloaded from MAST (cached locally). The full
validation campaign (17,372 injections, null tests, ROC, bootstrap, validity
domain) lives in validation/ in this repository.

CLI
---
    tess-decomb ensemble --sector 18
    tess-decomb decomb --lc <path> --sector 18 --k 5 --out <path>
    tess-decomb validate --sector 18
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import warnings

import numpy as np
import pandas as pd
import requests

from .lightcurve import ORBIT_H, detrend, fourier_amp, load_clean  # noqa: F401

DEFAULT_CACHE_DIR = os.environ.get("TESS_DECOMB_CACHE",
                                   os.path.join(os.getcwd(), "sysrem_cache"))

# ---- quality bits (SPOC/Kepler convention) --------------------------------
# We mask only cadences that are genuinely unusable (bad pointing / safe
# mode). We deliberately KEEP reaction-wheel-desaturation (bit 32),
# argabrightening (bit 16) and impulsive-outlier (bit 256) cadences in the
# reference-star ensemble -- those ARE the systematics we want the eigenbasis
# to learn.
QUAL_ATTITUDE_TWEAK = 1
QUAL_SAFE_MODE = 2
QUAL_COARSE_POINT = 4
QUAL_EARTH_POINT = 8
QUAL_MANUAL_EXCLUDE = 64
BAD_QUALITY_BITS = (QUAL_ATTITUDE_TWEAK | QUAL_SAFE_MODE | QUAL_COARSE_POINT
                     | QUAL_EARTH_POINT | QUAL_MANUAL_EXCLUDE)  # = 79

DEFAULT_N_STARS = 200
DEFAULT_N_HOLDOUT = 20
DEFAULT_BIN_MIN = 30.0
DEFAULT_K_MAX = 8
DEFAULT_SEED = 42
MIN_STAR_POINTS = 200          # minimum unbinned cadences to keep a reference star
BIN_COVERAGE_MIN = 0.30        # keep a time bin only if >=30% of stars have data there
STAR_OUTLIER_MAD = 8.0         # loose clip on raw per-star mag (cosmic rays only)
DOWNLOAD_WORKERS = 8
MAST_DL_URL = "https://mast.stsci.edu/api/v0.1/Download/file?uri={uri}"
FNAME_TIC_RE = re.compile(r"-s\d{4}-(\d+)-\d+-s_lc\.fits$")


# ============================================================================
# Sector target index (cached MAST query) + download
# ============================================================================

def _obs_index_path(cache_dir: str, sector: int) -> str:
    return os.path.join(cache_dir, f"sector{sector:04d}_obs_index.csv")


def _lc_dir(cache_dir: str, sector: int) -> str:
    d = os.path.join(cache_dir, f"sector{sector:04d}_lc")
    os.makedirs(d, exist_ok=True)
    return d


def get_sector_obs_index(sector: int, cache_dir: str = DEFAULT_CACHE_DIR) -> pd.DataFrame:
    """Cached table of (ticid, obsid) for every SPOC 2-min timeseries obs in a sector."""
    path = _obs_index_path(cache_dir, sector)
    if os.path.exists(path):
        return pd.read_csv(path, dtype={"ticid": str, "obsid": str})
    os.makedirs(cache_dir, exist_ok=True)
    warnings.filterwarnings("ignore")
    from astroquery.mast import Observations
    t0 = time.time()
    obs = Observations.query_criteria(obs_collection="TESS", dataproduct_type="timeseries",
                                       sequence_number=sector, provenance_name="SPOC")
    print(f"[sector {sector}] MAST query returned {len(obs)} obs in {time.time()-t0:.1f}s")
    df = pd.DataFrame({"ticid": [str(x) for x in obs["target_name"]],
                        "obsid": [str(x) for x in obs["obsid"]]})
    df = df.drop_duplicates(subset="ticid").reset_index(drop=True)
    df.to_csv(path, index=False)
    return df


def _download_one(uri: str, dest: str) -> bool:
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return True
    try:
        r = requests.get(MAST_DL_URL.format(uri=uri), timeout=60)
        r.raise_for_status()
        tmp = dest + ".part"
        with open(tmp, "wb") as f:
            f.write(r.content)
        os.replace(tmp, dest)
        return True
    except Exception as e:
        sys.stderr.write(f"  download failed {uri}: {e}\n")
        return False


def download_star_lightcurves(sector: int, ticids: list[str], cache_dir: str = DEFAULT_CACHE_DIR
                               ) -> dict[str, str]:
    """Download SAP LC fits for the given TIC ids (sector-wide SPOC), flat-cached
    as sysrem_cache/sector<NNNN>_lc/<ticid>.fits. Returns {ticid: fits_path} for
    every star that is now on disk (already-cached or freshly downloaded)."""
    warnings.filterwarnings("ignore")
    from astroquery.mast import Observations
    lc_dir = _lc_dir(cache_dir, sector)
    have = {t: os.path.join(lc_dir, f"{t}.fits") for t in ticids
            if os.path.exists(os.path.join(lc_dir, f"{t}.fits"))}
    missing = [t for t in ticids if t not in have]
    if not missing:
        return have

    idx = get_sector_obs_index(sector, cache_dir)
    idx = idx[idx.ticid.isin(missing)]
    obsids = list(idx.obsid)
    if not obsids:
        return have
    prods = Observations.get_product_list(obsids)
    lc = prods[prods["productSubGroupDescription"] == "LC"]
    jobs = []
    for uri, fname in zip(lc["dataURI"], lc["productFilename"]):
        m = FNAME_TIC_RE.search(str(fname))
        if not m:
            continue
        ticid = str(int(m.group(1)))  # strip zero-padding
        if ticid not in missing:
            continue
        jobs.append((str(uri), os.path.join(lc_dir, f"{ticid}.fits"), ticid))

    from concurrent.futures import ThreadPoolExecutor, as_completed
    t0 = time.time()
    ok = 0
    with ThreadPoolExecutor(max_workers=DOWNLOAD_WORKERS) as ex:
        futs = {ex.submit(_download_one, uri, dest): ticid for uri, dest, ticid in jobs}
        for fut in as_completed(futs):
            ticid = futs[fut]
            if fut.result():
                have[ticid] = os.path.join(lc_dir, f"{ticid}.fits")
                ok += 1
    print(f"[sector {sector}] downloaded {ok}/{len(jobs)} new LC files in {time.time()-t0:.1f}s")
    return have


# ============================================================================
# Per-star SAP-flux -> mag loading
# ============================================================================

def load_star_sap_mag(fits_path: str) -> tuple[np.ndarray, np.ndarray]:
    """Return (t_days BTJD, mag) for a SPOC LC fits file using SAP_FLUX
    (NOT PDCSAP -- PDC already removes the common-mode systematics we want to
    model), minimal quality masking (BAD_QUALITY_BITS only)."""
    from astropy.io import fits
    with fits.open(fits_path) as f:
        d = f[1].data
        t = np.asarray(d["TIME"], dtype=float)
        flux = np.asarray(d["SAP_FLUX"], dtype=float)
        qual = np.asarray(d["QUALITY"], dtype=int)
    good = np.isfinite(t) & np.isfinite(flux) & (flux > 0) & ((qual & BAD_QUALITY_BITS) == 0)
    t, flux = t[good], flux[good]
    if len(t) < MIN_STAR_POINTS:
        return t, np.array([])
    mag = -2.5 * np.log10(flux)
    med = np.median(mag)
    mad = 1.4826 * np.median(np.abs(mag - med))
    if mad > 0:
        keep = np.abs(mag - med) < STAR_OUTLIER_MAD * mad
        t, mag = t[keep], mag[keep]
    o = np.argsort(t)
    return t[o], mag[o]


def bin_to_grid(t_days: np.ndarray, y: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Median-bin (t,y) onto pre-defined bin edges (days); NaN where empty."""
    nb = len(edges) - 1
    out = np.full(nb, np.nan)
    if len(t_days) == 0:
        return out
    idx = np.clip(np.searchsorted(edges, t_days, side="right") - 1, 0, nb - 1)
    for b in np.unique(idx):
        out[b] = np.median(y[idx == b])
    return out


# ============================================================================
# Eigenbasis (reference ensemble -> SVD)
# ============================================================================

def _ensemble_cache_path(cache_dir: str, sector: int, n_stars: int, n_holdout: int,
                          bin_min: float, k_max: int, seed: int) -> str:
    return os.path.join(cache_dir,
                         f"ensemble_s{sector:04d}_n{n_stars}_h{n_holdout}"
                         f"_bin{int(bin_min)}_k{k_max}_seed{seed}.npz")


def build_eigenbasis(sector: int, n_stars: int = DEFAULT_N_STARS,
                      n_holdout: int = DEFAULT_N_HOLDOUT, bin_min: float = DEFAULT_BIN_MIN,
                      k_max: int = DEFAULT_K_MAX, cache_dir: str = DEFAULT_CACHE_DIR,
                      seed: int = DEFAULT_SEED, force: bool = False) -> dict:
    """Build (or load cached) top-K_max eigen-systematics for a sector from an
    ensemble of n_stars field stars (used for the SVD) plus n_holdout extra
    stars reserved untouched for the injection-recovery test (leave-one-out:
    never enter the SVD)."""
    cache_path = _ensemble_cache_path(cache_dir, sector, n_stars, n_holdout, bin_min, k_max, seed)
    if os.path.exists(cache_path) and not force:
        z = np.load(cache_path, allow_pickle=False)
        return {"sector": sector, "bin_min": bin_min, "k_max": k_max,
                "grid_edges": z["grid_edges"], "grid_centers": z["grid_centers"],
                "eigvectors": z["eigvectors"], "explained_var_ratio": z["explained_var_ratio"],
                "basis_ticids": [str(x) for x in z["basis_ticids"]],
                "holdout_ticids": [str(x) for x in z["holdout_ticids"]],
                "cache_dir": cache_dir}

    t0 = time.time()
    idx = get_sector_obs_index(sector, cache_dir)
    rng = np.random.default_rng(seed)
    pool = idx.ticid.sample(frac=1.0, random_state=seed).tolist()
    target_total = n_stars + n_holdout

    basis_ticids: list[str] = []
    holdout_ticids: list[str] = []
    star_curves: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    pool_i = 0
    # Draw candidates in growing batches, download, keep the ones with enough
    # good points, until we fill n_stars (basis) + n_holdout (holdout) or the
    # pool (capped at 3x target) is exhausted.
    max_draw = min(len(pool), 3 * target_total)
    while len(basis_ticids) + len(holdout_ticids) < target_total and pool_i < max_draw:
        batch = pool[pool_i: pool_i + max(40, target_total)]
        pool_i += len(batch)
        paths = download_star_lightcurves(sector, batch, cache_dir)
        for tic in batch:
            if len(basis_ticids) + len(holdout_ticids) >= target_total:
                break
            fp = paths.get(tic)
            if not fp:
                continue
            t_days, mag = load_star_sap_mag(fp)
            if len(t_days) < MIN_STAR_POINTS:
                continue
            star_curves[tic] = (t_days, mag)
            if len(basis_ticids) < n_stars:
                basis_ticids.append(tic)
            else:
                holdout_ticids.append(tic)

    if len(basis_ticids) < 0.5 * n_stars:
        raise RuntimeError(f"[sector {sector}] only {len(basis_ticids)} usable reference "
                            f"stars (< 50% of requested {n_stars}); SPOC coverage too thin")

    all_t = np.concatenate([star_curves[t][0] for t in basis_ticids])
    t_min, t_max = float(all_t.min()), float(all_t.max())
    step = bin_min / (24.0 * 60.0)
    edges = np.arange(t_min, t_max + step, step)
    centers = 0.5 * (edges[:-1] + edges[1:])

    matrix = np.full((len(basis_ticids), len(centers)), np.nan)
    for i, tic in enumerate(basis_ticids):
        t_days, mag = star_curves[tic]
        row = bin_to_grid(t_days, mag, edges)
        med = np.nanmedian(row)
        matrix[i] = row - med  # median-normalized mag deviation

    coverage = np.mean(np.isfinite(matrix), axis=0)
    keep_bins = coverage >= BIN_COVERAGE_MIN
    matrix = matrix[:, keep_bins]
    centers = centers[keep_bins]
    edges_kept = np.concatenate([edges[:-1][keep_bins], edges[1:][keep_bins][-1:]])
    matrix = np.nan_to_num(matrix, nan=0.0)

    U, S, Vt = np.linalg.svd(matrix, full_matrices=False)
    kk = min(k_max, Vt.shape[0])
    eigvectors = Vt[:kk]
    evr = (S ** 2 / np.sum(S ** 2))[:kk]

    os.makedirs(cache_dir, exist_ok=True)
    np.savez(cache_path, grid_edges=edges_kept, grid_centers=centers,
             eigvectors=eigvectors, explained_var_ratio=evr,
             basis_ticids=np.array(basis_ticids), holdout_ticids=np.array(holdout_ticids))
    print(f"[sector {sector}] ensemble built: {len(basis_ticids)} basis + "
          f"{len(holdout_ticids)} holdout stars, {matrix.shape[1]} time bins, "
          f"K={kk} (EVR top-3={np.round(evr[:3], 3).tolist()}) in {time.time()-t0:.1f}s")

    return {"sector": sector, "bin_min": bin_min, "k_max": kk,
            "grid_edges": edges_kept, "grid_centers": centers,
            "eigvectors": eigvectors, "explained_var_ratio": evr,
            "basis_ticids": basis_ticids, "holdout_ticids": holdout_ticids,
            "cache_dir": cache_dir, "build_seconds": time.time() - t0,
            "n_basis": len(basis_ticids), "n_holdout": len(holdout_ticids)}


# ============================================================================
# Projection fit / subtract
# ============================================================================

def get_eigencomponents_at(ensemble: dict, t_query_days: np.ndarray, K: int) -> np.ndarray:
    """Interpolate the first K eigen-systematics onto arbitrary timestamps,
    de-meaned over the query span so the fit's constant column fully absorbs
    the mean level (subtracting the systematics never shifts absolute mag)."""
    centers = ensemble["grid_centers"]
    eig = ensemble["eigvectors"][:K]
    out = np.empty((K, len(t_query_days)))
    for k in range(K):
        c = np.interp(t_query_days, centers, eig[k], left=eig[k, 0], right=eig[k, -1])
        out[k] = c - c.mean()
    return out


def fit_and_subtract(t_days: np.ndarray, mag: np.ndarray, ensemble: dict, K: int) -> dict:
    """LSQ-fit mag(t) = offset + sum_k coeff_k * eig_k(t); return the mag curve
    with only the systematics part removed (offset/mean preserved)."""
    comps = get_eigencomponents_at(ensemble, t_days, K)  # (K, N)
    design = np.column_stack([np.ones_like(t_days), *comps])
    coeffs, *_ = np.linalg.lstsq(design, mag, rcond=None)
    systematics = design[:, 1:] @ coeffs[1:]
    corrected = mag - systematics
    resid = mag - design @ coeffs
    ss_res = float(np.sum((mag - mag.mean()) ** 2))
    r2 = 1.0 - float(np.sum(resid ** 2)) / ss_res if ss_res > 0 else 0.0
    return {"corrected": corrected, "systematics": systematics, "coeffs": coeffs,
            "r2_systematics": r2}


# ============================================================================
# Lomb-Scargle helpers
# ============================================================================

def ls_power_at(t_days: np.ndarray, y: np.ndarray, period_h: float) -> float:
    from astropy.timeseries import LombScargle
    ls = LombScargle(t_days, y)
    return float(ls.power(np.array([24.0 / period_h]))[0])


def ls_scan(t_days: np.ndarray, y: np.ndarray, pmin_h: float, pmax_h: float, n: int = 20000):
    from astropy.timeseries import LombScargle
    per = np.linspace(pmin_h, pmax_h, n)
    ls = LombScargle(t_days, y)
    pw = ls.power(24.0 / per)
    i = int(np.argmax(pw))
    return per, pw, float(per[i]), float(pw[i])


def decomb_asteroid(lc_path: str, sector: int, K: int, cache_dir: str = DEFAULT_CACHE_DIR,
                     n_stars: int = DEFAULT_N_STARS, n_holdout: int = DEFAULT_N_HOLDOUT,
                     bin_min: float = DEFAULT_BIN_MIN) -> dict:
    ensemble = build_eigenbasis(sector, n_stars=n_stars, n_holdout=n_holdout, bin_min=bin_min,
                                 k_max=max(K, DEFAULT_K_MAX), cache_dir=cache_dir)
    t, m, e, n_clip = load_clean(lc_path)
    fit = fit_and_subtract(t, m, ensemble, K)
    return {"t": t, "m_raw": m, "err": e, "m_decomb": fit["corrected"],
            "systematics": fit["systematics"], "coeffs": fit["coeffs"],
            "r2_systematics": fit["r2_systematics"], "ensemble": ensemble}


# ============================================================================
# Validation tests
# ============================================================================

def test_injection(K: int = 5, cache_dir: str = DEFAULT_CACHE_DIR, P_inj_h: float = 90.0,
                    amp_pp_mag: float = 0.3, sector: int = 18) -> dict:
    ensemble = build_eigenbasis(sector, cache_dir=cache_dir, k_max=max(K, DEFAULT_K_MAX))
    lc_dir = _lc_dir(cache_dir, sector)
    # pick the holdout star with the smallest post-projection residual scatter ("quiet")
    best_tic, best_mad, best_curve = None, np.inf, None
    for tic in ensemble["holdout_ticids"]:
        fp = os.path.join(lc_dir, f"{tic}.fits")
        if not os.path.exists(fp):
            continue
        t_days, mag = load_star_sap_mag(fp)
        if len(t_days) < MIN_STAR_POINTS:
            continue
        fit = fit_and_subtract(t_days, mag, ensemble, K)
        mad = 1.4826 * float(np.median(np.abs(fit["corrected"] - np.median(fit["corrected"]))))
        if mad < best_mad:
            best_mad, best_tic, best_curve = mad, tic, (t_days, mag)
    if best_tic is None:
        raise RuntimeError(f"[sector {sector}] no usable holdout star for injection test")

    t_days, mag = best_curve
    amp_semi = amp_pp_mag / 2.0
    injected = mag + amp_semi * np.sin(2 * np.pi * t_days / (P_inj_h / 24.0))

    fit_raw = fit_and_subtract(t_days, injected, ensemble, K)  # decomb applied to injected curve
    y_raw_detrend = detrend(t_days, injected)
    y_decomb_detrend = detrend(t_days, fit_raw["corrected"])

    _, _, peak_P_before, peak_pw_before = ls_scan(t_days, y_raw_detrend, 20.0, 200.0)
    _, _, peak_P_after, peak_pw_after = ls_scan(t_days, y_decomb_detrend, 20.0, 200.0)

    amp_before = fourier_amp(t_days, y_raw_detrend, P_inj_h, K=1)[1]   # A1 = semi-amplitude
    amp_after = fourier_amp(t_days, y_decomb_detrend, P_inj_h, K=1)[1]
    amp_recovery_frac = amp_after / amp_semi if amp_semi else float("nan")
    period_err_frac = abs(peak_P_after - P_inj_h) / P_inj_h

    verdict = "PASS" if (period_err_frac <= 0.05 and amp_recovery_frac >= 0.90) else "FAIL"
    return dict(test="synthetic 90h/0.3mag injection", K=K, sector=sector,
                quiet_star_ticid=best_tic, quiet_star_holdout_mad=best_mad,
                P_inj_h=P_inj_h, amp_pp_injected_mag=amp_pp_mag, amp_semi_injected_mag=amp_semi,
                peak_period_before_h=peak_P_before, peak_period_after_h=peak_P_after,
                period_err_frac=period_err_frac,
                fourier_semiamp_before=amp_before, fourier_semiamp_after=amp_after,
                amp_recovery_frac=amp_recovery_frac, verdict=verdict)


def run_all_validation(K: int = 5, cache_dir: str = DEFAULT_CACHE_DIR,
                        out_json: str | None = None, sector: int = 18) -> dict:
    """Self-contained smoke validation: synthetic injection/recovery on a quiet
    held-out reference star (never entered the SVD). The survey-scale campaign
    (V1-V7) is archived under validation/ with its artifacts."""
    t0 = time.time()
    result = {"K": K, "generated": time.strftime("%Y-%m-%d %H:%M:%S")}
    result["injection"] = test_injection(K, cache_dir, sector=sector)
    result["total_seconds"] = time.time() - t0
    if out_json:
        with open(out_json, "w") as f:
            json.dump(result, f, indent=2, default=float)
    return result


# ============================================================================
# CLI
# ============================================================================

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("ensemble", help="build/cache the eigen-systematics for a sector")
    p1.add_argument("--sector", type=int, required=True)
    p1.add_argument("--n-stars", type=int, default=DEFAULT_N_STARS)
    p1.add_argument("--n-holdout", type=int, default=DEFAULT_N_HOLDOUT)
    p1.add_argument("--bin-min", type=float, default=DEFAULT_BIN_MIN)
    p1.add_argument("--k-max", type=int, default=DEFAULT_K_MAX)
    p1.add_argument("--cache-dir", default=DEFAULT_CACHE_DIR)
    p1.add_argument("--force", action="store_true")

    p2 = sub.add_parser("decomb", help="de-comb a single asteroid light curve")
    p2.add_argument("--lc", required=True)
    p2.add_argument("--sector", type=int, required=True)
    p2.add_argument("--k", type=int, default=5)
    p2.add_argument("--out", required=True)
    p2.add_argument("--cache-dir", default=DEFAULT_CACHE_DIR)
    p2.add_argument("--n-stars", type=int, default=DEFAULT_N_STARS)
    p2.add_argument("--n-holdout", type=int, default=DEFAULT_N_HOLDOUT)
    p2.add_argument("--bin-min", type=float, default=DEFAULT_BIN_MIN)

    p3 = sub.add_parser("validate", help="injection/recovery smoke test on a holdout star")
    p3.add_argument("--k", type=int, default=5)
    p3.add_argument("--sector", type=int, default=18)
    p3.add_argument("--cache-dir", default=DEFAULT_CACHE_DIR)
    p3.add_argument("--out-json", default=None)

    a = ap.parse_args()
    if a.cmd == "ensemble":
        build_eigenbasis(a.sector, n_stars=a.n_stars, n_holdout=a.n_holdout, bin_min=a.bin_min,
                          k_max=a.k_max, cache_dir=a.cache_dir, force=a.force)
    elif a.cmd == "decomb":
        run = decomb_asteroid(a.lc, a.sector, a.k, cache_dir=a.cache_dir,
                               n_stars=a.n_stars, n_holdout=a.n_holdout, bin_min=a.bin_min)
        out = pd.DataFrame({"time": run["t"], "mag_raw": run["m_raw"],
                             "mag_decomb": run["m_decomb"], "err": run["err"],
                             "systematics_model": run["systematics"]})
        out.to_csv(a.out, index=False)
        print(f"wrote {a.out}  (R^2 systematics={run['r2_systematics']:.3f})")
    elif a.cmd == "validate":
        result = run_all_validation(K=a.k, cache_dir=a.cache_dir, out_json=a.out_json,
                                    sector=a.sector)
        print(json.dumps(result, indent=2, default=float))
        if a.out_json:
            print(f"\nwrote {a.out_json}")


if __name__ == "__main__":
    main()
