#!/usr/bin/env python3
"""V4 + V5 (VALIDATION_PLAN.md).

V4 -- threshold calibration: place the survive (<=15% drop) and kill (>=50% drop &
R2>=0.3) thresholds on the V2 (injected real signals) and V3 (pure comb) drop
distributions; report error rates at the operating point and the ROC sweep.

V5 -- ensemble robustness: (a) bootstrap the eigen-basis by resampling the ~220 cached
field stars per sector (B=6 resamples of 200 with replacement, re-SVD, re-verdict all
triggered objects): verdict stability under ensemble composition. (b) split-ensemble
test: two disjoint 110-star halves -> two independent bases -> verdict agreement.
Fully offline. Outputs: V4_V5_RESULTS.md, v5_bootstrap.csv.
"""
import csv, glob, os, re, sys
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")
NG = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(NG, "tess-tools"))
sys.path.insert(0, os.path.join(NG, "tess-themis"))
import sysrem_decomb as sd                       # noqa: E402
from themis_census import load_clean, detrend    # noqa: E402
from astropy.timeseries import LombScargle       # noqa: E402

# ---------------- V4: threshold calibration ----------------
print("=== V4: threshold calibration ===", flush=True)
v2 = [r for r in csv.DictReader(open(os.path.join(NG, "v2_injections.csv"))) if r["verdict"]]
v3 = list(csv.DictReader(open(os.path.join(NG, "v3_null_test.csv"))))
drops_real = np.array([float(r["drop"]) for r in v2 if r["drop"]])
r2_real = np.array([float(r["r2"]) for r in v2 if r["r2"]])
drops_null = np.array([float(r["drop"]) for r in v3])
r2_null = np.array([float(r["r2"]) for r in v3])

lines = ["# V4 + V5 results\n", "## V4: verdict-threshold calibration (ROC on V2 real vs V3 null)\n"]
lines.append("| kill-drop threshold | false-kill% (real, R2>=0.3 arm) | true-kill% (null) |")
lines.append("|--|--|--|")
for kd in (0.30, 0.40, 0.50, 0.60, 0.70):
    fk = 100 * np.mean((drops_real >= kd) & (r2_real >= 0.30))
    tk = 100 * np.mean((drops_null >= kd) & (r2_null >= 0.30))
    mark = " <- operating point" if abs(kd - 0.50) < 1e-9 else ""
    lines.append(f"| {kd:.2f} | {fk:.2f} | {tk:.1f} |{mark}")
lines.append("")
lines.append("| survive-drop threshold | real signals surviving% | null surviving% (false-survive) |")
lines.append("|--|--|--|")
for sv in (0.05, 0.10, 0.15, 0.20, 0.25):
    rs = 100 * np.mean(drops_real <= sv)
    ns = 100 * np.mean(drops_null <= sv)
    mark = " <- operating point" if abs(sv - 0.15) < 1e-9 else ""
    lines.append(f"| {sv:.2f} | {rs:.1f} | {ns:.1f} |{mark}")
lines.append("\nReal-signal drop distribution: median "
             f"{np.median(drops_real):+.3f}, p95 {np.percentile(drops_real,95):+.3f}; "
             f"null drops: median {np.median(drops_null):+.3f}, p95 {np.percentile(drops_null,95):+.3f} (n={len(drops_null)}).")
print("\n".join(lines[-3:]), flush=True)

# ---------------- V5: ensemble robustness ----------------
print("=== V5: ensemble bootstrap ===", flush=True)
KILL_DROP, KILL_R2, SURVIVE_DROP = 0.50, 0.30, 0.15

def load_sector_stars(sector, npz):
    lcdir = sd._lc_dir(sd.DEFAULT_CACHE_DIR, sector)
    edges = npz["grid_edges"]
    curves = []
    for f in glob.glob(os.path.join(lcdir, "*")):
        try:
            t, y = sd.load_star_sap_mag(f)
            b = sd.bin_to_grid(t, y, edges)
            if np.isfinite(b).sum() >= 100:
                curves.append(b - np.nanmedian(b))
        except Exception:
            continue
    return np.array(curves), edges

def make_basis(curves, idx, K=5):
    M = curves[idx]
    M = np.where(np.isfinite(M), M, 0.0)
    M = M - M.mean(axis=1, keepdims=True)
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    return Vt[:K]

def verdict_for(t, m, eig, edges, P):
    # project eig (K,ngrid) sampled at the asteroid's binned grid points
    y0 = detrend(t, m)
    pw0 = float(LombScargle(t, y0).power(24.0 / P))
    if pw0 < 0.03:
        return None
    centers = (edges[:-1] + edges[1:]) / 2.0
    comps = np.vstack([np.interp(t, centers, e, left=np.nan, right=np.nan) for e in eig])
    ok = np.all(np.isfinite(comps), axis=0)
    if ok.sum() < 50:
        return None
    A = comps[:, ok].T
    coef, *_ = np.linalg.lstsq(A, m[ok] - np.median(m[ok]), rcond=None)
    model = np.zeros_like(m); model[ok] = A @ coef
    resid = m - model
    ss_tot = np.var(m[ok]); r2 = 1 - np.var(resid[ok]) / ss_tot if ss_tot > 0 else 0.0
    y1 = detrend(t[ok], resid[ok])
    pw1 = float(LombScargle(t[ok], y1).power(24.0 / P))
    drop = (pw0 - pw1) / pw0
    v = ("KILLED" if (drop >= KILL_DROP and r2 >= KILL_R2)
         else "survived" if drop <= SURVIVE_DROP else "inconclusive")
    return v

# triggered objects grouped by their (single strongest) sector for efficiency
import decomb_check as dc
DIRS = [d for d in glob.glob(os.path.join(NG, "tess-*")) if os.path.exists(os.path.join(d, "census_refined.csv"))]
jobs = {}   # sector -> list of (file, P)
for d in DIRS:
    for r in csv.DictReader(open(os.path.join(d, "census_refined.csv"))):
        fl = str(r.get("flags", ""))
        slow = float(r.get("P_rot_new") or 0) >= 30
        if not (dc.COMB_RE.search(fl) or slow or dc.CONTAM_RE.search(fl)):
            continue
        for f in glob.glob(os.path.join(d, f"lc_{r['num']}_s*.csv")):
            s = int(re.search(r"_s(\d+)\.csv", f).group(1))
            jobs.setdefault(s, []).append((f, float(r["P_phot_h"])))

B = 6
rows = []
rng = np.random.default_rng(7)
ens_files = {int(re.search(r"_s0*(\d+)_", f).group(1)): f
             for f in glob.glob(os.path.join(NG, "tess-tools/sysrem_cache/ensemble_s*.npz"))}
for si, (sector, items) in enumerate(sorted(jobs.items())):
    if sector not in ens_files:
        continue
    npz = np.load(ens_files[sector], allow_pickle=True)
    curves, edges = load_sector_stars(sector, npz)
    if len(curves) < 120:
        continue
    n = len(curves)
    # B bootstrap bases + 2 disjoint halves
    bases = [make_basis(curves, rng.integers(0, n, 200)) for _ in range(B)]
    perm = rng.permutation(n)
    half1 = make_basis(curves, perm[:n // 2]); half2 = make_basis(curves, perm[n // 2:])
    for f, P in items:
        try:
            t, m, e, _ = load_clean(f)
            if len(t) < 50 or (m.max() - m.min()) > 5.0:
                continue
            vb = [verdict_for(t, m, b, edges, P) for b in bases]
            vb = [v for v in vb if v]
            v1 = verdict_for(t, m, half1, edges, P); v2_ = verdict_for(t, m, half2, edges, P)
            if not vb:
                continue
            rows.append(dict(sector=sector, file=os.path.basename(f), P=P,
                             boot_verdicts="|".join(vb), boot_stable=int(len(set(vb)) == 1),
                             half1=v1 or "", half2=v2_ or "",
                             halves_agree=int(bool(v1) and v1 == v2_)))
        except Exception:
            continue
    if (si + 1) % 10 == 0:
        print(f"  sectors {si+1}/{len(jobs)}", flush=True)

with open(os.path.join(NG, "v5_bootstrap.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
nb = len(rows); st = sum(r["boot_stable"] for r in rows); ha = sum(r["halves_agree"] for r in rows)
h_both = sum(1 for r in rows if r["half1"] and r["half2"])
lines.append(f"\n## V5: ensemble-composition robustness ({nb} object-sector tests, B={B} bootstraps)")
lines.append(f"- Bootstrap-stable verdicts (identical across all {B} resampled 200-star bases): {st}/{nb} ({100*st/max(nb,1):.1f}%)")
lines.append(f"- Disjoint-half test (two independent ~110-star bases agree): {ha}/{h_both} ({100*ha/max(h_both,1):.1f}%)")
open(os.path.join(NG, "V4_V5_RESULTS.md"), "w").write("\n".join(lines) + "\n")
print(f"V5: bootstrap stable {st}/{nb}; halves agree {ha}/{h_both}", flush=True)
print("V4_V5_RESULTS.md written", flush=True)
