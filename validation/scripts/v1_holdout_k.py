#!/usr/bin/env python3
"""V1 (VALIDATION_PLAN.md): data-driven K via holdout-star cross-validation.

Phase A -- per sector: project the top-K eigenvectors (K=1..8) out of each of the 20
holdout field stars (never used in the SVD). Normalized holdout residual vs K gives
K*_min (argmin) and K*_1pct (first K whose marginal improvement < 1%). Overfit shows
up as the holdout residual bottoming out / rising while in-basis residual keeps falling.

Phase B -- verdict stability: re-run the catalog decomb verdicts at K in {3,4,5,6,7}
for every comb/slow/contamination-triggered object; report the fraction of verdicts
that change vs K=5.

All cached; no network. Outputs: v1_holdout_k_sectors.csv, v1_verdict_stability.csv,
V1_RESULTS.md.
"""
import csv, glob, os, re, sys
import numpy as np

NG = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(NG, "tess-tools"))
sys.path.insert(0, os.path.join(NG, "tess-themis"))
import sysrem_decomb as sd  # noqa: E402

# ---------- Phase A: holdout cross-validation ----------
def holdout_curves(sector, npz):
    lcdir = sd._lc_dir(sd.DEFAULT_CACHE_DIR, sector)
    edges = npz["grid_edges"]
    out = []
    for tic in npz["holdout_ticids"]:
        hits = glob.glob(os.path.join(lcdir, f"*{tic}*"))
        if not hits:
            continue
        try:
            t, y = sd.load_star_sap_mag(hits[0])
            b = sd.bin_to_grid(t, y, edges)
            if np.isfinite(b).sum() >= 50:
                out.append(b - np.nanmedian(b))
        except Exception:
            continue
    return out

def sector_kstar(sector, npz):
    eig = npz["eigvectors"]          # (8, ngrid)
    curves = holdout_curves(sector, npz)
    if len(curves) < 5:
        return None
    kmax = eig.shape[0]
    frac = np.full(kmax + 1, np.nan)  # frac[K] = mean normalized residual rms at rank K
    resid0 = [np.nanstd(c) for c in curves]
    frac[0] = 1.0
    for K in range(1, kmax + 1):
        vals = []
        for c, r0 in zip(curves, resid0):
            ok = np.isfinite(c)
            A = eig[:K, ok].T
            coef, *_ = np.linalg.lstsq(A, c[ok], rcond=None)
            r = c[ok] - A @ coef
            if r0 > 0:
                vals.append(np.std(r) / r0)
        frac[K] = np.mean(vals)
    kmin = int(np.nanargmin(frac[1:]) + 1)
    k1 = kmax
    for K in range(1, kmax):
        if frac[K] - frac[K + 1] < 0.01:   # <1% marginal improvement
            k1 = K
            break
    return dict(sector=sector, n_holdout=len(curves), k_min=kmin, k_1pct=k1,
                **{f"resid_k{K}": round(float(frac[K]), 4) for K in range(1, kmax + 1)})

print("=== Phase A: holdout K* per sector ===", flush=True)
rowsA = []
for f in sorted(glob.glob(os.path.join(NG, "tess-tools/sysrem_cache/ensemble_s*.npz"))):
    s = int(re.search(r"_s0*(\d+)_", f).group(1))
    npz = np.load(f, allow_pickle=True)
    r = sector_kstar(s, npz)
    if r:
        rowsA.append(r)
        print(f"  s{s}: K*_min={r['k_min']} K*_1pct={r['k_1pct']} (n={r['n_holdout']})", flush=True)
with open(os.path.join(NG, "v1_holdout_k_sectors.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rowsA[0].keys())); w.writeheader(); w.writerows(rowsA)

kmins = [r["k_min"] for r in rowsA]; k1s = [r["k_1pct"] for r in rowsA]
print(f"Phase A done: {len(rowsA)} sectors | K*_min median={np.median(kmins)} (range {min(kmins)}-{max(kmins)}) "
      f"| K*_1pct median={np.median(k1s)} (range {min(k1s)}-{max(k1s)})", flush=True)

# ---------- Phase B: verdict stability ----------
print("=== Phase B: verdict stability across K in {3,4,5,6,7} ===", flush=True)
import importlib
import decomb_check as dc
DIRS = {"tess-themis": "Themis", "tess-phocaea": "Phocaea", "tess-koronis": "Koronis",
        "tess-agnia": "Agnia", "tess-ursula": "Ursula", "tess-dora": "Dora",
        "tess-veritas": "Veritas", "tess-baptistina": "Baptistina",
        "tess-hoffmeister": "Hoffmeister", "tess-cameron": "Cameron", "tess-padua": "Padua",
        "tess-beltwide": "Beltwide", "tess-beltwide2": "Beltwide2", "tess-beltwide3": "Beltwide3",
        "tess-beltwide4": "Beltwide4", "tess-beltwide5": "Beltwide5", "tess-beltwide6": "Beltwide6",
        "tess-beltwide7": "Beltwide7", "tess-beltwide8": "Beltwide8", "tess-beltwide9": "Beltwide9",
        "tess-beltwide10": "Beltwide10"}
KS = [3, 4, 5, 6, 7]
rowsB = []
for d in DIRS:
    p = os.path.join(NG, d, "census_refined.csv")
    if not os.path.exists(p):
        continue
    for r in csv.DictReader(open(p)):
        fl = str(r.get("flags", ""))
        slow = float(r.get("P_rot_new") or 0) >= 30
        trig = bool(dc.COMB_RE.search(fl)) or slow or bool(dc.CONTAM_RE.search(fl))
        if not trig:
            continue
        num = str(r["num"]); Pp = float(r["P_phot_h"])
        verd = {}
        for K in KS:
            v = dc.check_object(num, Pp, os.path.join(NG, d), K)
            verd[K] = v.split("(")[0]
        rowsB.append(dict(dir=d, num=num, P_phot=Pp,
                          **{f"K{K}": verd[K] for K in KS},
                          stable=int(len(set(verd.values())) == 1),
                          same_as_k5=int(all(verd[K] == verd[5] for K in KS))))
        print(f"  {d}/{num}: " + " ".join(f"K{K}={verd[K]}" for K in KS), flush=True)
with open(os.path.join(NG, "v1_verdict_stability.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rowsB[0].keys())); w.writeheader(); w.writerows(rowsB)
n = len(rowsB); stable = sum(r["stable"] for r in rowsB)
print(f"Phase B done: {n} triggered objects | fully K-stable verdicts: {stable}/{n} ({100*stable/max(n,1):.1f}%)", flush=True)

# ---------- summary ----------
with open(os.path.join(NG, "V1_RESULTS.md"), "w") as f:
    f.write(f"""# V1 results: holdout-K cross-validation ({len(rowsA)} sectors) + verdict stability

## Phase A (holdout stars, never in the SVD basis)
- K*_min (argmin holdout residual): median {np.median(kmins):.0f}, range {min(kmins)}-{max(kmins)}
- K*_1pct (first K with <1% marginal improvement): median {np.median(k1s):.0f}, range {min(k1s)}-{max(k1s)}
- Per-sector table: v1_holdout_k_sectors.csv (residual-vs-K curves per sector)

## Phase B (catalog decomb verdicts re-run at K=3,4,5,6,7)
- Triggered objects tested: {n}
- Verdicts identical across ALL K in [3,7]: {stable}/{n} ({100*stable/max(n,1):.1f}%)
- Rows: v1_verdict_stability.csv (per-object verdict at each K)

Generated by v1_holdout_k.py (fully offline, cached ensembles + star curves).
""")
print("V1_RESULTS.md written", flush=True)
