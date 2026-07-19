#!/usr/bin/env python3
"""V6 (VALIDATION_PLAN.md): validity domain -- signal loss vs P/baseline.

Part A: bin the existing 17k V2 injections by P/baseline; recovery + false-kill vs ratio.
Part B: supplementary long-period injections (P up to 0.85 x baseline, the regime V2's
0.45 cap excluded) on ~60 tracks to map the fall-off.
Output: v6_validity.csv, V6_RESULTS.md, v6_validity_curve.png.
"""
import csv, glob, os, re, sys
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")
NG = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(NG, "tess-tools")); sys.path.insert(0, os.path.join(NG, "tess-themis"))
import sysrem_decomb as sd                     # noqa: E402
from themis_census import load_clean, detrend  # noqa: E402
from astropy.timeseries import LombScargle     # noqa: E402

# ---- Part A: P/baseline from existing V2 data ----
print("Part A: existing V2 injections vs P/baseline", flush=True)
# per-file baseline cache
base = {}
def baseline_of(fname):
    if fname in base:
        return base[fname]
    hits = glob.glob(os.path.join(NG, "tess-*", fname))
    if not hits:
        base[fname] = np.nan; return np.nan
    t, m, e, _ = load_clean(hits[0])
    b = (t.max() - t.min()) * 24.0
    base[fname] = b
    return b

rows = [r for r in csv.DictReader(open(os.path.join(NG, "v2_injections.csv"))) if r["verdict"]]
recs = []
for r in rows:
    b = baseline_of(r["file"])
    if not np.isfinite(b) or b <= 0:
        continue
    recs.append((float(r["P"]) / b, float(r["pw1"]) / max(float(r["pw0"]), 1e-9),
                 float(r["drop"]), r["verdict"], float(r["amp_in"]), float(r["amp_out"] or 0)))

# ---- Part B: long-P supplementary injections ----
print("Part B: long-P supplementary injections", flush=True)
KILL_DROP, KILL_R2, SURVIVE_DROP = 0.50, 0.30, 0.15
CACHED = {int(re.search(r"_s0*(\d+)_", f).group(1)) for f in glob.glob(os.path.join(NG, "tess-tools/sysrem_cache/ensemble_s*.npz"))}
subs = []
for cens in glob.glob(os.path.join(NG, "tess-*/census_objects.csv")):
    d = os.path.dirname(cens)
    try:
        cr = list(csv.DictReader(open(cens)))
    except Exception:
        continue
    for r in cr:
        if (r.get("tier") or "") != "NONE":
            continue
        for f in glob.glob(os.path.join(d, f"lc_{r['num']}_s*.csv")):
            s = int(re.search(r"_s(\d+)\.csv", f).group(1))
            if s in CACHED:
                subs.append((f, s)); break
seen = set(); subs = [x for x in subs if not (x[0] in seen or seen.add(x[0]))][:60]
RATIOS = [0.50, 0.60, 0.70, 0.85]
AMPS = [0.10, 0.30]
rng = np.random.default_rng(11)
nB = 0
for f, s in subs:
    try:
        ens = sd.build_eigenbasis(s, cache_dir=sd.DEFAULT_CACHE_DIR)
        t, m, e, _ = load_clean(f)
        if len(t) < 300:
            continue
        bh = (t.max() - t.min()) * 24.0
        for rat in RATIOS:
            P = rat * bh
            for A in AMPS:
                m_inj = m + (A / 2) * np.sin(2 * np.pi * t / (P / 24) + rng.uniform(0, 2 * np.pi))
                y0 = detrend(t, m_inj)
                p0 = float(LombScargle(t, y0).power(24 / P))
                fit = sd.fit_and_subtract(t, m_inj, ens, 5)
                y1 = detrend(t, fit["corrected"])
                p1 = float(LombScargle(t, y1).power(24 / P))
                drop = (p0 - p1) / max(p0, 1e-9); r2 = float(fit["r2_systematics"])
                v = ("killed" if (drop >= KILL_DROP and r2 >= KILL_R2)
                     else "survived" if drop <= SURVIVE_DROP else "inconclusive")
                recs.append((rat, p1 / max(p0, 1e-9), drop, v, A, np.nan))
                nB += 1
    except Exception:
        continue
print(f"Part B injections: {nB}", flush=True)

# ---- binning + outputs ----
recs = np.array(recs, dtype=object)
ratios = np.array([float(x[0]) for x in recs]); recov = np.array([float(x[1]) for x in recs])
verd = np.array([x[3] for x in recs])
bins = [0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.45, 0.60, 0.90]
lines = ["# V6 results: validity domain (signal retention vs P/baseline)\n",
         "| P/baseline bin | n | median recovery | false-kill% | survived% |", "|--|--|--|--|--|"]
xs, ys = [], []
with open(os.path.join(NG, "v6_validity.csv"), "w", newline="") as fo:
    w = csv.writer(fo); w.writerow(["ratio_lo", "ratio_hi", "n", "median_recovery", "false_kill_pct", "survived_pct"])
    for lo, hi in zip(bins[:-1], bins[1:]):
        sel = (ratios >= lo) & (ratios < hi)
        if sel.sum() < 20:
            continue
        mr = np.median(recov[sel]); fk = 100 * np.mean(verd[sel] == "killed"); sv = 100 * np.mean(verd[sel] == "survived")
        lines.append(f"| {lo:.2f}-{hi:.2f} | {sel.sum()} | {mr:.3f} | {fk:.2f} | {sv:.1f} |")
        w.writerow([lo, hi, int(sel.sum()), round(mr, 4), round(fk, 2), round(sv, 1)])
        xs.append((lo + hi) / 2); ys.append(mr)
open(os.path.join(NG, "V6_RESULTS.md"), "w").write("\n".join(lines) + "\n")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.figure(figsize=(7, 4.2))
plt.plot(xs, ys, "o-", color="crimson")
plt.axhline(0.85, ls="--", lw=0.8, color="gray"); plt.xlabel("P / sector baseline"); plt.ylabel("median power recovery after de-comb")
plt.title("De-comb signal retention vs period/baseline"); plt.tight_layout()
plt.savefig(os.path.join(NG, "v6_validity_curve.png"), dpi=130)
print("V6_RESULTS.md + v6_validity_curve.png written", flush=True)
