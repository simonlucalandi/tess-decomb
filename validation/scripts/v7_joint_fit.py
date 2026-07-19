#!/usr/bin/env python3
"""V7 (VALIDATION_PLAN.md): joint-fit prototype -- systematics eigen-basis and a Fourier
series at the trial period fitted SIMULTANEOUSLY, so the basis cannot absorb the signal.

Model: m(t) ~ c0 + poly3(t) + sum_k [a_k E_k(t)] + sum_j [b_j cos(j w t) + c_j sin(j w t)]
Signal significance: delta-BIC of the full model vs the no-Fourier model (collinearity
between basis and a comb-frequency signal is then handled naturally: Fourier terms only
win BIC if they explain variance the basis cannot).

Tests:
(A) Long-period injections (P/baseline 0.3-0.85): amplitude recovery joint vs sequential.
(B) The 46 V3 nulls: does the joint fit refuse to credit pure comb detections? (false-clear rate)
(C) Real slow rotators where sequential projection over-subtracted (3396 sectors, 25-41% drops).

Output: V7_RESULTS.md, v7_joint_injections.csv.
"""
import csv, glob, os, re, sys
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")
NG = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(NG, "tess-tools")); sys.path.insert(0, os.path.join(NG, "tess-themis"))
import sysrem_decomb as sd                     # noqa: E402
from themis_census import load_clean, detrend  # noqa: E402
from astropy.timeseries import LombScargle     # noqa: E402

K_HARM = 2

def design(t, ens, P=None):
    centers = (ens["grid_edges"][:-1] + ens["grid_edges"][1:]) / 2.0
    eig = ens["eigvectors"][:5]
    comps = np.vstack([np.interp(t, centers, e, left=np.nan, right=np.nan) for e in eig])
    ok = np.all(np.isfinite(comps), axis=0)
    x = (t - t.mean()) / max(t.max() - t.min(), 1e-9)
    cols = [np.ones_like(t), x, x**2, x**3] + [c for c in comps]
    if P is not None:
        w = 2 * np.pi * 24.0 / P
        for j in range(1, K_HARM + 1):
            cols += [np.cos(j * w * t), np.sin(j * w * t)]
    A = np.vstack(cols).T
    return A, ok

def fit_bic(A, y, ok):
    coef, *_ = np.linalg.lstsq(A[ok], y[ok], rcond=None)
    r = y[ok] - A[ok] @ coef
    n, k = ok.sum(), A.shape[1]
    rss = float(np.sum(r**2))
    bic = n * np.log(max(rss / n, 1e-12)) + k * np.log(n)
    return coef, bic

def joint_amp(t, m, ens, P):
    """Return (fitted peak-to-peak amp of the Fourier component, dBIC vs no-signal)."""
    A1, ok = design(t, ens, None)
    A2, _ = design(t, ens, P)
    if ok.sum() < 100:
        return np.nan, np.nan
    _, bic1 = fit_bic(A1, m, ok)
    coef2, bic2 = fit_bic(A2, m, ok)
    w = 2 * np.pi * 24.0 / P
    sig = np.zeros_like(t)
    for j in range(1, K_HARM + 1):
        b, c = coef2[-2 * (K_HARM - j + 1)], coef2[-2 * (K_HARM - j + 1) + 1]
        sig += b * np.cos(j * w * t) + c * np.sin(j * w * t)
    return float(sig.max() - sig.min()), float(bic1 - bic2)   # dBIC>0 favors signal

# ---------- (A) long-P injections: joint vs sequential ----------
print("(A) long-P injections: joint vs sequential amplitude recovery", flush=True)
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
seen = set(); subs = [x for x in subs if not (x[0] in seen or seen.add(x[0]))][:50]
rng = np.random.default_rng(23)
out = []
for f, s in subs:
    try:
        ens = sd.build_eigenbasis(s, cache_dir=sd.DEFAULT_CACHE_DIR)
        t, m, e, _ = load_clean(f)
        if len(t) < 300:
            continue
        bh = (t.max() - t.min()) * 24.0
        for rat in (0.30, 0.50, 0.70, 0.85):
            P = rat * bh
            for A0 in (0.10, 0.30):
                m_inj = m + (A0 / 2) * np.sin(2 * np.pi * t / (P / 24) + rng.uniform(0, 2 * np.pi))
                # sequential: project basis, then measure folded amp
                fit = sd.fit_and_subtract(t, m_inj, ens, 5)
                y1 = detrend(t, fit["corrected"])
                ph = (t / (P / 24)) % 1; nb = 25
                idx = np.clip((ph * nb).astype(int), 0, nb - 1)
                md = np.array([np.median(y1[idx == b]) if (idx == b).sum() > 2 else np.nan for b in range(nb)])
                amp_seq = float(np.nanpercentile(md, 95) - np.nanpercentile(md, 5))
                # joint
                amp_j, dbic = joint_amp(t, m_inj, ens, P)
                out.append(dict(file=os.path.basename(f), sector=s, ratio=rat, amp_in=A0,
                                amp_seq=round(amp_seq, 4), amp_joint=round(amp_j, 4), dbic=round(dbic, 1)))
    except Exception:
        continue
with open(os.path.join(NG, "v7_joint_injections.csv"), "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
print(f"(A) {len(out)} comparisons", flush=True)

lines = ["# V7 results: joint-fit prototype vs sequential projection\n",
         "## (A) Injected-amplitude recovery vs P/baseline (amp_out/amp_in, median)\n",
         "| P/baseline | amp_in | n | sequential | joint |", "|--|--|--|--|--|"]
for rat in (0.30, 0.50, 0.70, 0.85):
    for A0 in (0.10, 0.30):
        sel = [o for o in out if o["ratio"] == rat and o["amp_in"] == A0 and np.isfinite(o["amp_joint"])]
        if not sel:
            continue
        rs = np.median([o["amp_seq"] / A0 for o in sel]); rj = np.median([o["amp_joint"] / A0 for o in sel])
        lines.append(f"| {rat:.2f} | {A0} | {len(sel)} | {rs:.2f} | {rj:.2f} |")

# ---------- (B) nulls: false-clear rate ----------
print("(B) nulls", flush=True)
n_clear = 0; n_tot = 0
for r in csv.DictReader(open(os.path.join(NG, "v3_null_test.csv"))):
    hits = glob.glob(os.path.join(NG, "tess-*", r["file"]))
    if not hits:
        continue
    try:
        s = int(r["sector"]); P = float(r["P"])
        ens = sd.build_eigenbasis(s, cache_dir=sd.DEFAULT_CACHE_DIR)
        t, m, e, _ = load_clean(hits[0])
        amp_j, dbic = joint_amp(t, m, ens, P)
        n_tot += 1
        if np.isfinite(dbic) and dbic > 10:
            n_clear += 1
    except Exception:
        continue
lines.append(f"\n## (B) Null test: joint fit crediting pure comb detections (dBIC>10): {n_clear}/{n_tot} "
             f"({100*n_clear/max(n_tot,1):.0f}% false-clear)")

# ---------- (C) real slow rotators ----------
print("(C) real slow rotators", flush=True)
lines.append("\n## (C) Real slow rotators: sequential drop vs joint dBIC")
lines.append("| object | sector | P (h) | sequential drop | joint amp | joint dBIC |")
lines.append("|--|--|--|--|--|--|")
CASES = [("tess-beltwide6/lc_3396_s27.csv", 27, 61.963), ("tess-beltwide6/lc_3396_s42.csv", 42, 61.963),
         ("tess-beltwide6/lc_3396_s44.csv", 44, 61.963), ("tess-beltwide/lc_2869_s18.csv", 18, 155.234)]
for f, s, P in CASES:
    fp = os.path.join(NG, f)
    if not os.path.exists(fp) or s not in CACHED:
        continue
    try:
        ens = sd.build_eigenbasis(s, cache_dir=sd.DEFAULT_CACHE_DIR)
        t, m, e, _ = load_clean(fp)
        y0 = detrend(t, m)
        p0 = float(LombScargle(t, y0).power(24 / P))
        fit = sd.fit_and_subtract(t, m, ens, 5)
        y1 = detrend(t, fit["corrected"])
        p1 = float(LombScargle(t, y1).power(24 / P))
        drop = (p0 - p1) / max(p0, 1e-9)
        amp_j, dbic = joint_amp(t, m, ens, P)
        lines.append(f"| {os.path.basename(f)} | s{s} | {P} | {drop*100:+.0f}% | {amp_j:.3f} | {dbic:.0f} |")
    except Exception as ex:
        lines.append(f"| {os.path.basename(f)} | s{s} | {P} | ERROR {type(ex).__name__} | | |")
open(os.path.join(NG, "V7_RESULTS.md"), "w").write("\n".join(lines) + "\n")
print("V7_RESULTS.md written", flush=True)
