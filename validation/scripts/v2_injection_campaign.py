#!/usr/bin/env python3
"""V2 + V3 (VALIDATION_PLAN.md): injection-recovery campaign + null test.

V2: inject sinusoids into real NONE-tier moving-target tracks (real noise, no detection):
periods ON each momentum-dump comb tooth (328.8/n, n=2..15), BESIDE each tooth (+5%),
and off-comb controls; amplitudes 0.05-0.4 mag. Measure post-decomb recovery fraction,
amplitude bias, and the false-kill rate of the verdict thresholds.

V3 (same tracks, pre-injection): at every comb tooth where the RAW track shows a real
detection (Baluev FAP < 1e-3), apply the verdict machinery: how often does pure
systematics SURVIVE? (false-save rate).

Fully offline (cached ensembles + lc CSVs). Parallel over tracks.
Outputs: v2_injections.csv, v3_null_test.csv, V2_RESULTS.md.
"""
import csv, glob, os, re, sys
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

NG = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(NG, "tess-tools"))
sys.path.insert(0, os.path.join(NG, "tess-themis"))
import sysrem_decomb as sd            # noqa: E402
from themis_census import load_clean, detrend   # noqa: E402
from astropy.timeseries import LombScargle      # noqa: E402

TEETH = [328.8 / n for n in range(2, 16)]
BESIDE = [p * 1.05 for p in TEETH]
CONTROLS = [2.5, 4.0, 7.0, 12.0, 18.0, 35.0, 50.0, 75.0, 120.0, 200.0]
AMPS = [0.05, 0.10, 0.20, 0.40]
KILL_DROP, KILL_R2, SURVIVE_DROP = 0.50, 0.30, 0.15
CACHED = {int(re.search(r"_s0*(\d+)_", f).group(1))
          for f in glob.glob(os.path.join(NG, "tess-tools/sysrem_cache/ensemble_s*.npz"))}


def pick_substrates(max_tracks=120):
    """NONE-tier objects' longest track per object, sector must have a cached ensemble."""
    subs = []
    for cens in glob.glob(os.path.join(NG, "tess-*/census_objects.csv")):
        d = os.path.dirname(cens)
        try:
            rows = list(csv.DictReader(open(cens)))
        except Exception:
            continue
        for r in rows:
            if (r.get("tier") or "") != "NONE":
                continue
            best = None
            for f in glob.glob(os.path.join(d, f"lc_{r['num']}_s*.csv")):
                s = int(re.search(r"_s(\d+)\.csv", f).group(1))
                if s not in CACHED:
                    continue
                n = sum(1 for _ in open(f)) - 1
                if best is None or n > best[2]:
                    best = (f, s, n)
            if best and best[2] >= 300:
                subs.append(best[:2])
    # de-dup by file, spread across sectors
    seen = set(); out = []
    for f, s in subs:
        if f in seen:
            continue
        seen.add(f); out.append((f, s))
    out.sort(key=lambda x: x[1])
    return out[:max_tracks]


def pw_fap(t, y, P):
    ls = LombScargle(t, y)
    p = float(ls.power(24.0 / P))
    try:
        fap = float(ls.false_alarm_probability(p, method="baluev"))
    except Exception:
        fap = np.nan
    return p, fap


def folded_amp(t, y, P, nb=25):
    ph = (t / (P / 24.0)) % 1.0
    idx = np.clip((ph * nb).astype(int), 0, nb - 1)
    md = np.array([np.median(y[idx == b]) if (idx == b).sum() > 2 else np.nan for b in range(nb)])
    return float(np.nanpercentile(md, 95) - np.nanpercentile(md, 5))


def one_track(args):
    f, sector = args
    rowsV2, rowsV3 = [], []
    try:
        ens = sd.build_eigenbasis(sector, cache_dir=sd.DEFAULT_CACHE_DIR)
        t, m, e, _ = load_clean(f)
        if len(t) < 300:
            return rowsV2, rowsV3
        base_h = (t.max() - t.min()) * 24.0
        rng = np.random.default_rng(abs(hash(os.path.basename(f))) % (2**32))
        # ---- V3 null test (pre-injection) ----
        y0 = detrend(t, m)
        for P in TEETH:
            if P > 0.45 * base_h:
                continue
            p0, fap0 = pw_fap(t, y0, P)
            if not (np.isfinite(fap0) and fap0 < 1e-3 and p0 >= 0.03):
                continue   # no real comb detection here
            fit = sd.fit_and_subtract(t, m, ens, 5)
            ydc = detrend(t, fit["corrected"])
            p1, _ = pw_fap(t, ydc, P)
            drop = (p0 - p1) / p0 if p0 > 0 else np.nan
            r2 = float(fit["r2_systematics"])
            verdict = ("killed" if (drop >= KILL_DROP and r2 >= KILL_R2)
                       else "survived" if drop <= SURVIVE_DROP else "inconclusive")
            rowsV3.append(dict(file=os.path.basename(f), sector=sector, P=round(P, 3),
                               pw0=round(p0, 4), pw1=round(p1, 4), drop=round(drop, 3),
                               r2=round(r2, 3), verdict=verdict))
        # ---- V2 injections ----
        for P in TEETH + BESIDE + CONTROLS:
            if P > 0.45 * base_h or P < 1.6:
                continue
            kind = ("tooth" if P in TEETH else "beside" if P in BESIDE else "control")
            for A in AMPS:
                phi = rng.uniform(0, 2 * np.pi)
                m_inj = m + (A / 2.0) * np.sin(2 * np.pi * t / (P / 24.0) + phi)
                y_in = detrend(t, m_inj)
                p0, fap0 = pw_fap(t, y_in, P)
                fit = sd.fit_and_subtract(t, m_inj, ens, 5)
                ydc = detrend(t, fit["corrected"])
                p1, _ = pw_fap(t, ydc, P)
                r2 = float(fit["r2_systematics"])
                drop = (p0 - p1) / p0 if p0 > 0 else np.nan
                verdict = ("killed" if (drop >= KILL_DROP and r2 >= KILL_R2)
                           else "survived" if drop <= SURVIVE_DROP else "inconclusive")
                a_out = folded_amp(t, ydc, P)
                rowsV2.append(dict(file=os.path.basename(f), sector=sector, kind=kind,
                                   P=round(P, 3), amp_in=A, pw0=round(p0, 4), pw1=round(p1, 4),
                                   drop=round(drop, 3), r2=round(r2, 3), verdict=verdict,
                                   amp_out=round(a_out, 4), fap0=fap0))
    except Exception as ex:
        rowsV2.append(dict(file=os.path.basename(f), sector=sector, kind=f"ERROR:{type(ex).__name__}",
                           P="", amp_in="", pw0="", pw1="", drop="", r2="", verdict="", amp_out="", fap0=""))
    return rowsV2, rowsV3


def main():
    subs = pick_substrates()
    print(f"substrates: {len(subs)} NONE-tier tracks across {len({s for _, s in subs})} sectors", flush=True)
    from multiprocessing import Pool
    allV2, allV3 = [], []
    with Pool(20) as pool:
        for i, (r2, r3) in enumerate(pool.imap_unordered(one_track, subs)):
            allV2.extend(r2); allV3.extend(r3)
            if (i + 1) % 10 == 0:
                print(f"  {i+1}/{len(subs)} tracks", flush=True)
    with open(os.path.join(NG, "v2_injections.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(allV2[0].keys())); w.writeheader(); w.writerows(allV2)
    if allV3:
        with open(os.path.join(NG, "v3_null_test.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(allV3[0].keys())); w.writeheader(); w.writerows(allV3)

    # ---- summary ----
    ok = [r for r in allV2 if r["verdict"]]
    import collections
    def rate(rows, v): return 100.0 * sum(1 for r in rows if r["verdict"] == v) / max(len(rows), 1)
    lines = [f"# V2/V3 results ({len(subs)} tracks, {len(ok)} injections, {len(allV3)} null tests)\n"]
    lines.append("## V2: injected REAL signals (false-kill = injected signal wrongly killed)\n")
    lines.append("| kind | amp | n | survived% | inconclusive% | FALSE-KILL% | median recovery pw1/pw0 |")
    lines.append("|--|--|--|--|--|--|--|")
    for kind in ("tooth", "beside", "control"):
        for A in AMPS:
            sel = [r for r in ok if r["kind"] == kind and r["amp_in"] == A]
            if not sel:
                continue
            rec = np.median([r["pw1"] / r["pw0"] for r in sel if r["pw0"] and r["pw0"] > 0])
            lines.append(f"| {kind} | {A} | {len(sel)} | {rate(sel,'survived'):.1f} | "
                         f"{rate(sel,'inconclusive'):.1f} | {rate(sel,'killed'):.2f} | {rec:.3f} |")
    if allV3:
        k3 = rate(allV3, "killed"); s3 = rate(allV3, "survived"); i3 = rate(allV3, "inconclusive")
        lines.append(f"\n## V3: pure comb systematics (real pre-injection detections at teeth, n={len(allV3)})")
        lines.append(f"- killed {k3:.1f}% | inconclusive {i3:.1f}% | FALSE-SURVIVE {s3:.1f}%")
    open(os.path.join(NG, "V2_RESULTS.md"), "w").write("\n".join(lines) + "\n")
    print("V2_RESULTS.md written", flush=True)


if __name__ == "__main__":
    main()
