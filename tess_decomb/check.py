#!/usr/bin/env python3
"""Harmonic-aware de-comb verdict for a candidate asteroid rotation period.

For each supplied (light curve, sector) pair, every harmonic in {P/2, P, 2P}
that is a REAL per-sector detection (Baluev false-alarm probability below
DETECT_FAP) is de-combed with the K=5 eigen-systematics model and scored:

  survives      power drop <= 15%
  killed        power drop >= 50% AND systematics-model R^2 >= 0.3
  inconclusive  in between

The verdict is anchored on the object's STRONGEST detection (max pre-de-comb
power), survival-first across sectors (a clean sector's surviving period is
not condemned by one instrumental sector):

  survived(...)           primary detection preserved -> astrophysical
  KILLED(...)             all real detections collapse -> instrumental
  review-degraded(...)    the strongest detection is systematics; only a
                          weaker one survives -> human review
  inconclusive(...)       primary drop in the 15-50% gray zone

IMPORTANT (validated limits; see validation/ in this repository):
- The KILL side is the calibrated, trustworthy side (1.1% false-kill at 6.5%
  of injections killed, ROC knee). SURVIVAL of an on-tooth period is WEAK
  evidence: 61% of pure comb detections survive projection. Survival gains
  credibility; it does not certify.
- Validity domain: signal retention >= 88% and false-kill <= 3% for periods
  up to 0.45x the sector baseline; degraded beyond.

Usage:
  tess-decomb-check --period 25.876 --lc path/to/lc_s94.csv:94 [--lc f.csv:36]
"""
from __future__ import annotations

import argparse

import numpy as np

from . import sysrem as sd
from .lightcurve import load_clean, detrend

KILL_DROP, KILL_R2, SURVIVE_DROP = 0.50, 0.30, 0.15
RANGE_SKIP_MAG = 5.0
HARMONICS = (0.5, 1.0, 2.0)
DETECT_FAP = 1e-3   # Baluev FAP gate: normalises for N, unlike a raw-power floor


def ls_power_at(t, y, P_h):
    from astropy.timeseries import LombScargle
    return float(LombScargle(t, y).power(np.array([24.0 / P_h]))[0])


def ls_power_fap(t, y, P_h):
    from astropy.timeseries import LombScargle
    ls = LombScargle(t, y)
    pw = float(ls.power(24.0 / P_h))
    try:
        fap = float(ls.false_alarm_probability(pw, method="baluev"))
    except Exception:
        fap = np.nan
    return pw, fap


def check_period(pairs, P_phot_h, k=5, cache_dir=None):
    """pairs: list of (lc_path, sector). Returns the verdict string."""
    cache_dir = cache_dir or sd.DEFAULT_CACHE_DIR
    dets = []   # (pw0, sector, P_test, drop, r2, state)
    any_ensemble = False
    for f, sector in pairs:
        t, m, e, _ = load_clean(f)
        if len(t) < 50 or (m.max() - m.min()) > RANGE_SKIP_MAG:
            continue
        y = detrend(t, m)
        real = []
        for h in HARMONICS:
            pw0, fap = ls_power_fap(t, y, h * P_phot_h)
            if np.isfinite(fap) and fap < DETECT_FAP:
                real.append((h * P_phot_h, pw0))
        if not real:
            continue
        try:
            run = sd.decomb_asteroid(f, sector, k, cache_dir=cache_dir)
            y_dc = detrend(run["t"], run["m_decomb"])
            r2 = float(run["r2_systematics"])
            any_ensemble = True
        except Exception:
            continue
        for P_test, pw0 in real:
            drop = (pw0 - ls_power_at(run["t"], y_dc, P_test)) / pw0
            state = ("KILLED" if (drop >= KILL_DROP and r2 >= KILL_R2)
                     else "survived" if drop <= SURVIVE_DROP else "inconc")
            dets.append((pw0, sector, P_test, drop, r2, state))
    if not dets:
        return "ensemble-unavailable" if not any_ensemble else "no-data"
    dets.sort(reverse=True)                      # strongest detection first
    _, s, P, drop, r2, state = dets[0]           # the primary
    tag = f"dPW={drop*100:+.0f}%,R2={r2:.2f},s{s}@{P:.3f}h"
    if state == "survived":
        return f"survived({tag},{sum(1 for d in dets if d[5]=='survived')}sec)"
    if state == "KILLED":
        if any(d[5] == "survived" for d in dets):
            return f"review-degraded(primary-{tag}; weaker signal survives)"
        return f"KILLED({tag})"
    return f"inconclusive({tag})"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--period", type=float, required=True,
                    help="candidate photometric period (hours)")
    ap.add_argument("--lc", action="append", required=True, metavar="PATH:SECTOR",
                    help="light-curve CSV and its TESS sector, e.g. lc.csv:94 (repeatable)")
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--cache-dir", default=None)
    a = ap.parse_args()
    pairs = []
    for item in a.lc:
        path, _, sec = item.rpartition(":")
        pairs.append((path, int(sec)))
    print(check_period(pairs, a.period, k=a.k, cache_dir=a.cache_dir))


if __name__ == "__main__":
    main()
