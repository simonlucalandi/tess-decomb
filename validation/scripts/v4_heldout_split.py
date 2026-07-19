#!/usr/bin/env python3
"""V4 addendum: track-disjoint held-out check of the FIXED pre-registered verdict
thresholds (KILL drop>=0.50 & R2>=0.30; SURVIVE drop<=0.15; fixed 2026-07-17,
before the V2/V3 campaign). Tracks are split deterministically (md5 of filename)
into two disjoint halves; each half is evaluated independently at the fixed
thresholds. Nothing is re-tuned. Results in V4_V5_RESULTS.md (addendum table).
Runs against v2_injections.csv and v3_null_test.csv in the parent directory."""
import csv, hashlib, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
V = os.path.join(HERE, "..")

def half(f):
    return 'A' if int(hashlib.md5(f.encode()).hexdigest(), 16) % 2 == 0 else 'B'

rows = [r for r in csv.DictReader(open(os.path.join(V, 'v2_injections.csv'))) if r['verdict']]
for h in 'AB':
    sel = [r for r in rows if half(r['file']) == h]
    fk = 100 * sum(r['verdict'] == 'killed' for r in sel) / len(sel)
    rec = np.median([float(r['pw1']) / max(float(r['pw0']), 1e-9) for r in sel])
    byamp = {}
    for r in sel:
        byamp.setdefault(float(r['amp_in']), []).append(r['verdict'] == 'killed')
    amps = ', '.join(f"{a}:{100*np.mean(v):.2f}%" for a, v in sorted(byamp.items()))
    print(f"split {h}: n={len(sel)} tracks={len(set(r['file'] for r in sel))} "
          f"false-kill={fk:.2f}% median-recovery={rec:.3f} | by-amp: {amps}")
nulls = list(csv.DictReader(open(os.path.join(V, 'v3_null_test.csv'))))
for h in 'AB':
    sel = [r for r in nulls if half(r['file']) == h]
    if sel:
        sv = 100 * sum(r['verdict'] == 'survived' for r in sel) / len(sel)
        print(f"nulls split {h}: n={len(sel)} survive={sv:.0f}%")
