# V2/V3 results (120 tracks, 17372 injections, 46 null tests)

## V2: injected REAL signals (false-kill = injected signal wrongly killed)

| kind | amp | n | survived% | inconclusive% | FALSE-KILL% | median recovery pw1/pw0 |
|--|--|--|--|--|--|--|
| tooth | 0.05 | 1607 | 65.5 | 31.5 | 2.99 | 0.962 |
| tooth | 0.1 | 1607 | 71.3 | 27.0 | 1.74 | 0.955 |
| tooth | 0.2 | 1607 | 77.3 | 22.2 | 0.50 | 0.957 |
| tooth | 0.4 | 1607 | 85.4 | 14.3 | 0.31 | 0.961 |
| beside | 0.05 | 1602 | 64.6 | 33.1 | 2.25 | 0.942 |
| beside | 0.1 | 1602 | 71.4 | 27.3 | 1.31 | 0.954 |
| beside | 0.2 | 1602 | 79.3 | 20.1 | 0.56 | 0.955 |
| beside | 0.4 | 1602 | 86.8 | 13.2 | 0.00 | 0.961 |
| control | 0.05 | 1134 | 70.4 | 27.8 | 1.85 | 0.977 |
| control | 0.1 | 1134 | 75.5 | 23.5 | 1.06 | 0.972 |
| control | 0.2 | 1134 | 81.4 | 18.2 | 0.44 | 0.971 |
| control | 0.4 | 1134 | 88.5 | 11.4 | 0.09 | 0.977 |

## V3: pure comb systematics (real pre-injection detections at teeth, n=46)
- killed 6.5% | inconclusive 32.6% | FALSE-SURVIVE 60.9%

## Interpretation (the asymmetry finding)
- V2: the KILL verdict is trustworthy: false-kill <=3.0% at amp 0.05 (weakest), <=0.5% at amp >=0.2,
  ~0% at amp 0.4; median power recovery 94-98% everywhere, INCLUDING injections directly on comb teeth.
- V3: SURVIVAL is weak evidence for on-tooth periods: 61% of tooth-frequency detections in period-less
  tracks survive projection (6.5% killed). The eigen-basis removes the ensemble-SHARED comb component;
  track-specific comb power (the moving aperture's individual path through the scattered-light pattern)
  often remains. Caveats: n=46, and ground truth is impure (some "null" detections may be weak real
  signals in NONE-tier objects), so 61% is an upper bound on the false-survive rate.
- NET: the method is ASYMMETRIC: a high-confidence killer and a conservative clearer. Catalog practice
  already conforms (no object was promoted on decomb-survival alone; on-tooth claims require amplitude
  clearing + multi-sector/multi-instrument evidence). The V5 local-ensemble test is the natural upgrade
  path for the false-survive rate.
