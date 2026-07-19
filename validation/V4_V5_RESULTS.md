# V4 + V5 results

## V4: verdict-threshold calibration (ROC on V2 real vs V3 null)

| kill-drop threshold | false-kill% (real, R2>=0.3 arm) | true-kill% (null) |
|--|--|--|
| 0.30 | 2.95 | 8.7 |
| 0.40 | 1.91 | 8.7 |
| 0.50 | 1.13 | 6.5 | <- operating point
| 0.60 | 0.71 | 4.3 |
| 0.70 | 0.47 | 4.3 |

| survive-drop threshold | real signals surviving% | null surviving% (false-survive) |
|--|--|--|
| 0.05 | 54.5 | 45.7 |
| 0.10 | 67.9 | 52.2 |
| 0.15 | 76.2 | 60.9 | <- operating point
| 0.20 | 82.3 | 65.2 |
| 0.25 | 86.5 | 69.6 |

Real-signal drop distribution: median +0.038, p95 +0.441; null drops: median +0.071, p95 +0.735 (n=46).

## V5: ensemble-composition robustness (384 object-sector tests, B=6 bootstraps)
- Bootstrap-stable verdicts (identical across all 6 resampled 200-star bases): 297/384 (77.3%)
- Disjoint-half test (two independent ~110-star bases agree): 340/384 (88.5%)
