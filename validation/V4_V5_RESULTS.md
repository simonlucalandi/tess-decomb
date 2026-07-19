# V4 + V5 results

## V4: verdict-threshold VERIFICATION (ROC on V2 real vs V3 null). The thresholds (KILL drop>=0.50 & R2>=0.30; SURVIVE drop<=0.15) were fixed 2026-07-17 in decomb_check.py, BEFORE the V2/V3 campaign existed: the ROC verifies the pre-registered operating point sits near the knee; it did not tune it.

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

## V4 addendum (2026-07-19, post adversarial review): track-disjoint held-out check
To close the train-on-test objection: tracks were split deterministically (md5 of
filename) into two disjoint halves and the FIXED pre-registered thresholds evaluated
on each half independently. No quantity was re-tuned.

| split | injections | tracks | false-kill | median recovery | by-amp false-kill (0.05/0.1/0.2/0.4) |
|--|--|--|--|--|--|
| A | 9184 | 63 | 0.63% | 0.966 | 1.74% / 0.39% / 0.26% / 0.13% |
| B | 8188 | 57 | 1.66% | 0.957 | 3.18% / 2.54% / 0.78% / 0.15% |

V3 nulls by split: A 12/20 survive (60%), B 16/26 (62%) vs 61% on the full set.
All headline rates reproduce on tracks the thresholds never saw; the amplitude trend
is intact in both halves. Conservative note wording: false-kill "falling below ~1%
for amplitudes >= 0.2 mag" (split B shows 0.78% at 0.2 mag; the earlier "<0.5%"
held only for the 0.2+0.4 pooled subset).
