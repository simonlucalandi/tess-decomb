# V7 results: joint-fit prototype vs sequential projection

## (A) Injected-amplitude recovery vs P/baseline (amp_out/amp_in, median)

| P/baseline | amp_in | n | sequential | joint |
|--|--|--|--|--|
| 0.30 | 0.1 | 47 | 1.25 | 1.12 |
| 0.30 | 0.3 | 47 | 0.90 | 1.03 |
| 0.50 | 0.1 | 47 | 1.24 | 1.29 |
| 0.50 | 0.3 | 47 | 0.78 | 1.00 |
| 0.70 | 0.1 | 47 | 1.16 | 1.61 |
| 0.70 | 0.3 | 47 | 0.63 | 1.20 |
| 0.85 | 0.1 | 47 | 1.33 | 3.83 |
| 0.85 | 0.3 | 47 | 0.49 | 1.57 |

## (B) Null test: joint fit crediting pure comb detections (dBIC>10): 44/46 (96% false-clear)

## (C) Real slow rotators: sequential drop vs joint dBIC
| object | sector | P (h) | sequential drop | joint amp | joint dBIC |
|--|--|--|--|--|--|
| lc_3396_s27.csv | s27 | 61.963 | +28% | 0.053 | 443 |
| lc_3396_s42.csv | s42 | 61.963 | +25% | 0.078 | 411 |
| lc_3396_s44.csv | s44 | 61.963 | +41% | 0.061 | 223 |
| lc_2869_s18.csv | s18 | 155.234 | +21% | 0.582 | 1298 |

## Interpretation (honest, mixed)
- **(A) The joint fit ELIMINATES the slow-rotator over-subtraction in-domain:** at real
  amplitudes (0.3 mag), sequential recovery degrades 0.90 -> 0.49 across P/baseline
  0.3 -> 0.85, while joint recovery stays 1.00-1.20 up to ~0.7. BUT it overshoots in the
  few-cycle / low-amplitude limit (3.8x at amp 0.1, ratio 0.85): with < ~1.5 cycles the
  Fourier terms go quasi-degenerate with the trend + basis and absorb noise as "signal".
- **(B) The joint fit is NOT an artifact discriminator: 96% false-clear on pure comb
  nulls.** Expected in hindsight: the basis does not span track-specific comb power
  (V3), so Fourier terms legitimately explain that residual variance and dBIC rewards
  them. The joint fit cannot arbitrate on-tooth claims.
- **(C) On the real over-subtracted slow rotators (3396 sectors, 2869 s18) the joint fit
  decisively recovers the signal** (dBIC 223-1298) with amplitudes consistent with the
  adopted folds.
- **NET / production architecture:** the two tools are complementary, not substitutes.
  KEEP sequential projection as the KILL test (V4-calibrated, 1.1% false-kill); ADD the
  joint fit as the amplitude/significance ESTIMATOR for slow rotators within
  P/baseline <= ~0.7 (replacing raw "inconclusive" with a measured amplitude + dBIC).
  "Sequential kills; joint measures." This is the PASP paper's central methods result.
