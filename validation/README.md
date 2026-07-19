# Validation campaign (V1-V7)

Full validation of the de-comb method, run 2026-07-19 against the survey's
extracted light-curve base (~600 TESS sectors x objects). Results documents
are the V*_RESULTS.md files; raw artifacts are the CSVs.

| Stage | Question | Files |
|--|--|--|
| V1 | Is K = 5 justified? | V1_RESULTS.md, v1_holdout_k_sectors.csv, v1_verdict_stability.csv |
| V2 | Injection/recovery at scale (17,372 injections) | V2_RESULTS.md, v2_injections.csv |
| V3 | Null test: pure comb detections | V2_RESULTS.md (section), v3_null_test.csv |
| V4 | ROC: where do the thresholds sit? | V4_V5_RESULTS.md |
| V5 | Ensemble-composition bootstrap | V4_V5_RESULTS.md, v5_bootstrap.csv |
| V6 | Validity domain vs period/baseline | V6_RESULTS.md, v6_validity.csv, v6_validity_curve.png |
| V7 | Joint (simultaneous) fit prototype | V7_RESULTS.md, v7_joint_injections.csv |

`scripts/` contains the exact campaign scripts for transparency. They expect
the survey's internal directory layout (`tess-<family>/lc_<num>_s<sector>.csv`
plus per-family census tables) and private loader modules, so they are NOT
runnable from this repository alone; the algorithms they implement are all
present in the released `tess_decomb` package. Column meanings are documented
in the results files.

Headline numbers and the two design consequences (kill/survive asymmetry,
0.45x-baseline validity domain) are summarized in the top-level README.
