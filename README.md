# tess-decomb

Field-star eigen-systematics "de-comb" for TESS moving-target (asteroid)
photometry: separate real rotation signals from the TESS momentum-dump alias
comb (periods at 328.8/n h) and scattered-light systematics, by projecting
each light curve onto an eigen-basis learned from ~200 ordinary field stars
observed in the same sector.

Companion code to the RNAAS research note (Landi 2026, DOI to be added) and
the associated rotation-period survey papers.

## Method in one paragraph

SPOC 2-min SAP light curves of ~200 field stars from the target's sector are
median-binned onto a common 30-min grid; the SVD of the (stars x time) matrix
yields the sector's shared instrumental time series ("eigen-systematics").
An asteroid's unbinned light curve is then fitted as offset + a linear
combination of the top K = 5 eigen-vectors interpolated to its own
timestamps, and the systematics part is subtracted. A candidate period is
judged harmonic-aware ({P/2, P, 2P}, each gated by Baluev FAP < 1e-3):
a detection whose Lomb-Scargle power drops >= 50% under projection (with
systematics-model R^2 >= 0.3) is instrumental; a drop <= 15% survives.
Verdicts are anchored on the strongest detection and survival-first across
sectors.

## Validated performance (full campaign in `validation/`)

| Test | Result |
|--|--|
| V1 holdout-K | verdicts plateau for K in 2-8; ground-truth cases bracket K in [4, 5] |
| V2 17,372 injections | false-kill 3% -> 0% with amplitude; recovery 94-98%, incl. on-tooth periods |
| V3 null test | 61% of pure comb detections SURVIVE -> survival is weak evidence; the kill side is the calibrated side |
| V4 ROC | kill thresholds sit at the knee: 1.13% false-kill at 6.5% of injections killed; no survive threshold separates the populations |
| V5 bootstrap | 77% verdict stability across 6 resampled ensembles, 89% disjoint-halves agreement, zero confirmed-object flips |
| V6 validity domain | retention >= 88% and false-kill <= 3% up to P = 0.45x sector baseline; ~70% beyond |
| V7 joint fit | simultaneous basis+Fourier fit removes slow-rotator over-subtraction in-domain, but 96% false-clear on comb nulls: architecture is "sequential kills; joint measures" |

Two design consequences worth internalizing before use:
1. **Asymmetry.** Use this tool to KILL comb artifacts, not to certify on-tooth
   survivors. Survival gains credibility; it does not prove a period.
2. **Validity domain.** Trust verdicts for periods up to ~0.45x the sector
   baseline (~250 h for a full sector); beyond that, prefer multi-sector
   phase coherence.

## Install

```bash
pip install numpy pandas astropy astroquery requests
pip install -e .
```

## Usage

```bash
# one-time per sector: build the eigen-basis from MAST (public SPOC data)
tess-decomb ensemble --sector 18

# de-comb a light curve (CSV: time [BTJD days], mag, err)
tess-decomb decomb --lc my_asteroid_s18.csv --sector 18 --k 5 --out decombed.csv

# harmonic-aware verdict for a candidate period across sectors
tess-decomb-check --period 25.876 --lc lc_s94.csv:94 --lc lc_s36.csv:36

# self-contained smoke test (synthetic injection on a held-out field star)
tess-decomb validate --sector 18
```

Python API: `build_eigenbasis`, `fit_and_subtract`, `decomb_asteroid`
(`tess_decomb.sysrem`), `check_period` (`tess_decomb.check`).

The eigen-basis cache defaults to `./sysrem_cache` (override with the
`TESS_DECOMB_CACHE` environment variable). Building one sector downloads
~220 SPOC light curves (~8 min) once; everything after is offline.

## Repository layout

- `tess_decomb/` -- the module (eigen-basis, projection, verdict CLI).
- `validation/` -- the V1-V7 validation campaign: results documents, raw
  artifacts (17,372-injection table, null tests, bootstrap, validity curve),
  and the campaign scripts (archival: they ran against the survey's internal
  directory layout and are included for transparency, not as runnable tools).

## Citing

If you use this code, cite the RNAAS note (DOI to be added on publication)
and this repository (see `CITATION.cff`). The method builds on SysRem
(Tamuz, Mazeh & Zucker 2005) and the TFA/CBV family of ensemble-systematics
approaches (Kovacs et al. 2005; Smith et al. 2012), applied to moving-target
photometry extracted with `tess-asteroids` (Tuson et al. 2025).

## License

MIT (see LICENSE).
