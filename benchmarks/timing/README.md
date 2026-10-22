# Expanded validation runtime after optimization

Measurement on 2026-10-22: Python 3.12.12, 240 families, 9,999 permutations, seven scenarios and bins 2/4/8. Each version ran 21 analyses, with OR=1 and replicate 0.

The sequential estimate for the full configuration of 63,000 analyses fell from **73.0 to 34.7 hours**, a **2.10×** speedup. Cumulative analysis time fell from 87.6 to 41.7 seconds.

Contributions from fully selected strata are now precomputed; only partially selected strata use the random number generator. Sampling without replacement and shared draws across terms are preserved. The sequence for a given seed may change relative to the previous implementation. Code hashes distinguish runs; previous checkpoints require a separate directory.

| Scenario | Bins | Before (s) | After (s) | Fraction of informative strata |
|---|---:|---:|---:|---:|
| exact | 2 | 1.74 | 1.42 | 1.000 |
| exact | 4 | 1.47 | 1.25 | 1.000 |
| exact | 8 | 1.46 | 1.25 | 1.000 |
| continuous | 2 | 1.46 | 1.22 | 1.000 |
| continuous | 4 | 5.82 | 4.17 | 0.603 |
| continuous | 8 | 10.82 | 2.16 | 0.093 |
| missing_size | 2 | 1.43 | 1.20 | 1.000 |
| missing_size | 4 | 5.83 | 4.26 | 0.619 |
| missing_size | 8 | 10.24 | 3.28 | 0.102 |
| dependent | 2 | 1.51 | 1.30 | 1.000 |
| dependent | 4 | 1.55 | 1.30 | 1.000 |
| dependent | 8 | 1.60 | 1.30 | 1.000 |
| sparse | 2 | 1.49 | 1.08 | 0.750 |
| sparse | 4 | 6.20 | 2.53 | 0.333 |
| sparse | 8 | 13.53 | 1.09 | 0.037 |
| derived_degree | 2 | 1.44 | 1.24 | 1.000 |
| derived_degree | 4 | 5.81 | 4.25 | 0.655 |
| derived_degree | 8 | 10.00 | 3.77 | 0.219 |
| unconfounded | 2 | 1.41 | 1.21 | 1.000 |
| unconfounded | 4 | 1.41 | 1.21 | 1.000 |
| unconfounded | 8 | 1.41 | 1.21 | 1.000 |

These are point estimates without uncertainty intervals. OR=2/4, other replicates and machine load may change runtimes. Data generation, checkpoints, aggregation and reports are excluded from the extrapolation. Tests ran partly in parallel with the measurement. This comparison measures computational runtime; the full statistical validation (63,000 analyses) was subsequently executed and consolidated in `benchmark_results/validation/` (see the [benchmark protocol](../../docs/benchmark-protocol.md)).

With eight bins, some scenarios have few informative strata, which should be considered when assessing power and support.

Verification: 44 tests passed, including the exact distribution of a small universe with six possible pairs, shared draws across duplicate terms, and strata with fixed and variable selection.

Reproduce from the repository root:

```bash
PYTHONPATH=. .venv/bin/python benchmarks/timing/measure.py
```

`runtime.json` contains the new timings, configuration, versions and hashes; `baseline/runtime.json` preserves the earlier measurement. `measurements.json` is updated after each analysis. Rerunning replaces the current measurement.
