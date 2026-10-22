# Audit of expanded validation

Audit performed on 2026-10-22 without modifying the original checkpoints or tables in `benchmark_results/validation`.

The results passed structural checks and sampled numerical reproduction. They may be used for exploratory analysis with the provenance caveat below. This audit does not retrospectively certify the origin of every historical task.

## Integrity and reconstruction

- All 4,200 expected batches exist, with correct names and complete, unique coverage of scenario, OR, replicate, bins, method and term.
- Batch copies in the root directory and scenario directories have identical SHA256 hashes.
- The audit checked 126,000 metric rows, corresponding to 63,000 analyses, and 450,000 term-level rows.
- `replicates.tsv`, `terms.tsv`, `summary.tsv` and `calibration.tsv` match reconstruction from the checkpoints, including intervals calculated by the runner. Numerical comparisons use `rtol=1e-8` and `atol=1e-10`.
- The audit repeated 21 analyses with 9,999 permutations: seven scenarios × three OR values, always replicate 0 and two bins. All term-level fields matched within tolerance. The sample does not cover other replicates or four/eight bins.

Slurm accounting confirms all 84 tasks of `846705` and report job `846789` as `COMPLETED`, with exit code `0:0`. The array started at 11:58:29 and ended at 12:37:44, São Paulo time: 39 minutes and 15 seconds. The report finished at 12:38:36. The previous seven-task array, `846696`, was cancelled and supplied checkpoints that were later reused; the final array's runtime does not include all earlier work.

## Historical provenance

The recorded NumPy (2.4.3), pandas (3.0.1) and SciPy (1.17.1) versions match the audited environment. Hashes of `benchmarks/expanded.py` and the numerical modules match the current sources.

Two files differ from the manifest: `scbenrich/__init__.py` and `scbenrich/utils/io.py`. Their historical versions were found in commit `a98e4638a4ec061c09425aac675e51869d159f72`, with SHA256 hashes matching the manifest. The current files incorporate input-preparation functionality added after the record was created. Sampled reproduction used the current package.

The old run has no `runner_metadata.json`: historical orchestration module hashes were not recorded. The seven scenario manifests match the central manifest and each scenario's configuration, but this does not prove that every historical task checked them before execution. No retrospective manifest was created in the original directory.

For a run with full traceability from the outset, use a new directory, prepare manifests with the corrected runner and repeat validation. The current audit supports exploratory use of the existing results but is not equivalent to that new run.

## Observed statistical limitations

The estimates below use OR=1, 1,000 replicates and the stratified method. Intervals are the replicate bootstrap intervals produced by the benchmark, without simultaneous adjustment across scenarios.

| Scenario | Bins | Empirical FDR | 95% bootstrap interval | Mean fraction of informative strata |
|---|---:|---:|---:|---:|
| Exact | 2 | 3.9% | 2.8–5.1% | 99.7% |
| Continuous | 2 | 6.0% | 4.6–7.5% | 100.0% |
| Sparse | 2 | 26.5% | 23.8–29.3% | 80.8% |
| Sparse | 4 | 4.2% | 3.0–5.5% | 28.6% |
| Sparse | 8 | 0.1% | 0.0–0.3% | 3.4% |

Reference truth in these continuous scenarios comes from the generator conditional on continuous covariates. The implemented test conditions on bins. The sparse result with two bins demonstrates a limitation of this approximation and does not support a general guarantee of FDR at 5%. Low rejection with eight bins accompanies a large loss of informative strata and must be interpreted alongside power. Binning must not be chosen retrospectively based only on the lowest observed FDR.

The `derived_degree` scenario retains unknown conditional truth. Its missing FDR cells are intentional.

## Evidence and reproduction

- `audit.json`: checks, sampled reproduction results and caveats.
- `checkpoint_hashes.json`: SHA256 hashes of the 4,200 audited batches.
- `observed-manifest.json`: copy of the original manifest.
- `slurm-accounting.psv`: accounting history retrieved from Slurm.
- `summary.tsv`: copy of the verified summary.

From the repository root, with the original results present:

```bash
PYTHONPATH=. .venv/bin/python docs/validation-audit/reproduce.py
```

The script reconstructs tables in memory and records the audit in this directory. It does not modify the original results. Reproduction depends on the versions and files listed in the manifest; Slurm accounting and the historical source locations were checked separately.
