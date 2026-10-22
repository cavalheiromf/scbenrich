# Benchmark: hypotheses, metrics and interpretation

## Reference truth — schema 2

The generator defines truth from the population distribution, without consulting observed p-values. `term_truth.tsv` exports `marginal_null` and `conditional_null` for each term and scenario. The conditional target is continuous size; the test uses bins and may retain residual confounding within them.

| Scenario | `null_*` terms: marginal null | `null_*` terms: conditional null | `effect` term |
|---|---|---|---|
| simple_null | true | true | absent |
| confounded_null | false | true | absent |
| simple_positive | true | true | marginal and conditional alternative |
| confounded_positive | false | true | marginal and conditional alternative |

The historical name `confounded_null` refers to a conditional null. Selection and annotation increase with size, producing positive marginal covariance. A Fisher discovery for these terms may represent a real marginal association.

## Metrics

- `target`: marginal hypothesis for Fisher; size-conditional hypothesis for stratification.
- `false_discoveries`: discoveries of terms that are null for the method's target.
- `n_null_terms`: conceptual denominator identifying scenarios with no null hypotheses for that target. Zero FDR when there are no nulls is trivial and does not demonstrate calibration.
- `conditional_target_mismatch`: number of conditionally null terms discovered by Fisher in confounded scenarios. This is a difference in targets, not marginal type I error.
- `probability_any_false_discovery`: frequency of replicates with at least one false discovery against the corresponding target. It equals FDR only under the global null.
- `empirical_fdr`: mean false discovery proportion per replicate, assigning zero when there are no discoveries.
- `power`: detection frequency for the `effect` term, when present. Direct power comparisons between methods require matching targets; in confounded scenarios, interpret this as detection of the predefined term under different hypotheses.
- `any_false_discovery_ci_low/high`: binomial confidence interval for the event of at least one false discovery; it is not an FDR interval in mixtures of null and alternative terms.
- `ci_coverage`: coverage of the known marginal OR only in the simple_positive scenario, without prior selection by significance.

## Current benchmark protocol

Default configuration: seed 42, 100 replicates, 999 permutations, bins 2/4/8, 120 families, eight `null_*` terms and an effect log-OR of 1.5. Scenarios are evaluated using the same term matrix for both methods. BH and a strict 0.05 threshold are applied separately. Binning is not chosen retrospectively based on the best result.

This experiment does not yet model three covariates, incomplete GO annotations or phylogenetic structure. It also does not provide FDR confidence intervals in a mixture scenario. These components belong to expanded validation. Resolution and uncertainty intervals must accompany any conclusion.

## Historical pilot

`benchmarks/pilot` is preserved without changes to its numerical files. Its previous labels assessed discoveries against the conditional null for both methods. In confounded scenarios, Fisher's false-discovery columns cannot be interpreted as errors of its marginal test. See `benchmarks/pilot/INTERPRETATION.md`.

## Expanded validation

`benchmarks.expanded` adds seven scenarios and three covariates. In discrete scenarios, values are balanced so that quantile bins separate each binary covariate exactly. In the continuous scenario, bins approximate conditioning. `missing_size` applies size-dependent annotation loss; `dependent` creates two children and an ancestor defined by their union; `sparse` uses nearly deterministic selection; `derived_degree` computes degree from the matrix itself. `unconfounded` provides a known marginal OR for measuring coverage.

The conditional truth for `derived_degree` is unknown: conditioning on a degree calculated from the outcome can change the hypothesis. The marginal truth for `missing_size` is also unknown, because differential missingness can reverse or cancel the association. FDR is not assigned to these cells. OR=1 represents the null; OR>1 represents a positive effect. Effect detection rates should not be compared between methods with different targets. Zero FDR in scenarios without marginally null terms is trivial.

Short run to verify functionality:

```bash
python -m benchmarks.expanded --config benchmarks/configs/smoke.json --out-dir benchmark_results/smoke
python -m benchmarks.report benchmark_results/smoke --plot
```

For the expanded configuration:

```bash
python -m benchmarks.expanded --config benchmarks/configs/validation.json --out-dir benchmark_results/validation
```

The expanded configuration specifies 1,000 replicates, 9,999 permutations, 240 families, bins 2/4/8 and OR 1/2/4. This configuration is computationally larger and had not been executed at this delivery stage. The short run used five replicates, 99 permutations, 120 families, bins 2/4 and OR 1/4; it demonstrates execution, not calibration.

Each batch is written atomically. Repeating the command reuses completed batches if the configuration, versions and code hashes match. Changes require another directory. `SeedSequence` receives scenario, effect, replicate, purpose and bins: reordering scenarios does not change the data. On resumption, an incomplete batch is rerun deterministically. Scenario IDs are stable and must not be renumbered.

`terms.tsv` preserves p-values, adjustments, reference truth, support and exceedances. `replicates.tsv` preserves discoveries, false discovery proportions, effect detection, coverage and support. `summary.tsv` calculates FDR using replicate bootstrap intervals and binomial intervals for detection/coverage. Cells with fewer than two replicates receive no intervals. Bootstrap intervals that collapse to zero with few replicates do not demonstrate the absence of error; the short run is not suitable for nominal calibration validation.

For each threshold, `calibration.tsv` records the mean across replicates of the rejection fractions among known null terms. SVG figures are descriptive calibration curves without simultaneous intervals. Discrete p-values may produce conservatism. Paired bootstrap intervals for differences between bins, simultaneous curve uncertainty, phylogenetic simulation and evaluation of alternative IEA evidence policies are not yet available.

## Executed external comparisons

The `benchmarks.comparators.run` runner harmonizes family-level inputs and audits counts before comparing p-values. Two-sided GOATOOLS is compared with two-sided Fisher; topGO classic and clusterProfiler are compared with one-sided Fisher. topGO elim is identified as a different algorithm. Versions, omissions and common BH adjustments are documented in `benchmarks/comparators/README.md`. The example run verifies integration; it does not yet constitute a comparator calibration benchmark or a real biological study.

## Measured cost of the full run

A measurement of 21 analyses with 9,999 permutations estimated approximately 35 hours of sequential execution after optimizing deterministic strata (73 hours before) for the full configuration of 63,000 analyses. The extrapolation uses only OR=1 and one replicate; it excludes part of the cost of generation, writing and aggregation. See [timings, support and reproduction](../benchmarks/timing/README.md). Results from the subsequent full run are described below.

## Full run completed on 2026-10-22

Full validation used the statistical configuration in `benchmarks/configs/validation_checkpointed.json` (1,000 replicates, 9,999 permutations, 240 families, bins 2/4/8, odds ratios 1/2/4 and seven scenarios).

To fit the execution within the Slurm cluster resources:

- Analyses were distributed across a job array of 84 simultaneous tasks (7 scenarios × 3 odds ratios × 4 replicate chunks of 250 replicates / 50 batches).
- Total wall-clock time was approximately **39 minutes**, compared with the sequential estimate of ~35 hours.
- All **4,200 batches** (`batch_*.json`) were generated and aggregated successfully, totaling **63,000 analyses**.

Final files are consolidated in `benchmark_results/validation/` (ignored by Git):

- `summary.tsv`: summary of 126 combinations with empirical FDR, power, coverage and bootstrap/binomial intervals.
- `calibration.tsv`: observed null rejection rates at nominal thresholds for calibration assessment.
- `calibration_*.svg`: descriptive vector calibration curves by scenario.
- `replicates.tsv` and `terms.tsv`: raw data by replicate and term.

## Slurm resumption contracts

Run `prepare` once before submission. Each task checks the configuration, NumPy/pandas/SciPy versions, statistical code hashes and scenario manifest without modifying the manifests. `runner_metadata.json` additionally records hashes of the orchestration modules (`parallel.py`, `slurm.py` and `checkpoints.py`). Changes to these files after preparation prevent execution and collection; use a new directory for changed code.

```bash
.venv/bin/python -m benchmarks.slurm prepare --config benchmarks/configs/validation_checkpointed.json --out-dir benchmark_results/validation
```

Tasks partition global batch indices, preserving their names and compatibility with sequential checkpoints. When there are more tasks than batches, some tasks are empty. The default array contains 84 tasks for seven scenarios, three OR values and four chunks per effect; adjust the array range if this configuration changes.

Resumption and collection reject unreadable checkpoints, missing or duplicate rows and mismatched permutation counts. Collection requires exact batch names and complete coverage of scenario, OR, replicate, bins, method and term before consolidating results. These checks verify structure and declared provenance; they cannot retrospectively reconstruct the origin of results produced before manifest validation was implemented.

## Audit of historical results

All 4,200 batches and four consolidated tables were checked, and 21 analyses were repeated without numerical discrepancies. Historical provenance remains partial, and the sparse scenario with two bins had an empirical FDR of 26.5% against the generator's null. Consult the [audit report](validation-audit/README.md) before interpreting or reusing the results.
