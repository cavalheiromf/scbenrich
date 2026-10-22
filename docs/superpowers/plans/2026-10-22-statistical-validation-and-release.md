# SCB-Enrich Statistical Validation and Release Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. Execute inline unless delegation is explicitly authorized.

**Goal:** Address the inferential and operational issues identified during review, produce validation evidence proportional to the claims, and prepare a reproducible exploratory release.

**Architecture:** Preserve family-level tests and the basic input interface. Separate statistical decisions, randomization support and covariate diagnostics; centralize validation and provenance in small modules. Validation must distinguish marginal and conditional hypotheses before comparing methods.

**Tech Stack:** Python >=3.8, NumPy, pandas, SciPy >=1.10, pytest; existing `.venv` environment. Matplotlib as an optional benchmark dependency for figures.

## Global Constraints

- Primary hypothesis: overrepresentation conditional on the specified strata, under exchangeability within them.
- Fisher's test and the OR conditional on the margins remain complementary marginal analyses.
- Do not interpret loss of significance as an artifact, causality or absence of an effect.
- Preserve existing local changes; `.ll` is outside the scope of this work.
- Do not run `git add .`; select the files for each task.
- Do not use benchmark results to retrospectively choose the most favorable default.
- Do not treat the current 16 tests as evidence of general calibration.
- Do not promise equivalent results across dependency versions solely because the same seed is used.
- Changes to columns and filenames require migration notes.
- Commits, publication and global installation are not part of creating this plan.

## Initial State and Scope

BH adjustment for stratified p-values, partial handling of undefined SMD, stratum reports, OR through SciPy, simple GO normalization, reduction maps, partial manifests and a pilot benchmark already exist. These components should be corrected and extended without unnecessary reimplementation.

The plan is divided into independently reviewable deliveries:

1. Inference contracts, diagnostics and input safety.
2. GO, outputs and reproducibility.
3. Methodological validation and release preparation.

## File Map

| File | Responsibility |
|---|---|
| `scbenrich/core/inference.py` (new) | Independent decision and status fields. |
| `scbenrich/utils/validation.py` (new) | Shared CLI and API validation. |
| `scbenrich/core/stratified.py` | Sampling, term-level support and null-sample diagnostics. |
| `scbenrich/core/stats.py` | Multiple-testing adjustments and numerical OR interface. |
| `scbenrich/core/obo.py` | Parsing, ID normalization and cycle detection. |
| `scbenrich/core/dag.py` | Presentation summaries and removal maps. |
| `scbenrich/utils/provenance.py` (new) | Manifests, hashes and code identification. |
| `scbenrich/cli.py` | Orchestration and export without duplicated validation. |
| `benchmarks/simulate.py` | Generators with explicit marginal and conditional truths. |
| `benchmarks/report.py` (new) | Summaries, uncertainty and figures from raw results. |
| `tests/test_inference.py`, `tests/test_validation.py`, `tests/test_ontology.py`, `tests/test_provenance.py` (new) | Contracts and regression tests for new components. |
| `tests/test_benchmark.py` | Simulated truth, metrics and reproduction. |
| `README.md`, `docs/methods.md`, `docs/migration.md` | Method, assumptions, limitations and migration. |

## Task 1 — Define Interpretation and the Decision Contract

**Dependencies:** none.

**Files:** create `scbenrich/core/inference.py`, `tests/test_inference.py`, `docs/methods.md`, `docs/migration.md`; modify `scbenrich/cli.py`, `scbenrich/__init__.py`, `README.md`.

**Interface:** `decision_fields(q_marginal, q_stratified, analysis_status, balance_status, alpha=0.05) -> dict`.

Return `marginal_significant`, `stratified_significant`, `analysis_status`, `balance_status`. Analysis states: `ok`, `no_randomization`, `constant_statistic`, `not_run`. Balance states: `within_threshold`, `above_threshold`, `undefined`; these fields do not filter the stratified discovery list.

- [ ] Write a test for independence between significance and diagnostics:

```python
result = decision_fields(0.8, 0.01, "ok", "undefined")
assert result["stratified_significant"] is True
assert result["marginal_significant"] is False
assert result["balance_status"] == "undefined"
assert not decision_fields(0.01, 0.01, "no_randomization", "undefined")["stratified_significant"]
```

- [ ] Run `.venv/bin/python -m pytest tests/test_inference.py -q`; confirm failure because the interface is missing.
- [ ] Implement decisions using full-precision comparisons: `q < alpha`; a missing p-value never produces a discovery. Validate `0 < alpha < 1`.
- [ ] Replace decision-making uses of `Evidence_Class`; document removal of this column in `docs/migration.md`.
- [ ] Rename the OR description to `Marginal_OR` and keep confidence intervals explicitly marginal; record migration of old columns.
- [ ] Remove docstring promises of explicit control of all covariates and rigorous null models without assumptions.
- [ ] Run task tests and existing CLI tests; update expectations only for intentional changes.

**Acceptance:** a stratified discovery may coexist with an undefined diagnostic without contradiction between the full table and the exported list; no confirmatory or artifact label remains.

## Task 2 — Shared Validation and API Contracts

**Dependency:** task 1 for exported states.

**Files:** create `scbenrich/utils/validation.py`, `tests/test_validation.py`; modify `scbenrich/utils/io.py`, `scbenrich/core/stratified.py`, `scbenrich/cli.py`.

**Interfaces:** `validate_family_inputs(covariates, foreground_ids, term_matrix) -> None`; `validate_parameters(permutations, bins, alpha) -> None`.

- [ ] Create parameterized tests for duplicate, null or whitespace-only IDs, nonbinary matrices, extra/missing matrix IDs, empty/full foregrounds, nonnumeric/nonfinite covariates, bins <2, and negative or fractional permutation counts.

```python
with pytest.raises(ValueError, match="binary"):
    validate_family_inputs(cov, {"A"}, pd.DataFrame({"T": [1, 2]}, index=["A", "B"]))
```

- [ ] Run `.venv/bin/python -m pytest tests/test_validation.py -q` and verify the expected failures.
- [ ] Read IDs as strings, preserving leading zeros; reject values that are empty after checking whitespace. Do not silently merge distinct IDs.
- [ ] In the API, align matrix order by index only after validating equality of ID sets.
- [ ] For column-based foregrounds, reject missing and unknown tokens; explicitly define accepted false values (`0`, `false`, `no`, `not_expanded`).
- [ ] Reject annotations outside the universe by default, with counts and examples; do not silently discard them.
- [ ] When no terms are eligible, export empty tables with headers and `no_testable_terms` status, without indexing errors.
- [ ] Run new tests and the full suite.

**Acceptance:** CLI and API reject the same invalid inputs, preserve IDs and produce readable results when there are no testable hypotheses.

## Task 3 — Term-Level Support, Bins and Monte Carlo Uncertainty

**Dependencies:** tasks 1–2.

**Files:** modify `scbenrich/core/stratified.py`, `scbenrich/cli.py`, `tests/test_delivery_one.py`; create `tests/test_randomization.py`.

**New interface:** `run_analysis(foreground_ids, term_matrix, n_permutations=10000) -> dict`, containing `terms` (DataFrame), `strata` (DataFrame), `balance` (DataFrame). Keep `run_permutations` as an adapter returning the previous dictionary.

- [ ] Test a universe of four families with two selected against the six enumerable selections; verify the observed statistic, possible bounds and samples without duplicates.
- [ ] Test a term that is constant across all subsets despite global randomization being possible.
- [ ] For each stratum with `N` families, `K` annotated and `k` selected, sum the bounds `max(0, k-(N-K))` and `min(k,K)` to obtain support for each term count.
- [ ] Mark `constant_statistic` when bounds coincide; return p=1 and no discovery. Record `no_randomization` when all selections are fixed.
- [ ] Export exceedances, actual B, `(b+1)/(B+1)` and resolution `1/(B+1)`.
- [ ] Calculate a binomial interval for the Monte Carlo exceedance probability using `scipy.stats.binomtest`; describe it as uncertainty in the tail estimate, not a biological effect interval.
- [ ] Expose `--bins` in the CLI and preserve the current default when omitted. Record effective bin boundaries, including collapsed repeated boundaries.
- [ ] Export means and quantiles for each covariate in the observed foreground and the distribution of resampled means; keep term-level SMD as a separate description.
- [ ] Test determinism for the same seed and state advancement when reusing an instance.
- [ ] Run `.venv/bin/python -m pytest tests/test_randomization.py tests/test_delivery_one.py -q`.

**Acceptance:** users can distinguish insufficient support, low Monte Carlo precision and lack of significance; the null definition and stopping rule remain fixed.

## Task 4 — GO Without Ambiguous Replacements or Cyclic Recursion

**Dependency:** task 2.

**Files:** modify `scbenrich/core/obo.py`, `scbenrich/cli.py`; create `tests/test_ontology.py`.

**Interface:** `normalize_go_id(term_id, terms, alt_map, replacements) -> str`; `replacements` conceptually becomes `dict[str, set[str]]`, using annotations compatible with Python 3.8.

- [ ] Create small synthetic OBO files with two `replaced_by` entries, a replacement chain, nonexistent destination, cycle and multiple valid parents.

```python
with pytest.raises(ValueError, match="ambiguous"):
    normalize_go_id("GO:0000001", terms, {}, {"GO:0000001": {"GO:0000002", "GO:0000003"}})
```

- [ ] Run `.venv/bin/python -m pytest tests/test_ontology.py -q` and confirm failures.
- [ ] Preserve all replacements in the parser. Resolve chains iteratively using a visited set; reject cycles, multiple destinations and missing destinations.
- [ ] For ancestors, use unvisited/processing/completed states; an edge to a node being processed must raise an error identifying the cycle.
- [ ] Distinguish unknown GO terms from non-GO terms; with OBO supplied, reject unknown GO terms and preserve non-GO terms without propagation.
- [ ] Export the complete normalization path and reason for each transformation. Update names and categories using the canonical ID.
- [ ] Keep only `is_a` and `part_of`; explicitly test that `has_part` and `regulates` do not propagate.
- [ ] Run ontology and CLI integration tests.

**Acceptance:** no ambiguity is resolved by OBO line order; cyclic graphs fail with an identifiable message.

## Task 5 — Summaries and Provenance Without Mixing Runs

**Dependencies:** tasks 1–4.

**Files:** create `scbenrich/utils/provenance.py`, `tests/test_provenance.py`; modify `scbenrich/cli.py`, `scbenrich/core/dag.py`, `README.md`.

**Interfaces:** `build_manifest(args, inputs, outputs, status) -> dict`; `write_manifest(path, manifest) -> None`.

- [ ] Test that two different foreground columns produce different manifests even with identical input hashes.
- [ ] Test that a directory containing artifacts from another run is rejected, preserving all existing files.
- [ ] Create a dedicated directory for each run or require an empty directory. Never silently clear a user-supplied directory.
- [ ] Record normalized arguments, foreground file or column, hashes, versions, RNG generator, seed, B, bins, adjustment method, commit, local-change indicator and hashes of modules used.
- [ ] Treat missing Git information as an explicit state; code version/hash remains mandatory.
- [ ] Write the initial manifest with `running`, followed by `completed` or `failed`; on failure, record error type and message without replacing the original error.
- [ ] Export marginal and stratified GO summaries in separate files. Reduction must use exactly the list declared in the filename.
- [ ] Keep separate maps for each summary. Report the rule, threshold and representative for each removal.
- [ ] Test a term significant only in the stratified analysis: it must appear in the stratified summary and be absent from the marginal summary.
- [ ] Run `.venv/bin/python -m pytest tests/test_provenance.py tests/test_delivery_two.py -q`.

**Acceptance:** each file belongs to an identifiable run; the selection used can be reconstructed; summaries do not silently omit discoveries from the analysis they represent.

## Task 6 — Correct Benchmark Reference Truth

**Dependencies:** task 1; executable before the remaining operational tasks.

**Files:** modify `benchmarks/simulate.py`, `tests/test_benchmark.py`; create `docs/benchmark-protocol.md`; update `README.md`.

**Contract:** each simulated term has `marginal_null` and `conditional_null`. Neither truth is inferred from the resulting p-value.

- [ ] Create a deterministic test for a scenario in which size increases selection and annotation: expected marginal covariance is positive even though additional conditional association is zero.
- [ ] Label Fisher results against the conditional null as `conditional_target_mismatch`, rather than Fisher type I error.
- [ ] Keep each method's error metrics relative to its own target; compare power directly only in scenarios with matching targets.
- [ ] Preserve the original pilot as a historical artifact with an interpretation note; do not overwrite old results with new definitions.
- [ ] Rename `rejection_rate` to `probability_any_false_discovery`. Under the global null it equals FDR; in a mixture scenario it does not.
- [ ] Specify scenarios, seeds, bins, metrics and analysis rules in `docs/benchmark-protocol.md` before the expanded run.
- [ ] Run `.venv/bin/python -m pytest tests/test_benchmark.py -q`.

**Acceptance:** no superiority claim depends on labeling a true marginal association as a marginal false positive.

## Task 7 — Expand Simulations and Uncertainty

**Dependencies:** tasks 3 and 6.

**Files:** modify `benchmarks/simulate.py`, `tests/test_benchmark.py`; create `benchmarks/report.py`, `benchmarks/configs/validation.json`, `tests/test_benchmark_report.py`.

**Outputs:** term-level results, replicate-level results, scenario summaries, configuration, code hashes and versions. All units used in error estimates must be identified.

- [ ] Add size, reference homolog count and annotation degree. Include variants where degree is external and where it is calculated from the tested matrix.
- [ ] Include five scenario families: exact exchangeability; residual continuous confounding; size-dependent incomplete annotation; terms dependent through ancestor structure; sparse strata and nearly deterministic selection.
- [ ] Fix initial grids of bins `[2, 4, 8]`, OR effects `[1, 2, 4]` and missingness rates `[0, 0.3, 0.6]`; first run predefined combinations that isolate one factor at a time, avoiding purposeless combinatorial expansion.
- [ ] Run a smoke test with 5 replicates and 99 permutations for each scenario, checking files, states and reproduction.
- [ ] Run validation with 1,000 replicates per selected scenario and 9,999 permutations, in resumable batches of 50 replicates; record runtime and failures for each batch.
- [ ] Record independent seeds by scenario, replicate and method via `SeedSequence`; adding a scenario must not change previous simulations.
- [ ] Report FDR as the mean false discovery proportion per replicate; use replicate bootstrap intervals for FDR and paired differences. For binary power and coverage, use binomial intervals when the unit is an independent replicate.
- [ ] Generate null-term calibration curves, power/FDR plots by bins and confidence interval coverage plots. Do not treat thousands of terms in the same replicate as independent observations.
- [ ] Calculate OR coverage only when its marginal target is known; do not compare marginal intervals with an OR conditional on size.
- [ ] Do not change defaults after observing which configuration performed best; results outside nominal levels motivate investigation and new validation under a separate protocol.

**Acceptance:** a reproducible report with uncertainty, assumptions and failure conditions; absence of significant differences is not proof of equivalence.

## Task 8 — External Comparators and Case Studies

**Dependency:** task 6; expanded analyses depend on task 7.

**Files:** create `benchmarks/comparators/README.md` and specific adapters only for tools that are actually installed and version-pinned.

- [ ] Define two comparisons: gene-level and family-level, without mixing them in the same error estimate.
- [ ] Export identical matrices and universes for Fisher, GOATOOLS, clusterProfiler and explicitly identified topGO algorithms.
- [ ] Record propagation, filters, GO relationships, test alternatives and multiple-testing adjustment for each comparator.
- [ ] Confirm that two Fisher implementations recover the same results for the same problem; investigate differences before comparing distinct methods.
- [ ] Add kinase, P450 and NLR case studies only with recoverable data and versions; treat them as biological plausibility, not automatic calibration truth.
- [ ] If real data are unavailable, mark the section as unexecuted and complete simulated validation separately; do not invent controls.

**Acceptance:** comparisons make each tool's hypothesis explicit and allow its universe and annotations to be reproduced.

## Task 9 — Integrated Review and Exploratory Release

**Dependencies:** tasks 1–6 for an exploratory release; tasks 7–8 for broader methodological claims.

**Files:** `README.md`, `docs/methods.md`, `docs/migration.md`, `pyproject.toml`, `.github/workflows/tests.yml` (new).

- [ ] Run `.venv/bin/python -m pytest -q` and `git diff --check`; require success.
- [ ] Run the CLI with examples, without permutations, with permutations, without eligible terms and with an intentional GO error; check states and files for each run.
- [ ] Test clean installation and the advertised minimum Python version as well as a modern version; update advertised support if the matrix cannot be met.
- [ ] Review all text containing `rigorous`, `robust`, `confirmatory`, `artifact` or equivalents; each claim must point to evidence or be rephrased.
- [ ] Review benchmark results separately from software test success.
- [ ] Prepare migration notes for columns, files, input rules and GO replacements.
- [ ] Create commits for each reviewed delivery, explicitly selecting files. Package publication is a later, separate action.

**Acceptance:** the release is described as exploratory, without general FDR control claims; known limitations are not hidden by passing tests.

## Execution Order and Review Checkpoints

1. Tasks 1–3: finalize inference contracts and diagnostics.
2. Tasks 4–5: finalize GO handling and traceability.
3. Task 6: correct pilot interpretation before expanding its conclusions.
4. Task 7: run validation in batches with progress updates.
5. Task 8: run comparisons using available data and environments.
6. Task 9: review the complete work and prepare the release.

Creating this document does not implement these tasks. Checkboxes must be marked only after the corresponding code, tests and documentation have been verified.

## Execution checkpoint

First implementation block verified: 23 tests passed. Decision fields, marginal OR names, shared API validation, ID string preservation, foreground token validation, term support, Monte Carlo intervals, null-mean diagnostics and CLI bins are implemented. Tasks 1–3 remain partially open: empty-term output, exact bin-boundary export, broader parameterized validation and distribution diagnostics require completion. Tasks 4–9 have not been executed under this plan. No release or calibration claim is complete.

## Second execution checkpoint

Implemented: empty-term output, actual bin-edge export, GO replacement sets/chains/ambiguity rejection, cycle detection, canonical GO names, complete mapping paths, empty output-directory enforcement, running/failed manifests, normalized arguments, Git state and source hashes, and separate marginal/stratified summaries. Tests cover these contracts; remaining checklist items are not automatically marked complete. Tasks 6–9 (benchmark interpretation, expanded validation, external comparators and release) remain open. Distribution diagnostics currently summarize means; quantiles of foreground values and wider validation tests remain to expand.

## Benchmark interpretation checkpoint

Task 6 implemented: population truth exported per term for marginal and continuous-size-conditional targets; method-specific false-discovery metrics, explicit conditional_target_mismatch, renamed probability_any_false_discovery, historical pilot interpretation preserved, and benchmark protocol documented. Tests verify positive population covariance under conditional null and method-specific error counting. Expanded simulations and external comparator work remain open.

## Expanded validation checkpoint

Task 7 infrastructure implemented in benchmarks/expanded.py and benchmarks/report.py: seven scenarios, three covariates, unknown truth preserved, independent SeedSequence streams, atomic restartable batches, source/version validation, term-level outputs, replicate bootstrap FDR intervals, binomial power/coverage intervals and calibration figures. Executed smoke configuration: five replicates, 99 permutations, OR 1/4, bins 2/4. Summary artifacts are preserved in benchmarks/expanded_smoke. Full 1,000-replicate/9,999-permutation configuration is supplied but not executed. Coverage assertions use the unconfounded marginal target only. Paired difference intervals and fuller distribution diagnostics remain open; these tasks are not marked as fully complete. Suite: 35 tests passed.

## External comparator checkpoint

Task 8 family-level adapters implemented and executed with installed GOATOOLS 1.6.4, topGO 2.62.0 classic/elim and clusterProfiler 4.18.4. Exported shared matrix/universe and audited counts, denominators, alternatives, p-values and harmonized BH. GOATOOLS two-sided, topGO classic greater and clusterProfiler greater agreed with matched Fisher references on mock data; elim differed as expected. clusterProfiler zero-overlap omissions are reconstructed as upper-tail p=1 with explicit source marking. Results and versions are preserved in benchmarks/comparators/example_results. Gene-level comparators and biological case studies remain open. Suite: 37 tests passed. No release or completed large-scale validation is claimed.

## Exploratory release checkpoint

Task 9 preparation implemented: Python 3.8/3.12 CI workflow, wheel build and installation outside the checkout on Python 3.12, four subprocess CLI cases, corrected output names and exploratory wording, CHANGELOG and release-check limits. The local wheel environment inherited scientific dependencies; fully isolated resolution and remote CI remain unverified. Python 3.8 syntax compilation passed; numerical execution is pending dependencies/CI. Suite: 41 tests passed; git diff --check passed. No release publication or version bump occurred. Full simulation configuration, gene-level comparisons, biological studies and other previously listed validation gaps remain open.
