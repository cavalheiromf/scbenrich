# SCB-Enrich: Stratified Covariate-Balanced Enrichment for Gene Families

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**SCB-Enrich** is a Python library and command-line tool for exploratory functional overrepresentation analysis using gene families as the units of analysis. It combines family-level Fisher tests with stratified resampling based on user-supplied covariates. Tables can be prepared from CAFE, OrthoFinder, or other family analyses; automatic conversion from their native outputs is not currently provided.

**Status:** the current implementation provides exploratory statistics and adjusts stratified permutation p-values with Benjamini–Hochberg. Statistical calibration, FDR control under the full data-generating process, and confidence-interval coverage across realistic family datasets have not been established by benchmarks presented here.

## Biological question and scope

Gene-level enrichment tools, including topGO, clusterProfiler, and GOATOOLS, can accept lists containing paralogs; they do not require single-copy genes. The relevant issue is whether the unit of analysis, universe, and null model match the biological question and the process used to select the foreground.

Larger families may have more opportunities to acquire annotations, while foreground selection may also depend on family size or annotation availability. These associations can affect enrichment results when omitted from the null model. However, a large family does not necessarily have high annotation density: a poorly characterized expansion may contain many members but few annotated functions. Total annotation counts, distinct terms per family, and annotations per gene are different quantities and must not be used interchangeably.

SCB-Enrich compares two questions:

- **Marginal overrepresentation:** is a term more frequent among foreground families than among the remaining eligible families?
- **Stratified overrepresentation:** is the observed foreground term count unusually high relative to random family sets with the same numbers of families in the selected covariate strata?

Conditioning changes the question. It can remove biologically relevant signal associated with expansion itself, and does not establish causality or identify a uniquely “true” enrichment. Covariates should be justified biologically; annotation degree calculated from the tested annotation matrix may itself contain information about the tested function.

## Installation

```bash
git clone https://github.com/mfcaval/scbenrich.git
cd scbenrich
pip install .
```

For development:

```bash
pip install -e ".[dev]"
pytest
```

## Quick start

### Command-line interface

```bash
scb-enrich \
  --master-table examples/mock_master_table.tsv \
  --annotations examples/mock_annotations.tsv \
  --foreground is_expanded \
  --covariates Total_Count,Count_Ath,Annotation_Degree \
  --obo examples/mock_go.obo \
  --permutations 10000 \
  --min-universe 2 \
  --max-universe-fraction 0.95 \
  --seed 42 \
  --out-dir results/
```

The frequency thresholds above are for the small mock dataset, not recommended biological defaults. The CLI defaults are at least 5 annotated families per term and at most 90% of the universe. Use `scb-enrich --help` for options. Setting `--permutations 0` skips stratified resampling.

### Python API

```python
import pandas as pd
from scbenrich import (
    StratifiedPermuter, fisher_exact_test, exact_conditional_or, fdr_bh,
)

# Marginal family-level comparison; these ORs are not covariate-adjusted.
sample_or, p_fisher = fisher_exact_test(a=15, b=5, c=2, d=18)
conditional_mle, ci_low, ci_high = exact_conditional_or(a=15, b=5, c=2, d=18)

# term_matrix: binary family-by-term DataFrame, aligned to the eligible universe.
# foreground_ids: IDs of selected families within that universe.
cov_df = pd.read_csv("master_table.tsv", sep="\t", index_col="Family_ID")
permuter = StratifiedPermuter(
    cov_df,
    covariates=["Total_Count", "Count_Ath", "Annotation_Degree"],
    seed=42,
)
p_strat = permuter.run_permutations(
    foreground_ids, term_matrix, n_permutations=10000,
)
# An explicit API-level BH adjustment over all terms in this analysis:
terms = list(p_strat)
p_strat_bh = dict(zip(terms, fdr_bh([p_strat[t] for t in terms])))
```

The CLI also applies BH to stratified p-values. Applying BH does not by itself validate the null model or its dependence assumptions.

## Inputs and universe

### Master table (`--master-table`)

Provide one row per unique `Family_ID`, with numeric, finite covariates and a foreground indicator if selection is not supplied as a separate file. TSV (`.tsv`, `.tab`, `.txt`) and CSV files are supported.

You can construct a validated master table programmatically using `build_master_table`:

```python
from scbenrich import build_master_table

# Automatically merges counts, calculates Annotation_Degree,
# extracts GFF3 structural covariates from one or multiple genomes, and adds foreground flag
master_df = build_master_table(
    gene_counts="Orthogroups.GeneCount.tsv",
    annotations="annotations.tsv",
    gff3="annotations_dir/",                 # Single file, list of GFFs, directory, or dict {species: path}
    gene_to_family="Orthogroups.tsv",        # OrthoFinder gene-to-family mapping (TSV or txt)
    by_species=False,                        # True to also include species-specific structural covariates
    foreground={"OG0000001", "OG0000002", "OG0000007"},
    foreground_col="is_expanded",
    output_path="master_table.tsv",
)
```

| Family_ID | Total_Count | Count_Ath | Annotation_Degree | is_expanded |
| :--- | :--- | :--- | :--- | :--- |
| OG0000001 | 25 | 4 | 12 | 1 |
| OG0000002 | 18 | 3 | 8 | 1 |
| OG0000003 | 5 | 1 | 2 | 0 |

All master-table rows form the eligible universe. The comparison background is the universe minus the foreground. Families without annotations remain in the universe and contribute zeros to the annotation matrix; zero means no supplied annotation, not demonstrated absence of biological function. Choose eligibility filters before testing and document their relationship to foreground selection.

`--foreground` accepts a column name or a text file with one family ID per line. Boolean columns use `True`; other columns select values matching `1`, `true`, `yes`, or `expanded`, case-insensitively. File IDs outside the universe raise an error. Foreground and comparison background must both be nonempty. Input validation covers required columns, family IDs, numeric finite covariates, and key CLI parameters; other biological consistency checks remain the user's responsibility.

The example covariate names are conventions, not automatically computed metrics:

- `Total_Count`: document which organisms and member types contribute to the family count.
- `Count_Ath`: document the reference organism, homolog definition, and counting procedure; another numeric reference column can be supplied instead.
- `Annotation_Degree`: define whether this counts distinct direct terms, propagated terms, or another quantity, and which evidence sources are included.

Resolve duplicate IDs, isoforms, missing values, and many-to-many mappings upstream. Organism-specific annotation retrieval, orthology inference, and family construction are outside the current pipeline.

### Annotations (`--annotations`)

Provide `Family_ID` and `Term_ID`; `Term_Name` and `Category` are optional.

| Family_ID | Term_ID | Term_Name | Category |
| :--- | :--- | :--- | :--- |
| OG0000001 | GO:0006950 | response to stress | BP |
| OG0000001 | GO:0009408 | response to heat | BP |

Repeated family–term pairs are collapsed to binary presence. Each family contributes at most once to each term, regardless of its number of genes. Document how member-level annotations were aggregated into family annotations, including transfer by homology and any consensus rule. Only families in the master table contribute to term counts.

## Statistical procedure and assumptions

### Marginal test and effect estimate

For each retained term, the CLI builds a 2×2 table of foreground/background families by annotated/not annotated status. It runs a one-sided Fisher test (`greater`) and applies Benjamini–Hochberg to all retained terms together within that run, across categories. The output `q_value` is this adjusted Fisher p-value, not a local FDR estimate and not an adjusted permutation p-value.

`Marginal_OR` estimates the marginal family-level association using a conditional likelihood based on Fisher's noncentral hypergeometric distribution. “Conditional” refers to the fixed table margins, **not adjustment for the stratification covariates**. SciPy computes the conditional MLE and equal-tailed interval endpoints. Individual 95% intervals are not simultaneous or selection-adjusted intervals.

The conditional OR and equal-tailed interval now use SciPy's `odds_ratio(..., kind="conditional")`. The CLI exports `OR_Status`; a table without an identifiable OR receives `not_identifiable` and missing estimates. Legitimate boundary estimates of zero or infinity remain possible. These individual intervals are not adjusted for term selection or family covariates.

### Stratified resampling

For each covariate, the current bin-count heuristic is `min(10, max(2, n_unique // 5))`. Quantile bins are constructed with `pandas.qcut`, dropping duplicate boundaries; on an exception, the implementation falls back to `pandas.cut`. Joint strata combine the resulting labels. Bin granularity is configurable with `--bins`.

Each replicate samples families without replacement within each stratum from the entire eligible universe, including observed foreground families. It preserves the observed foreground count in each stratum. Families can recur across replicates. The same sampled family set is used for all terms, preserving their observed coannotation structure.

The statistic is the number of sampled families annotated with a term. With `b` replicates whose count is at least the observed count and `B` replicates in total, `p_strat = (b + 1) / (B + 1)`. The default 10,000 replicates yield a minimum p-value of approximately 0.0001; they do not guarantee adequate precision for large multiple-testing problems. Monte Carlo tail-probability intervals are exported.

Validity requires foreground labels to be exchangeable within strata under the null. Exact stratum counts do not imply exact balance of continuous covariate values or control of unmeasured factors. The CLI exports stratum counts and the number of possible within-stratum selections, and rejects a run if no stratum permits randomization. Sparse informative strata can still limit power; automatic merging is not implemented.

### Balance diagnostic

`Max_SMD` is the maximum absolute standardized mean difference across supplied covariates, comparing **foreground and background families annotated with the particular term**. Standardization uses the pooled within-group sample variance. It is not a comparison between the observed foreground and the resampled null sets, nor a demonstration of residual balance within strata.

Empty comparison groups, insufficient variance data, or constant groups with unequal means now produce an undefined diagnostic. Constant groups with equal means produce zero. A threshold of 0.10 is a heuristic; inspect group sizes and covariate distributions before interpreting it.

## Decisions and diagnostics

The output contains independent `marginal_significant` and `stratified_significant` decisions, plus `analysis_status` and `balance_status`. Undefined balance does not hide a statistical discovery. `Marginal_OR`, `Marginal_CI_Low_95`, and `Marginal_CI_High_95` describe the marginal comparison.

The CLI accepts `--bins` and exports term-level support, exceedance counts, permutation resolution, and binomial tail-probability intervals. `balance_diagnostics.tsv` describes the observed foreground means and the distribution of null sample means. These diagnostics do not establish exchangeability or absence of confounding.

See [method contracts](docs/methods.md) and [migration notes](docs/migration.md).

## GO propagation and redundancy summary

With `--obo`, annotations propagate through `is_a` and `part_of` ancestor paths, including multiple parents. `regulates` and `has_part` are not traversed. The three GO root terms are excluded during propagation. Without an OBO file, supplied annotations are used directly.

The CLI normalizes alternative IDs and explicit obsolete-term replacements, writing changes to `go_id_mappings.tsv`. Replacement chains are resolved completely; ambiguous replacements, missing GO IDs and cycles raise explicit errors. The ontology version is recorded in `run_metadata.json`. Use an acyclic ontology for the supported relations. Evidence-code and annotation-qualifier filtering, including IEA and NOT handling, must occur upstream; the input mapping does not encode these semantics.

Nonredundant outputs provide separate greedy presentation summaries for marginal and stratified discovery lists. For CLI results, representatives are ordered by marginal q-value and term ID; it removes later ancestor terms and terms with family-set Jaccard overlap ≥ 0.75. Previously selected ancestors are not retrospectively removed, so not all ancestor–descendant pairs necessarily disappear. It is not a validated semantic equivalence procedure or a multiple-testing correction. Broad terms may carry useful information: retain and inspect the full table and `term_reduction_map_marginal.tsv` / `term_reduction_map_stratified.tsv`.

## Outputs and reproducibility

| File | Contents |
| :--- | :--- |
| `enrichment_results_all.tsv` | All tested terms, counts, marginal OR and interval, Fisher p-value, Fisher BH adjustment, raw and BH-adjusted p_strat, Max_SMD, independent decisions and diagnostic states. |
| `enrichment_results_significant.tsv` | Terms with Fisher `q_value < 0.05`; this is the marginal discovery list. |
| `enrichment_results_stratified.tsv` | Terms with `p_strat_adjusted < 0.05`, when permutations are enabled. |
| `balance_diagnostics.tsv` | Observed covariate means and the distribution of null sample means. |
| `bin_edges.json` | Actual quantile-bin edges for each covariate. |
| `strata_diagnostics.tsv` | Eligible, foreground and background counts, subset count, and informative flag for each stratum. |
| `go_id_mappings.tsv` | Applied GO identifier changes when an OBO file is supplied. |
| `term_reduction_map_marginal.tsv` / `term_reduction_map_stratified.tsv` | Removed term, retained representative, and removal reason when redundancy reduction runs. |
| `run_metadata.json` | Parameters, input hashes and paths, ontology version, seed, and runtime versions. |
| `enrichment_results_nonredundant_marginal.tsv` / `enrichment_results_nonredundant_stratified.tsv` | Separate summaries of the marginal and stratified discovery lists when ontology ancestors are available. |

The CLI reports progress to the console and writes `run_metadata.json`. The manifest includes normalized arguments (including foreground), repository commit and dirty state, and source-module hashes. Execution status is recorded as running, completed, no_testable_terms or failed. Output directories must be empty to avoid mixing runs. The implementation uses NumPy's `default_rng`; the same seed alone does not guarantee identical results across changed inputs, ordering, software versions, or generator state. Reusing a permuter instance advances its generator state.

## Validation priorities

An executable pilot benchmark is provided in `benchmarks/simulate.py`:

```bash
python -m benchmarks.simulate --replicates 30 --permutations 199 \
  --bins 2 4 8 --seed 42 --out-dir benchmark_results
```

It simulates 120 families and eight null terms in four scenarios: an unconfounded null, a null where family size influences both foreground selection and term annotation, an unconfounded known effect, and an effect with the same size-related confounding. It compares one-sided family-level Fisher/BH with stratified permutation/BH using the same term matrix. The known effect has an odds ratio of `exp(1.5)` conditional on size. Marginal OR interval coverage is evaluated only in the unconfounded effect scenario, where that parameter matches the fitted comparison.

The checked-in [pilot summary](benchmarks/pilot/summary.tsv), [replicate results](benchmarks/pilot/replicates.tsv), and [configuration](benchmarks/pilot/config.json) use 30 replicates and 199 permutations and retain their historical metric schema. In the confounded scenarios, Fisher can detect real marginal associations among terms with no additional size-conditional effect: the historical 60% rate is not a marginal false-positive rate. See the [pilot interpretation](benchmarks/pilot/INTERPRETATION.md).

New runs export `term_truth.tsv` with separate marginal and conditional null indicators. Error metrics use each method's target; `conditional_target_mismatch` separately counts Fisher discoveries that do not match the conditional objective. `probability_any_false_discovery` names the event measured by the binomial intervals. The [benchmark protocol](docs/benchmark-protocol.md) explains why these intervals are not FDR intervals in mixed scenarios and why FDR zero with no marginal nulls is uninformative about calibration. Direct power comparisons require matching targets.

An expanded validation runner is available in `benchmarks/expanded.py`, with three covariates, missing annotations, dependent terms, sparse strata and outcome-derived annotation degree. It provides independent random streams, resumable batches, method-specific truth, bootstrap FDR intervals, binomial detection/coverage intervals and term-level calibration tables. See the [expanded protocol](docs/benchmark-protocol.md) and [functional smoke results](benchmarks/expanded_smoke/README.md). The five-replicate smoke run verifies execution; the 1,000-replicate validation configuration has not been run. Install plotting dependencies with `pip install -e ".[benchmark]"`.

Before making confirmatory claims, the project still needs:

- Larger null and known-effect simulations measuring rejection rates, empirical FDR, power, effect bias, interval coverage, and Monte Carlo uncertainty, including annotation-dependent selection and sparse strata.
- Sensitivity analyses for binning, covariate choice, reference organism, annotation evidence, and missing annotations; explicit assessment of biological signal removed by conditioning.
- Comparisons with family-level Fisher and established enrichment tools using harmonized universes, annotations, and hypotheses. Kinases, P450s, and NLRs are useful case studies, not automatic ground-truth controls.
- Exhaustive small-universe checks, numerical edge-case validation, input validation, and regression tests for ontology handling and the full CLI.
- Richer balance diagnostics, numerical status reporting, and machine-readable provenance before relying on automated evidence labels.

The existing test suite is a starting point and is not evidence of population-level statistical calibration. Native input converters, workflow wrappers, and versioned containers are future integration opportunities.

## Methodological references

- [topGO documentation](https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO_manual.html) and [GOATOOLS](https://github.com/tanghaibao/goatools): enrichment inputs and available methods.
- [Phipson and Smyth (2010)](https://pubmed.ncbi.nlm.nih.gov/21044043/): nonzero Monte Carlo permutation p-values.
- [SciPy conditional odds ratio documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.contingency.odds_ratio.html): conditioning on table margins.
- [GO relations](https://geneontology.org/docs/ontology-relations/) and [annotation principles](https://www.geneontology.org/docs/go-annotations/): relation-specific propagation and interpretation of missing annotation.

## Citation

An associated manuscript is in preparation:

> Cavalheiro, M. F., et al. (2026). Macroevolutionary dynamics and comparative genomics in *Vellozia* (*campos rupestres*). *In preparation*.

This is not a published validation reference. For reproducible use, also record the repository URL and the exact software version or commit.

## License

MIT License. Copyright (c) 2026 Mariana Ferreira Cavalheiro.

Runs with no eligible terms export empty tables with headers and metadata status `no_testable_terms`. Bin boundaries are exported to `bin_edges.json` when permutations run.

External family-level baseline adapters are available for GOATOOLS, topGO classic/elim and clusterProfiler. The [comparator protocol and executed example](benchmarks/comparators/README.md) document shared inputs, alternatives, effective universes and harmonized BH audits. Agreement on the mock data checks integration; it is not evidence of general calibration or biological superiority.

Release preparation and validation limits are recorded in [release checks](docs/release-checks.md) and [changes](CHANGELOG.md). CI is configured for Python 3.8 and 3.12; remote CI execution is pending.
