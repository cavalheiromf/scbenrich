# SCB-Enrich: Stratified Covariate-Balanced Enrichment for Gene Families

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**SCB-Enrich** is a Python library and command-line pipeline for exploratory functional overrepresentation analysis of **gene families (orthogroups)**. It reports marginal family-level associations and stratified permutation results as independent statistical decisions.

Each family contributes at most once per functional term. The permutation test preserves the observed number of selected families within each joint covariate bin. It does not exactly condition on continuous covariate values, and residual confounding can remain within bins. **The available validation does not establish general FDR control at 5%.**

## Quick Start

Run a complete example using the small dataset included in the repository. From a terminal:

```bash
git clone https://github.com/mfcaval/scbenrich.git
cd scbenrich
python -m venv .venv
source .venv/bin/activate
python -m pip install .

scb-enrich \
  --master-table examples/mock_master_table.tsv \
  --annotations examples/mock_annotations.tsv \
  --foreground is_expanded \
  --min-universe 2 \
  --out-dir quickstart_results
```

If the repository and environment are already set up, run only the `scb-enrich` command from the repository root. The output directory must be absent or empty.

This example selects families marked `is_expanded=1` as the foreground and uses all other families in the master table as the background. It uses the default covariates, 10,000 permutations, automatic bins and seed 42. `--min-universe 2` allows terms annotated in at least two families to be tested in this small example; the normal default is five. No OBO file is needed.

Open `quickstart_results/enrichment_results_all.tsv` for the complete results. The `enrichment_results_significant.tsv` and `enrichment_results_stratified.tsv` files contain the separate marginal and stratified discovery lists; either list may be empty. `run_metadata.json` records the execution status and provenance.

This example demonstrates the workflow, not biological validation. For your own data, replace the two input files and the foreground column or ID list, then review the [CLI arguments and defaults](#cli-arguments-and-defaults) and [input formats](#input-formats).

## Biological Motivation and Scope

A family-level analysis matches questions about family expansion and avoids counting multiple annotated paralogs as separate family observations. It does not establish independence between families or account for shared phylogenetic history.

Family size, reference homolog count and annotation density can influence annotation opportunity. A large family need not be densely annotated: a poorly characterized expansion can contain many members and few known functions. These quantities should be examined separately when selecting covariates.

GO hierarchy reduction and covariate adjustment address different aspects of an analysis. The comparator adapters use harmonized family-level inputs to check baseline calculations. Their example results do not establish superiority over other tools or show that additional discoveries from a gene-level analysis are false.

## Analysis Workflow

1. Define the eligible family universe, foreground and binary family-to-term annotations.
2. Compute marginal Fisher tests, marginal odds ratios and BH-adjusted p-values.
3. Resample foreground families within joint quantile bins, preserving stratum counts, and apply BH to the stratified p-values.
4. Report independent significance decisions, randomization support and descriptive covariate balance diagnostics.
5. Optionally produce GO redundancy summaries with removal-to-representative audit maps.

Permutation inference requires exchangeability within the selected strata. Coarse bins can leave important covariate differences, while finer bins can leave few possible exchanges and reduce power. Binning should be specified before examining significance results; a lower rejection rate alone does not establish a better configuration.

Covariate choice also changes the target hypothesis. In particular, annotation degree calculated from the tested term matrix is outcome-derived. The benchmark treats that setting separately and leaves its conditional truth unknown.

## Validation Results and Limitations

The expanded historical run contains 1,000 replicates, 9,999 permutations per analysis and seven scenarios. An audit checked all 4,200 checkpoints, reconstructed the four consolidated tables and reproduced 21 sampled analyses without numerical discrepancies. Historical orchestration hashes were not recorded, so provenance remains partial. See the [audit report and evidence](docs/validation-audit/README.md).

Selected results for the stratified method under OR=1 are:

| Scenario | Bins per covariate | Empirical FDR | 95% bootstrap interval | Mean fraction of informative strata |
| :--- | ---: | ---: | ---: | ---: |
| Exact | 2 | 3.9% | 2.8–5.1% | 99.7% |
| Continuous | 2 | 6.0% | 4.6–7.5% | 100.0% |
| Sparse | 2 | 26.5% | 23.8–29.3% | 80.8% |
| Sparse | 4 | 4.2% | 3.0–5.5% | 28.6% |
| Sparse | 8 | 0.1% | 0.0–0.3% | 3.4% |

These estimates use 1,000 independent replicates per configuration. The intervals are not simultaneous across configurations. In continuous scenarios, reference truth is conditional on the generator's continuous covariates, whereas the implemented test conditions on their bins. The sparse scenario with two bins shows that this approximation can be inadequate. The low rejection rate with eight bins accompanies a substantial loss of informative strata and must be interpreted alongside power.

BH adjustment is applied separately to each test family; it cannot repair p-values affected by residual confounding. Statistical significance does not establish biological specialization, and loss of significance after stratification does not establish that a marginal discovery is a size artifact.

The historical results support exploratory assessment under these limitations. A run with the corrected Slurm runner and complete provenance from the outset remains necessary for a fully traceable validation. Consult the [full summary](docs/validation-audit/summary.tsv), [benchmark protocol](docs/benchmark-protocol.md), [comparator audits](benchmarks/comparators/README.md) and [release checks](docs/release-checks.md).

---

## Installation

### From Source

```bash
git clone https://github.com/mfcaval/scbenrich.git
cd scbenrich
pip install .
```

A standard `pip install .` installs the required runtime dependencies: NumPy, SciPy, pandas and PyYAML. The package requires Python 3.8 or later. Optional extras are installed only when requested:

| Extra | Adds | Needed for |
| :--- | :--- | :--- |
| `dev` | pytest | Running the test suite. |
| `benchmark` | Matplotlib | Generating benchmark plots. |

### For Development & Benchmarks

```bash
pip install -e ".[dev,benchmark]"
pytest
```

---

## Detailed Usage

### 1. Constructing the Master Table Programmatically

Use `build_master_table` to merge OrthoFinder counts, extract GFF3 gene architecture covariates, compute annotation degrees, and assign foreground flags:

```python
from scbenrich import build_master_table

# Build the master table; inspect the covariates before analysis
master_df = build_master_table(
    gene_counts="Orthogroups.GeneCount.tsv",          # OrthoFinder gene count matrix
    annotations="annotations.tsv",                   # Functional annotations (Family_ID, Term_ID)
    gff3="annotations_dir/",                         # Single GFF3, list of GFFs, or directory
    gene_to_family="Orthogroups.tsv",                # OrthoFinder gene-to-family mapping
    by_species=False,                                # True to also include species-specific covariates
    foreground={"OG0000001", "OG0000002"},           # Set, list, or Series of foreground IDs
    foreground_col="is_expanded",
    output_path="master_table.tsv",
)
```

### 2. Command-Line Interface (CLI)

The CLI requires only `--master-table`, `--annotations` and `--foreground`. A minimal invocation is:

```bash
scb-enrich \
  --master-table master_table.tsv \
  --annotations annotations.tsv \
  --foreground is_expanded
```

This uses the default covariate columns, runs 10,000 permutations, writes to `scb_enrich_results/` and does not load an ontology. The default covariate columns must exist and contain finite numeric values. To use different columns, pass `--covariates` explicitly.

The following example supplies optional settings explicitly and enables GO ancestor propagation by providing `--obo`:


```bash
scb-enrich \
  --master-table master_table.tsv \
  --annotations annotations.tsv \
  --foreground is_expanded \
  --covariates Total_Count,Count_Ath,Annotation_Degree \
  --obo go-basic.obo \
  --permutations 10000 \
  --min-universe 5 \
  --max-universe-fraction 0.90 \
  --seed 42 \
  --out-dir results/
```

#### CLI Arguments and Defaults

| Argument | Required? | Default when omitted | Purpose |
| :--- | :--- | :--- | :--- |
| `--master-table`, `-m` | Yes | None | TSV/CSV table with one row per family, a `Family_ID` column and numeric covariates. |
| `--annotations`, `-a` | Yes | None | TSV/CSV family-to-term mapping with `Family_ID` and `Term_ID`; `Term_Name` and `Category` are optional columns. |
| `--foreground`, `-f` | Yes | None | Foreground column name in the master table, or a text file with one family ID per line. |
| `--covariates`, `-c` | No | `Total_Count,Count_Ath,Annotation_Degree` | Comma-separated numeric columns used for joint stratification. Omitting the argument still requires these default columns in the input. |
| `--obo` | No | `None` (no ontology loaded) | OBO file, such as `go-basic.obo`, for GO ID normalization, ancestor propagation and redundancy summaries. |
| `--permutations`, `-p` | No | `10000` | Number of stratified resampling replicates. Set to `0` for marginal tests only. |
| `--min-universe` | No | `5` | Minimum number of eligible families annotated with a term for that term to be tested. |
| `--max-universe-fraction` | No | `0.90` | Maximum fraction of eligible families annotated with a term for that term to be tested. |
| `--bins` | No | `None` (automatic heuristic) | Requested quantile bins per covariate; an explicit value must be an integer ≥2. |
| `--seed` | No | `42` | Random seed for resampling. |
| `--out-dir`, `-o` | No | `scb_enrich_results` | Output directory. It must be absent or empty; failed runs retain failure metadata. |
| `--help`, `-h` | No | Not requested | Print the CLI help and exit. |

Without `--obo`, the analysis uses the terms supplied directly in the annotation table. Marginal and stratified tests still run, but GO ID normalization, ancestor propagation and hierarchy-based redundancy summaries are disabled.

When `--bins` is omitted, the requested number for each covariate is `min(10, max(2, n_unique // 5))`, where `n_unique` is its number of distinct values. Repeated quantile boundaries can reduce the actual number of bins; constant covariates form one group. The automatic setting is a heuristic, not a guarantee of calibration or adequate randomization support.

`--min-universe` and `--max-universe-fraction` filter terms, not families. All rows in the master table define the eligible family universe; the background consists of all eligible families outside the foreground.

### 3. Python API

```python
import pandas as pd
from scbenrich import StratifiedPermuter, fisher_exact_test, exact_conditional_or, fdr_bh

# 1. Marginal 2x2 comparison
sample_or, p_fisher = fisher_exact_test(a=15, b=5, c=2, d=18)
conditional_mle, ci_low, ci_high = exact_conditional_or(a=15, b=5, c=2, d=18)

# 2. Stratified Permutations
cov_df = pd.read_csv("master_table.tsv", sep="\t", index_col="Family_ID")
permuter = StratifiedPermuter(
    cov_df,
    covariates=["Total_Count", "Count_Ath", "Annotation_Degree"],
    seed=42,
)

# term_matrix: binary family-by-term DataFrame
# foreground_ids: list/set of selected family IDs
p_strat = permuter.run_permutations(
    foreground_ids=foreground_families,
    term_matrix=binary_matrix,
    n_permutations=10000,
)

# Adjust stratified p-values using Benjamini-Hochberg
terms = list(p_strat.keys())
p_strat_adj = dict(zip(terms, fdr_bh([p_strat[t] for t in terms])))
```

---

## Input Formats

### 1. Master Table (`--master-table`)

A tab-separated file with one row per unique `Family_ID`, containing all numeric covariates and the foreground indicator:

| Family_ID | Total_Count | Count_Ath | Annotation_Degree | is_expanded |
| :--- | :---: | :---: | :---: | :---: |
| OG0000001 | 25 | 4 | 12 | 1 |
| OG0000002 | 18 | 3 | 8 | 1 |
| OG0000003 | 5 | 1 | 2 | 0 |

*Universe Definition:* All rows in the master table define the eligible universe. Families without annotations remain in the universe and contribute zeros to the annotation matrix, ensuring that background size is not distorted.

### 2. Annotation Table (`--annotations`)

A tab-separated file linking families to functional terms. Multiple genes in the same family are collapsed to binary presence (each family contributes at most once per term):

| Family_ID | Term_ID | Term_Name | Category |
| :--- | :--- | :--- | :--- |
| OG0000001 | GO:0006950 | response to stress | BP |
| OG0000001 | GO:0009408 | response to heat | BP |
| OG0000002 | GO:0006950 | response to stress | BP |

---

## Statistical Procedure and Decision Rules

### 1. Marginal Overrepresentation
For each term, a 2×2 contingency table is constructed across eligible families:

- Evaluates one-sided Fisher's exact test (`greater`).
- Calculates the exact conditional maximum likelihood odds ratio (`Marginal_OR`) with 95% confidence intervals based on Fisher's noncentral hypergeometric distribution.
- Computes Benjamini–Hochberg adjusted q-values (`q_value`).

### 2. Stratified Permutations (`p_strat`)

- Families are jointly partitioned into strata using multi-dimensional quantile binning across the selected covariates.
- In each of the replicates (default: 10,000), random foreground sets are drawn without replacement within each stratum, strictly preserving the observed foreground size and stratum proportions.
- The empirical stratified p-value is computed as:

```text
p_strat = (b + 1) / (B + 1)
```

where `b` is the number of permutations where the sampled foreground count equals or exceeds the observed count, and `B` is the total number of permutations.

Stratified p-values are adjusted across the tested terms using Benjamini–Hochberg (`p_strat_adjusted`). The marginal and stratified tests target different hypotheses. Their significance decisions are reported independently.

### 3. Covariate Balance Diagnostic (`Max_SMD`)

For each term, the code compares the covariates of annotated foreground and annotated background families using a sample-size-weighted pooled standard deviation:

```text
s_pool = sqrt(((n_fg - 1) * var(X_fg) + (n_bg - 1) * var(X_bg))
              / (n_fg + n_bg - 2))
SMD = (mean(X_fg) - mean(X_bg)) / s_pool
Max_SMD = maximum absolute SMD across the selected covariates
```

Variances use sample estimates when at least two observations are available. Empty groups, insufficient variance information or zero pooled variance with unequal means produce an undefined diagnostic. Equal constant groups have SMD zero.

`balance_status` is `within_threshold` when `Max_SMD < 0.10`, `above_threshold` otherwise, and `undefined` when it cannot be calculated. This is a descriptive comparison within annotated families. It is not proof of exchangeability, removed confounding or FDR control, and it does not select discoveries. `balance_diagnostics.tsv` separately compares observed foreground covariate means with their distributions under resampling.

### 4. Independent Decisions and Support States

| Field | Meaning |
| :--- | :--- |
| `marginal_significant` | Finite marginal `q_value < 0.05`. |
| `stratified_significant` | Finite `p_strat_adjusted < 0.05` and `analysis_status == "ok"`. |
| `analysis_status` | `ok`, `constant_statistic`, `no_randomization` or `not_run`. |
| `balance_status` | Descriptive balance state, separate from significance decisions. |

The Python API reports `no_randomization` when no stratum permits exchanges. The CLI rejects that configuration and retains failure metadata. Constant term statistics receive `constant_statistic`; disabling permutations produces `not_run` in CLI decisions.

`Marginal_OR` and its interval describe the marginal family-level comparison, without covariate adjustment. The Monte Carlo `tail_ci_low` and `tail_ci_high` fields describe uncertainty in the exceedance probability; they are not intervals for an effect size or an adjusted p-value.

See [methods](docs/methods.md) and [migration notes](docs/migration.md) for details.

---

## Output Files and Diagnostics

The CLI outputs a structured set of tables and audit diagnostics:

| Output File | Description |
| :--- | :--- |
| `enrichment_results_all.tsv` | Complete results for all tested terms: counts, marginal OR, 95% CI, Fisher p-value, BH q-value, raw and adjusted p_strat, Max_SMD, independent significance decisions and diagnostic states. |
| `enrichment_results_significant.tsv` | Subset of terms meeting marginal significance (`q_value < 0.05`). |
| `enrichment_results_stratified.tsv` | Subset of terms meeting stratified permutation significance (`p_strat_adjusted < 0.05` and `analysis_status == "ok"`). |
| `enrichment_results_nonredundant_*.tsv` | Nonredundant term summaries after greedy DAG ancestor reduction (separate for marginal and stratified lists). |
| `term_reduction_map_*.tsv` | Removed terms, retained representatives and reasons for greedy ancestor or overlap reduction. These are descriptive summaries; FDR control is not guaranteed for the reduced lists. |
| `balance_diagnostics.tsv` | Observed covariate means vs. empirical distribution of null sample means. |
| `strata_diagnostics.tsv` | Stratum-by-stratum counts of eligible universe, foreground, background, possible selections and informative-stratum indicators. |
| `bin_edges.json` | Exact quantile-bin boundaries used for each continuous covariate. |
| `go_id_mappings.tsv` | Record of any alternative or obsolete GO IDs normalized via the OBO ontology. |
| `run_metadata.json` | Execution manifest with available parameters, input file SHA256 hashes, software versions, and random seed. |

---

## Validation & Benchmarks

The repository includes a simulation suite in `benchmarks/` to evaluate calibration and statistical power:

- **Pilot Simulation (`benchmarks/simulate.py`)**: Tests unconfounded nulls, size-confounded nulls, and true effect scenarios across multi-bin grids.
- **Expanded Validation (`benchmarks/expanded.py` & `benchmarks/slurm.py`)**: Sequential and Slurm runners for calibration, power and coverage assessment. The audited historical run contains 63,000 analyses; interpretation and provenance limitations are described above. The corrected Slurm runner verifies prepared manifests and checkpoint coverage and partitions work along complete batch boundaries.
- **External Comparator Adapters (`benchmarks/comparators/`)**: Baseline wrappers comparing SCB-Enrich directly with `clusterProfiler`, `topGO`, and `GOATOOLS` on harmonized datasets.

To run the pilot benchmark:

```bash
python -m benchmarks.simulate --replicates 30 --permutations 199 --bins 2 4 8 --seed 42 --out-dir benchmark_results/pilot
```

---

## Use of AI Assistance

AI assistance, including OpenAI Codex, was used during software development, test creation, documentation drafting and translation, and code and statistical-method review. AI-assisted review is not independent scientific validation. The project authors remain responsible for the code, methodological choices, interpretation of results and final documentation. Validation evidence and known limitations are documented in the sections above and the linked audit report.

---

## Citation

If you use SCB-Enrich in your research, please cite:

> Cavalheiro, M. F., et al. (2026). Macroevolutionary dynamics and comparative genomics in *Vellozia* (*campos rupestres*). *In preparation*.

For computational reproducibility, please also report the software version and repository commit hash:

```bash
python -c "import scbenrich; print(scbenrich.__version__)"
git rev-parse HEAD
```

---

## License

This project is licensed under the [MIT License](LICENSE).  
Copyright (c) 2026 Mariana Ferreira Cavalheiro.
