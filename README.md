# SCB-Enrich: Stratified Covariate-Balanced Enrichment for Gene Families

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Tests: Passing](https://img.shields.io/badge/tests-52%20passed-brightgreen.svg)](#installation)

**SCB-Enrich** is a Python library and command-line pipeline designed for robust functional overrepresentation analysis of **gene families (orthogroups)**. 

By modeling gene families rather than individual genes and evaluating candidates against a **stratified covariate-balanced permutation null model**, SCB-Enrich eliminates the false discoveries caused by **paralog pseudoreplication** and **extrinsic family-size confounding** in comparative genomics and macroevolutionary studies.

---

## The Problem: Why Conventional Gene-Level Tools Fail for Gene Families

Standard functional enrichment workflows—such as `clusterProfiler::enricher`, `topGO`, or `GOATOOLS`—were designed for differential expression, where single genes serve as independent observations. When applied to macroevolutionary gene-family expansions (e.g., from CAFE or OrthoFinder), these tools suffer from three fundamental biases:

1. **Paralog Pseudoreplication**: Lineage-specific expansions typically arise via tandem or proximal duplications, yielding dozens of co-annotated paralogous genes derived from a single evolutionary event. Treating each paralog as an independent draw severely violates statistical independence, artificially inflating test statistics ($p < 10^{-50}$ to $10^{-180}$) and generating hundreds of spurious functional discoveries.
2. **Extrinsic Family-Size Confounding**: Expanded families are systematically larger and more thoroughly annotated across reference species. Standard hypergeometric or Fisher tests confound genuine biological specialization with overall family size and annotation density.
3. **The Inadequacy of Topological Graph Decorrelation (`topGO elim`)**: While graph-decorrelation algorithms (e.g., `topGO elim`) resolve vertical parent–child redundancies within the Gene Ontology directed acyclic graph (DAG), they operate entirely along ontology edges. They remain blind to horizontal paralog clusters and family-size disparities, continuing to treat dozens of tandem duplicates as independent events.

### The SCB-Enrich Solution

SCB-Enrich resolves these limitations through a principled two-tier statistical framework:

- **Family as the Unit of Analysis**: Orthogroups serve as the primary, independent sampling units. Paralogous genes within a family contribute exactly once per functional term.
- **Stratified Permutation Null Model**: To distinguish genuine functional selection from family-size artifacts, candidate terms are evaluated against 10,000 stratified permutations that preserve the exact joint distribution of key covariates:
  - **Family Size** (total gene count across the comparative clade)
  - **Reference Homolog Count** (mapped homologs in a model species such as *Arabidopsis thaliana*)
  - **Annotation Degree** (total number of distinct functional terms assigned)
- **Covariate Balance Auditing (SMD)**: For each category, balance between foreground and background families is validated using Standardized Mean Differences, requiring $|SMD| < 0.10$.

```
                        Input Data (OrthoFinder / CAFE / GFF3)
                                          │
                                          ▼
                      ┌───────────────────────────────────────┐
                      │    build_master_table() / CLI Input   │
                      │  - Gene counts & family size          │
                      │  - Reference homologs (Ath count)     │
                      │  - GFF3 gene architecture covariates  │
                      │  - Functional annotation degree       │
                      └───────────────────┬───────────────────┘
                                          │
                                          ▼
                      ┌───────────────────────────────────────┐
                      │    Step 1: Marginal Family ORA        │
                      │    Fisher's exact test & Odds Ratio   │
                      │    Benjamini-Hochberg FDR (q_local)   │
                      └───────────────────┬───────────────────┘
                                          │
                                          ▼
                      ┌───────────────────────────────────────┐
                      │    Step 2: Stratified Null Model      │
                      │    10,000 resamples matched across    │
                      │    joint size & annotation strata     │
                      │    Empirical p_strat calculation      │
                      └───────────────────┬───────────────────┘
                                          │
                                          ▼
                      ┌───────────────────────────────────────┐
                      │    Step 3: Balance Audit & Reduction  │
                      │    Balance check: |SMD| < 0.10        │
                      │    DAG ancestor redundancy reduction  │
                      │    Four-tier evidence classification  │
                      └───────────────────────────────────────┘
```

---

## Empirical Benchmark Comparison

Evaluating lineage-specific expansions in *Vellozia peripherica* (217 expanded families, 2,344 genes, 125 families annotated with `GO_BP`) against all modeled monocot families demonstrates how SCB-Enrich prevents false discoveries:

| Tool / Framework | Sampling Unit | Graph Treatment | Covariate Control (Size / Degree) | Significant GO_BP Terms | Lowest Observed p-value |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`clusterProfiler::enricher`** | Gene | None | None | **329** ($q < 0.05$) | $4.14 \times 10^{-188}$ |
| **`topGO classic` (Fisher)** | Gene | None | None | **380** ($p < 0.05$) | $< 10^{-30}$ |
| **`topGO elim` (Fisher)** | Gene | DAG Decorrelation | None | **143** ($p < 0.05$) | $< 10^{-30}$ |
| **Marginal Family ORA** | Family | None | None | **28** ($q < 0.05$) | $1.07 \times 10^{-4}$ |
| **`scbenrich` (SCB-Enrich)** | Family | DAG Reduction | **Stratified Permutations + SMD < 0.10** | **3** (2 Robust + 1 Suggestive) | $0.0013$ ($q_{local}$) / $0.017$ ($p_{strat}$) |

*Takeaway:* Gene-level tools inflate discoveries by over an order of magnitude due to tandem paralog duplication. Even at the family level, 25 of the 28 marginal terms (89%) fail under covariate-balanced permutations because they are driven by large family size. SCB-Enrich isolates the 2 genuinely specialized pathways (*response to salicylic acid* and *recognition of pollen*).

---

## Installation

### From Source

```bash
git clone https://github.com/mfcaval/scbenrich.git
cd scbenrich
pip install .
```

### For Development & Benchmarks

```bash
pip install -e ".[dev,benchmark]"
pytest
```

---

## Quick Start

### 1. Constructing the Master Table Programmatically

Use `build_master_table` to merge OrthoFinder counts, extract GFF3 gene architecture covariates, compute annotation degrees, and assign foreground flags:

```python
from scbenrich import build_master_table

# Automatically builds a validated master table
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

Run full stratified enrichment with GO DAG propagation:

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

Key CLI options:
- `--master-table`: Path to TSV/CSV with one row per family and numeric covariates.
- `--annotations`: Path to family-to-term mapping (`Family_ID`, `Term_ID`, optional `Term_Name`, `Category`).
- `--foreground`: Column name in master table (boolean or 1/0) or path to a text file with one family ID per line.
- `--covariates`: Comma-separated list of numeric covariates to balance (e.g. `Total_Count,Count_Ath,Annotation_Degree`).
- `--obo`: Path to `go-basic.obo` for hierarchical ancestor propagation and redundancy reduction.
- `--permutations`: Number of stratified resampling replicates (default: `10000`; set to `0` for marginal tests only).
- `--bins`: Number of quantile strata per covariate (default: dynamic quantile heuristic).
- `--out-dir`: Destination folder for all output tables, diagnostics, and run metadata.

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
- Computes Benjamini–Hochberg adjusted $q$-values (`q_value`).

### 2. Stratified Permutations ($p_{strat}$)
- Families are jointly partitioned into strata using multi-dimensional quantile binning across the selected covariates.
- In each of the $B$ replicates (default: 10,000), random foreground sets are drawn without replacement within each stratum, strictly preserving the observed foreground size and stratum proportions.
- The empirical stratified p-value is computed as:
  $$p_{strat} = \frac{b + 1}{B + 1}$$
  where $b$ is the number of permutations where the sampled foreground count equals or exceeds the observed count.

### 3. Covariate Balance Diagnostic (`Max_SMD`)
To verify that significance is not driven by residual covariate imbalance between foreground and background families within a term, SCB-Enrich computes the absolute Standardized Mean Difference:
$$SMD = \frac{\bar{X}_{fg} - \bar{X}_{bg}}{\sqrt{\frac{s_{fg}^2 + s_{bg}^2}{2}}}$$
Values of $|SMD| < 0.10$ indicate that foreground and background families for that category have balanced opportunity distributions.

### 4. Four-Tier Evidence Classification

Every tested term is classified according to statistical stringency:

| Classification | Marginal FDR ($q_{local}$) | Stratified Permutation ($p_{strat}$) | Covariate Balance ($|SMD|$) | Interpretation |
| :--- | :---: | :---: | :---: | :--- |
| **Confirmatory (Robust)** | $< 0.05$ | $< 0.05$ | $< 0.10$ | Highly supported candidate; functional signal is independent of family size and annotation density. |
| **Suggestive** | $< 0.05$ | $< 0.05$ | $\ge 0.10$ (or marginal $q_{global}$) | Candidate with empirical permutation support but residual covariate imbalance. |
| **Size Artifact** | $< 0.05$ | $\ge 0.05$ | Any | False discovery in standard ORA; significance is explained away by family size or annotation bias. |
| **Not Significant** | $\ge 0.05$ | Any | Any | No marginal overrepresentation detected. |

---

## Output Files and Diagnostics

The CLI outputs a structured set of tables and audit diagnostics:

| Output File | Description |
| :--- | :--- |
| `enrichment_results_all.tsv` | Complete results for all tested terms: counts, marginal OR, 95% CI, Fisher p-value, BH $q$-value, raw and adjusted $p_{strat}$, Max_SMD, and evidence classification. |
| `enrichment_results_significant.tsv` | Subset of terms meeting marginal significance ($q < 0.05$). |
| `enrichment_results_stratified.tsv` | Subset of terms meeting stratified permutation significance ($p_{strat\_adj} < 0.05$). |
| `enrichment_results_nonredundant_*.tsv` | Nonredundant term summaries after greedy DAG ancestor reduction (separate for marginal and stratified lists). |
| `term_reduction_map_*.tsv` | Audit mapping showing which descendant terms pruned which ancestor terms during redundancy reduction. |
| `balance_diagnostics.tsv` | Observed covariate means vs. empirical distribution of null sample means. |
| `strata_diagnostics.tsv` | Stratum-by-stratum counts of eligible universe, foreground, background, and permutation degrees of freedom. |
| `bin_edges.json` | Exact quantile-bin boundaries used for each continuous covariate. |
| `go_id_mappings.tsv` | Record of any alternative or obsolete GO IDs normalized via the OBO ontology. |
| `run_metadata.json` | Complete execution manifest: parameters, input file SHA256 hashes, software versions, and random seed. |

---

## Validation & Benchmarks

The repository includes a simulation suite in `benchmarks/` to evaluate calibration and statistical power:

- **Pilot Simulation (`benchmarks/simulate.py`)**: Tests unconfounded nulls, size-confounded nulls, and true effect scenarios across multi-bin grids.
- **Large-Scale Validation (`benchmarks/expanded.py` & `benchmarks/slurm.py`)**: Runs high-throughput calibration simulations (63,000 analyses across 1,000 replicates in 7 scenarios) on Slurm clusters to quantify empirical FDR control, power, and interval coverage.
- **External Comparator Adapters (`benchmarks/comparators/`)**: Baseline wrappers comparing SCB-Enrich directly with `clusterProfiler`, `topGO`, and `GOATOOLS` on harmonized datasets.

To run the pilot benchmark:

```bash
python -m benchmarks.simulate --replicates 30 --permutations 199 --bins 2 4 8 --seed 42 --out-dir benchmark_results
```

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
