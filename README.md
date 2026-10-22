# SCB-Enrich: Stratified Covariate-Balanced Enrichment for Gene Families

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Tests: Passing](https://img.shields.io/badge/tests-52%20passed-brightgreen.svg)](#installation)

**SCB-Enrich** tests functional overrepresentation in **gene families (orthogroups)** rather than individual genes. 

Standard gene-level enrichment tools (`clusterProfiler`, `topGO`) assume each gene is an independent evolutionary event. In gene family expansions, tandem duplications create dozens of co-annotated paralogs, inflating test statistics (`p < 1e-50`) and generating hundreds of false positives. Expanded families are also systematically larger and more annotated. 

SCB-Enrich solves this by treating the orthogroup as the unit of analysis and testing candidates against a **stratified permutation null model** that balances family size, reference homolog count, and annotation density.

---

## Quick Start

### 1. Installation

```bash
pip install scbenrich
# Or from source:
git clone https://github.com/mfcaval/scbenrich.git
cd scbenrich
pip install .
```

### 2. Run the CLI

```bash
scb-enrich \
  --master-table master_table.tsv \
  --annotations annotations.tsv \
  --foreground is_expanded \
  --obo go-basic.obo \
  --out-dir results/
```

This runs a marginal Fisher's exact test and 10,000 stratified permutations, saving complete results to `results/enrichment_results_all.tsv`.

### 3. Or in Python

```python
import pandas as pd
from scbenrich import StratifiedPermuter, build_master_table

# Build a validated master table from OrthoFinder and annotations
master_df = build_master_table(
    gene_counts="Orthogroups.GeneCount.tsv",
    annotations="annotations.tsv",
    foreground={"OG0000001", "OG0000002"},
    foreground_col="is_expanded",
    output_path="master_table.tsv",
)

# Run stratified permutations
cov_df = pd.read_csv("master_table.tsv", sep="\t", index_col="Family_ID")
permuter = StratifiedPermuter(
    cov_df,
    covariates=["Total_Count", "Count_Ath", "Annotation_Degree"],
    seed=42,
)

# term_matrix is a binary DataFrame (Family_ID x Term_ID)
p_strat = permuter.run_permutations(
    foreground_ids=foreground_families,
    term_matrix=binary_matrix,
    n_permutations=10000,
)
```

---

## Why SCB-Enrich? (Simulation Benchmark)

When large gene families expand, standard enrichment tests easily mistake overall family size for functional specialization. The controlled simulation suite in `benchmarks/simulate.py` demonstrates this behavior:

| Method | Sampling Unit | Unconfounded Null (FDR) | Size-Confounded Null (False Discovery Rate) | True Effect Recovery |
| :--- | :---: | :---: | :---: | :---: |
| **Gene-level ORA** (`clusterProfiler`, `topGO`) | Gene | Inflated by tandem paralogs | Highly inflated (spurious `p < 1e-50`) | Confounded with family scale |
| **Marginal Family ORA** (Standard Fisher) | Family | ~5% | **~60%** (size artifacts mistaken for enrichment) | Detects signal, but retains size noise |
| **SCB-Enrich** (Stratified Permutations) | Family | **~5%** | **~4–5%** (strictly controlled) | **High power**, isolates true functional signal |

- **Gene-level tools** (`clusterProfiler`, `topGO`) treat each paralog as an independent observation, inflating test statistics when gene families expand via tandem duplication.
- **Topological tools** (`topGO elim`) decorrelate parent–child terms in the Gene Ontology DAG, but remain blind to horizontal paralog clusters and family-size disparities.
- **SCB-Enrich** controls false discoveries by evaluating candidates against an empirical null model that preserves the joint distribution of family size, reference homologs, and annotation degree.

You can reproduce this benchmark directly from the repository:

```bash
python -m benchmarks.simulate --replicates 100 --permutations 1000 --bins 2 4 --seed 42 --out-dir benchmark_results
```

---

## Input File Formats

### 1. Master Table (`--master-table`)

A tab-separated file with one row per family, numeric covariates, and a foreground indicator:

| Family_ID | Total_Count | Count_Ath | Annotation_Degree | is_expanded |
| :--- | :---: | :---: | :---: | :---: |
| OG0000001 | 25 | 4 | 12 | 1 |
| OG0000002 | 18 | 3 | 8 | 1 |
| OG0000003 | 5 | 1 | 2 | 0 |

- `Family_ID`: Unique orthogroup ID.
- `Total_Count`: Total genes in the family across all analyzed species.
- `Count_Ath`: Homologs in a reference model (e.g., *Arabidopsis*).
- `Annotation_Degree`: Number of distinct functional terms mapped to this family.
- `is_expanded`: `1` (or `True`) for foreground families; `0` for background.

You can generate this table automatically with `build_master_table()` from OrthoFinder counts, GFF3 files, and annotations.

### 2. Annotation Table (`--annotations`)

A tab-separated file mapping families to functional terms (Gene Ontology, KEGG, Pfam):

| Family_ID | Term_ID | Term_Name | Category |
| :--- | :--- | :--- | :--- |
| OG0000001 | GO:0006950 | response to stress | BP |
| OG0000001 | GO:0009408 | response to heat | BP |
| OG0000002 | GO:0006950 | response to stress | BP |

Each family is counted at most once per term, regardless of its gene count.

---

## CLI Options

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--master-table`, `-m` | *Required* | Path to master table TSV/CSV. |
| `--annotations`, `-a` | *Required* | Path to annotations TSV/CSV. |
| `--foreground`, `-f` | *Required* | Foreground column name or path to a text file with foreground IDs. |
| `--covariates`, `-c` | `Total_Count,Count_Ath,Annotation_Degree` | Comma-separated columns used for stratified resampling. |
| `--obo` | `None` | Path to `go-basic.obo` for ancestor term propagation and redundancy reduction. |
| `--permutations`, `-p` | `10000` | Number of stratified resampling permutations (set to `0` for marginal test only). |
| `--min-universe` | `5` | Minimum number of annotated families required to test a term. |
| `--max-universe-fraction` | `0.90` | Maximum fraction of universe annotated with a term (filters general terms). |
| `--bins` | `None` (auto) | Number of quantile bins per covariate. |
| `--seed` | `42` | Random seed for reproducible permutations. |
| `--out-dir`, `-o` | `scb_enrich_results` | Output directory (must be empty or will be created). |

---

## How to Read the Outputs

The pipeline generates structured tables in `--out-dir`:

| Output File | Contents |
| :--- | :--- |
| **`enrichment_results_all.tsv`** | Main results table for all tested terms, including odds ratios, Fisher `q_value`, permutation `p_strat`, and `Max_SMD`. |
| **`enrichment_results_stratified.tsv`** | Candidate terms that remain significant after stratified permutations (`p_strat_adjusted < 0.05`). |
| **`enrichment_results_significant.tsv`** | Terms meeting standard marginal significance (`q_value < 0.05`). |
| **`enrichment_results_nonredundant_*.tsv`** | Reduced term lists after removing ancestor terms in the GO DAG. |
| **`balance_diagnostics.tsv`** | Covariate means for foreground vs. background across permutations. |
| **`run_metadata.json`** | Run parameters, input file hashes, and software version for reproducibility. |

### Decision Criteria

| Result | `q_value` | `p_strat` | `Max_SMD` | Meaning |
| :--- | :---: | :---: | :---: | :--- |
| **Confirmatory (Robust)** | `< 0.05` | `< 0.05` | `< 0.10` | True functional enrichment, independent of family size or annotation bias. |
| **Suggestive** | `< 0.05` | `< 0.05` | `≥ 0.10` | Enriched under permutations, but with residual covariate imbalance. |
| **Size Artifact** | `< 0.05` | `≥ 0.05` | Any | Significant in marginal tests only because large families expanded; not specialized. |
| **Not Significant** | `≥ 0.05` | Any | Any | No evidence of enrichment. |

---

## Methodological Summary

1. **Marginal Overrepresentation**: Builds a 2x2 contingency table (family in foreground vs. background x family has term vs. does not). Tests one-sided Fisher's exact test and computes the conditional odds ratio. Controls false discovery with Benjamini-Hochberg (`q_value`).
2. **Stratified Permutations**: Groups families into joint strata based on their covariates (family size, reference homologs, annotation degree). Resamples foreground sets of the same size with matching stratum proportions. Computes empirical `p_strat = (b + 1) / (B + 1)`.
3. **Covariate Balance**: Computes the Standardized Mean Difference (`Max_SMD`) between foreground and background families annotated with each term. A threshold of `|SMD| < 0.10` confirms that opportunities are balanced.

For simulation protocols, validation audits, and mathematical derivations, see [docs/methods.md](docs/methods.md) and [docs/validation-audit/README.md](docs/validation-audit/README.md).

---

## Use of AI Assistance

Generative AI tools (including Claude and OpenAI Codex) were used to assist with code drafting, unit test development, documentation, and refactoring. The authors designed the study, verified the statistical models, interpreted all biological results, and maintain full responsibility for the code and conclusions.

---

## Citation

If you use SCB-Enrich in your research, please cite:

> Cavalheiro, M. F., et al. (2026). Macroevolutionary dynamics and comparative genomics in *Vellozia* (*campos rupestres*). *In preparation*.

```bash
python -c "import scbenrich; print(scbenrich.__version__)"
git rev-parse HEAD
```

---

## License

MIT License. Copyright (c) 2026 Mariana Ferreira Cavalheiro.
