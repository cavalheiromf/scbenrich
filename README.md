# SCB-Enrich: Stratified Covariate-Balanced Enrichment for Gene Families

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**SCB-Enrich** is a statistically rigorous Python framework and command-line tool designed for functional overrepresentation analysis of multi-gene families (e.g. from CAFE, OrthoFinder, or custom copy-number variation datasets).

Standard gene set enrichment tools (e.g., topGO, clusterProfiler, GOATOOLS) were designed for discrete single-copy gene lists (such as RNA-seq differential expression). When applied to **multigene families**, standard hypergeometric and Fisher tests suffer from severe **family-size and annotation-density biases**: large gene families (such as protein kinases, cytochrome P450s, or NLRs) have more paralogs, more cross-species homologs, and consequently much higher annotation density. They are frequently misidentified as "significantly enriched" purely due to copy-number abundance rather than true biological overrepresentation.

SCB-Enrich solves this problem by integrating:
1. **Multi-dimensional Stratification:** Bins gene families into joint strata based on total family size, reference organism homologs, and annotation degree.
2. **Stratified Null Resampling (10,000 Permutations):** Generates an empirical null distribution preserving the exact covariate composition of the foreground.
3. **Exact Conditional Odds Ratios:** Uses Brent root-finding on the noncentral hypergeometric distribution to obtain exact conditional Maximum Likelihood Estimates and 95% confidence intervals.
4. **Covariate Balance Metrics:** Evaluates the Standardized Mean Difference (SMD < 0.10) to guard against residual covariate imbalance.
5. **Stringency Evidence Classification:** Classifies terms as **Robust** (confirmatory), **Suggestive**, or **Size Artifacts**.
6. **GO DAG Semantic Redundancy Reduction:** Collapses redundant ancestor-descendant terms along the Gene Ontology directed acyclic graph.

---

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

---

## Quick Start

### 1. Command-Line Interface (CLI)

```bash
scb-enrich run \
  --master-table examples/mock_master_table.tsv \
  --annotations examples/mock_annotations.tsv \
  --foreground is_expanded \
  --covariates Total_Count,Count_Ath,Annotation_Degree \
  --obo examples/mock_go.obo \
  --permutations 10000 \
  --out-dir results/
```

### 2. Python API

```python
import pandas as pd
from scbenrich import StratifiedPermuter, fisher_exact_test, exact_conditional_or, fdr_bh

# Parametric Fisher & Exact OR
odds_ratio, p_val = fisher_exact_test(a=15, b=5, c=2, d=18)
mle_or, ci_low, ci_high = exact_conditional_or(a=15, b=5, c=2, d=18)

# Stratified resampling
cov_df = pd.read_csv("master_table.tsv", sep="\t", index_col="Family_ID")
permuter = StratifiedPermuter(cov_df, covariates=["Total_Count", "Count_Ath", "Annotation_Degree"])
```

---

## Inputs and Formats

### Master Table (`--master-table`)
A TSV or CSV table indexed by `Family_ID` containing family counts and covariate metrics:

| Family_ID | Total_Count | Count_Ath | Annotation_Degree | is_expanded |
| :--- | :--- | :--- | :--- | :--- |
| `OG0000001` | 25 | 4 | 12 | 1 |
| `OG0000002` | 18 | 3 | 8 | 1 |
| `OG0000003` | 5 | 1 | 2 | 0 |

### Annotations (`--annotations`)
A TSV or CSV table mapping families to functional terms:

| Family_ID | Term_ID | Term_Name | Category |
| :--- | :--- | :--- | :--- |
| `OG0000001` | `GO:0006950` | response to stress | BP |
| `OG0000001` | `GO:0009408` | response to heat | BP |

---

## Evidence Classification

| Classification | Definition |
| :--- | :--- |
| **Robust** | $q_{local} < 0.05$, $p_{strat} < 0.05$, and residual $\text{SMD} < 0.10$. Fully supported across both parametric and covariate-controlled non-parametric tests. |
| **Suggestive** | $q_{local} < 0.05$ and $p_{strat} < 0.05$, but with marginal residual covariate imbalance ($\text{SMD} \ge 0.10$). |
| **Size Artifact** | $q_{local} < 0.05$, but $p_{strat} \ge 0.05$. Significant under naive Fisher test, but rejected when conditioned on family size / annotation density. |
| **Not Significant** | $q_{local} \ge 0.05$. |

---

## Citation

If you use SCB-Enrich in your research, please cite:
> Cavalheiro, M. F., et al. (2026). Macroevolutionary dynamics and comparative genomics in *Vellozia* (*campos rupestres*). *In preparation*.

---

## License
MIT License. Copyright (c) 2026 Mariana Ferreira Cavalheiro.

