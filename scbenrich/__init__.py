"""
SCB-Enrich: Stratified Covariate-Balanced Enrichment for Gene Families.

A robust framework for functional overrepresentation testing in multi-gene families
that explicitly controls for family size, reference homolog count, and annotation degree
through multi-dimensional stratification, null resampling, and semantic DAG reduction.
"""

__version__ = "0.1.0"
__author__ = "Mariana Ferreira Cavalheiro"

from .core.stats import fdr_bh, fisher_exact_test, exact_conditional_or
from .core.stratified import StratifiedPermuter, compute_smd
from .core.obo import parse_obo, propagate_go_annotations
from .core.dag import reduce_go_redundancy

__all__ = [
    "__version__",
    "fdr_bh",
    "fisher_exact_test",
    "exact_conditional_or",
    "StratifiedPermuter",
    "compute_smd",
    "parse_obo",
    "propagate_go_annotations",
    "reduce_go_redundancy",
]

