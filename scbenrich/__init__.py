"""
SCB-Enrich: Stratified Covariate-Balanced Enrichment for Gene Families.

Exploratory family-level enrichment using covariate strata. Statistical validity
requires exchangeability within the selected strata.
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

