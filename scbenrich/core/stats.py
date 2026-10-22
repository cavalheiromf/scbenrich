"""
Statistical tests, exact odds ratios, and false discovery rate corrections.
"""

import numpy as np
import scipy.stats as stats
from scipy.optimize import brentq
from scipy.special import gammaln


def fdr_bh(pvals):
    """
    Benjamini-Hochberg False Discovery Rate (FDR) adjustment.
    
    Parameters
    ----------
    pvals : array-like
        Array of nominal p-values.
        
    Returns
    -------
    np.ndarray
        Array of adjusted q-values maintaining original order.
    """
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ranked_p = pvals[order]
    qvals = np.zeros(n)
    accum = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked_p[i] * n / rank
        accum = min(accum, val)
        qvals[i] = accum
    orig_qvals = np.zeros(n)
    orig_qvals[order] = qvals
    return orig_qvals


def log_comb(n, k):
    """Logarithm of binomial coefficient C(n, k)."""
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


def exact_conditional_or(a, b, c, d, alpha=0.05):
    """Conditional 2x2 odds ratio and equal-tailed exact confidence interval.

    Returns NaN estimates for tables with no identifiable odds ratio. This is
    conditional on table margins, not on family covariates.
    """
    cells = (a, b, c, d)
    if any(isinstance(v, bool) or not isinstance(v, (int, np.integer)) or v < 0 for v in cells):
        raise ValueError("table counts must be nonnegative integers")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be in (0, 1)")
    from scipy.stats.contingency import odds_ratio

    result = odds_ratio([[a, b], [c, d]], kind="conditional")
    if np.isnan(result.statistic):
        return float("nan"), float("nan"), float("nan")
    interval = result.confidence_interval(confidence_level=1 - alpha)
    return float(result.statistic), float(interval.low), float(interval.high)


def fisher_exact_test(a, b, c, d, alternative="greater"):
    """
    Computes Fisher's exact test on a 2x2 contingency table.
    
    Parameters
    ----------
    a, b, c, d : int
        Contingency table counts.
    alternative : str, optional
        'greater' (overrepresentation, default), 'less', or 'two-sided'.
        
    Returns
    -------
    tuple (float, float)
        (odds_ratio, p_value)
    """
    table = [[a, b], [c, d]]
    res = stats.fisher_exact(table, alternative=alternative)
    return float(res[0]), float(res[1])

