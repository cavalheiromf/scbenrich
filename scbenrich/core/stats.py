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
    """
    Exact Conditional Maximum Likelihood Estimation of Odds Ratio and 95% Confidence Interval.
    Solves noncentral hypergeometric distribution via Brent root finding.
    
    Contingency table:
                Term+    Term-
    Foreground    a        b     | n1 = a + b
    Background    c        d     | n2 = c + d
    ----------------------------
    Total         m1       m2    | N = a + b + c + d
    
    Parameters
    ----------
    a, b, c, d : int
        Cell counts in the 2x2 contingency table.
    alpha : float, optional
        Significance level for confidence interval (default: 0.05 for 95% CI).
        
    Returns
    -------
    tuple (float, float, float)
        (point_estimate, ci_low, ci_high)
    """
    n1 = a + b
    n2 = c + d
    m1 = a + c
    x = a
    
    x_min = max(0, m1 - n2)
    x_max = min(n1, m1)
    support = np.arange(x_min, x_max + 1)
    
    if len(support) == 1:
        return 1.0, 0.0, float("inf")
        
    log_weights = np.array([log_comb(n1, k) + log_comb(n2, m1 - k) for k in support])
    max_lw = np.max(log_weights)
    weights = np.exp(log_weights - max_lw)
    
    # Boundary cases: degenerate tables
    if x == x_min:
        mle_or = 0.0
    elif x == x_max:
        mle_or = float("inf")
    else:
        def mean_diff(log_theta):
            theta = np.exp(log_theta)
            pw = weights * (theta ** (support - x_min))
            p = pw / np.sum(pw)
            return np.sum(support * p) - x
            
        try:
            log_theta_mle = brentq(mean_diff, -20.0, 20.0)
            mle_or = float(np.exp(log_theta_mle))
        except Exception:
            mle_or = float((a * d) / (b * c)) if (b * c) > 0 else float("inf")
            
    # Exact conditional confidence interval
    def p_ge(theta):
        pw = weights * (theta ** (support - x_min))
        s = np.sum(pw)
        if s == 0:
            return 1.0
        return np.sum(pw[support >= x]) / s
        
    def p_le(theta):
        pw = weights * (theta ** (support - x_min))
        s = np.sum(pw)
        if s == 0:
            return 1.0
        return np.sum(pw[support <= x]) / s
        
    target_low = alpha / 2.0
    target_high = 1.0 - alpha / 2.0
    
    # Lower bound
    if x == x_min:
        ci_low = 0.0
    else:
        try:
            ci_low = float(brentq(lambda t: p_ge(t) - target_low, 1e-6, 1e4))
        except Exception:
            ci_low = 0.0
            
    # Upper bound
    if x == x_max:
        ci_high = float("inf")
    else:
        try:
            ci_high = float(brentq(lambda t: p_le(t) - target_low, 1e-4, 1e6))
        except Exception:
            ci_high = float("inf")
            
    return mle_or, ci_low, ci_high


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

