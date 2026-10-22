"""
Covariate-adjusted logistic regression model for sensitivity audits of enriched terms.
"""

import numpy as np
from scipy.optimize import minimize
import scipy.stats as stats


def fit_covariate_model(y, fg_vec, covariate_matrix):
    """
    Fits a covariate-adjusted logistic regression model:
        logit(P(Y=1)) = beta_0 + beta_fg * FG + sum(beta_k * log1p(covariate_k))
        
    Parameters
    ----------
    y : np.ndarray
        Binary vector indicating whether each family carries the functional term.
    fg_vec : np.ndarray
        Binary indicator of foreground membership (1 if foreground, 0 if background).
    covariate_matrix : np.ndarray
        2D array of continuous covariates (e.g. Total_Count, Ath_Count, Annotation_Degree).
        
    Returns
    -------
    tuple (float, float)
        (adjusted_odds_ratio, p_likelihood_ratio_test)
    """
    y = np.asarray(y, dtype=float)
    fg_vec = np.asarray(fg_vec, dtype=float)
    covariates = np.asarray(covariate_matrix, dtype=float)
    n = len(y)
    
    if y.sum() == 0 or y.sum() == n:
        return 1.0, 1.0

    # Log1p transformation and standardization of continuous covariates
    log_covs = np.log1p(np.maximum(0, covariates))
    standardized_covs = np.zeros_like(log_covs)
    for col in range(log_covs.shape[1]):
        std = np.std(log_covs[:, col])
        if std > 0:
            standardized_covs[:, col] = (log_covs[:, col] - np.mean(log_covs[:, col])) / std

    X_full = np.column_stack([np.ones(n), fg_vec, standardized_covs])
    X_null = np.column_stack([np.ones(n), standardized_covs])

    def nll(w, X):
        z = X.dot(w)
        p = 1.0 / (1.0 + np.exp(-np.clip(z, -25, 25)))
        eps = 1e-12
        return -np.sum(y * np.log(p + eps) + (1.0 - y) * np.log(1.0 - p + eps))

    try:
        res_full = minimize(nll, np.zeros(X_full.shape[1]), args=(X_full,), method="BFGS")
        res_null = minimize(nll, np.zeros(X_null.shape[1]), args=(X_null,), method="BFGS")
        
        beta_fg = res_full.x[1]
        adj_or = float(np.exp(beta_fg))
        
        # Likelihood ratio test statistic
        lrt_stat = max(0.0, 2.0 * (res_null.fun - res_full.fun))
        p_lrt = float(stats.chi2.sf(lrt_stat, df=1))
        
        # One-sided adjustment if beta_fg < 0
        if beta_fg < 0:
            p_lrt = 1.0 - p_lrt / 2.0
        else:
            p_lrt = p_lrt / 2.0
            
        return adj_or, p_lrt
    except Exception:
        return 1.0, 1.0

