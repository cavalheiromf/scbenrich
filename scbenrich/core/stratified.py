"""
Multi-dimensional covariate stratification, balance metrics (SMD), and stratified null resampling.
"""

import numpy as np
import pandas as pd
from collections import defaultdict


def compute_smd(x1, x0, s_pool=None):
    """
    Computes Standardized Mean Difference (Cohen's d) between two groups.
    
    Parameters
    ----------
    x1 : array-like
        Foreground group values.
    x0 : array-like
        Background group values.
    s_pool : float, optional
        Pre-computed pooled standard deviation. If None, calculated directly.
        
    Returns
    -------
    float
        Standardized mean difference.
    """
    x1 = np.asarray(x1, dtype=float)
    x0 = np.asarray(x0, dtype=float)
    if len(x1) == 0 or len(x0) == 0:
        return 0.0
    if s_pool is None:
        n1 = len(x1)
        n0 = len(x0)
        v1 = np.var(x1, ddof=1) if n1 > 1 else 0.0
        v0 = np.var(x0, ddof=1) if n0 > 1 else 0.0
        denom = (n1 - 1) + (n0 - 1)
        if denom <= 0:
            return 0.0
        s_pool = np.sqrt(((n1 - 1) * v1 + (n0 - 1) * v0) / denom)
    if s_pool == 0 or np.isnan(s_pool):
        return 0.0
    return float((np.mean(x1) - np.mean(x0)) / s_pool)


class StratifiedPermuter:
    """
    Manages multi-dimensional covariate binning and stratified resampling
    to generate rigorous null distributions for gene-family enrichment testing.
    """

    def __init__(self, covariate_df, covariates=None, seed=42):
        """
        Parameters
        ----------
        covariate_df : pd.DataFrame
            DataFrame indexed by family_id with covariate columns.
        covariates : list of str, optional
            List of column names to stratify on. Defaults to all numeric columns.
        seed : int, optional
            Random number generator seed for reproducible resampling.
        """
        self.df = covariate_df.copy()
        self.covariates = covariates or list(self.df.select_dtypes(include=[np.number]).columns)
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.strata_ = None
        self._build_strata()

    def _build_strata(self, n_bins_dict=None):
        """
        Partitions gene families into joint multi-dimensional strata bins.
        """
        if n_bins_dict is None:
            # Default heuristic: up to 10 quantiles per covariate, adjusted for unique values
            n_bins_dict = {}
            for col in self.covariates:
                nunique = self.df[col].nunique()
                n_bins_dict[col] = min(10, max(2, nunique // 5))

        stratum_labels = []
        for col in self.covariates:
            n_b = n_bins_dict[col]
            # Use qcut with duplicates='drop'
            try:
                binned = pd.qcut(self.df[col], q=n_b, duplicates="drop")
            except Exception:
                binned = pd.cut(self.df[col], bins=n_b)
            stratum_labels.append(binned.astype(str))

        joint_strata = stratum_labels[0]
        for sl in stratum_labels[1:]:
            joint_strata = joint_strata + "|" + sl
        self.df["_stratum"] = joint_strata
        self.strata_ = self.df["_stratum"]

    def run_permutations(self, foreground_ids, term_matrix, n_permutations=10000):
        """
        Performs stratified resampling preserving exact stratum counts of the foreground.
        
        Parameters
        ----------
        foreground_ids : set or list
            Set of family IDs in the foreground.
        term_matrix : pd.DataFrame
            Binary presence/absence DataFrame (families x terms).
        n_permutations : int, optional
            Number of permutations (default: 10,000).
            
        Returns
        -------
        dict
            Mapping of term -> empirical p-value (p_strat).
        """
        all_families = self.df.index.values
        fg_mask = np.isin(all_families, list(foreground_ids))
        fg_strata = self.df["_stratum"].values[fg_mask]
        
        # Count foreground families per stratum
        stratum_fg_counts = pd.Series(fg_strata).value_counts().to_dict()
        
        # Index families per stratum across the whole universe
        stratum_to_indices = defaultdict(list)
        for idx, s in enumerate(self.df["_stratum"].values):
            stratum_to_indices[s].append(idx)
        for s in stratum_to_indices:
            stratum_to_indices[s] = np.array(stratum_to_indices[s])

        # Observed counts for terms in foreground
        common_fams = [f for f in all_families if f in term_matrix.index]
        tm = term_matrix.reindex(all_families, fill_value=0).values
        n_terms = tm.shape[1]
        
        obs_counts = np.sum(tm[fg_mask, :], axis=0)
        
        # Resampling simulation
        exceed_counts = np.zeros(n_terms, dtype=int)
        
        for _ in range(n_permutations):
            sample_indices = []
            for s, req_count in stratum_fg_counts.items():
                pool = stratum_to_indices[s]
                # If required count exceeds pool size (should rarely happen), take all
                c = min(req_count, len(pool))
                chosen = self.rng.choice(pool, size=c, replace=False)
                sample_indices.extend(chosen)
                
            sample_indices = np.array(sample_indices)
            perm_counts = np.sum(tm[sample_indices, :], axis=0)
            exceed_counts += (perm_counts >= obs_counts)

        empirical_p = (exceed_counts + 1.0) / (n_permutations + 1.0)
        
        results = {}
        for idx, term in enumerate(term_matrix.columns):
            results[term] = float(empirical_p[idx])
            
        return results
