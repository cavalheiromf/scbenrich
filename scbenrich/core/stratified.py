"""
Multi-dimensional covariate stratification, balance metrics (SMD), and stratified null resampling.
"""

import numpy as np
import pandas as pd
from collections import defaultdict
from math import comb
from scbenrich.utils.validation import validate_family_inputs, validate_parameters


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
        return float("nan")
    if s_pool is None:
        n1 = len(x1)
        n0 = len(x0)
        v1 = np.var(x1, ddof=1) if n1 > 1 else 0.0
        v0 = np.var(x0, ddof=1) if n0 > 1 else 0.0
        denom = (n1 - 1) + (n0 - 1)
        if denom <= 0:
            return float("nan")
        s_pool = np.sqrt(((n1 - 1) * v1 + (n0 - 1) * v0) / denom)
    if np.isnan(s_pool):
        return float("nan")
    if s_pool == 0:
        return 0.0 if np.mean(x1) == np.mean(x0) else float("nan")
    return float((np.mean(x1) - np.mean(x0)) / s_pool)


class StratifiedPermuter:
    """
    Manages multi-dimensional covariate binning and stratified resampling
    under an assumed exchangeable selection mechanism within strata.
    """

    def __init__(self, covariate_df, covariates=None, seed=42, bins=None):
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
        validate_parameters(0, bins)
        self.df = covariate_df.copy()
        self.covariates = covariates or list(self.df.select_dtypes(include=[np.number]).columns)
        validate_family_inputs(self.df[self.covariates])
        self.seed = seed
        self.bins = bins
        self.rng = np.random.default_rng(seed)
        self.strata_ = None
        self._build_strata(n_bins_dict={c: bins for c in self.covariates} if bins is not None else None)

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
        if any(not isinstance(v, int) or v < 2 for v in n_bins_dict.values()):
            raise ValueError("bins must be an integer >= 2")

        self.bin_edges_ = {}
        stratum_labels = []
        for col in self.covariates:
            n_b = n_bins_dict[col]
            if self.df[col].nunique() == 1:
                self.bin_edges_[col] = [float(self.df[col].iloc[0])]
                stratum_labels.append(pd.Series("constant", index=self.df.index))
                continue
            # Use qcut with duplicates='drop'
            try:
                binned, edges = pd.qcut(self.df[col], q=n_b, duplicates="drop", retbins=True)
            except Exception:
                binned, edges = pd.cut(self.df[col], bins=n_b, retbins=True)
            self.bin_edges_[col] = edges.tolist()
            stratum_labels.append(binned.astype(str))

        joint_strata = stratum_labels[0]
        for sl in stratum_labels[1:]:
            joint_strata = joint_strata + "|" + sl
        self.df["_stratum"] = joint_strata
        self.strata_ = self.df["_stratum"]

    def support_diagnostics(self, foreground_ids):
        selected = self.df.index.isin(foreground_ids)
        rows = []
        for label, group in self.df.groupby("_stratum", sort=True):
            n = len(group)
            k = int(selected[self.df.index.get_indexer(group.index)].sum())
            rows.append({"stratum": label, "eligible": n, "foreground": k,
                         "background": n - k, "combinations": comb(n, k),
                         "informative": 0 < k < n})
        return pd.DataFrame(rows)

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
        validate_parameters(n_permutations, self.bins)
        validate_family_inputs(self.df[self.covariates], foreground_ids, term_matrix)
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
        
        lower = np.zeros(n_terms, dtype=int)
        upper = np.zeros(n_terms, dtype=int)
        for label, pool in stratum_to_indices.items():
            k = stratum_fg_counts.get(label, 0)
            annotated = tm[pool].sum(axis=0)
            lower += np.maximum(0, k - (len(pool) - annotated)).astype(int)
            upper += np.minimum(k, annotated).astype(int)
        cov_values = self.df[self.covariates].to_numpy(dtype=float)
        null_means = []
        # Fully selected strata contribute fixed counts and covariate sums.
        # Only informative strata need random choices. Each choice is shared by
        # all terms, preserving their joint randomization distribution.
        fixed_indices = []
        random_strata = []
        for label, count in stratum_fg_counts.items():
            pool = stratum_to_indices[label]
            if count == len(pool):
                fixed_indices.extend(pool)
            elif count:
                random_strata.append((pool, count))
        fixed_indices = np.asarray(fixed_indices, dtype=int)
        fixed_counts = tm[fixed_indices].sum(axis=0)
        fixed_cov_sum = cov_values[fixed_indices].sum(axis=0)
        foreground_count = int(fg_mask.sum())
        exceed_counts = np.zeros(n_terms, dtype=int)

        for _ in range(n_permutations):
            samples = [self.rng.choice(pool, size=count, replace=False)
                       for pool, count in random_strata]
            sample_indices = np.concatenate(samples) if samples else np.empty(0, dtype=int)
            null_means.append((fixed_cov_sum + cov_values[sample_indices].sum(axis=0))
                              / foreground_count)
            perm_counts = fixed_counts + tm[sample_indices].sum(axis=0)
            exceed_counts += (perm_counts >= obs_counts)

        empirical_p = (exceed_counts + 1.0) / (n_permutations + 1.0)
        
        from scipy.stats import binomtest
        support = self.support_diagnostics(foreground_ids)
        statuses = ["no_randomization" if not support.informative.any() else
                    ("constant_statistic" if lo == hi else "ok") for lo, hi in zip(lower, upper)]
        intervals = [binomtest(int(b), n_permutations).proportion_ci() for b in exceed_counts] if n_permutations else []
        self.term_diagnostics_ = pd.DataFrame({
            "Term_ID": term_matrix.columns, "observed_count": obs_counts,
            "support_min": lower, "support_max": upper, "analysis_status": statuses,
            "exceedances": exceed_counts, "permutations": n_permutations,
            "p_strat": empirical_p, "mc_resolution": 1 / (n_permutations + 1),
            "tail_ci_low": [ci.low for ci in intervals] if intervals else np.nan,
            "tail_ci_high": [ci.high for ci in intervals] if intervals else np.nan,
        })
        self.balance_diagnostics_ = pd.DataFrame([
            {"covariate": col, "observed_mean": cov_values[fg_mask, i].mean(),
             "null_mean": np.mean(np.asarray(null_means)[:, i]) if null_means else np.nan,
             "null_mean_q025": np.quantile(np.asarray(null_means)[:, i], .025) if null_means else np.nan,
             "null_mean_q975": np.quantile(np.asarray(null_means)[:, i], .975) if null_means else np.nan}
            for i, col in enumerate(self.covariates)])
        results = {}
        for idx, term in enumerate(term_matrix.columns):
            results[term] = float(empirical_p[idx])
            
        return results


    def run_analysis(self, foreground_ids, term_matrix, n_permutations=10000):
        """Return term support, tail uncertainty, strata and null-mean diagnostics."""
        self.run_permutations(foreground_ids, term_matrix, n_permutations)
        return {"terms": self.term_diagnostics_.copy(),
                "strata": self.support_diagnostics(foreground_ids),
                "balance": self.balance_diagnostics_.copy()}
