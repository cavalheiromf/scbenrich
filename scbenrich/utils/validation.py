"""Shared input contracts for family randomization."""
from numbers import Integral
import numpy as np
import pandas as pd


def validate_parameters(permutations, bins=None, alpha=.05):
    if isinstance(permutations, bool) or not isinstance(permutations, Integral) or permutations < 0:
        raise ValueError('permutations must be a nonnegative integer')
    if bins is not None and (isinstance(bins, bool) or not isinstance(bins, Integral) or bins < 2):
        raise ValueError('bins must be an integer >= 2')
    if not 0 < alpha < 1:
        raise ValueError('alpha must be in (0, 1)')


def validate_family_inputs(covariates, foreground_ids=None, term_matrix=None):
    ids = covariates.index
    if ids.has_duplicates or ids.isna().any() or any(not str(x).strip() for x in ids):
        raise ValueError('Family IDs must be unique and nonempty')
    if covariates.empty or not len(covariates.columns):
        raise ValueError('at least one family and covariate are required')
    for col in covariates:
        if not pd.api.types.is_numeric_dtype(covariates[col]) or not np.isfinite(covariates[col].to_numpy(dtype=float)).all():
            raise ValueError('covariates must contain finite numeric values')
    if foreground_ids is not None:
        foreground = set(foreground_ids)
        if not foreground or not foreground < set(ids):
            raise ValueError('foreground must be a nonempty proper subset of family IDs')
    if term_matrix is not None:
        if term_matrix.index.has_duplicates or set(term_matrix.index) != set(ids):
            raise ValueError('matrix family IDs must match covariate IDs exactly')
        if term_matrix.columns.has_duplicates:
            raise ValueError('term IDs must be unique')
        if not term_matrix.isin([0, 1]).all().all():
            raise ValueError('term matrix must be binary without missing values')
