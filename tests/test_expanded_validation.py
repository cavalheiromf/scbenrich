import json
import numpy as np
import pandas as pd
import pytest
from benchmarks.expanded import generate, run_validation


def test_generator_truth_and_dependency():
    data = generate(np.random.default_rng(4), 'dependent', 80, 4)
    assert data['covariates'].shape == (80, 3)
    assert (data['terms']['parent'] >= data['terms']['child_0']).all()
    assert data['truth'].loc['effect', 'conditional_null'] == False
    derived = generate(np.random.default_rng(4), 'derived_degree', 80, 1)
    assert derived['truth'].conditional_null.isna().all()


def test_resume_and_scenario_order_independence(tmp_path):
    config = dict(seed=17, replicates=2, permutations=19, n_families=60,
                  bins=[2], odds_ratios=[1, 4], scenarios=['exact', 'missing_size'], batch_size=1)
    raw, summary = run_validation(config, tmp_path)
    times = {p.name: p.stat().st_mtime_ns for p in tmp_path.glob('batch_*.json')}
    again, _ = run_validation(config, tmp_path)
    pd.testing.assert_frame_equal(raw, again)
    assert times == {p.name: p.stat().st_mtime_ns for p in tmp_path.glob('batch_*.json')}
    reverse = dict(config, scenarios=list(reversed(config['scenarios'])))
    other, _ = run_validation(reverse, tmp_path / 'other')
    columns = ['scenario', 'odds_ratio', 'replicate', 'bins', 'method', 'fdp', 'effect_detected']
    pd.testing.assert_frame_equal(raw[columns].sort_values(columns[:5]).reset_index(drop=True), other[columns].sort_values(columns[:5]).reset_index(drop=True))
    assert {'fdr_ci_low', 'fdr_ci_high', 'power_ci_low', 'power_ci_high'} <= set(summary)
    with pytest.raises(ValueError, match='configuration'):
        run_validation(dict(config, permutations=29), tmp_path)


def test_discrete_exact_strata_and_thinning_truth():
    from scbenrich.core.stratified import StratifiedPermuter
    data = generate(np.random.default_rng(1), 'exact', 80, 1)
    permuter = StratifiedPermuter(data['covariates'], bins=2)
    for _, group in permuter.df.groupby('_stratum'):
        assert (group[data['covariates'].columns].nunique() == 1).all()
    missing = generate(np.random.default_rng(1), 'missing_size', 80, 1)
    assert missing['truth'].marginal_null.isna().all()
    assert missing['truth'].conditional_null.all()


def test_calibration_does_not_treat_terms_as_independent_replicates():
    from benchmarks.report import calibration
    terms = pd.DataFrame({'scenario': ['exact'] * 3, 'odds_ratio': [1] * 3,
                          'bins': [2] * 3, 'replicate': [0, 0, 1],
                          'p_fisher': [.01, .9, .01], 'p_strat': [.01, .9, .01],
                          'marginal_null': [True] * 3, 'conditional_null': [True] * 3})
    result = calibration(terms)
    assert result.loc[result.threshold == .05, 'observed_rate'].eq(.75).all()
    assert result.replicates.eq(2).all()
