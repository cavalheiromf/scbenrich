import numpy as np
import pytest

from benchmarks.simulate import evaluate_one, run_benchmark, simulate_one


def test_simulation_has_known_term_only_in_positive_scenarios():
    rng = np.random.default_rng(9)
    for scenario in ("simple_null", "confounded_null", "simple_positive", "confounded_positive"):
        cov, foreground, terms, effect = simulate_one(rng, scenario, n_families=40, n_null_terms=3)
        assert len(cov) == len(terms) == 40
        assert 0 < len(foreground) < 40
        assert effect == ("effect" if "positive" in scenario else None)


def test_benchmark_is_reproducible_and_summary_matches_replicates():
    raw, summary = run_benchmark(replicates=3, permutations=29, bins=(2, 4), seed=11,
                                 n_families=40, n_null_terms=3)
    raw_again, summary_again = run_benchmark(replicates=3, permutations=29, bins=(2, 4), seed=11,
                                             n_families=40, n_null_terms=3)
    assert raw.equals(raw_again)
    assert summary.equals(summary_again)
    assert len(raw) == 4 * 3 * 2 * 2
    assert summary.replicates.eq(3).all()
    assert summary.probability_any_false_discovery.between(0, 1).all()
    assert summary.empirical_fdr.between(0, 1).all()
    assert summary.loc[summary.scenario == "simple_positive", "ci_coverage"].notna().any()
    assert summary.loc[summary.scenario == "confounded_positive", "ci_coverage"].isna().all()


def test_bin_parameter_validation():
    with pytest.raises(ValueError, match="bin counts"):
        run_benchmark(replicates=1, permutations=3, bins=(1,))


def test_reference_truth_distinguishes_marginal_and_conditional_nulls():
    from benchmarks.simulate import term_truth, _expit
    truth = term_truth('confounded_null', ['null_0'])
    assert not truth.loc[0, 'marginal_null']
    assert truth.loc[0, 'conditional_null']
    # Population covariance is positive because both probabilities increase with size.
    x = (np.arange(10000) + .5) / 10000
    selection = _expit(-2 + 3 * x)
    annotation = _expit(-2.5 + 3 * x)
    assert np.mean(selection * annotation) - selection.mean() * annotation.mean() > 0


def test_fisher_confounded_discoveries_are_not_marginal_false_positives():
    raw, summary = run_benchmark(replicates=2, permutations=19, bins=(2,), seed=7,
                                 n_families=40, n_null_terms=3)
    fisher = raw[(raw.scenario == 'confounded_null') & (raw.method == 'fisher')]
    assert fisher.false_discoveries.eq(0).all()
    assert fisher.n_null_terms.eq(0).all()
    assert fisher.conditional_target_mismatch.equals(fisher.discoveries)
    assert 'probability_any_false_discovery' in summary
    assert 'rejection_rate' not in summary
