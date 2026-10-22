import pandas as pd
from scbenrich.core.stratified import StratifiedPermuter


def test_constant_term_has_support_and_tail_interval():
    cov = pd.DataFrame({'x': [1, 1, 1, 1]}, index=list('ABCD'))
    matrix = pd.DataFrame({'constant': [1, 1, 1, 1], 'variable': [1, 1, 0, 0]}, index=cov.index)
    result = StratifiedPermuter(cov, ['x'], seed=1).run_analysis({'A', 'C'}, matrix, 20)
    term = result['terms'].set_index('Term_ID').loc['constant']
    assert term.analysis_status == 'constant_statistic'
    assert term.support_min == term.support_max == 2
    assert term.p_strat == 1
    assert term.exceedances == 20
    assert len(result['balance']) == 1


def test_deterministic_strata_do_not_consume_random_draws():
    import numpy as np

    class InformativeOnlyRng:
        def __init__(self):
            self.generator = np.random.default_rng(42)

        def choice(self, pool, size, replace):
            assert 0 < size < len(pool), 'deterministic stratum consumed randomness'
            return self.generator.choice(pool, size=size, replace=replace)

    cov = pd.DataFrame({'x': [0, 0, 1, 1, 2, 2]}, index=list('ABCDEF'))
    matrix = pd.DataFrame({'term': [1, 0, 1, 0, 1, 0]}, index=cov.index)
    permuter = StratifiedPermuter(cov, ['x'], bins=3)
    permuter.rng = InformativeOnlyRng()
    result = permuter.run_analysis({'A', 'B', 'C'}, matrix, 999)
    term = result['terms'].iloc[0]
    assert term.support_min == 1
    assert term.support_max == 2
    assert abs(term.p_strat - .5) < .06


def test_joint_sampling_preserves_exact_tiny_universe_distribution():
    cov = pd.DataFrame({'x': [0, 0, 0, 0]}, index=list('ABCD'))
    matrix = pd.DataFrame({'first': [1, 1, 0, 0], 'duplicate': [1, 1, 0, 0],
                           'complement': [0, 0, 1, 1]}, index=cov.index)
    result = StratifiedPermuter(cov, ['x'], seed=7).run_analysis({'A', 'B'}, matrix, 9999)
    terms = result['terms'].set_index('Term_ID')
    # Of the six equiprobable pairs, only AB reaches the observed first count.
    assert abs(terms.loc['first', 'p_strat'] - 1 / 6) < .025
    assert terms.loc['first', 'exceedances'] == terms.loc['duplicate', 'exceedances']
    assert terms.loc['complement', 'p_strat'] == 1
