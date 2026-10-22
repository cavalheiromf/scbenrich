import pytest
from scbenrich.core.inference import decision_fields


def test_significance_independent_of_balance():
    result = decision_fields(0.8, 0.01, 'ok', 'undefined')
    assert result['stratified_significant'] is True
    assert result['marginal_significant'] is False
    assert result['balance_status'] == 'undefined'
    assert not decision_fields(.01, .01, 'no_randomization', 'undefined')['stratified_significant']


def test_invalid_alpha():
    with pytest.raises(ValueError):
        decision_fields(.01, .01, 'ok', 'undefined', alpha=0)
