import numpy as np
import pandas as pd
import pytest
from scbenrich.utils.validation import validate_family_inputs, validate_parameters


def test_matrix_must_be_binary_and_aligned():
    cov = pd.DataFrame({'size': [1, 2]}, index=['A', 'B'])
    with pytest.raises(ValueError, match='binary'):
        validate_family_inputs(cov, {'A'}, pd.DataFrame({'T': [1, 2]}, index=cov.index))
    with pytest.raises(ValueError, match='IDs'):
        validate_family_inputs(cov, {'A'}, pd.DataFrame({'T': [1, 0]}, index=['A', 'C']))


@pytest.mark.parametrize('permutations', [-1, 1.5, True])
def test_invalid_permutations(permutations):
    with pytest.raises(ValueError):
        validate_parameters(permutations, None, .05)
