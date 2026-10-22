import pandas as pd
import pytest
from benchmarks.comparators.run import export_inputs, compare_tables


def test_export_and_mismatched_counts(tmp_path):
    matrix = pd.DataFrame({'GO:0006950': [1, 1, 0, 0]}, index=['A', 'B', 'C', 'D'])
    reference = export_inputs(matrix, {'A', 'C'}, tmp_path, min_count=1, max_fraction=1)
    assert reference.loc[0, 'foreground_count'] == 1
    assert reference.loc[0, 'universe_count'] == 2
    assert (tmp_path / 'universe.txt').read_text().splitlines() == list('ABCD')
    external = reference[['Term_ID', 'foreground_count', 'universe_count']].copy()
    external['p_value'] = reference.p_greater
    result = compare_tables(reference, external, 'greater')
    assert result.counts_match.all() and result.p_match.all()
    external.loc[0, 'universe_count'] = 1
    result = compare_tables(reference, external, 'greater')
    assert not result.counts_match.any()
    assert not result.p_match.any()


def test_goatools_matches_two_sided_with_unannotated_families(tmp_path):
    pytest.importorskip('goatools')
    from pathlib import Path
    from benchmarks.comparators.run import run_goatools
    matrix = pd.DataFrame({'GO:0006950': [1, 1, 0, 0]}, index=['A', 'B', 'C', 'D'])
    reference = export_inputs(matrix, {'A'}, tmp_path, min_count=1, max_fraction=1)
    external = run_goatools(tmp_path, Path(__file__).resolve().parents[1] / 'examples/mock_go.obo')
    audit = compare_tables(reference, external, 'two-sided')
    assert audit.counts_match.all()
    assert audit.p_match.all()
    assert audit.q_match.all()
