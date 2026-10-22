"""Subprocess tests of the installed CLI contract."""
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('permutations,min_universe,status', [(0, 2, 'completed'), (29, 2, 'completed'), (29, 30, 'no_testable_terms')])
def test_cli_subprocess(tmp_path, permutations, min_universe, status):
    output = tmp_path / 'results'
    command = [sys.executable, '-m', 'scbenrich.cli', '--master-table', str(ROOT / 'examples/mock_master_table.tsv'),
               '--annotations', str(ROOT / 'examples/mock_annotations.tsv'), '--foreground', 'is_expanded',
               '--obo', str(ROOT / 'examples/mock_go.obo'), '--permutations', str(permutations),
               '--min-universe', str(min_universe), '--max-universe-fraction', '.95', '--out-dir', str(output)]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    manifest = json.loads((output / 'run_metadata.json').read_text())
    assert manifest['status'] == status
    assert manifest['arguments']['foreground'] == 'is_expanded'
    frame = pd.read_csv(output / 'enrichment_results_all.tsv', sep='\t')
    assert 'stratified_significant' in frame
    if status == 'no_testable_terms':
        assert frame.empty
    if permutations == 0:
        assert not frame.stratified_significant.any()


def test_cli_go_error_preserves_failed_metadata(tmp_path):
    ann = tmp_path / 'annotations.tsv'
    pd.DataFrame({'Family_ID': ['OG0000001'], 'Term_ID': ['GO:9999999']}).to_csv(ann, sep='\t', index=False)
    output = tmp_path / 'results'
    command = [sys.executable, '-m', 'scbenrich.cli', '--master-table', str(ROOT / 'examples/mock_master_table.tsv'),
               '--annotations', str(ann), '--foreground', 'is_expanded', '--obo', str(ROOT / 'examples/mock_go.obo'),
               '--out-dir', str(output)]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert result.returncode != 0
    assert 'unknown GO term' in result.stderr
    manifest = json.loads((output / 'run_metadata.json').read_text())
    assert manifest['status'] == 'failed'
    assert manifest['error']['type'] == 'ValueError'
