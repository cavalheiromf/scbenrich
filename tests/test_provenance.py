from argparse import Namespace
import json
import pytest
from scbenrich.cli import run_pipeline


def test_nonempty_output_is_preserved(tmp_path):
    keep = tmp_path / 'keep.txt'
    keep.write_text('original')
    args = Namespace(permutations=10, out_dir=str(tmp_path))
    with pytest.raises(ValueError, match='empty'):
        run_pipeline(args)
    assert keep.read_text() == 'original'


def test_failed_execution_writes_manifest(tmp_path):
    out = tmp_path / 'out'
    args = Namespace(permutations=10, out_dir=str(out), master_table=str(tmp_path / 'missing'),
                     annotations=str(tmp_path / 'missing2'), foreground='selected', min_universe=1,
                     max_universe_fraction=.9, covariates='size', obo=None, seed=1)
    with pytest.raises(FileNotFoundError):
        run_pipeline(args)
    metadata = json.loads((out / 'run_metadata.json').read_text())
    assert metadata['status'] == 'failed'
    assert metadata['arguments']['foreground'] == 'selected'
    assert metadata['error']['type'] == 'FileNotFoundError'
    assert metadata['code_hashes']
