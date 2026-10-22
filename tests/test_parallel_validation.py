import pandas as pd
from benchmarks.expanded import run_validation
from benchmarks.parallel import run_parallel


def test_parallel_matches_sequential_and_reuses_checkpoints(tmp_path):
    config = dict(seed=42, replicates=2, permutations=19, n_families=24,
                  bins=[2], odds_ratios=[1], scenarios=['exact', 'unconfounded'], batch_size=1)
    expected, summary = run_validation(config, tmp_path / 'sequential')
    actual, parallel_summary = run_parallel(config, tmp_path / 'parallel', workers=2)
    pd.testing.assert_frame_equal(expected, actual)
    pd.testing.assert_frame_equal(summary, parallel_summary)
    checkpoint = next((tmp_path / 'parallel').glob('batch_*.json'))
    before = checkpoint.stat().st_mtime_ns
    resumed, _ = run_parallel(config, tmp_path / 'parallel', workers=2)
    pd.testing.assert_frame_equal(actual, resumed)
    assert checkpoint.stat().st_mtime_ns == before
