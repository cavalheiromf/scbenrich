import json
import pandas as pd
import pytest
from benchmarks.expanded import run_validation
from benchmarks.parallel import prepare, collect
from benchmarks.slurm import run_slice, task_parameters


def config(**changes):
    return dict(dict(seed=42, replicates=2, permutations=9, n_families=24,
                     bins=[2], odds_ratios=[1], scenarios=['exact'], batch_size=1), **changes)


def test_task_rejects_different_configuration(tmp_path):
    c = config()
    prepare(c, tmp_path)
    with pytest.raises(ValueError, match='configuration|manifest'):
        run_slice(dict(c, permutations=19), tmp_path, 'exact', 1, 0, 2)
    assert not list(tmp_path.rglob('batch_*.json'))


def test_task_requires_prepared_manifest(tmp_path):
    with pytest.raises(ValueError, match='manifest|prepare'):
        run_slice(config(), tmp_path, 'exact', 1, 0, 2)


@pytest.mark.parametrize('damage', ['versions', 'source_hashes'])
def test_task_rejects_changed_environment(tmp_path, damage):
    c = config()
    prepare(c, tmp_path)
    path = tmp_path / 'config.json'
    manifest = json.loads(path.read_text())
    manifest[damage] = {}
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='configuration|manifest'):
        run_slice(c, tmp_path, 'exact', 1, 0, 2)


@pytest.mark.parametrize('damage', ['duplicate_replica', 'missing_term', 'duplicate_term', 'wrong_permutations'])
def test_rejects_invalid_checkpoint_on_collect_and_resume(tmp_path, damage):
    c = config()
    prepare(c, tmp_path)
    run_slice(c, tmp_path, 'exact', 1, 0, 2)
    paths = sorted((tmp_path / 'shards/exact').glob('batch_*.json'))
    payload = json.loads(paths[1].read_text())
    if damage == 'duplicate_replica':
        payload = json.loads(paths[0].read_text())
    elif damage == 'missing_term':
        payload['terms'].pop()
    elif damage == 'duplicate_term':
        payload['terms'].append(payload['terms'][0])
    else:
        payload['terms'][0]['permutations'] = 19
    paths[1].write_text(json.dumps(payload))
    with pytest.raises(ValueError, match='checkpoint'):
        collect(c, tmp_path)
    with pytest.raises(ValueError, match='checkpoint'):
        run_slice(c, tmp_path, 'exact', 1, 0, 2)


def test_chunks_align_with_serial_batches_and_match_results(tmp_path):
    c = config(replicates=10, batch_size=5)
    expected, summary = run_validation(c, tmp_path / 'serial')
    out = tmp_path / 'slurm'
    prepare(c, out)
    for index in range(4):
        scenario, effect, start, end = task_parameters(c, index, 4)
        assert start % 5 == 0
        assert end % 5 == 0 or end == c['replicates']
        run_slice(c, out, scenario, effect, start, end)
    actual, other_summary = collect(c, out)
    pd.testing.assert_frame_equal(expected, actual)
    pd.testing.assert_frame_equal(summary, other_summary)
    before = {p.name: p.stat().st_mtime_ns for p in (out/'shards/exact').glob('batch_*.json')}
    run_slice(c, out, 'exact', 1, 0, 10)
    assert before == {p.name: p.stat().st_mtime_ns for p in (out/'shards/exact').glob('batch_*.json')}


def test_unaligned_manual_slice_rejected(tmp_path):
    c = config(replicates=10, batch_size=5)
    prepare(c, tmp_path)
    with pytest.raises(ValueError, match='align|batch'):
        run_slice(c, tmp_path, 'exact', 1, 0, 3)


@pytest.mark.parametrize('replicates,batch_size,chunks', [(11, 3, 4), (3, 5, 8), (1000, 5, 4)])
def test_partition_covers_each_global_batch_once(replicates, batch_size, chunks):
    c = config(replicates=replicates, batch_size=batch_size)
    starts = []
    for index in range(chunks):
        _, _, start, end = task_parameters(c, index, chunks)
        starts.extend(range(start, end, batch_size))
    assert starts == list(range(0, replicates, batch_size))


def test_runner_fingerprint_cannot_be_silently_replaced(tmp_path):
    c = config()
    prepare(c, tmp_path)
    manifest = tmp_path / 'runner_metadata.json'
    manifest.write_text('{}')
    with pytest.raises(ValueError, match='runner'):
        run_slice(c, tmp_path, 'exact', 1, 0, 2)
    with pytest.raises(ValueError, match='runner'):
        prepare(c, tmp_path)


def test_shard_manifest_is_checked_by_task(tmp_path):
    c = config()
    prepare(c, tmp_path)
    (tmp_path / 'shards/exact/config.json').write_text('{}')
    with pytest.raises(ValueError, match='manifest'):
        run_slice(c, tmp_path, 'exact', 1, 0, 2)
