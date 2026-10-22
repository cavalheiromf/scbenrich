"""Run independent scenarios in separate processes, preserving serial checkpoints."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
from importlib.metadata import version
import json
import multiprocessing
from pathlib import Path
import shutil

import pandas as pd

from benchmarks import expanded


def _scenario(config, directory):
    expanded.run_validation(config, directory)
    return str(directory)


def prepare(config, out_dir):
    scenarios = config['scenarios']
    if not scenarios or len(set(scenarios)) != len(scenarios) or set(scenarios) - set(expanded.SCENARIOS):
        raise ValueError('invalid scenarios')
    out = Path(out_dir)
    root = Path(expanded.__file__).resolve().parents[1]
    sources = [Path(expanded.__file__)] + list((root / 'scbenrich').rglob('*.py'))
    identity = dict(config=config,
                    source_hashes={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources)},
                    versions={name: version(name) for name in ('numpy', 'pandas', 'scipy')})
    manifest = out / 'config.json'
    if manifest.exists():
        if json.loads(manifest.read_text()) != identity:
            raise ValueError('configuration or code differs from existing run')
    else:
        out.mkdir(parents=True, exist_ok=True)
        if any(out.iterdir()):
            raise ValueError('output directory contains unidentified results')
        manifest.write_text(json.dumps(identity, indent=2, sort_keys=True) + '\n')

    jobs = []
    for scenario in scenarios:
        subset = dict(config, scenarios=[scenario])
        directory = out / 'shards' / scenario
        directory.mkdir(parents=True, exist_ok=True)
        shard_identity = dict(identity, config=subset)
        shard_manifest = directory / 'config.json'
        if shard_manifest.exists():
            if json.loads(shard_manifest.read_text()) != shard_identity:
                raise ValueError('incompatible shard: ' + scenario)
        else:
            shard_manifest.write_text(json.dumps(shard_identity, indent=2, sort_keys=True) + '\n')
        for checkpoint in out.glob(f'batch_{scenario}_*.json'):
            target = directory / checkpoint.name
            if not target.exists():
                shutil.copy2(checkpoint, target)
        jobs.append((subset, directory))

    return jobs


def collect(config, out_dir):
    out = Path(out_dir)
    jobs = prepare(config, out_dir)
    records, terms = [], []
    for _, directory in jobs:
        expected = len(config['odds_ratios']) * ((config['replicates'] + config['batch_size'] - 1) // config['batch_size'])
        if len(list(directory.glob('batch_*.json'))) != expected:
            raise ValueError('incomplete scenario: ' + directory.name)
        for checkpoint in sorted(directory.glob('batch_*.json')):
            payload = json.loads(checkpoint.read_text())
            records.extend(payload['replicates'])
            terms.extend(payload['terms'])
            target = out / checkpoint.name
            if not target.exists():
                temporary = target.with_suffix('.tmp')
                shutil.copy2(checkpoint, temporary)
                temporary.replace(target)
    raw = pd.DataFrame(records).sort_values(['scenario', 'odds_ratio', 'replicate', 'bins', 'method']).reset_index(drop=True)
    summary = expanded.summarize(raw, config['seed'])
    raw.to_csv(out / 'replicates.tsv', sep='\t', index=False)
    pd.DataFrame(terms).sort_values(['scenario', 'odds_ratio', 'replicate', 'bins', 'Term_ID']).to_csv(out / 'terms.tsv', sep='\t', index=False)
    summary.to_csv(out / 'summary.tsv', sep='\t', index=False)
    return raw, summary


def run_parallel(config, out_dir, workers=7):
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 1:
        raise ValueError('workers must be a positive integer')
    jobs = prepare(config, out_dir)
    with ProcessPoolExecutor(max_workers=min(workers, len(jobs)),
                             mp_context=multiprocessing.get_context('spawn')) as executor:
        futures = [executor.submit(_scenario, subset, directory) for subset, directory in jobs]
        for future in as_completed(futures):
            print('Completed scenario:', Path(future.result()).name, flush=True)
    return collect(config, out_dir)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=7)
    args = parser.parse_args()
    run_parallel(json.loads(args.config.read_text()), args.out_dir, args.workers)


if __name__ == '__main__':
    main()
