"""Prepare, execute one Slurm array task, or collect validation scenarios."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

from benchmarks.parallel import prepare, collect, validate_prepared
from benchmarks.checkpoints import read_checkpoint, validate_checkpoint
from benchmarks.expanded import run_validation, generate, evaluate, _seed


def task_parameters(config, task_index, chunks_per_effect=4):
    if isinstance(chunks_per_effect, bool) or not isinstance(chunks_per_effect, int) or chunks_per_effect < 1:
        raise ValueError('chunks_per_effect must be a positive integer')
    scenarios, effects = config['scenarios'], config['odds_ratios']
    total = len(scenarios) * len(effects) * chunks_per_effect
    if isinstance(task_index, bool) or not isinstance(task_index, int) or not 0 <= task_index < total:
        raise IndexError('task index out of range')
    chunk = task_index % chunks_per_effect
    effect = effects[(task_index // chunks_per_effect) % len(effects)]
    scenario = scenarios[task_index // (chunks_per_effect * len(effects))]
    # Partition the global batch indices, never individual replicas. Empty chunks
    # are valid when more tasks than batches are requested.
    batch_size = config['batch_size']
    batches = (config['replicates'] + batch_size - 1) // batch_size
    start_rep = min(config['replicates'], (chunk * batches // chunks_per_effect) * batch_size)
    end_rep = min(config['replicates'], ((chunk + 1) * batches // chunks_per_effect) * batch_size)

    return scenario, effect, start_rep, end_rep


def run_slice(config, out_dir, scenario, effect, start_rep, end_rep):
    validate_prepared(config, out_dir, scenario)
    if scenario not in config['scenarios'] or effect not in config['odds_ratios']:
        raise ValueError('scenario or effect outside prepared configuration')
    batch_size = config['batch_size']
    if (not 0 <= start_rep <= end_rep <= config['replicates'] or
            (start_rep != config['replicates'] and start_rep % batch_size) or
            (end_rep != config['replicates'] and end_rep % batch_size)):
        raise ValueError('slice boundaries must align with complete batches')
    directory = Path(out_dir) / 'shards' / scenario
    directory.mkdir(parents=True, exist_ok=True)
    tag = hashlib.sha256(str(effect).encode()).hexdigest()[:12]

    for start in range(start_rep, end_rep, batch_size):
        path = directory / f'batch_{scenario}_{tag}_{start:06d}.json'
        if path.exists():
            read_checkpoint(path, config, scenario, effect, start)
            continue
        records, term_records = [], []
        rep_end = min(start + batch_size, end_rep)
        for rep in range(start, rep_end):
            data = generate(np.random.default_rng(_seed(config['seed'], scenario, effect, rep, 0)),
                            scenario, config['n_families'], effect)
            for bins in config['bins']:
                table, fraction = evaluate(data, bins, config['permutations'],
                                           _seed(config['seed'], scenario, effect, rep, 1, bins))
                common = dict(scenario=scenario, odds_ratio=effect, replicate=rep, bins=bins)
                for term, row in table.iterrows():
                    term_records.append(dict(common, Term_ID=term,
                                             **{k: (v.item() if isinstance(v, np.generic) else v) for k, v in row.items()}))
                for method, q in [('fisher', 'q_fisher'), ('stratified', 'q_stratified')]:
                    truth = table['marginal_null' if method == 'fisher' else 'conditional_null']
                    discovered = table[q] < .05
                    if method == 'stratified':
                        discovered &= table.analysis_status.eq('ok')
                    known = truth.notna().all()
                    false = int((discovered & truth.fillna(False).astype(bool)).sum()) if known else None
                    intervals = table.loc['effect']
                    records.append(dict(
                        common, method=method, informative_strata_fraction=fraction,
                        no_randomization=bool(table.analysis_status.eq('no_randomization').all()),
                        discoveries=int(discovered.sum()),
                        n_null_terms=int(truth.fillna(False).sum()) if known else None,
                        fdp=false / max(1, int(discovered.sum())) if known else None,
                        effect_detected=int(discovered['effect']) if effect > 1 and known else None,
                        coverage=int(intervals.ci_low <= effect <= intervals.ci_high) if scenario == 'unconfounded' and method == 'fisher' else None
                    ))
        payload = json.loads(pd.Series({'replicates': records, 'terms': term_records}).to_json())
        validate_checkpoint(payload, config, scenario, effect, start)
        temporary = path.with_suffix('.tmp')
        temporary.write_text(json.dumps(payload, allow_nan=False))
        temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'task', 'collect'])
    parser.add_argument('--config', required=True)
    parser.add_argument('--out-dir', required=True)
    parser.add_argument('--index', type=int)
    parser.add_argument('--chunks-per-effect', type=int, default=4)
    args = parser.parse_args()
    with open(args.config) as handle:
        config = json.load(handle)

    if args.action == 'prepare':
        prepare(config, args.out_dir)
    elif args.action == 'task':
        if args.index is None:
            parser.error('task requires --index')
        if args.chunks_per_effect <= 0:
            if not 0 <= args.index < len(config['scenarios']):
                parser.error('task requires a valid --index for scenarios')
            scenario = config['scenarios'][args.index]
            directory = Path(args.out_dir) / 'shards' / scenario
            validate_prepared(config, args.out_dir, scenario)
            subset = dict(config, scenarios=[scenario])
            run_validation(subset, directory)
        else:
            total_tasks = len(config['scenarios']) * len(config['odds_ratios']) * args.chunks_per_effect
            if not 0 <= args.index < total_tasks:
                parser.error(f'task index must be between 0 and {total_tasks - 1}')
            scenario, effect, start_rep, end_rep = task_parameters(config, args.index, args.chunks_per_effect)
            run_slice(config, args.out_dir, scenario, effect, start_rep, end_rep)
    else:
        collect(config, args.out_dir)


if __name__ == '__main__':
    main()
