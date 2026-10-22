"""Prepare, execute one Slurm array task, or collect validation scenarios."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

from benchmarks.parallel import prepare, collect
from benchmarks.expanded import run_validation, generate, evaluate, _seed


def task_parameters(config, task_index, chunks_per_effect=4):
    scenarios = config['scenarios']
    odds_ratios = config['odds_ratios']
    n_chunks = chunks_per_effect

    chunk_idx = task_index % n_chunks
    rem = task_index // n_chunks
    or_idx = rem % len(odds_ratios)
    scen_idx = rem // len(odds_ratios)

    if not (0 <= scen_idx < len(scenarios)):
        raise IndexError(f"Task index {task_index} out of range")

    scenario = scenarios[scen_idx]
    effect = odds_ratios[or_idx]

    total_reps = config['replicates']
    reps_per_chunk = (total_reps + n_chunks - 1) // n_chunks
    start_rep = chunk_idx * reps_per_chunk
    end_rep = min(total_reps, (chunk_idx + 1) * reps_per_chunk)

    return scenario, effect, start_rep, end_rep


def run_slice(config, out_dir, scenario, effect, start_rep, end_rep):
    directory = Path(out_dir) / 'shards' / scenario
    directory.mkdir(parents=True, exist_ok=True)
    batch_size = config['batch_size']
    tag = hashlib.sha256(str(effect).encode()).hexdigest()[:12]

    for start in range(start_rep, end_rep, batch_size):
        path = directory / f'batch_{scenario}_{tag}_{start:06d}.json'
        if path.exists():
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
