"""Checkpoint contracts shared by parallel and Slurm validation runners."""
import hashlib
import itertools
import json
from pathlib import Path


def batch_path(directory, scenario, effect, start):
    tag = hashlib.sha256(str(effect).encode()).hexdigest()[:12]
    return Path(directory) / f'batch_{scenario}_{tag}_{start:06d}.json'


def validate_checkpoint(payload, config, scenario, effect, start):
    """Reject missing/duplicate rows and batches inconsistent with their filename."""
    message = f'invalid checkpoint: {scenario}, effect={effect}, start={start}'
    try:
        reps = range(start, min(start + config['batch_size'], config['replicates']))
        base = list(itertools.product([scenario], [effect], reps, config['bins']))
        expected_metrics = {key + (method,) for key in base for method in ['fisher', 'stratified']}
        terms = {'effect'} | {f'null_{i}' for i in range(6)}
        if scenario == 'dependent':
            terms = (terms - {'null_0', 'null_1'}) | {'child_0', 'child_1', 'parent'}
        expected_terms = {key + (term,) for key in base for term in terms}
        common = ['scenario', 'odds_ratio', 'replicate', 'bins']
        metrics = [tuple(row[k] for k in common + ['method']) for row in payload['replicates']]
        annotations = [tuple(row[k] for k in common + ['Term_ID']) for row in payload['terms']]
        if (len(metrics) != len(expected_metrics) or set(metrics) != expected_metrics or
                len(annotations) != len(expected_terms) or set(annotations) != expected_terms or
                any(row['permutations'] != config['permutations'] for row in payload['terms'])):
            raise ValueError(message)
    except (KeyError, TypeError) as error:
        raise ValueError(message) from error
    return payload


def read_checkpoint(path, config, scenario, effect, start):
    try:
        payload = json.loads(Path(path).read_text())
    except (OSError, ValueError) as error:
        raise ValueError('unreadable checkpoint: ' + str(path)) from error
    return validate_checkpoint(payload, config, scenario, effect, start)
