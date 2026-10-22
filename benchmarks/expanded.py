"""Expanded simulation validation with stable random streams and resumable batches."""
import argparse
import hashlib
import json
from pathlib import Path
from importlib.metadata import version

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import binomtest

from scbenrich.core.stats import fdr_bh, fisher_exact_test, exact_conditional_or
from scbenrich.core.stratified import StratifiedPermuter

SCENARIOS = {'exact': 1, 'continuous': 2, 'missing_size': 3, 'dependent': 4,
             'sparse': 5, 'derived_degree': 6, 'unconfounded': 7}


def generate(rng, scenario, n_families, odds_ratio):
    """Generate observed annotations and population truth; unknown truth stays missing."""
    if scenario not in SCENARIOS:
        raise ValueError('unknown scenario')
    if n_families % 2:
        raise ValueError('n_families must be even for balanced discrete covariates')
    x = np.column_stack([rng.permutation(np.tile([0., 1.], n_families // 2)) for _ in range(3)])
    if scenario in ('continuous', 'missing_size', 'sparse', 'derived_degree'):
        x = rng.uniform(0, 1, (n_families, 3))
    score = x.sum(axis=1)
    propensity = expit(-2 + score)
    if scenario == 'sparse':
        propensity = expit(10 * (score - 1.5))
    if scenario == 'unconfounded':
        propensity = np.full(n_families, .3)
    selected = rng.random(n_families) < propensity
    if not 0 < selected.sum() < n_families:
        raise ValueError('degenerate foreground draw')
    baseline = -2.5 + score
    if scenario == 'unconfounded':
        baseline = np.full(n_families, -1.4)
    terms = {f'null_{i}': (rng.random(n_families) < expit(baseline)).astype(int) for i in range(6)}
    terms['effect'] = (rng.random(n_families) < expit(baseline + np.log(odds_ratio) * selected)).astype(int)
    if scenario == 'dependent':
        terms['child_0'] = terms.pop('null_0')
        terms['child_1'] = terms.pop('null_1')
        terms['parent'] = np.maximum(terms['child_0'], terms['child_1'])
    if scenario == 'missing_size':
        for term in terms:
            terms[term] *= (rng.random(n_families) < expit(2 - 3 * x[:, 0])).astype(int)
    ids = [f'F{i:05d}' for i in range(n_families)]
    matrix = pd.DataFrame(terms, index=ids)
    cov = pd.DataFrame(x, index=ids, columns=['size', 'reference_homologs', 'external_annotation_degree'])
    if scenario == 'derived_degree':
        cov['external_annotation_degree'] = matrix.sum(axis=1)
        cov = cov.rename(columns={'external_annotation_degree': 'derived_annotation_degree'})
    truth = pd.DataFrame({'Term_ID': matrix.columns,
                          'conditional_null': [not (t == 'effect' and odds_ratio > 1) for t in matrix],
                          'marginal_null': [scenario == 'unconfounded' and not (t == 'effect' and odds_ratio > 1) for t in matrix]}).set_index('Term_ID')
    if scenario == 'missing_size':
        # Size-dependent thinning can reverse or cancel marginal associations.
        truth['marginal_null'] = None
    if scenario == 'derived_degree':
        # Conditioning on an outcome-derived degree changes the null; do not assert it.
        truth['conditional_null'] = None
    return {'covariates': cov, 'foreground': set(np.array(ids)[selected]), 'terms': matrix, 'truth': truth}


def _seed(seed, scenario, odds_ratio, replicate, purpose, bins=0):
    digest = hashlib.sha256(str(odds_ratio).encode()).digest()
    key = int.from_bytes(digest[:4], 'little')
    return np.random.SeedSequence([seed, SCENARIOS[scenario], key, replicate, purpose, bins])


def evaluate(data, bins, permutations, seed):
    cov, matrix, fg = data['covariates'], data['terms'], data['foreground']
    mask = matrix.index.isin(fg)
    fisher, intervals = [], []
    for term in matrix:
        values = matrix[term].to_numpy(dtype=bool)
        cells = [int((values & mask).sum()), int((~values & mask).sum()),
                 int((values & ~mask).sum()), int((~values & ~mask).sum())]
        fisher.append(fisher_exact_test(*cells)[1])
        intervals.append(exact_conditional_or(*cells)[1:])
    permuter = StratifiedPermuter(cov, list(cov), seed=seed, bins=bins)
    result = permuter.run_analysis(fg, matrix, permutations)
    table = result['terms'].set_index('Term_ID')
    table['p_fisher'] = fisher
    table['q_fisher'] = fdr_bh(fisher)
    table['q_stratified'] = fdr_bh(table.p_strat.to_numpy())
    table['ci_low'] = [ci[0] for ci in intervals]
    table['ci_high'] = [ci[1] for ci in intervals]
    table = table.join(data['truth'])
    fraction = float(result['strata'].informative.mean())
    return table, fraction


def summarize(raw, seed=42):
    """Uncertainty is calculated across independent replicates, not pooled terms."""
    rng = np.random.default_rng(seed)
    rows = []
    for key, group in raw.groupby(['scenario', 'odds_ratio', 'bins', 'method'], sort=True):
        row = dict(zip(['scenario', 'odds_ratio', 'bins', 'method'], key))
        row['replicates'] = len(group)
        row['informative_strata_fraction'] = group.informative_strata_fraction.mean()
        for metric, source in [('fdr', 'fdp'), ('power', 'effect_detected'), ('coverage', 'coverage')]:
            values = group[source].dropna().to_numpy(dtype=float)
            row[metric] = float(values.mean()) if len(values) else np.nan
            row[metric + '_n'] = len(values)
            low = high = np.nan
            if len(values) >= 2:
                if metric == 'fdr':
                    means = rng.choice(values, (2000, len(values)), replace=True).mean(axis=1)
                    low, high = np.quantile(means, [.025, .975])
                else:
                    ci = binomtest(int(values.sum()), len(values)).proportion_ci()
                    low, high = ci.low, ci.high
            row[metric + '_ci_low'], row[metric + '_ci_high'] = low, high
        row['no_randomization_rate'] = group.no_randomization.mean()
        rows.append(row)
    return pd.DataFrame(rows)


def run_validation(config, out_dir):
    """Resume only identical configurations and source fingerprints."""
    required = ['seed', 'replicates', 'permutations', 'n_families', 'bins', 'odds_ratios', 'scenarios', 'batch_size']
    if any(k not in config for k in required):
        raise ValueError('missing configuration field')
    for name in ['seed', 'replicates', 'permutations', 'n_families', 'batch_size']:
        if isinstance(config[name], bool) or not isinstance(config[name], int) or config[name] < (0 if name == 'seed' else 1):
            raise ValueError('invalid configuration: ' + name)
    if config['n_families'] < 4 or not config['bins'] or not config['odds_ratios'] or not config['scenarios']:
        raise ValueError('configuration requires families, bins, effects and scenarios')
    if any(not isinstance(b, int) or b < 2 for b in config['bins']):
        raise ValueError('invalid bins')
    if any(not np.isfinite(o) or o < 1 for o in config['odds_ratios']):
        raise ValueError('odds ratios must be >= 1')
    if set(config['scenarios']) - set(SCENARIOS):
        raise ValueError('unknown scenario')
    for name in ['bins', 'odds_ratios', 'scenarios']:
        if len(set(config[name])) != len(config[name]):
            raise ValueError('duplicate configuration values')
    root = Path(__file__).resolve().parents[1]
    sources = [Path(__file__)] + list((root / 'scbenrich').rglob('*.py'))
    identity = {'config': config, 'source_hashes': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources)},
                'versions': {name: version(name) for name in ('numpy', 'pandas', 'scipy')}}
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    manifest = out / 'config.json'
    if manifest.exists():
        if json.loads(manifest.read_text()) != identity:
            raise ValueError('configuration or code differs from existing run')
    elif any(out.iterdir()):
        raise ValueError('output directory contains unidentified results')
    else:
        manifest.write_text(json.dumps(identity, indent=2, sort_keys=True) + '\n')
    all_records, all_terms = [], []
    for scenario in config['scenarios']:
        for effect in config['odds_ratios']:
            for start in range(0, config['replicates'], config['batch_size']):
                tag = hashlib.sha256(str(effect).encode()).hexdigest()[:12]
                path = out / f'batch_{scenario}_{tag}_{start:06d}.json'
                if path.exists():
                    payload = json.loads(path.read_text())
                else:
                    records, term_records = [], []
                    for rep in range(start, min(start + config['batch_size'], config['replicates'])):
                        data = generate(np.random.default_rng(_seed(config['seed'], scenario, effect, rep, 0)), scenario, config['n_families'], effect)
                        for bins in config['bins']:
                            table, fraction = evaluate(data, bins, config['permutations'], _seed(config['seed'], scenario, effect, rep, 1, bins))
                            common = dict(scenario=scenario, odds_ratio=effect, replicate=rep, bins=bins)
                            for term, row in table.iterrows():
                                term_records.append(dict(common, Term_ID=term, **{k: (v.item() if isinstance(v, np.generic) else v) for k, v in row.items()}))
                            for method, q in [('fisher', 'q_fisher'), ('stratified', 'q_stratified')]:
                                truth = table['marginal_null' if method == 'fisher' else 'conditional_null']
                                discovered = table[q] < .05
                                if method == 'stratified':
                                    discovered &= table.analysis_status.eq('ok')
                                known = truth.notna().all()
                                false = int((discovered & truth.fillna(False).astype(bool)).sum()) if known else None
                                intervals = table.loc['effect']
                                records.append(dict(common, method=method, informative_strata_fraction=fraction,
                                    no_randomization=bool(table.analysis_status.eq('no_randomization').all()),
                                    discoveries=int(discovered.sum()), n_null_terms=int(truth.fillna(False).sum()) if known else None,
                                    fdp=false / max(1, int(discovered.sum())) if known else None,
                                    effect_detected=int(discovered['effect']) if effect > 1 and known else None,
                                    coverage=int(intervals.ci_low <= effect <= intervals.ci_high) if scenario == 'unconfounded' and method == 'fisher' else None))
                    # JSON is also the restart checkpoint; NaNs in term data become null.
                    payload = json.loads(pd.Series({'replicates': records, 'terms': term_records}).to_json())
                    temporary = path.with_suffix('.tmp')
                    temporary.write_text(json.dumps(payload, allow_nan=False))
                    temporary.replace(path)
                all_records.extend(payload['replicates'])
                all_terms.extend(payload['terms'])
    raw = pd.DataFrame(all_records).sort_values(['scenario', 'odds_ratio', 'replicate', 'bins', 'method']).reset_index(drop=True)
    summary = summarize(raw, config['seed'])
    raw.to_csv(out / 'replicates.tsv', sep='\t', index=False)
    pd.DataFrame(all_terms).to_csv(out / 'terms.tsv', sep='\t', index=False)
    summary.to_csv(out / 'summary.tsv', sep='\t', index=False)
    return raw, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    _, summary = run_validation(json.loads(args.config.read_text()), args.out_dir)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
