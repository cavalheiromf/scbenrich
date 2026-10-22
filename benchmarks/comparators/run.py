"""Audited family-level baseline comparisons on a shared binary annotation matrix."""
import argparse
import hashlib
import io
import json
from importlib.metadata import version
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pandas as pd

from scbenrich.core.stats import fisher_exact_test, fdr_bh
from scbenrich.core.obo import parse_obo, normalize_go_id, propagate_go_annotations
from scbenrich.utils.io import load_master_table, load_annotations


def export_inputs(matrix, foreground, out, min_count=1, max_fraction=.9):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    ids = list(matrix.index)
    fg = matrix.index.isin(foreground)
    if not 0 < fg.sum() < len(ids):
        raise ValueError('foreground must be a nonempty proper subset')
    retained = matrix.loc[:, (matrix.sum() >= min_count) & (matrix.sum() <= int(len(ids) * max_fraction))]
    if retained.empty:
        raise ValueError('no eligible terms')
    (out / 'universe.txt').write_text('\n'.join(ids) + '\n')
    (out / 'foreground.txt').write_text('\n'.join(matrix.index[fg]) + '\n')
    pairs = [{'Family_ID': fam, 'Term_ID': term} for term in retained for fam in retained.index[retained[term] == 1]]
    pd.DataFrame(pairs).to_csv(out / 'term2family.tsv', sep='\t', index=False)
    retained.to_csv(out / 'matrix.tsv', sep='\t', index_label='Family_ID')
    records = []
    for term in retained:
        a = int(retained.loc[fg, term].sum())
        c = int(retained.loc[~fg, term].sum())
        cells = (a, int(fg.sum()) - a, c, len(ids) - int(fg.sum()) - c)
        records.append(dict(Term_ID=term, foreground_count=a, universe_count=a+c,
                            foreground_total=int(fg.sum()), universe_total=len(ids),
                            p_greater=fisher_exact_test(*cells, alternative='greater')[1],
                            p_two_sided=fisher_exact_test(*cells, alternative='two-sided')[1]))
    reference = pd.DataFrame(records)
    reference['q_greater'] = fdr_bh(reference.p_greater.to_numpy())
    reference['q_two_sided'] = fdr_bh(reference.p_two_sided.to_numpy())
    reference.to_csv(out / 'reference.tsv', sep='\t', index=False)
    return reference


def compare_tables(reference, external, alternative):
    result = reference.merge(external, on='Term_ID', how='outer', suffixes=('_reference', '_external'), indicator=True)
    match = result['_merge'].eq('both')
    for col in ['foreground_count', 'universe_count', 'foreground_total', 'universe_total']:
        if col + '_external' in result:
            match &= result[col + '_reference'].eq(result[col + '_external'])
    result['counts_match'] = match
    pcol = 'p_greater' if alternative == 'greater' else 'p_two_sided'
    result['p_match'] = match & np.isclose(result[pcol], result.p_value, rtol=1e-9, atol=1e-12, equal_nan=False)
    result['q_harmonized'] = fdr_bh(result.p_value.to_numpy()) if result.p_value.notna().all() else np.nan
    qcol = 'q_greater' if alternative == 'greater' else 'q_two_sided'
    result['q_match'] = match & np.isclose(result[qcol], result.q_harmonized, rtol=1e-9, atol=1e-12, equal_nan=False)
    return result.drop(columns='_merge')


def run_goatools(out, obo):
    from goatools.obo_parser import GODag
    from goatools.go_enrichment import GOEnrichmentStudy
    population = (out / 'universe.txt').read_text().splitlines()
    foreground = (out / 'foreground.txt').read_text().splitlines()
    pairs = pd.read_csv(out / 'term2family.tsv', sep='\t')
    association = pairs.groupby('Family_ID').Term_ID.apply(set).to_dict()
    log = io.StringIO()
    dag = GODag(str(obo), prt=log)
    study = GOEnrichmentStudy(population, association, dag, propagate_counts=False,
                              methods=['fdr_bh'], log=log)
    results = study.run_study(foreground, prt=log)
    frame = pd.DataFrame([dict(Term_ID=r.GO, foreground_count=r.study_count,
                               universe_count=r.pop_count, foreground_total=r.study_n,
                               universe_total=r.pop_n, p_value=r.p_uncorrected) for r in results])
    (out / 'goatools.log').write_text(log.getvalue())
    frame.to_csv(out / 'goatools.tsv', sep='\t', index=False)
    return frame


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--master-table', required=True)
    parser.add_argument('--annotations', required=True)
    parser.add_argument('--foreground', required=True, help='Master-table column')
    parser.add_argument('--obo', required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--min-count', type=int, default=2)
    parser.add_argument('--max-fraction', type=float, default=.95)
    args = parser.parse_args()
    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise ValueError('output directory must be empty')
    master = load_master_table(args.master_table)
    fg = set(master.index[master[args.foreground].astype(str).str.lower().isin(['1', 'true', 'yes', 'expanded'])])
    ann = load_annotations(args.annotations)
    terms, aliases, replacements, ancestors, go_version = parse_obo(args.obo)
    ann['Term_ID'] = ann.Term_ID.map(lambda t: normalize_go_id(t, terms, aliases, replacements))
    mapping = propagate_go_annotations(ann.groupby('Family_ID').Term_ID.apply(set).to_dict(), ancestors)
    matrix = pd.DataFrame({t: [int(t in mapping.get(f, set())) for f in master.index]
                           for t in sorted(set().union(*mapping.values()))}, index=master.index)
    reference = export_inputs(matrix, fg, args.out_dir, args.min_count, args.max_fraction)
    manifest = {'unit': 'family', 'propagation': 'precomputed is_a/part_of; no GOATOOLS propagation',
                'go_version': go_version, 'bh_scope': 'common retained terms across categories',
                'foreground_column': args.foreground, 'min_count': args.min_count,
                'max_fraction': args.max_fraction,
                'python_versions': {name: version(name) for name in ['scbenrich', 'numpy', 'pandas', 'scipy']},
                'adapter_hashes': {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                   for name in ['run.py', 'run_r.R']},
                'tools': {}, 'inputs': {name: hashlib.sha256(Path(path).read_bytes()).hexdigest()
                 for name, path in [('master', args.master_table), ('annotations', args.annotations), ('obo', args.obo)]}}
    try:
        frame = run_goatools(args.out_dir, args.obo)
        audit = compare_tables(reference, frame, 'two-sided')
        audit.to_csv(args.out_dir / 'audit_goatools.tsv', sep='\t', index=False)
        manifest['tools']['goatools'] = {'status': 'completed', 'version': version('goatools'),
                    'alternative': 'two-sided', 'all_counts_match': bool(audit.counts_match.all()),
                    'all_p_match': bool(audit.p_match.all())}
    except ImportError as exc:
        manifest['tools']['goatools'] = {'status': 'unavailable', 'reason': str(exc)}
    if shutil.which('Rscript'):
        script = Path(__file__).with_name('run_r.R')
        run = subprocess.run(['Rscript', str(script), str(args.out_dir.resolve())], capture_output=True, text=True, timeout=240)
        (args.out_dir / 'R.log').write_text(run.stdout + run.stderr)
        if run.returncode:
            manifest['tools']['R'] = {'status': 'failed', 'returncode': run.returncode}
        r_meta = args.out_dir / 'r_versions.tsv'
        if r_meta.exists():
            manifest['r_versions'] = pd.read_csv(r_meta, sep='\t').to_dict('records')
        for tool in ['topgo_classic', 'topgo_elim', 'clusterprofiler']:
            path = args.out_dir / (tool + '.tsv')
            if path.exists():
                external = pd.read_csv(path, sep='\t')
                audit = compare_tables(reference, external, 'greater')
                audit.to_csv(args.out_dir / ('audit_' + tool + '.tsv'), sep='\t', index=False)
                manifest['tools'][tool] = {'status': 'completed', 'alternative': 'greater',
                    'all_counts_match': bool(audit.counts_match.all()), 'all_p_match': bool(audit.p_match.all()),
                    'p_equality_expected': tool != 'topgo_elim'}
            else:
                manifest['tools'][tool] = {'status': 'not_completed', 'reason': 'See R.log'}
    else:
        manifest['tools']['R'] = {'status': 'unavailable'}
    manifest['export_hashes'] = {name: hashlib.sha256((args.out_dir / name).read_bytes()).hexdigest()
                               for name in ['universe.txt', 'foreground.txt', 'matrix.tsv', 'term2family.tsv']}
    (args.out_dir / 'comparison_metadata.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest['tools'], indent=2))


if __name__ == '__main__':
    main()
