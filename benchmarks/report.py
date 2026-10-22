"""Calibration tables and optional figures from expanded term-level results."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def calibration(terms):
    rows = []
    for method, p_col, truth_col in [('fisher', 'p_fisher', 'marginal_null'),
                                      ('stratified', 'p_strat', 'conditional_null')]:
        nulls = terms.loc[terms[truth_col].eq(True)]
        for key, group in nulls.groupby(['scenario', 'odds_ratio', 'bins'], sort=True):
            for threshold in [.001, .005, .01, .025, .05, .1, .25, .5, .75, 1.]:
                # Average within each replica before summarizing across replicas.
                rates = group.assign(reject=group[p_col] <= threshold).groupby('replicate').reject.mean()
                rows.append(dict(zip(['scenario', 'odds_ratio', 'bins'], key), method=method,
                                 threshold=threshold, observed_rate=rates.mean(), replicates=len(rates)))
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    terms = pd.read_csv(args.directory / 'terms.tsv', sep='\t')
    result = calibration(terms)
    result.to_csv(args.directory / 'calibration.tsv', sep='\t', index=False)
    if args.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        for scenario, group in result.groupby('scenario'):
            fig, ax = plt.subplots(figsize=(6, 5))
            for key, curve in group.groupby(['odds_ratio', 'bins', 'method']):
                ax.plot(curve.threshold, curve.observed_rate, marker='.', label=str(key))
            ax.plot([0, 1], [0, 1], color='black', linestyle='--')
            ax.set(xlabel='Nominal threshold', ylabel='Mean null rejection fraction', title=scenario)
            ax.legend(fontsize=6)
            fig.tight_layout()
            fig.savefig(args.directory / f'calibration_{scenario}.svg')
            plt.close(fig)


if __name__ == '__main__':
    main()
