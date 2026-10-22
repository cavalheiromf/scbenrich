"""Reproducible pilot benchmark for family-level enrichment statistics.

Run: python benchmarks/simulate.py --replicates 100 --permutations 999 --out-dir results/benchmark
These simulations assess the chosen data-generating models, not universal FDR control.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from scbenrich.core.stats import exact_conditional_or, fdr_bh, fisher_exact_test
from scbenrich.core.stratified import StratifiedPermuter


def _expit(x):
    return 1 / (1 + np.exp(-x))


def term_truth(scenario, terms):
    """Population truth for marginal and continuous-size-conditional targets.

    All effects are positive in this benchmark. In confounded scenarios,
    increasing annotation and selection probabilities imply marginal association
    even for terms with zero direct conditional effect.
    """
    if scenario not in ("simple_null", "simple_positive", "confounded_null", "confounded_positive"):
        raise ValueError("unknown scenario")
    confounded = scenario.startswith("confounded")
    return pd.DataFrame([{"Term_ID": t,
                          "marginal_null": not confounded and t != "effect",
                          "conditional_null": t != "effect"}
                         for t in terms])


def simulate_one(rng, scenario, n_families=120, n_null_terms=8, log_or=1.5):
    """Return covariates, foreground IDs, binary terms, and the true term name."""
    size = rng.uniform(0, 1, n_families)
    ids = [f"F{i:04d}" for i in range(n_families)]
    if scenario in ("simple_null", "simple_positive"):
        foreground = rng.random(n_families) < 0.3
        null_probability = np.full(n_families, 0.2)
    elif scenario in ("confounded_null", "confounded_positive"):
        foreground = rng.random(n_families) < _expit(-2 + 3 * size)
        null_probability = _expit(-2.5 + 3 * size)
    else:
        raise ValueError("unknown scenario")
    if foreground.sum() == 0 or foreground.sum() == n_families:
        raise ValueError("degenerate foreground draw")
    null = rng.random((n_families, n_null_terms)) < null_probability[:, None]
    data = {f"null_{i}": null[:, i].astype(int) for i in range(n_null_terms)}
    true_term = None
    if scenario in ("simple_positive", "confounded_positive"):
        baseline = -1.4 if scenario == "simple_positive" else -2.5 + 3 * size
        probability = _expit(baseline + log_or * foreground)
        data["effect"] = (rng.random(n_families) < probability).astype(int)
        true_term = "effect"
    return (pd.DataFrame({"size": size}, index=ids), set(np.array(ids)[foreground]),
            pd.DataFrame(data, index=ids), true_term)


def evaluate_one(covariates, foreground_ids, terms, bins, permutations, seed):
    """Use the project's Fisher, BH and stratified permutation implementations."""
    foreground = terms.index.isin(foreground_ids)
    fisher = []
    intervals = {}
    for term in terms:
        present = terms[term].to_numpy(dtype=bool)
        a = int(np.sum(present & foreground))
        b = int(np.sum(~present & foreground))
        c = int(np.sum(present & ~foreground))
        d = int(np.sum(~present & ~foreground))
        fisher.append(fisher_exact_test(a, b, c, d, alternative="greater")[1])
        intervals[term] = exact_conditional_or(a, b, c, d)[1:]
    permuter = StratifiedPermuter(covariates, ["size"], seed=seed, bins=bins)
    p_strat = permuter.run_permutations(foreground_ids, terms, permutations)
    names = list(terms)
    return pd.DataFrame({"Term_ID": names, "p_fisher": fisher,
                         "q_fisher": fdr_bh(fisher),
                         "p_strat": [p_strat[t] for t in names],
                         "q_strat": fdr_bh([p_strat[t] for t in names]),
                         "ci_low": [intervals[t][0] for t in names],
                         "ci_high": [intervals[t][1] for t in names]})


def run_benchmark(replicates=100, permutations=999, bins=(2, 4, 8), seed=42,
                  n_families=120, n_null_terms=8, log_or=1.5):
    if log_or <= 0:
        raise ValueError("log_or must be positive for the known positive-effect scenarios")
    if replicates < 1 or permutations < 1 or n_families < 4 or n_null_terms < 1:
        raise ValueError("replicates, permutations, n_families and n_null_terms must be positive")
    if any(not isinstance(b, int) or b < 2 for b in bins):
        raise ValueError("all bin counts must be integers >= 2")
    rng = np.random.default_rng(seed)
    records = []
    for scenario in ("simple_null", "confounded_null", "simple_positive", "confounded_positive"):
        for replicate in range(replicates):
            cov, foreground, terms, true_term = simulate_one(
                rng, scenario, n_families, n_null_terms, log_or)
            for bin_count in bins:
                results = evaluate_one(cov, foreground, terms, bin_count, permutations,
                                       int(rng.integers(0, 2**32 - 1)))
                for method, column in (("fisher", "q_fisher"), ("stratified", "q_strat")):
                    discoveries = set(results.loc[results[column] < 0.05, "Term_ID"])
                    truth = term_truth(scenario, results.Term_ID)
                    target = "marginal_null" if method == "fisher" else "conditional_null"
                    null_terms = set(truth.loc[truth[target], "Term_ID"])
                    conditional_nulls = set(truth.loc[truth.conditional_null, "Term_ID"])
                    false_count = len(discoveries & null_terms)
                    mismatch = len(discoveries & conditional_nulls) if method == "fisher" and scenario.startswith("confounded") else 0
                    records.append({"scenario": scenario, "replicate": replicate,
                                    "bins": bin_count, "method": method,
                                    "discoveries": len(discoveries), "false_discoveries": false_count,
                                    "target": "marginal" if method == "fisher" else "conditional_on_size",
                                    "n_null_terms": len(null_terms),
                                    "conditional_target_mismatch": mismatch,
                                    "any_false_discovery": int(false_count > 0),
                                    "false_discovery_proportion": false_count / max(1, len(discoveries)),
                                    "power": int(true_term in discoveries) if true_term else np.nan,
                                    "ci_coverage": (int(results.loc[results.Term_ID == true_term, "ci_low"].iloc[0]
                                                        <= np.exp(log_or) <=
                                                        results.loc[results.Term_ID == true_term, "ci_high"].iloc[0])
                                                    if scenario == "simple_positive" and method == "fisher" else np.nan)})
    raw = pd.DataFrame(records)
    grouped = raw.groupby(["scenario", "bins", "method"], sort=True)
    summary = grouped.agg(replicates=("replicate", "size"),
                          probability_any_false_discovery=("any_false_discovery", "mean"),
                          empirical_fdr=("false_discovery_proportion", "mean"),
                          mean_conditional_target_mismatch=("conditional_target_mismatch", "mean"),
                          n_null_terms=("n_null_terms", "first"),
                          power=("power", "mean"), ci_coverage=("ci_coverage", "mean")).reset_index()
    summary["any_false_discovery_ci_low"] = [binomtest(round(row.probability_any_false_discovery * row.replicates),
                                               int(row.replicates)).proportion_ci().low
                                   for row in summary.itertuples()]
    summary["any_false_discovery_ci_high"] = [binomtest(round(row.probability_any_false_discovery * row.replicates),
                                                int(row.replicates)).proportion_ci().high
                                    for row in summary.itertuples()]
    return raw, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--bins", type=int, nargs="+", default=[2, 4, 8])
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out-dir", type=Path, default=Path("benchmark_results"))
    args = parser.parse_args()
    raw, summary = run_benchmark(args.replicates, args.permutations, tuple(args.bins), args.seed)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    raw.to_csv(args.out_dir / "replicates.tsv", sep="\t", index=False)
    summary.to_csv(args.out_dir / "summary.tsv", sep="\t", index=False)
    truth_tables = []
    for scenario in ("simple_null", "confounded_null", "simple_positive", "confounded_positive"):
        names = [f"null_{i}" for i in range(8)] + (["effect"] if "positive" in scenario else [])
        truth_tables.append(term_truth(scenario, names).assign(scenario=scenario))
    pd.concat(truth_tables).to_csv(args.out_dir / "term_truth.tsv", sep="\t", index=False)
    (args.out_dir / "config.json").write_text(json.dumps({"replicates": args.replicates,
        "permutations": args.permutations, "bins": args.bins, "seed": args.seed,
        "metric_schema": 2, "conditional_target": "continuous_size", "n_families": 120, "n_null_terms": 8, "log_or": 1.5}, indent=2) + "\n")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
