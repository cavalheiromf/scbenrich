from argparse import Namespace
from itertools import combinations

import numpy as np
import pandas as pd
import pytest

from scbenrich.cli import run_pipeline
from scbenrich.core.stratified import StratifiedPermuter, compute_smd


def test_small_universe_matches_exact_randomization():
    ids = list("ABCD")
    cov = pd.DataFrame({"size": [1, 1, 1, 1]}, index=ids)
    terms = pd.DataFrame({"T": [1, 1, 0, 0]}, index=ids)
    permuter = StratifiedPermuter(cov, ["size"], seed=1)
    assert permuter.support_diagnostics({"A", "C"})["combinations"].iloc[0] == 6
    counts = [sum(terms.loc[list(pair), "T"]) for pair in combinations(ids, 2)]
    assert sum(c >= 1 for c in counts) / len(counts) == pytest.approx(5 / 6)
    assert permuter.run_permutations({"A", "C"}, terms, 5000)["T"] == pytest.approx(5 / 6, abs=0.025)


def test_smd_undefined_is_not_balance():
    assert np.isnan(compute_smd([], [1, 2]))
    assert np.isnan(compute_smd([1, 1], [2, 2]))
    assert compute_smd([1, 1], [1, 1]) == 0


def test_pipeline_rejects_invalid_inputs_and_exports_support(tmp_path):
    master = tmp_path / "master.tsv"
    ann = tmp_path / "ann.tsv"
    pd.DataFrame({"Family_ID": list("ABCD"), "size": [1, 1, 1, 1], "selected": [1, 0, 1, 0]}).to_csv(master, sep="\t", index=False)
    pd.DataFrame({"Family_ID": list("AB"), "Term_ID": ["T", "T"]}).to_csv(ann, sep="\t", index=False)
    args = Namespace(master_table=str(master), annotations=str(ann), foreground="selected", covariates="size", obo=None,
                     permutations=100, min_universe=1, max_universe_fraction=0.9, seed=1, out_dir=str(tmp_path / "out"))
    run_pipeline(args)
    result = pd.read_csv(tmp_path / "out/enrichment_results_all.tsv", sep="\t")
    assert "p_strat_adjusted" in result
    assert result.loc[0, "p_strat_adjusted"] >= result.loc[0, "p_strat"]
    assert result.loc[0, "balance_status"] == "undefined"
    support = pd.read_csv(tmp_path / "out/strata_diagnostics.tsv", sep="\t")
    assert support.loc[0, "combinations"] == 6
    args.permutations = -1
    with pytest.raises(ValueError, match="permutations"):
        run_pipeline(args)


@pytest.mark.parametrize("bad_size, expected", [([1, np.nan, 2, 3], "finite numeric"), ([1, 1, 2, 2], "no randomization")])
def test_invalid_covariate_or_no_randomization(tmp_path, bad_size, expected):
    master = tmp_path / "master.tsv"
    ann = tmp_path / "ann.tsv"
    selected = [1, 0, 1, 0] if expected == "finite numeric" else [1, 1, 0, 0]
    pd.DataFrame({"Family_ID": list("ABCD"), "size": bad_size, "selected": selected}).to_csv(master, sep="\t", index=False)
    pd.DataFrame({"Family_ID": ["A", "B"], "Term_ID": ["T", "T"]}).to_csv(ann, sep="\t", index=False)
    args = Namespace(master_table=str(master), annotations=str(ann), foreground="selected", covariates="size", obo=None,
                     permutations=20, min_universe=1, max_universe_fraction=0.9, seed=1, out_dir=str(tmp_path / "out"))
    with pytest.raises(ValueError, match=expected):
        run_pipeline(args)


def test_empty_term_run_exports_headers(tmp_path):
    master = tmp_path / 'master.tsv'
    ann = tmp_path / 'annotations.tsv'
    pd.DataFrame({'Family_ID': list('ABCD'), 'size': [1, 1, 1, 1], 'selected': [1, 0, 1, 0]}).to_csv(master, sep='\t', index=False)
    pd.DataFrame({'Family_ID': ['A'], 'Term_ID': ['T']}).to_csv(ann, sep='\t', index=False)
    args = Namespace(master_table=str(master), annotations=str(ann), foreground='selected', covariates='size', obo=None,
                     permutations=10, min_universe=5, max_universe_fraction=.9, seed=1, out_dir=str(tmp_path / 'out'))
    run_pipeline(args)
    result = pd.read_csv(tmp_path / 'out/enrichment_results_all.tsv', sep='\t')
    assert result.empty
    assert 'analysis_status' in result.columns
