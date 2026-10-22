"""
Unit tests for SCB-Enrich core functions and pipeline execution.
"""

import os
import shutil
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from scbenrich.core.stats import fdr_bh, fisher_exact_test, exact_conditional_or
from scbenrich.core.stratified import compute_smd, StratifiedPermuter
from scbenrich.core.obo import parse_obo, propagate_go_annotations
from scbenrich.core.dag import reduce_go_redundancy
from scbenrich.cli import run_pipeline


def test_fdr_bh():
    pvals = [0.001, 0.01, 0.05, 0.5, 0.9]
    qvals = fdr_bh(pvals)
    assert len(qvals) == len(pvals)
    assert np.all(qvals >= pvals)
    assert qvals[0] == 0.005


def test_fisher_and_exact_or():
    # 2x2 table: strong enrichment
    a, b, c, d = 15, 5, 2, 18
    or_val, p_val = fisher_exact_test(a, b, c, d)
    assert p_val < 0.001
    assert or_val > 10.0
    
    mle_or, ci_low, ci_high = exact_conditional_or(a, b, c, d)
    assert mle_or > 5.0
    assert ci_low > 1.0
    assert ci_high > mle_or


def test_compute_smd():
    x1 = [10, 11, 12, 10, 11]
    x0 = [2, 3, 2, 4, 3]
    smd = compute_smd(x1, x0)
    assert smd > 5.0
    
    # Balanced
    smd_bal = compute_smd([5, 5, 5], [5, 5, 5])
    assert smd_bal == 0.0


def test_obo_and_propagation():
    test_dir = Path(__file__).parent.parent / "examples"
    obo_file = test_dir / "mock_go.obo"
    terms, alt, rep, anc, ver = parse_obo(obo_file)
    assert "GO:0009408" in terms
    assert "GO:0006950" in anc["GO:0009408"]
    assert "GO:0008150" in anc["GO:0009408"]
    
    direct = {"OG1": {"GO:0009408"}}
    prop = propagate_go_annotations(direct, anc)
    # Root term GO:0008150 should be excluded, parent GO:0006950 included
    assert "GO:0006950" in prop["OG1"]
    assert "GO:0008150" not in prop["OG1"]


def test_pipeline_execution():
    test_dir = Path(__file__).parent.parent / "examples"
    tmp_out = Path(tempfile.mkdtemp())
    
    class Args:
        master_table = str(test_dir / "mock_master_table.tsv")
        annotations = str(test_dir / "mock_annotations.tsv")
        foreground = "is_expanded"
        covariates = "Total_Count,Count_Ath,Annotation_Degree"
        obo = str(test_dir / "mock_go.obo")
        permutations = 200
        min_universe = 2
        max_universe_fraction = 0.95
        seed = 42
        out_dir = str(tmp_out)
        
    try:
        run_pipeline(Args())
        assert (tmp_out / "enrichment_results_all.tsv").exists()
        assert (tmp_out / "enrichment_results_significant.tsv").exists()
        
        df_sig = pd.read_csv(tmp_out / "enrichment_results_significant.tsv", sep="\t")
        assert not df_sig.empty
        # GO:0006950 should be significantly enriched
        assert "GO:0006950" in df_sig["Term_ID"].values
    finally:
        shutil.rmtree(tmp_out)

