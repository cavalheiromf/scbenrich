import json
from argparse import Namespace

import pandas as pd
import pytest

from scbenrich.core.stats import exact_conditional_or
from scbenrich.core.dag import reduce_go_redundancy
from scbenrich.cli import run_pipeline


def test_conditional_or_rejects_invalid_and_marks_degenerate():
    with pytest.raises(ValueError):
        exact_conditional_or(-1, 2, 3, 4)
    assert exact_conditional_or(0, 0, 2, 2) == (float("nan"), float("nan"), float("nan")) or all(
        pd.isna(x) for x in exact_conditional_or(0, 0, 2, 2)
    )


def test_reduction_map_records_reason():
    df = pd.DataFrame({"Term_ID": ["child", "parent"], "q_value": [0.001, 0.01]})
    reduced, mapping = reduce_go_redundancy(df, {"child": {"parent"}}, return_mapping=True)
    assert reduced.Term_ID.tolist() == ["child"]
    assert mapping.to_dict("records") == [{"Term_ID": "parent", "Representative": "child", "Reason": "ancestor"}]


def test_pipeline_exports_manifest_and_resolves_alt_go_id(tmp_path):
    master = tmp_path / "master.tsv"
    ann = tmp_path / "ann.tsv"
    obo = tmp_path / "go.obo"
    pd.DataFrame({"Family_ID": list("ABCD"), "size": [1, 1, 1, 1], "selected": [1, 0, 1, 0]}).to_csv(master, sep="\t", index=False)
    pd.DataFrame({"Family_ID": ["A", "B"], "Term_ID": ["GO:0000002", "GO:0000001"]}).to_csv(ann, sep="\t", index=False)
    obo.write_text("format-version: 1.2\ndata-version: test\n\n[Term]\nid: GO:0000001\nname: test\nnamespace: biological_process\nalt_id: GO:0000002\n")
    args = Namespace(master_table=str(master), annotations=str(ann), foreground="selected", covariates="size", obo=str(obo),
                     permutations=10, min_universe=1, max_universe_fraction=0.9, seed=1, out_dir=str(tmp_path / "out"))
    run_pipeline(args)
    result = pd.read_csv(tmp_path / "out/enrichment_results_all.tsv", sep="\t")
    assert result.Term_ID.tolist() == ["GO:0000001"]
    manifest = json.loads((tmp_path / "out/run_metadata.json").read_text())
    assert manifest["seed"] == 1 and manifest["go_version"] == "test"
    assert "sha256" in manifest["inputs"]["master_table"]
