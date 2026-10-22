import pytest
import numpy as np
import pandas as pd
from pathlib import Path
from scbenrich import build_master_table, load_master_table


def test_build_master_table_basic():
    counts_df = pd.DataFrame({
        "Orthogroup": ["OG1", "OG2", "OG3"],
        "Total_Count": [10, 20, 30],
        "Count_Sp1": [5, 10, 15],
    })
    master = build_master_table(counts_df)
    assert master.index.name == "Family_ID"
    assert list(master.index) == ["OG1", "OG2", "OG3"]
    assert "Total_Count" in master.columns
    assert "Count_Sp1" in master.columns
    assert len(master) == 3


def test_build_master_table_with_annotations(tmp_path):
    counts_df = pd.DataFrame({
        "Family_ID": ["OG1", "OG2", "OG3"],
        "Size": [10, 20, 30],
    })
    annot_df = pd.DataFrame({
        "Family_ID": ["OG1", "OG1", "OG2"],
        "Term_ID": ["GO:001", "GO:002", "GO:001"],
    })
    
    # Using DataFrame
    master = build_master_table(counts_df, annotations=annot_df)
    assert "Annotation_Degree" in master.columns
    assert master.loc["OG1", "Annotation_Degree"] == 2
    assert master.loc["OG2", "Annotation_Degree"] == 1
    assert master.loc["OG3", "Annotation_Degree"] == 0

    # Using File path
    annot_path = tmp_path / "annotations.tsv"
    annot_df.to_csv(annot_path, sep="\t", index=False)
    master_file = build_master_table(counts_df, annotations=annot_path)
    pd.testing.assert_frame_equal(master, master_file)


def test_build_master_table_with_foreground(tmp_path):
    counts_df = pd.DataFrame({
        "Family_ID": ["OG1", "OG2", "OG3"],
        "Size": [10, 20, 30],
    })
    
    # 1. Using a set/list
    master = build_master_table(counts_df, foreground={"OG1", "OG3"}, foreground_col="is_expanded")
    assert "is_expanded" in master.columns
    assert master.loc["OG1", "is_expanded"] == 1
    assert master.loc["OG2", "is_expanded"] == 0
    assert master.loc["OG3", "is_expanded"] == 1

    # 2. Using a text file
    fg_file = tmp_path / "foreground.txt"
    fg_file.write_text("OG1\nOG3\n")
    master_txt = build_master_table(counts_df, foreground=fg_file, foreground_col="is_expanded")
    pd.testing.assert_frame_equal(master, master_txt)

    # 3. Using existing column in counts
    counts_with_fg = counts_df.copy()
    counts_with_fg["expanded_status"] = ["yes", "no", "yes"]
    master_col = build_master_table(counts_with_fg, foreground="expanded_status", foreground_col="is_expanded")
    pd.testing.assert_frame_equal(master[["Size", "is_expanded"]], master_col[["Size", "is_expanded"]])


def test_build_master_table_save_and_reload(tmp_path):
    counts_path = tmp_path / "counts.tsv"
    counts_df = pd.DataFrame({
        "Orthogroup": ["OG1", "OG2", "OG3"],
        "Size": [10, 20, 30],
    })
    counts_df.to_csv(counts_path, sep="\t", index=False)

    out_path = tmp_path / "master_table.tsv"
    master = build_master_table(
        gene_counts=counts_path,
        foreground={"OG2"},
        foreground_col="is_expanded",
        output_path=out_path,
    )

    assert out_path.exists()
    reloaded = load_master_table(out_path)
    assert reloaded.index.name == "Family_ID"
    assert list(reloaded.index) == ["OG1", "OG2", "OG3"]
    assert reloaded.loc["OG2", "is_expanded"] == 1
    assert reloaded.loc["OG1", "is_expanded"] == 0


def test_build_master_table_errors(tmp_path):
    # Duplicate IDs
    dups = pd.DataFrame({"Family_ID": ["OG1", "OG1"], "Size": [1, 2]})
    with pytest.raises(ValueError, match="Duplicate Family_IDs"):
        build_master_table(dups)

    # Missing ID column
    no_id = pd.DataFrame({"ColA": [1, 2], "ColB": [3, 4]})
    with pytest.raises(ValueError, match="must contain a Family_ID"):
        build_master_table(no_id)

    # Non-existent foreground file
    valid = pd.DataFrame({"Family_ID": ["OG1", "OG2"], "Size": [1, 2]})
    with pytest.raises(FileNotFoundError):
        build_master_table(valid, foreground=tmp_path / "missing.txt")

    # Foreground has no overlap with universe
    with pytest.raises(ValueError, match="None of the foreground IDs were found"):
        build_master_table(valid, foreground={"UNKNOWN1", "UNKNOWN2"})

    # Non-numeric covariate
    bad_cov = pd.DataFrame({"Family_ID": ["OG1", "OG2"], "Size": ["abc", "def"]})
    with pytest.raises(ValueError, match="could not be converted to numeric"):
        build_master_table(bad_cov, covariates=["Size"])

    # Non-finite covariate
    nan_cov = pd.DataFrame({"Family_ID": ["OG1", "OG2"], "Size": [1.0, np.nan]})
    with pytest.raises(ValueError, match="contains non-finite"):
        build_master_table(nan_cov, covariates=["Size"])


def test_load_gene_to_family_mapping(tmp_path):
    from scbenrich import load_gene_to_family_mapping

    # 1. From dict
    d = {"g1": "OG1", "g2": "OG2"}
    assert load_gene_to_family_mapping(d) == d

    # 2. From Orthogroups.txt
    txt_path = tmp_path / "Orthogroups.txt"
    txt_path.write_text("OG0000001: geneA1 geneA2\nOG0000002: geneB1\n")
    m_txt = load_gene_to_family_mapping(txt_path)
    assert m_txt == {"geneA1": "OG0000001", "geneA2": "OG0000001", "geneB1": "OG0000002"}

    # 3. From Orthogroups.tsv
    tsv_path = tmp_path / "Orthogroups.tsv"
    tsv_path.write_text("Orthogroup\tSp1\tSp2\nOG1\tgA1, gA2\tgB1\nOG2\tgA3\t\n")
    m_tsv = load_gene_to_family_mapping(tsv_path)
    assert m_tsv == {"gA1": "OG1", "gA2": "OG1", "gB1": "OG1", "gA3": "OG2"}


def test_extract_gff3_covariates_and_integration(tmp_path):
    from scbenrich import extract_gff3_covariates

    gff_content = (
        "##gff-version 3\n"
        "chr1\tsrc\tgene\t1000\t2500\t.\t+\t.\tID=gene:g1\n"
        "chr1\tsrc\tmRNA\t1000\t2500\t.\t+\t.\tID=transcript:g1.1;Parent=gene:g1\n"
        "chr1\tsrc\texon\t1000\t1200\t.\t+\t.\tParent=transcript:g1.1\n"
        "chr1\tsrc\tCDS\t1050\t1200\t.\t+\t0\tParent=transcript:g1.1\n"
        "chr1\tsrc\texon\t1800\t2500\t.\t+\t.\tParent=transcript:g1.1\n"
        "chr1\tsrc\tCDS\t1800\t2400\t.\t+\t0\tParent=transcript:g1.1\n"
        "chr1\tsrc\tgene\t5000\t6000\t.\t-\t.\tID=g2\n"
        "chr1\tsrc\tmRNA\t5000\t6000\t.\t-\t.\tID=g2.1;Parent=g2\n"
        "chr1\tsrc\texon\t5000\t6000\t.\t-\t.\tParent=g2.1\n"
        "chr1\tsrc\tCDS\t5100\t5900\t.\t-\t0\tParent=g2.1\n"
    )
    gff_file = tmp_path / "sample.gff3"
    gff_file.write_text(gff_content)

    mapping = {"g1": "OG1", "g2": "OG1"}
    covs = extract_gff3_covariates(gff_file, mapping, agg_func="mean")
    assert "OG1" in covs.index
    # g1: gene_length=1501, cds_length=(151+601)=752, exon_count=2
    # g2: gene_length=1001, cds_length=801, exon_count=1
    # mean gene_length: (1501 + 1001)/2 = 1251.0
    assert covs.loc["OG1", "Mean_Gene_Length"] == 1251.0
    assert covs.loc["OG1", "Mean_Exon_Count"] == 1.5

    # Integration test with build_master_table
    counts_df = pd.DataFrame({
        "Orthogroup": ["OG1", "OG_MISSING"],
        "Total_Count": [2, 5],
    })
    master = build_master_table(
        gene_counts=counts_df,
        gff3=gff_file,
        gene_to_family=mapping,
        fill_gff_missing="median",
    )
    assert "Mean_Gene_Length" in master.columns
    assert "Mean_CDS_Length" in master.columns
    assert "Mean_Exon_Count" in master.columns
    assert master.loc["OG1", "Mean_Gene_Length"] == 1251.0
    # OG_MISSING gets filled with median
    assert master.loc["OG_MISSING", "Mean_Gene_Length"] == 1251.0

