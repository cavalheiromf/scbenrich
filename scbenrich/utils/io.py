"""
Data loaders, matrix formatters, and export helpers.
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd


def load_annotations(annotation_path):
    """
    Loads family-to-term mapping table.
    
    Expected format: TSV or CSV with columns:
    - Family_ID (or Orthogroup)
    - Term_ID (e.g. GO:0006950 or ath00940)
    - Term_Name (optional)
    - Category / Namespace (optional: BP, MF, CC, KEGG)
    """
    sep = "\t" if str(annotation_path).endswith((".tsv", ".tab", ".txt")) else ","
    df = pd.read_csv(annotation_path, sep=sep, dtype=str)
    
    # Standardize column names
    col_map = {}
    for c in df.columns:
        cl = c.lower().strip()
        if cl in ("family_id", "familyid", "orthogroup", "og", "family"):
            col_map[c] = "Family_ID"
        elif cl in ("term_id", "termid", "go_id", "goid", "term"):
            col_map[c] = "Term_ID"
        elif cl in ("term_name", "termname", "name", "description"):
            col_map[c] = "Term_Name"
        elif cl in ("category", "namespace", "aspect", "ontology"):
            col_map[c] = "Category"
            
    df = df.rename(columns=col_map)
    return df


def load_master_table(master_path):
    """
    Loads master table with family metrics, foreground indicators, and covariates.
    
    Expected format: TSV or CSV indexed by Family_ID.
    """
    sep = "\t" if str(master_path).endswith((".tsv", ".tab", ".txt")) else ","
    header = pd.read_csv(master_path, sep=sep, nrows=0)
    id_cols = [c for c in header if c.lower().strip() in ("family_id", "familyid", "orthogroup", "og", "family")]
    if not id_cols:
        raise ValueError("master table requires Family_ID column")
    df = pd.read_csv(master_path, sep=sep, dtype={id_cols[0]: str})
    
    # Standardize family id column
    for c in df.columns:
        if c.lower().strip() in ("family_id", "familyid", "orthogroup", "og", "family"):
            df = df.rename(columns={c: "Family_ID"})
            break
            
    if "Family_ID" in df.columns:
        df = df.set_index("Family_ID")
        
    return df


def build_master_table(
    gene_counts,
    annotations=None,
    foreground=None,
    foreground_col="is_foreground",
    covariates=None,
    output_path=None,
):
    """
    Constructs a validated master table for scb-enrich from gene counts, annotations, and foreground data.

    Parameters
    ----------
    gene_counts : str, Path, or pd.DataFrame
        Table of gene family counts (e.g. OrthoFinder Orthogroups.GeneCount.tsv).
        Must contain a family ID column (Family_ID, Orthogroup, OG, etc.).
    annotations : str, Path, or pd.DataFrame, optional
        Table of family-to-term annotations. If provided, computes 'Annotation_Degree'
        (count of unique terms per family). Families with no annotations receive 0.
    foreground : set, list, tuple, pd.Series, str, or Path, optional
        Foreground family IDs (as a collection or text file with one ID per line),
        or the name of an existing column in gene_counts.
    foreground_col : str, default 'is_foreground'
        Column name to use for the binary foreground indicator (1 for foreground, 0 for background).
    covariates : list of str, optional
        Specific covariate columns to retain. If None, all numeric columns from gene_counts
        (plus Annotation_Degree if computed) are kept.
    output_path : str or Path, optional
        If provided, writes the master table to TSV or CSV.

    Returns
    -------
    pd.DataFrame
        Validated master table indexed by Family_ID.
    """
    # 1. Parse and standardize gene_counts
    if isinstance(gene_counts, (str, Path)):
        sep = "\t" if str(gene_counts).endswith((".tsv", ".tab", ".txt")) else ","
        df = pd.read_csv(gene_counts, sep=sep)
    elif isinstance(gene_counts, pd.DataFrame):
        df = gene_counts.copy()
    else:
        raise TypeError("gene_counts must be a file path (str or Path) or pandas DataFrame")

    id_aliases = ("family_id", "familyid", "orthogroup", "og", "family")
    id_col = None
    if df.index.name and df.index.name.lower().strip() in id_aliases:
        df = df.reset_index()
        id_col = df.columns[0]
    else:
        for c in df.columns:
            if str(c).lower().strip() in id_aliases:
                id_col = c
                break

    if id_col is None:
        raise ValueError("gene_counts must contain a Family_ID or Orthogroup column")

    df = df.rename(columns={id_col: "Family_ID"})
    df["Family_ID"] = df["Family_ID"].astype(str).str.strip()

    if df["Family_ID"].duplicated().any():
        dups = df["Family_ID"][df["Family_ID"].duplicated()].tolist()[:5]
        raise ValueError(f"Duplicate Family_IDs found in gene_counts: {dups}")
    if df["Family_ID"].isna().any() or (df["Family_ID"] == "").any():
        raise ValueError("Empty or NaN Family_IDs found in gene_counts")

    # 2. Compute Annotation_Degree if annotations provided
    if annotations is not None:
        if isinstance(annotations, (str, Path)):
            df_annot = load_annotations(annotations)
        elif isinstance(annotations, pd.DataFrame):
            df_annot = annotations.copy()
            col_map = {}
            for c in df_annot.columns:
                cl = str(c).lower().strip()
                if cl in id_aliases:
                    col_map[c] = "Family_ID"
                elif cl in ("term_id", "termid", "go_id", "goid", "term"):
                    col_map[c] = "Term_ID"
            df_annot = df_annot.rename(columns=col_map)
        else:
            raise TypeError("annotations must be a file path or pandas DataFrame")

        if "Family_ID" not in df_annot.columns or "Term_ID" not in df_annot.columns:
            raise ValueError("annotations must contain Family_ID and Term_ID columns")

        df_annot["Family_ID"] = df_annot["Family_ID"].astype(str).str.strip()
        df_annot["Term_ID"] = df_annot["Term_ID"].astype(str).str.strip()

        degree = (
            df_annot.groupby("Family_ID")["Term_ID"]
            .nunique()
            .rename("Annotation_Degree")
        )

        if "Annotation_Degree" in df.columns:
            df = df.drop(columns=["Annotation_Degree"])
        df = df.merge(degree, on="Family_ID", how="left")
        df["Annotation_Degree"] = df["Annotation_Degree"].fillna(0).astype(int)

    # 3. Handle Foreground
    if foreground is not None:
        if isinstance(foreground, str) and foreground in df.columns:
            col = df[foreground]
            is_fg = col.astype(str).str.lower().isin(["1", "true", "yes", "expanded"])
            df[foreground_col] = is_fg.astype(int)
        else:
            if isinstance(foreground, (str, Path)):
                fg_path = Path(foreground)
                if not fg_path.exists():
                    raise FileNotFoundError(f"Foreground file not found: {foreground}")
                with open(fg_path) as f:
                    fg_set = set(line.strip() for line in f if line.strip())
            elif isinstance(foreground, (set, list, tuple, pd.Series)):
                fg_set = set(str(x).strip() for x in foreground)
            else:
                raise TypeError("foreground must be a collection, file path, or column name in gene_counts")

            overlap = fg_set.intersection(set(df["Family_ID"]))
            if not overlap:
                raise ValueError("None of the foreground IDs were found in the universe of gene_counts")
            df[foreground_col] = df["Family_ID"].isin(fg_set).astype(int)

    # 4. Select and validate Covariates
    if covariates is not None:
        selected_covs = list(covariates)
        for c in selected_covs:
            if c not in df.columns:
                raise ValueError(f"Covariate column '{c}' not found in data")
            try:
                s = pd.to_numeric(df[c], errors="raise")
            except (ValueError, TypeError):
                raise ValueError(f"Covariate column '{c}' could not be converted to numeric")
            if not np.isfinite(s).all():
                raise ValueError(f"Covariate column '{c}' contains non-finite values (NaN or Inf)")
            df[c] = s
    else:
        exclude = {"Family_ID"}
        if foreground is not None:
            exclude.add(foreground_col)
            if isinstance(foreground, str) and foreground in df.columns:
                exclude.add(foreground)
        numeric_cols = []
        for c in df.columns:
            if c not in exclude:
                try:
                    s = pd.to_numeric(df[c], errors="coerce")
                    if s.notna().all() and np.isfinite(s).all():
                        df[c] = s
                        numeric_cols.append(c)
                except Exception:
                    pass
        selected_covs = numeric_cols

    # 5. Format final columns
    final_cols = ["Family_ID"] + selected_covs
    if foreground is not None and foreground_col not in final_cols:
        final_cols.append(foreground_col)

    df_out = df[final_cols].copy()

    # 6. Save if output_path is provided
    if output_path is not None:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        sep = "\t" if str(out).endswith((".tsv", ".tab", ".txt")) else ","
        df_out.to_csv(out, sep=sep, index=False)

    return df_out.set_index("Family_ID")


def save_results(df, output_path):
    """Saves DataFrame as TSV or CSV based on extension."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sep = "\t" if output_path.suffix in (".tsv", ".tab") else ","
    df.to_csv(output_path, sep=sep, index=False)


