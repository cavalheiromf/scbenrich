"""
Data loaders, matrix formatters, and export helpers.
"""

import os
from pathlib import Path
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
    df = pd.read_csv(annotation_path, sep=sep)
    
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
    df = pd.read_csv(master_path, sep=sep)
    
    # Standardize family id column
    for c in df.columns:
        if c.lower().strip() in ("family_id", "familyid", "orthogroup", "og", "family"):
            df = df.rename(columns={c: "Family_ID"})
            break
            
    if "Family_ID" in df.columns:
        df = df.set_index("Family_ID")
        
    return df


def save_results(df, output_path):
    """Saves DataFrame as TSV or CSV based on extension."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sep = "\t" if output_path.suffix in (".tsv", ".tab") else ","
    df.to_csv(output_path, sep=sep, index=False)

