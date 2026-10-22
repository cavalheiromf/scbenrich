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


def load_gene_to_family_mapping(mapping_source):
    """
    Loads gene-to-family mapping from OrthoFinder results, files, or dictionaries.

    Supported formats:
    - dict or pd.Series: mapping gene_id -> family_id.
    - Orthogroups.txt: text file with lines like 'OG0000001: gene1 gene2 ...'
    - Orthogroups.tsv: table where col 0 is Orthogroup, and other cols are species with comma-separated gene IDs.
    - 2-column TSV/CSV: columns for Gene_ID and Family_ID.
    """
    if isinstance(mapping_source, dict):
        return {str(k).strip(): str(v).strip() for k, v in mapping_source.items()}
    elif isinstance(mapping_source, pd.Series):
        return {str(k).strip(): str(v).strip() for k, v in mapping_source.items()}
    elif isinstance(mapping_source, pd.DataFrame):
        df = mapping_source
    elif isinstance(mapping_source, (str, Path)):
        path = Path(mapping_source)
        if not path.exists():
            raise FileNotFoundError(f"Gene-to-family mapping file not found: {mapping_source}")

        if str(path).endswith(".txt"):
            mapping = {}
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                first_line = f.readline()
                f.seek(0)
                if ":" in first_line and not first_line.startswith("#"):
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("#") or ":" not in line:
                            continue
                        og, genes = line.split(":", 1)
                        og = og.strip()
                        for g in genes.split():
                            mapping[g.strip()] = og
                    return mapping

        sep = "\t" if str(path).endswith((".tsv", ".tab", ".txt")) else ","
        df = pd.read_csv(path, sep=sep, dtype=str)
    else:
        raise TypeError("mapping_source must be a dict, pd.Series, pd.DataFrame, or file path")

    mapping = {}
    id_aliases = ("family_id", "familyid", "orthogroup", "og", "family")

    # 2-column DataFrame
    if df.shape[1] == 2:
        cols = list(df.columns)
        if any(c.lower().strip() in id_aliases for c in cols):
            og_col = [c for c in cols if c.lower().strip() in id_aliases][0]
            gene_col = [c for c in cols if c != og_col][0]
        else:
            og_col, gene_col = cols[1], cols[0]
        for _, row in df.iterrows():
            g = str(row[gene_col]).strip()
            fam = str(row[og_col]).strip()
            if g and fam and g != "nan" and fam != "nan":
                mapping[g] = fam
        return mapping

    # Orthogroups.tsv format
    og_col = df.columns[0]
    for _, row in df.iterrows():
        fam = str(row[og_col]).strip()
        if not fam or fam == "nan":
            continue
        for col in df.columns[1:]:
            val = row[col]
            if pd.notna(val):
                for g in str(val).split(","):
                    g = g.strip()
                    if g and g != "nan":
                        mapping[g] = fam
    return mapping


def _resolve_gff_files(gff_input):
    """
    Resolves gff_input into a list of (label, Path) pairs.
    Supports:
    - single path (str or Path)
    - list/tuple/set of paths
    - directory path (finds all .gff, .gff3, .gtf, and .gz variants)
    - dict of {label/species: path}
    """
    files = []
    if isinstance(gff_input, dict):
        for label, p in gff_input.items():
            path = Path(p)
            if not path.exists():
                raise FileNotFoundError(f"GFF file not found: {path}")
            files.append((str(label), path))
    elif isinstance(gff_input, (list, tuple, set)):
        for p in gff_input:
            path = Path(p)
            if not path.exists():
                raise FileNotFoundError(f"GFF file not found: {path}")
            files.append((path.stem.replace(".gff3", "").replace(".gff", "").replace(".gtf", ""), path))
    elif isinstance(gff_input, (str, Path)):
        path = Path(gff_input)
        if not path.exists():
            raise FileNotFoundError(f"GFF path not found: {path}")
        if path.is_dir():
            patterns = ("*.gff", "*.gff3", "*.gtf", "*.gff.gz", "*.gff3.gz", "*.gtf.gz")
            found = []
            for pat in patterns:
                found.extend(path.glob(pat))
            found = sorted(set(found))
            if not found:
                raise FileNotFoundError(f"No GFF/GTF files found in directory: {path}")
            for p in found:
                label = p.name
                for ext in (".gff3.gz", ".gff.gz", ".gtf.gz", ".gff3", ".gff", ".gtf"):
                    if label.endswith(ext):
                        label = label[:-len(ext)]
                        break
                files.append((label, p))
        else:
            files.append((path.stem.replace(".gff3", "").replace(".gff", "").replace(".gtf", ""), path))
    else:
        raise TypeError("gff_input must be a file path, directory path, list of paths, or dict")
    return files


def _parse_gff_file(gff_path):
    def strip_prefix(s):
        if not s:
            return s
        for p in ("gene:", "gene-", "transcript:", "transcript-", "rna:", "rna-", "protein:"):
            if s.startswith(p):
                return s[len(p):]
        return s

    genes = {}
    from collections import defaultdict
    txs = defaultdict(lambda: {"gene": None, "cds_len": 0, "exon_cnt": 0, "len": 0})

    if str(gff_path).endswith(".gz"):
        import gzip
        opener = gzip.open(gff_path, "rt", encoding="utf-8", errors="replace")
    else:
        opener = open(gff_path, "r", encoding="utf-8", errors="replace")

    with opener as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            _, _, ftype, start_s, end_s, _, _, _, attr_str = parts
            try:
                start, end = int(start_s), int(end_s)
            except ValueError:
                continue
            length = max(1, end - start + 1)
            ftype_lower = ftype.lower()

            attrs = {}
            for item in attr_str.split(";"):
                item = item.strip()
                if not item:
                    continue
                if "=" in item:
                    k, v = item.split("=", 1)
                    attrs[k.strip()] = v.strip().strip('"').strip("'")
                elif " " in item:
                    subparts = item.split(None, 1)
                    if len(subparts) == 2:
                        attrs[subparts[0].strip()] = subparts[1].strip().strip('"').strip("'")

            fid = attrs.get("ID") or attrs.get("gene_id") or attrs.get("transcript_id")
            parent = attrs.get("Parent") or attrs.get("transcript_id") or attrs.get("gene_id")

            clean_fid = strip_prefix(fid) if fid else None
            clean_parent = strip_prefix(parent) if parent else None

            if ftype_lower == "gene":
                if clean_fid:
                    genes[clean_fid] = length
            elif ftype_lower in ("mrna", "transcript"):
                if clean_fid:
                    txs[clean_fid]["len"] = length
                    if clean_parent:
                        txs[clean_fid]["gene"] = clean_parent
            elif ftype_lower == "exon":
                if clean_parent:
                    for p in clean_parent.split(","):
                        txs[strip_prefix(p.strip())]["exon_cnt"] += 1
            elif ftype_lower == "cds":
                if clean_parent:
                    for p in clean_parent.split(","):
                        txs[strip_prefix(p.strip())]["cds_len"] += length

    gene_to_txs = defaultdict(list)
    for tid, tinfo in txs.items():
        if tinfo["gene"]:
            gene_to_txs[tinfo["gene"]].append(tinfo)

    features = {}
    for gid in set(genes.keys()).union(set(gene_to_txs.keys())):
        isoforms = gene_to_txs.get(gid, [])
        g_len = genes.get(gid, max([t["len"] for t in isoforms], default=1))
        cds_len = max([t["cds_len"] for t in isoforms], default=0)
        exon_cnt = max([t["exon_cnt"] for t in isoforms], default=1)
        features[gid] = {
            "gene_length": g_len if g_len > 0 else 1,
            "cds_length": cds_len,
            "exon_count": exon_cnt if exon_cnt > 0 else 1,
        }
    for tid, t in txs.items():
        if tid not in features:
            g_len = genes.get(t["gene"], t["len"])
            features[tid] = {
                "gene_length": g_len if g_len > 0 else 1,
                "cds_length": t["cds_len"],
                "exon_count": t["exon_cnt"] if t["exon_cnt"] > 0 else 1,
            }
    return features


def extract_gff3_covariates(gff, gene_to_family, agg_func="mean", by_species=False):
    """
    Extracts structural covariates (gene length, CDS length, exon count) from one or multiple
    GFF3/GTF files and aggregates them by gene family.

    Parameters
    ----------
    gff : str, Path, list, tuple, set, or dict
        Path to a GFF3/GTF file, list/tuple of file paths, directory containing GFF files,
        or dict mapping species names to file paths.
    gene_to_family : dict, pd.Series, pd.DataFrame, str, or Path
        Mapping from gene/transcript ID to Family_ID (e.g. from OrthoFinder).
    agg_func : {"mean", "median"}, default "mean"
        Aggregation function across genes in each family.
    by_species : bool, default False
        If True and multiple species/GFF files are provided, also adds species-specific
        covariates in addition to the pooled metrics.

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by Family_ID with columns:
        - Mean_Gene_Length (or Median_Gene_Length)
        - Mean_CDS_Length (or Median_CDS_Length)
        - Mean_Exon_Count (or Median_Exon_Count)
        (Plus species-specific columns if by_species=True).
    """
    files = _resolve_gff_files(gff)
    mapping = load_gene_to_family_mapping(gene_to_family)
    if not mapping:
        raise ValueError("gene_to_family mapping is empty")

    records = []
    seen_genes = set()
    for label, file_path in files:
        file_features = _parse_gff_file(file_path)
        for gid, metrics in file_features.items():
            fam = mapping.get(gid)
            if not fam:
                continue
            if (gid, fam) in seen_genes:
                continue
            seen_genes.add((gid, fam))
            records.append({
                "Gene_ID": gid,
                "Family_ID": fam,
                "Species": label,
                **metrics,
            })

    if not records:
        raise ValueError("None of the gene IDs in gene_to_family were found in the provided GFF file(s)")

    df_genes = pd.DataFrame(records)

    agg = "median" if agg_func == "median" else "mean"
    prefix = "Median" if agg == "median" else "Mean"

    fam_df = df_genes.groupby("Family_ID")[["gene_length", "cds_length", "exon_count"]].agg(agg).rename(columns={
        "gene_length": f"{prefix}_Gene_Length",
        "cds_length": f"{prefix}_CDS_Length",
        "exon_count": f"{prefix}_Exon_Count",
    })

    if by_species and df_genes["Species"].nunique() > 1:
        sp_df = df_genes.pivot_table(
            index="Family_ID",
            columns="Species",
            values=["gene_length", "cds_length", "exon_count"],
            aggfunc=agg
        )
        sp_cols = []
        for val_name, sp_name in sp_df.columns:
            clean_val = val_name.replace("gene_length", "Gene_Length").replace("cds_length", "CDS_Length").replace("exon_count", "Exon_Count")
            sp_cols.append(f"{prefix}_{clean_val}_{sp_name}")
        sp_df.columns = sp_cols
        fam_df = fam_df.join(sp_df, how="left")

    return fam_df


def build_master_table(
    gene_counts,
    annotations=None,
    foreground=None,
    foreground_col="is_foreground",
    covariates=None,
    gff3=None,
    gene_to_family=None,
    gff_agg="mean",
    by_species=False,
    fill_gff_missing="median",
    output_path=None,
):
    """
    Constructs a validated master table for scb-enrich from gene counts, annotations, foreground data,
    and optional GFF3 genomic features from one or multiple species.

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
        (plus Annotation_Degree and GFF3 features if computed) are kept.
    gff3 : str, Path, list, tuple, set, or dict, optional
        Path to a GFF3/GTF file, list of files, directory containing GFFs, or dict of {species: path}.
    gene_to_family : dict, pd.Series, pd.DataFrame, str, or Path, optional
        Mapping from gene IDs to Family_ID (e.g. OrthoFinder Orthogroups.tsv / Orthogroups.txt).
        Required if gff3 is provided.
    gff_agg : {"mean", "median"}, default "mean"
        Aggregation function across genes in each family for GFF features.
    by_species : bool, default False
        If True and multiple species GFFs are provided, adds species-specific structural covariates.
    fill_gff_missing : {"median", "mean", float}, default "median"
        Strategy to fill structural covariates for families not present in the GFF file.
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

    # 3. Extract GFF3 Covariates if provided
    if gff3 is not None:
        if gene_to_family is None:
            raise ValueError("Extracting covariates from GFF3 requires gene_to_family mapping")
        df_gff = extract_gff3_covariates(gff3, gene_to_family, agg_func=gff_agg, by_species=by_species)
        for c in df_gff.columns:
            if c in df.columns:
                df = df.drop(columns=[c])
        df = df.merge(df_gff, on="Family_ID", how="left")
        for c in df_gff.columns:
            if df[c].isna().any():
                if fill_gff_missing == "median":
                    fill_val = float(df_gff[c].median()) if len(df_gff) > 0 and pd.notna(df_gff[c].median()) else 0.0
                elif fill_gff_missing == "mean":
                    fill_val = float(df_gff[c].mean()) if len(df_gff) > 0 and pd.notna(df_gff[c].mean()) else 0.0
                else:
                    fill_val = float(fill_gff_missing)
                df[c] = df[c].fillna(fill_val)

    # 4. Handle Foreground
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


