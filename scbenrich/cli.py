"""
Command Line Interface (CLI) for SCB-Enrich.
"""

import sys
import argparse
import hashlib
import json
import platform
from importlib.metadata import version
from pathlib import Path
import numpy as np
import pandas as pd

from scbenrich.utils.validation import validate_parameters, validate_family_inputs
from scbenrich.core.inference import decision_fields
from scbenrich.core.stats import fdr_bh, fisher_exact_test, exact_conditional_or
from scbenrich.core.stratified import StratifiedPermuter, compute_smd
from scbenrich.core.obo import parse_obo, propagate_go_annotations, normalize_go_id
from scbenrich.core.dag import reduce_go_redundancy
from scbenrich.utils.io import load_annotations, load_master_table, save_results


def _run_pipeline(args):
    validate_parameters(args.permutations, getattr(args, "bins", None))
    if args.permutations < 0:
        raise ValueError("permutations must be nonnegative")
    if args.min_universe < 1:
        raise ValueError("min_universe must be positive")
    if not 0 < args.max_universe_fraction <= 1:
        raise ValueError("max_universe_fraction must be in (0, 1]")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"[*] Starting SCB-Enrich analysis...")
    print(f"    Master table: {args.master_table}")
    print(f"    Annotations:  {args.annotations}")
    print(f"    Output dir:   {out_dir}")
    
    # 1. Load Master Table & Covariates
    df_master = load_master_table(args.master_table)
    if df_master.index.has_duplicates or df_master.index.isna().any():
        raise ValueError("Family_ID values must be unique and nonempty")
    all_families = list(df_master.index)
    print(f"[*] Loaded {len(all_families)} gene families in universe.")
    
    # Identify foreground
    if Path(args.foreground).is_file():
        with open(args.foreground) as f:
            foreground_set = set(line.strip() for line in f if line.strip())
    elif args.foreground in df_master.columns:
        col = df_master[args.foreground]
        tokens = col.astype(str).str.lower()
        if col.isna().any() or not tokens.isin(["1", "true", "yes", "expanded", "0", "false", "no", "not_expanded"]).all():
            raise ValueError("foreground column contains missing or unknown values")
        if col.dtype == bool:
            foreground_set = set(df_master.index[col])
        else:
            foreground_set = set(df_master.index[col.astype(str).str.lower().isin(["1", "true", "yes", "expanded"])])
    else:
        raise ValueError(f"Foreground '{args.foreground}' is neither an existing file nor a column in master table.")
        
    unknown = foreground_set.difference(all_families)
    if unknown:
        raise ValueError(f"Foreground contains unknown Family_ID values: {sorted(unknown)[:5]}")
    background_set = set(all_families) - foreground_set
    if not foreground_set or not background_set:
        raise ValueError("foreground and background must both be nonempty")
    print(f"    Foreground size: {len(foreground_set)} families.")
    print(f"    Background size: {len(background_set)} families.")
    
    # Parse covariates
    cov_cols = [c.strip() for c in args.covariates.split(",") if c.strip()]
    if not cov_cols:
        raise ValueError("at least one covariate is required")
    for c in cov_cols:
        if c not in df_master.columns:
            raise ValueError(f"Covariate column '{c}' not found in master table.")
    cov_df = df_master[cov_cols].copy()
    for col in cov_cols:
        if not pd.api.types.is_numeric_dtype(cov_df[col]) or not np.isfinite(cov_df[col].to_numpy(dtype=float)).all():
            raise ValueError(f"covariate '{col}' must contain finite numeric values")
    
    # 2. Load Annotations & Propagate if OBO provided
    df_ann = load_annotations(args.annotations)
    if not {"Family_ID", "Term_ID"}.issubset(df_ann.columns):
        raise ValueError("annotations require Family_ID and Term_ID columns")
    if df_ann[["Family_ID", "Term_ID"]].isna().any().any():
        raise ValueError("annotations contain missing Family_ID or Term_ID")
    if not set(df_ann["Family_ID"]).issubset(all_families):
        raise ValueError("annotations contain unknown Family_ID values")
    if df_ann["Term_ID"].str.strip().eq("").any():
        raise ValueError("annotations contain empty Term_ID")
    validate_family_inputs(cov_df, foreground_set)
    go_version = None
    id_changes = []
    
    obo_ancestors = {}
    term_names = {}
    term_categories = {}
    
    if "Term_Name" in df_ann.columns:
        term_names = df_ann.drop_duplicates("Term_ID").set_index("Term_ID")["Term_Name"].to_dict()
    if "Category" in df_ann.columns:
        term_categories = df_ann.drop_duplicates("Term_ID").set_index("Term_ID")["Category"].to_dict()

    if args.obo:
        print(f"[*] Loading Gene Ontology OBO from: {args.obo}")
        obo_terms, alt_map, replaced_map, obo_ancestors, go_version = parse_obo(args.obo)
        def canonical_id(term):
            original = term
            term, path = normalize_go_id(term, obo_terms, alt_map, replaced_map, return_path=True)
            if term != original:
                id_changes.append({"Original_ID": original, "Canonical_ID": term,
                                   "Path": " -> ".join(path), "Reason": "alias_or_replacement"})
            return term
        df_ann["Term_ID"] = df_ann["Term_ID"].map(canonical_id)
        save_results(pd.DataFrame(id_changes, columns=["Original_ID", "Canonical_ID", "Path", "Reason"]), out_dir / "go_id_mappings.tsv")
        for tid, tinfo in obo_terms.items():
            if "name" in tinfo:
                term_names[tid] = tinfo["name"]
            if "namespace" in tinfo:
                term_categories[tid] = tinfo["namespace"]
    direct_map = df_ann.groupby("Family_ID")["Term_ID"].apply(set).to_dict()
    if args.obo:
        print(f"    Propagating terms along GO hierarchical ancestor closure...")
        prop_map = propagate_go_annotations(direct_map, obo_ancestors)
    else:
        prop_map = direct_map

    # 3. Build Term Presence/Absence Matrix
    print("[*] Assembling binary family-term matrix...")
    unique_terms = sorted(list(set(t for terms in prop_map.values() for t in terms)))
    term_universe_counts = {t: 0 for t in unique_terms}
    for fam in all_families:
        for t in prop_map.get(fam, []):
            if t in term_universe_counts:
                term_universe_counts[t] += 1
                
    # Filter universe frequency
    min_count = args.min_universe
    max_count = int(args.max_universe_fraction * len(all_families))
    tested_terms = [t for t, cnt in term_universe_counts.items() if min_count <= cnt <= max_count]
    print(f"    Testing {len(tested_terms)} terms meeting universe thresholds ({min_count} <= size <= {max_count}).")

    # Construct dense/sparse matrix for tested terms
    term_idx = {t: i for i, t in enumerate(tested_terms)}
    tm = np.zeros((len(all_families), len(tested_terms)), dtype=np.int8)
    for i, fam in enumerate(all_families):
        for t in prop_map.get(fam, []):
            if t in term_idx:
                tm[i, term_idx[t]] = 1
                
    term_matrix = pd.DataFrame(tm, index=all_families, columns=tested_terms)
    
    # 4. Parametric Fisher's Exact Test & Exact Conditional OR
    print("[*] Computing parametric Fisher exact tests and exact conditional odds ratios...")
    fg_mask = np.isin(all_families, list(foreground_set))
    bg_mask = ~fg_mask
    
    results = []
    pvals = []
    
    for term in tested_terms:
        col = term_matrix[term].values
        a = int(np.sum(col[fg_mask]))
        b = int(len(foreground_set) - a)
        c = int(np.sum(col[bg_mask]))
        d = int(len(background_set) - c)
        
        _, p_fisher = fisher_exact_test(a, b, c, d, alternative="greater")
        pvals.append(p_fisher)
        
        mle_or, ci_low, ci_high = exact_conditional_or(a, b, c, d)
        
        results.append({
            "Term_ID": term,
            "Term_Name": term_names.get(term, "Uncharacterized"),
            "Category": term_categories.get(term, "General"),
            "Foreground_Count": a,
            "Foreground_Total": len(foreground_set),
            "Background_Count": c,
            "Background_Total": len(background_set),
            "Marginal_OR": round(mle_or, 4),
            "Marginal_CI_Low_95": round(ci_low, 4),
            "Marginal_CI_High_95": round(ci_high, 4),
            "OR_Status": "not_identifiable" if np.isnan(mle_or) else "ok",
            "p_value": p_fisher,
        })
        
    df_res = pd.DataFrame(results, columns=["Term_ID", "Term_Name", "Category", "Foreground_Count",
        "Foreground_Total", "Background_Count", "Background_Total", "Marginal_OR", "Marginal_CI_Low_95",
        "Marginal_CI_High_95", "OR_Status", "p_value"])
    df_res["q_value"] = fdr_bh(pvals)
    
    # 5. Stratified Permutations (Null Resampling)
    n_perms = args.permutations
    if n_perms > 0 and not df_res.empty:
        print(f"[*] Executing {n_perms:,} stratified null permutations balancing covariates: {cov_cols}...")
        permuter = StratifiedPermuter(cov_df, covariates=cov_cols, seed=args.seed, bins=getattr(args, "bins", None))
        (out_dir / "bin_edges.json").write_text(json.dumps(permuter.bin_edges_, indent=2) + "\n")
        support = permuter.support_diagnostics(foreground_set)
        save_results(support, out_dir / "strata_diagnostics.tsv")
        if not support["informative"].any():
            raise ValueError("strata provide no randomization: every stratum is fully selected or unselected")
        p_strat_map = permuter.run_permutations(foreground_set, term_matrix, n_permutations=n_perms)
        df_res["p_strat"] = [p_strat_map.get(t, 1.0) for t in df_res["Term_ID"]]
        df_res["p_strat_adjusted"] = fdr_bh(df_res["p_strat"].to_numpy())
        diagnostics = permuter.term_diagnostics_.set_index("Term_ID")
        for col in diagnostics.columns:
            if col != "p_strat":
                df_res[col] = df_res.Term_ID.map(diagnostics[col])
        save_results(permuter.balance_diagnostics_, out_dir / "balance_diagnostics.tsv")
        
        # Calculate maximum standardized mean difference across covariates
        smds = []
        for _, row in df_res.iterrows():
            tid = row["Term_ID"]
            t_fams = set(term_matrix.index[term_matrix[tid] == 1])
            fg_t = list(t_fams.intersection(foreground_set))
            bg_t = list(t_fams.intersection(background_set))
            max_smd = float("nan")
            if fg_t and bg_t:
                values = []
                for col in cov_cols:
                    smd = abs(compute_smd(cov_df.loc[fg_t, col], cov_df.loc[bg_t, col]))
                    values.append(smd)
                if values and np.all(np.isfinite(values)):
                    max_smd = max(values)
            smds.append(max_smd)
        df_res["Max_SMD"] = smds
        
    else:
        df_res["p_strat"] = np.nan
        df_res["p_strat_adjusted"] = np.nan
        df_res["Max_SMD"] = np.nan

    decisions = []
    for _, row in df_res.iterrows():
        smd = row["Max_SMD"]
        balance = "undefined" if not np.isfinite(smd) else ("within_threshold" if smd < 0.10 else "above_threshold")
        decisions.append(decision_fields(row["q_value"], row["p_strat_adjusted"],
                                         row["analysis_status"] if n_perms > 0 else "not_run", balance))
    df_res = df_res.drop(columns=["analysis_status"], errors="ignore")
    df_res = pd.concat([df_res, pd.DataFrame(decisions, index=df_res.index, columns=["marginal_significant", "stratified_significant", "analysis_status", "balance_status"])], axis=1)

    # Sort results
    df_res = df_res.sort_values(by=["q_value", "p_value"]).reset_index(drop=True)
    
    # Save full results
    full_out = out_dir / "enrichment_results_all.tsv"
    save_results(df_res, full_out)
    print(f"    Saved full results: {full_out}")
    
    # Filter to significant
    df_sig = df_res[df_res["q_value"] < 0.05].copy()
    sig_out = out_dir / "enrichment_results_significant.tsv"
    save_results(df_sig, sig_out)
    print(f"    Saved significant results ({len(df_sig)} terms): {sig_out}")
    if n_perms > 0:
        strat_out = out_dir / "enrichment_results_stratified.tsv"
        save_results(df_res.loc[df_res["stratified_significant"].astype(bool)], strat_out)
        print(f"    Saved stratified discoveries: {strat_out}")
    
    # Presentation summaries follow the discovery list named in each file.
    if obo_ancestors:
        lists = {"marginal": df_sig}
        if n_perms > 0:
            lists["stratified"] = df_res.loc[df_res["stratified_significant"].astype(bool)]
        for label, discoveries in lists.items():
            reduced, mapping = reduce_go_redundancy(discoveries, obo_ancestors,
                term_matrix=term_matrix, return_mapping=True)
            save_results(reduced, out_dir / f"enrichment_results_nonredundant_{label}.tsv")
            save_results(mapping, out_dir / f"term_reduction_map_{label}.tsv")

    input_paths = {"master_table": args.master_table, "annotations": args.annotations}
    if args.obo:
        input_paths["obo"] = args.obo
    if Path(args.foreground).is_file():
        input_paths["foreground"] = args.foreground
    metadata = {
        "software_version": version("scbenrich"),
        "python_version": platform.python_version(),
        "dependencies": {name: version(name) for name in ("numpy", "pandas", "scipy")},
        "seed": args.seed, "permutations": args.permutations,
        "covariates": cov_cols, "min_universe": args.min_universe,
        "max_universe_fraction": args.max_universe_fraction,
        "go_version": go_version,
        "status": "no_testable_terms" if df_res.empty else "completed",
        "inputs": {name: {"path": str(Path(path).resolve()),
                           "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest()}
                   for name, path in input_paths.items()},
    }
    (out_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")

    print("[*] SCB-Enrich completed successfully!")


def run_pipeline(args):
    """Execute in a fresh output directory and preserve failure provenance."""
    from scbenrich.utils.provenance import build_manifest, write_manifest
    validate_parameters(args.permutations, getattr(args, "bins", None))
    out = Path(args.out_dir)
    if out.exists() and any(out.iterdir()):
        raise ValueError("output directory must be empty; existing files were preserved")
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = out / "run_metadata.json"
    manifest = build_manifest(args)
    write_manifest(manifest_path, manifest)
    try:
        _run_pipeline(args)
    except Exception as exc:
        manifest.update(status="failed", error={"type": type(exc).__name__, "message": str(exc)})
        write_manifest(manifest_path, manifest)
        raise
    completed = json.loads(manifest_path.read_text())
    manifest.update(completed)
    manifest['outputs'] = sorted(p.name for p in out.iterdir() if p.is_file())
    write_manifest(manifest_path, manifest)


def main():
    parser = argparse.ArgumentParser(
        prog="scb-enrich",
        description="SCB-Enrich: Stratified Covariate-Balanced Enrichment for Gene Families."
    )
    parser.add_argument("--master-table", "-m", required=True,
                        help="Path to TSV/CSV master table with gene families, counts, and covariates.")
    parser.add_argument("--annotations", "-a", required=True,
                        help="Path to TSV/CSV annotations mapping Family_ID to Term_ID.")
    parser.add_argument("--foreground", "-f", required=True,
                        help="Column name in master table OR path to text file with foreground Family IDs.")
    parser.add_argument("--covariates", "-c", default="Total_Count,Count_Ath,Annotation_Degree",
                        help="Comma-separated covariate columns to balance (default: Total_Count,Count_Ath,Annotation_Degree).")
    parser.add_argument("--obo", default=None,
                        help="Optional path to go.obo for ancestor propagation and DAG reduction.")
    parser.add_argument("--permutations", "-p", type=int, default=10000,
                        help="Number of stratified null permutations (default: 10,000; set 0 to skip).")
    parser.add_argument("--min-universe", type=int, default=5,
                        help="Minimum number of families with term in universe (default: 5).")
    parser.add_argument("--max-universe-fraction", type=float, default=0.90,
                        help="Maximum universe fraction for uninformative terms (default: 0.90).")
    parser.add_argument("--bins", type=int, default=None, help="Quantile bins per covariate (>=2).")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for permutations (default: 42).")
    parser.add_argument("--out-dir", "-o", default="scb_enrich_results",
                        help="Output directory (default: scb_enrich_results).")
                        
    args = parser.parse_args()
    run_pipeline(args)


if __name__ == "__main__":
    main()

