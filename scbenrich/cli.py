"""
Command Line Interface (CLI) for SCB-Enrich.
"""

import sys
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

from scbenrich.core.stats import fdr_bh, fisher_exact_test, exact_conditional_or
from scbenrich.core.stratified import StratifiedPermuter, compute_smd
from scbenrich.core.obo import parse_obo, propagate_go_annotations
from scbenrich.core.dag import reduce_go_redundancy
from scbenrich.utils.io import load_annotations, load_master_table, save_results


def run_pipeline(args):
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"[*] Starting SCB-Enrich analysis...")
    print(f"    Master table: {args.master_table}")
    print(f"    Annotations:  {args.annotations}")
    print(f"    Output dir:   {out_dir}")
    
    # 1. Load Master Table & Covariates
    df_master = load_master_table(args.master_table)
    all_families = list(df_master.index)
    print(f"[*] Loaded {len(all_families)} gene families in universe.")
    
    # Identify foreground
    if Path(args.foreground).is_file():
        with open(args.foreground) as f:
            foreground_set = set(line.strip() for line in f if line.strip())
    elif args.foreground in df_master.columns:
        col = df_master[args.foreground]
        if col.dtype == bool:
            foreground_set = set(df_master.index[col])
        else:
            foreground_set = set(df_master.index[col.astype(str).str.lower().isin(["1", "true", "yes", "expanded"])])
    else:
        raise ValueError(f"Foreground '{args.foreground}' is neither an existing file nor a column in master table.")
        
    foreground_set = foreground_set.intersection(all_families)
    background_set = set(all_families) - foreground_set
    print(f"    Foreground size: {len(foreground_set)} families.")
    print(f"    Background size: {len(background_set)} families.")
    
    # Parse covariates
    cov_cols = [c.strip() for c in args.covariates.split(",") if c.strip()]
    for c in cov_cols:
        if c not in df_master.columns:
            raise ValueError(f"Covariate column '{c}' not found in master table.")
    cov_df = df_master[cov_cols].copy()
    
    # 2. Load Annotations & Propagate if OBO provided
    df_ann = load_annotations(args.annotations)
    direct_map = df_ann.groupby("Family_ID")["Term_ID"].apply(set).to_dict()
    
    obo_ancestors = {}
    term_names = {}
    term_categories = {}
    
    if "Term_Name" in df_ann.columns:
        term_names = df_ann.drop_duplicates("Term_ID").set_index("Term_ID")["Term_Name"].to_dict()
    if "Category" in df_ann.columns:
        term_categories = df_ann.drop_duplicates("Term_ID").set_index("Term_ID")["Category"].to_dict()

    if args.obo:
        print(f"[*] Loading Gene Ontology OBO from: {args.obo}")
        obo_terms, _, _, obo_ancestors, _ = parse_obo(args.obo)
        for tid, tinfo in obo_terms.items():
            if tid not in term_names and "name" in tinfo:
                term_names[tid] = tinfo["name"]
            if tid not in term_categories and "namespace" in tinfo:
                term_categories[tid] = tinfo["namespace"]
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
            "Odds_Ratio": round(mle_or, 4),
            "CI_Low_95": round(ci_low, 4),
            "CI_High_95": round(ci_high, 4),
            "p_value": p_fisher,
        })
        
    df_res = pd.DataFrame(results)
    df_res["q_value"] = fdr_bh(pvals)
    
    # 5. Stratified Permutations (Null Resampling)
    n_perms = args.permutations
    if n_perms > 0:
        print(f"[*] Executing {n_perms:,} stratified null permutations balancing covariates: {cov_cols}...")
        permuter = StratifiedPermuter(cov_df, covariates=cov_cols, seed=args.seed)
        p_strat_map = permuter.run_permutations(foreground_set, term_matrix, n_permutations=n_perms)
        df_res["p_strat"] = [p_strat_map.get(t, 1.0) for t in df_res["Term_ID"]]
        
        # Calculate maximum standardized mean difference across covariates
        smds = []
        for _, row in df_res.iterrows():
            tid = row["Term_ID"]
            t_fams = set(term_matrix.index[term_matrix[tid] == 1])
            fg_t = list(t_fams.intersection(foreground_set))
            bg_t = list(t_fams.intersection(background_set))
            max_smd = 0.0
            if fg_t and bg_t:
                for col in cov_cols:
                    smd = abs(compute_smd(cov_df.loc[fg_t, col], cov_df.loc[bg_t, col]))
                    if smd > max_smd:
                        max_smd = smd
            smds.append(round(max_smd, 4))
        df_res["Max_SMD"] = smds
        
        # Evidence Classification
        evidence = []
        for _, row in df_res.iterrows():
            q = row["q_value"]
            p_s = row["p_strat"]
            smd = row["Max_SMD"]
            if q < 0.05 and p_s < 0.05 and smd < 0.10:
                evidence.append("Robust")
            elif q < 0.05 and p_s < 0.05:
                evidence.append("Suggestive")
            elif q < 0.05 and p_s >= 0.05:
                evidence.append("Size_Artifact")
            else:
                evidence.append("Not_Significant")
        df_res["Evidence_Class"] = evidence
    else:
        df_res["p_strat"] = np.nan
        df_res["Max_SMD"] = np.nan
        df_res["Evidence_Class"] = np.where(df_res["q_value"] < 0.05, "Significant", "Not_Significant")

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
    
    # 6. Redundancy Reduction
    if obo_ancestors and not df_sig.empty:
        print("[*] Reducing GO DAG semantic redundancy...")
        df_reduced = reduce_go_redundancy(df_sig, obo_ancestors, term_matrix=term_matrix)
        red_out = out_dir / "enrichment_results_nonredundant.tsv"
        save_results(df_reduced, red_out)
        print(f"    Saved non-redundant results ({len(df_reduced)} terms): {red_out}")

    print("[*] SCB-Enrich completed successfully!")


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
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for permutations (default: 42).")
    parser.add_argument("--out-dir", "-o", default="scb_enrich_results",
                        help="Output directory (default: scb_enrich_results).")
                        
    args = parser.parse_args()
    run_pipeline(args)


if __name__ == "__main__":
    main()

