"""
Semantic redundancy reduction along the Gene Ontology directed acyclic graph (DAG).
"""

from collections import defaultdict


def reduce_go_redundancy(df_enrichment, ancestors_map, term_matrix=None, max_jaccard=0.75):
    """
    Reduces redundant ancestor-descendant and highly overlapping GO terms,
    selecting the most informative representative based on stringency class,
    statistical significance, and specificity.
    
    Parameters
    ----------
    df_enrichment : pd.DataFrame
        DataFrame containing enriched terms. Must include columns:
        'Term_ID', 'q_value' (or 'q_local'), and optionally 'Evidence_Class' / 'Specificity'.
    ancestors_map : dict
        Mapping of GO term -> set of DAG ancestor terms.
    term_matrix : pd.DataFrame, optional
        Binary family x term DataFrame to compute Jaccard family overlap.
    max_jaccard : float, optional
        Threshold for clustering terms with high family overlap (default: 0.75).
        
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame containing only representative non-redundant terms.
    """
    if df_enrichment.empty:
        return df_enrichment.copy()

    df = df_enrichment.copy()
    
    # Ranking score: prioritize stringency class, then lowest q-value, then foreground size
    evidence_priority = {
        "Robust": 3,
        "Confirmatory": 3,
        "Suggestive": 2,
        "Qualified": 2,
        "Size_Artifact": 1,
        "Not_Significant": 0,
    }
    
    q_col = "q_local" if "q_local" in df.columns else ("q_value" if "q_value" in df.columns else None)
    
    def score_row(row):
        ev = str(row.get("Evidence_Class", "Suggestive"))
        ev_score = evidence_priority.get(ev, 1)
        q = row[q_col] if q_col else 1.0
        # Lower q is better -> higher negative log
        log_q = -1.0 * (0.0 if q <= 0 else q)
        return (ev_score, log_q)

    # Sort terms by priority
    df["_rank_score"] = [score_row(r) for _, r in df.iterrows()]
    df = df.sort_values(by=["_rank_score"], ascending=False).reset_index(drop=True)
    
    kept_terms = []
    redundant_terms = set()
    
    for idx, row in df.iterrows():
        tid = row["Term_ID"]
        if tid in redundant_terms:
            continue
            
        kept_terms.append(idx)
        
        # Mark ancestors as redundant if they are less specific
        anc = ancestors_map.get(tid, set())
        redundant_terms.update(anc)
        
        # Check family Jaccard overlap if matrix provided
        if term_matrix is not None and tid in term_matrix.columns:
            fams_tid = set(term_matrix.index[term_matrix[tid] == 1])
            if fams_tid:
                for other_idx in range(idx + 1, len(df)):
                    other_tid = df.loc[other_idx, "Term_ID"]
                    if other_tid in redundant_terms or other_tid not in term_matrix.columns:
                        continue
                    fams_other = set(term_matrix.index[term_matrix[other_tid] == 1])
                    if fams_other:
                        jaccard = len(fams_tid & fams_other) / len(fams_tid | fams_other)
                        if jaccard >= max_jaccard:
                            redundant_terms.add(other_tid)

    df_reduced = df.loc[kept_terms].drop(columns=["_rank_score"]).reset_index(drop=True)
    return df_reduced

