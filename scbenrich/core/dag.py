"""Greedy presentation summary of overlapping GO terms."""

import pandas as pd


def reduce_go_redundancy(df_enrichment, ancestors_map, term_matrix=None,
                         max_jaccard=0.75, return_mapping=False):
    """Select representatives and optionally report every removed term."""
    if not 0 <= max_jaccard <= 1:
        raise ValueError("max_jaccard must be in [0, 1]")
    mapping_columns = ["Term_ID", "Representative", "Reason"]
    if df_enrichment.empty:
        result = (df_enrichment.copy(), pd.DataFrame(columns=mapping_columns))
        return result if return_mapping else result[0]

    priority = {"Stratified_Significant": 2, "Stratified_Balance_Concern": 1}
    df = df_enrichment.copy()
    df["_priority"] = df.get("Evidence_Class", pd.Series(index=df.index, dtype=str)).map(priority).fillna(0)
    q_col = "q_local" if "q_local" in df else "q_value"
    df = df.sort_values(["_priority", q_col, "Term_ID"], ascending=[False, True, True]).drop(columns="_priority")
    kept = []
    removed = {}
    for _, row in df.iterrows():
        term = row["Term_ID"]
        if term in removed:
            continue
        kept.append(term)
        for other in df.Term_ID:
            if other == term or other in removed or other in kept:
                continue
            if other in ancestors_map.get(term, set()):
                removed[other] = (term, "ancestor")
            elif term_matrix is not None and term in term_matrix and other in term_matrix:
                left = set(term_matrix.index[term_matrix[term] == 1])
                right = set(term_matrix.index[term_matrix[other] == 1])
                union = left | right
                if union and len(left & right) / len(union) >= max_jaccard:
                    removed[other] = (term, "jaccard")
    reduced = df[df.Term_ID.isin(kept)].reset_index(drop=True)
    mapping = pd.DataFrame([{"Term_ID": term, "Representative": rep, "Reason": reason}
                            for term, (rep, reason) in removed.items()], columns=mapping_columns)
    return (reduced, mapping) if return_mapping else reduced
