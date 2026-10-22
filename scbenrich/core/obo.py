"""
Gene Ontology (OBO) parsing, relationship extraction, and hierarchical propagation.
"""

from collections import defaultdict


def parse_obo(obo_path):
    """
    Parses a Gene Ontology OBO file and computes full ancestor closures.
    
    Parameters
    ----------
    obo_path : str or Path
        Path to go.obo / go-basic.obo file.
        
    Returns
    -------
    tuple
        (terms_dict, alt_map, replaced_map, ancestors_map, data_version)
    """
    terms = {}
    alt_map = {}
    replaced_map = {}
    parents_map = defaultdict(set)
    children_map = defaultdict(set)
    
    current_term = None
    obo_version = "unknown"
    
    with open(obo_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith("data-version:"):
                obo_version = line.split("data-version:")[1].strip()
            elif line == "[Term]":
                if current_term and "id" in current_term:
                    tid = current_term["id"]
                    terms[tid] = current_term
                current_term = {"parents": set(), "is_obsolete": False}
            elif line == "[Typedef]":
                if current_term and "id" in current_term:
                    tid = current_term["id"]
                    terms[tid] = current_term
                current_term = None
            elif current_term is not None:
                if line.startswith("id: GO:"):
                    current_term["id"] = line.split("id: ")[1].strip()
                elif line.startswith("name: "):
                    current_term["name"] = line.split("name: ")[1].strip()
                elif line.startswith("namespace: "):
                    ns = line.split("namespace: ")[1].strip()
                    if ns == "biological_process":
                        current_term["namespace"] = "BP"
                    elif ns == "molecular_function":
                        current_term["namespace"] = "MF"
                    elif ns == "cellular_component":
                        current_term["namespace"] = "CC"
                    else:
                        current_term["namespace"] = ns
                elif line.startswith("alt_id: GO:"):
                    alt = line.split("alt_id: ")[1].strip()
                    alt_map[alt] = current_term["id"]
                elif line.startswith("is_obsolete: true"):
                    current_term["is_obsolete"] = True
                elif line.startswith("replaced_by: GO:"):
                    rep = line.split("replaced_by: ")[1].strip()
                    replaced_map[current_term["id"]] = rep
                elif line.startswith("is_a: GO:"):
                    p = line.split("is_a: ")[1].split()[0].strip()
                    parents_map[current_term["id"]].add(p)
                    children_map[p].add(current_term["id"])
                elif line.startswith("relationship: part_of GO:"):
                    p = line.split("relationship: part_of ")[1].split()[0].strip()
                    parents_map[current_term["id"]].add(p)
                    children_map[p].add(current_term["id"])
                    
        if current_term and "id" in current_term:
            tid = current_term["id"]
            terms[tid] = current_term
            
    # Compute recursive ancestor map
    ancestors_map = {}
    
    def get_ancestors(t):
        if t in ancestors_map:
            return ancestors_map[t]
        anc = set()
        for p in parents_map.get(t, []):
            anc.add(p)
            anc.update(get_ancestors(p))
        ancestors_map[t] = anc
        return anc

    for t in terms:
        get_ancestors(t)

    return terms, alt_map, replaced_map, ancestors_map, obo_version


def propagate_go_annotations(direct_annotations, ancestors_map, root_terms=None):
    """
    Propagates direct GO term annotations to all their DAG ancestors.
    
    Parameters
    ----------
    direct_annotations : dict
        Mapping of family_id -> set of direct GO IDs.
    ancestors_map : dict
        Mapping of GO ID -> set of ancestor GO IDs.
    root_terms : set, optional
        Set of uninformative root terms to exclude (e.g. GO:0008150, GO:0003674, GO:0005575).
        
    Returns
    -------
    dict
        Mapping of family_id -> set of propagated GO IDs.
    """
    if root_terms is None:
        root_terms = {"GO:0008150", "GO:0003674", "GO:0005575"}
        
    propagated = {}
    for fam, terms in direct_annotations.items():
        all_terms = set(terms)
        for t in terms:
            if t in ancestors_map:
                all_terms.update(ancestors_map[t])
        all_terms.difference_update(root_terms)
        propagated[fam] = all_terms
        
    return propagated

