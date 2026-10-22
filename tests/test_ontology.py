import pytest
from scbenrich.core.obo import normalize_go_id, parse_obo


def test_chain_and_ambiguity():
    terms = {'GO:1': {'is_obsolete': True}, 'GO:2': {'is_obsolete': True}, 'GO:3': {}}
    assert normalize_go_id('GO:alt', terms, {'GO:alt': 'GO:1'}, {'GO:1': {'GO:2'}, 'GO:2': {'GO:3'}}) == 'GO:3'
    with pytest.raises(ValueError, match='ambiguous'):
        normalize_go_id('GO:1', terms, {}, {'GO:1': {'GO:2', 'GO:3'}})
    with pytest.raises(ValueError, match='cycle'):
        normalize_go_id('GO:1', terms, {}, {'GO:1': {'GO:2'}, 'GO:2': {'GO:1'}})
    with pytest.raises(ValueError, match='unknown'):
        normalize_go_id('GO:missing', terms, {}, {})


def test_parser_detects_cycle(tmp_path):
    path = tmp_path / 'cycle.obo'
    path.write_text('[Term]\nid: GO:1\nis_a: GO:2\n\n[Term]\nid: GO:2\nis_a: GO:1\n')
    with pytest.raises(ValueError, match='cycle'):
        parse_obo(path)


def test_parser_preserves_multiple_replacements_and_relation_semantics(tmp_path):
    path = tmp_path / 'terms.obo'
    path.write_text('[Term]\nid: GO:1\nreplaced_by: GO:2\nreplaced_by: GO:3\nrelationship: regulates GO:2\nrelationship: has_part GO:3\n\n[Term]\nid: GO:2\n\n[Term]\nid: GO:3\nrelationship: part_of GO:2\n')
    _, _, replacements, ancestors, _ = parse_obo(path)
    assert replacements['GO:1'] == {'GO:2', 'GO:3'}
    assert ancestors['GO:1'] == set()
    assert ancestors['GO:3'] == {'GO:2'}
