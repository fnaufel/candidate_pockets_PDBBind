from pathlib import Path

from biosensia_pocket_library.rcsb import enrich_from_mmcif


def test_mmcif_mapping_citation_uniprot_and_ligand_candidates(tmp_path: Path):
    path = tmp_path / "1abc.cif"
    path.write_text("""data_1abc
loop_
_entity_poly.entity_id
_entity_poly.type
1 'polypeptide(L)'
loop_
_entity.id
_entity.pdbx_description
1 'Example protein'
loop_
_atom_site.group_PDB
_atom_site.auth_asym_id
_atom_site.label_asym_id
_atom_site.label_entity_id
ATOM A X 1
ATOM A X 1
loop_
_struct_ref.id
_struct_ref.db_name
_struct_ref.pdbx_db_accession
1 UNP P12345
loop_
_struct_ref_seq.align_id
_struct_ref_seq.ref_id
_struct_ref_seq.pdbx_strand_id
_struct_ref_seq.seq_align_beg
_struct_ref_seq.seq_align_end
_struct_ref_seq.db_align_beg
_struct_ref_seq.db_align_end
1 1 A 1 10 5 14
loop_
_pdbx_nonpoly_scheme.asym_id
_pdbx_nonpoly_scheme.entity_id
_pdbx_nonpoly_scheme.mon_id
_pdbx_nonpoly_scheme.pdb_strand_id
_pdbx_nonpoly_scheme.pdb_seq_num
L 2 LIG B 101
loop_
_citation.id
_citation.title
_citation.journal_abbrev
_citation.year
_citation.journal_volume
_citation.page_first
_citation.page_last
_citation.pdbx_database_id_DOI
_citation.pdbx_database_id_PubMed
primary 'Example structure' JTEST 2020 1 1 9 10.1/example 1234
loop_
_citation_author.citation_id
_citation_author.name
_citation_author.ordinal
primary 'Doe, J.' 1
""")
    chains = [{"pocket_instance_id": "pocket", "pdbbind_auth_chain_id": "A", "selected_atom_count": 2}]
    ligands = [({"ligand_instance_id": "ligand"}, "LIG")]
    result = enrich_from_mmcif("1abc", path, False, chains, ligands)
    assert result["chain_mapping_candidates"][0]["selected"] is True
    assert result["chain_mapping_candidates"][0]["mapping_status"] == "exact_identifier_match"
    assert result["chain_uniprot_mappings"][0]["uniprot_accession"] == "P12345"
    assert result["rcsb_ligand_mapping_candidates"][0]["selected"] is True
    assert result["citations"][0]["doi"] == "10.1/example"
    assert result["citation_authors"][0]["author_name"] == "Doe, J."
