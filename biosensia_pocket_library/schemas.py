"""Explicit versioned Arrow schemas and relational constraints."""

from __future__ import annotations

from dataclasses import dataclass, replace

import pyarrow as pa

S = pa.string()
I = pa.int64()
F = pa.float64()
B = pa.bool_()
LS = pa.list_(S)
LI = pa.list_(I)
TS = pa.timestamp("us", tz="UTC")


@dataclass(frozen=True)
class TableSpec:
    schema: pa.Schema
    primary_key: tuple[str, ...]
    foreign_keys: tuple[tuple[str, str, str], ...] = ()
    sort_by: tuple[str, ...] = ()
    allowed_enums: dict[str, frozenset[str]] | None = None
    volatile_columns: tuple[str, ...] = ()


def _schema(fields: list[tuple[str, pa.DataType]], name: str) -> pa.Schema:
    return pa.schema(fields, metadata={b"schema_name": name.encode(), b"semantic_version": b"1.0.0"})


TABLES: dict[str, TableSpec] = {}


def _add(name: str, fields, pk, fk=(), sort=()):
    TABLES[name] = TableSpec(_schema(fields, name), tuple(pk), tuple(fk), tuple(sort or pk))


_add("source_files", [
    ("source_file_id", S), ("source_kind", S), ("pdb_id", S), ("path", S), ("size_bytes", I),
    ("sha256", S), ("modified_time_utc", TS), ("download_url", S), ("downloaded_at_utc", TS),
    ("http_etag", S), ("http_last_modified", S), ("validation_status", S), ("warning_codes", LS),
], ["source_file_id"])
_add("index_record_occurrences", [("complex_id", S), ("source_index_line_number", I), ("source_line_sha256", S)],
     ["complex_id", "source_index_line_number"], sort=["complex_id", "source_index_line_number"])
_add("binding_measurements", [
    ("measurement_id", S), ("complex_id", S), ("pdb_id", S), ("measurement_type_raw", S),
    ("measurement_type_normalized", S), ("relation_raw", S), ("relation_normalized", S),
    ("value_raw", S), ("value_numeric", F), ("unit_raw", S), ("unit_normalized", S),
    ("value_molar", F), ("value_inverse_molar", F), ("p_measurement_name", S), ("p_relation", S),
    ("p_value", F), ("normalization_kind", S), ("measurement_raw", S), ("parse_status", S),
    ("parse_warning_codes", LS), ("source_index_line_number", I),
], ["measurement_id"], (("complex_id", "complexes", "complex_id"),))
_add("complexes", [(name, kind) for name, kind in [
    ("complex_id", S), ("pdb_id", S), ("distribution_id", S), ("nominal_complex_set_version", S),
    ("structure_processing_version", S), ("index_revision_date", S), ("primary_index_line_number", I),
    ("index_line_redacted", S), ("source_line_sha256", S), ("release_year", I), ("resolution_raw", S),
    ("resolution_angstrom", F), ("experimental_method_hint", S), ("ligand_label", S), ("index_comment", S),
    ("complex_directory", S), ("protein_file_id", S), ("ligand_sdf_file_id", S),
    ("ligand_mol2_file_id", S), ("pdbbind_pocket_file_id", S), ("rcsb_entry_status", S),
    ("processing_status", S), ("geometry_quality_tier", S), ("pocket_comparison_quality", S),
    ("structure_mapping_quality", S), ("bibliography_quality", S), ("warning_count", I), ("error_count", I),
]], ["complex_id"])
_add("ligand_instances", [(name, kind) for name, kind in [
    ("ligand_instance_id", S), ("complex_id", S), ("pdb_id", S), ("selected_source_format", S),
    ("selected_source_file_id", S), ("ligand_geometry_content_hash", S), ("ligand_derivation_hash", S),
    ("rdkit_parse_status", S), ("rdkit_sanitization_status", S), ("canonical_smiles", S),
    ("isomeric_smiles", S), ("inchi", S), ("inchikey", S), ("molecular_formula", S),
    ("formal_charge", I), ("molecular_weight", F), ("atom_count", I), ("heavy_atom_count", I),
    ("component_count", I), ("element_counts_json", S), ("stereochemistry_status", S),
    ("sdf_mol2_comparison_status", S), ("sdf_mol2_coordinate_rmsd", F),
    ("rcsb_ligand_match_overall_status", S), ("warnings", LS),
]], ["ligand_instance_id"], (("complex_id", "complexes", "complex_id"),))
_add("ligand_components", [("ligand_instance_id", S), ("component_index", I), ("atom_indices", LI),
    ("atom_count", I), ("heavy_atom_count", I), ("element_counts_json", S), ("formal_charge", I),
    ("centroid_x", F), ("centroid_y", F), ("centroid_z", F),
    ("minimum_other_component_separation", F), ("is_pocket_defining", B)],
    ["ligand_instance_id", "component_index"], (("ligand_instance_id", "ligand_instances", "ligand_instance_id"),))
_add("pockets", [(name, kind) for name, kind in [
    ("pocket_instance_id", S), ("complex_id", S), ("ligand_instance_id", S), ("pdb_id", S),
    ("pocket_geometry_content_hash", S), ("pocket_derivation_hash", S), ("extraction_schema_version", S),
    ("distance_cutoff_angstrom", F), ("selected_model_id", I), ("model_count", I), ("altloc_policy", S),
    ("hydrogen_policy", S), ("contact_atom_count", I), ("residue_expanded_atom_count", I),
    ("exported_atom_count", I), ("contact_residue_count", I), ("drugclip_export_view", S),
    ("contributing_chain_count", I), ("contributing_auth_chain_ids", LS),
    ("minimum_ligand_distance_min", F), ("minimum_ligand_distance_mean", F),
    ("minimum_ligand_distance_median", F), ("minimum_ligand_distance_max", F), ("crop_applied", B),
    ("crop_max_atoms", I), ("maximum_retained_ligand_distance", F), ("minimum_discarded_ligand_distance", F),
    ("all_elements_supported", B), ("processing_status", S), ("geometry_quality_tier", S),
    ("pocket_comparison_quality", S), ("structure_mapping_quality", S), ("bibliography_quality", S),
    ("warning_codes", LS), ("error_codes", LS), ("lmdb_profile_memberships", LS),
]], ["pocket_instance_id"], (("complex_id", "complexes", "complex_id"), ("ligand_instance_id", "ligand_instances", "ligand_instance_id")))
_add("pocket_residues", [("pocket_instance_id", S), ("pdb_id", S), ("model_id", I),
    ("auth_chain_id", S), ("auth_residue_number", I), ("insertion_code", S), ("residue_name", S),
    ("minimum_ligand_distance", F), ("selected_atom_count", I), ("total_heavy_atom_count", I),
    ("rcsb_mapping_status", S), ("rcsb_label_asym_id", S), ("rcsb_label_seq_id", I),
    ("rcsb_polymer_entity_id", S)],
    ["pocket_instance_id", "model_id", "auth_chain_id", "auth_residue_number", "insertion_code", "residue_name"],
    (("pocket_instance_id", "pockets", "pocket_instance_id"),))
_add("pocket_atoms", [(name, kind) for name, kind in [
    ("pocket_instance_id", S), ("pdbbind_atom_key", S), ("source_order", I), ("model_id", I),
    ("record_type", S), ("auth_chain_id", S), ("auth_residue_number", I), ("insertion_code", S),
    ("residue_name", S), ("atom_name", S), ("altloc", S), ("element", S), ("occupancy", F),
    ("b_factor", F), ("x", F), ("y", F), ("z", F), ("minimum_ligand_distance", F),
    ("in_contact_atom_view", B), ("in_residue_expanded_atom_view", B), ("retained_after_crop", B),
    ("export_order", I), ("element_supported_by_drugclip", B), ("rcsb_atom_mapping_status", S),
    ("rcsb_label_asym_id", S), ("rcsb_label_seq_id", I), ("rcsb_atom_id", S),
    ("rcsb_polymer_entity_id", S),
]], ["pocket_instance_id", "pdbbind_atom_key"], (("pocket_instance_id", "pockets", "pocket_instance_id"),),
    ["pocket_instance_id", "source_order"])

_comparison_fields = [
    ("pocket_instance_id", S), ("comparison_view", S), ("pdbbind_pocket_file_id", S), ("comparison_status", S),
    ("reextracted_atom_count", I), ("pdbbind_atom_count", I), ("reextracted_heavy_atom_count", I),
    ("pdbbind_heavy_atom_count", I), ("common_atom_exact_count", I), ("common_atom_fallback_count", I),
    ("only_reextracted_atom_count", I), ("only_pdbbind_atom_count", I), ("atom_jaccard", F),
    ("reextracted_residue_count", I), ("pdbbind_residue_count", I), ("common_residue_count", I),
    ("only_reextracted_residue_count", I), ("only_pdbbind_residue_count", I), ("residue_jaccard", F),
    ("reextracted_chain_ids", LS), ("pdbbind_chain_ids", LS), ("chain_sets_equal", B),
    ("reextracted_subset_of_pdbbind", B), ("pdbbind_subset_of_reextracted", B),
    ("common_atom_coordinate_rmsd", F), ("common_atom_max_coordinate_difference", F),
    ("reextracted_maximum_ligand_distance", F), ("pdbbind_maximum_ligand_distance", F),
    ("reextracted_mean_ligand_distance", F), ("pdbbind_mean_ligand_distance", F),
    ("reextracted_median_ligand_distance", F), ("pdbbind_median_ligand_distance", F),
    ("reextracted_p95_ligand_distance", F), ("pdbbind_p95_ligand_distance", F), ("warning_codes", LS),
]
_add("pocket_comparisons", _comparison_fields, ["pocket_instance_id", "comparison_view"],
     (("pocket_instance_id", "pockets", "pocket_instance_id"),))
_add("pocket_atom_differences", [("pocket_instance_id", S), ("comparison_view", S),
    ("comparison_class", S), ("reextracted_atom_key", S), ("pdbbind_atom_key", S), ("match_method", S),
    ("coordinate_distance", F), ("auth_chain_id", S), ("auth_residue_number", I), ("insertion_code", S),
    ("residue_name", S), ("atom_name", S), ("element", S)],
    ["pocket_instance_id", "comparison_view", "comparison_class", "reextracted_atom_key", "pdbbind_atom_key"])
_add("protein_chains", [("pocket_instance_id", S), ("pdb_id", S), ("pdbbind_auth_chain_id", S),
    ("selected_atom_count", I), ("selected_residue_count", I), ("rcsb_mapping_status", S), ("warnings", LS)],
    ["pocket_instance_id", "pdbbind_auth_chain_id"])

# Enrichment tables are normalized and always emitted, even when offline/empty.
_add("chain_mapping_candidates", [("chain_mapping_candidate_id", S), ("pocket_instance_id", S), ("pdb_id", S),
    ("pdbbind_auth_chain_id", S), ("rcsb_label_asym_id", S), ("rcsb_auth_asym_id", S),
    ("rcsb_polymer_entity_id", S), ("entity_type", S), ("entity_description", S), ("organism_name", S),
    ("taxonomy_id", I), ("sequence_length", I), ("mapping_method", S), ("mapping_status", S),
    ("candidate_rank", I), ("selected", B), ("identifier_match_count", I), ("atom_match_count", I),
    ("residue_coverage", F), ("alignment_atom_count", I), ("alignment_rmsd", F), ("transform_json", S),
    ("evidence_codes", LS), ("warning_codes", LS)], ["chain_mapping_candidate_id"])
_add("chain_uniprot_mappings", [("pocket_instance_id", S), ("pdb_id", S), ("pdbbind_auth_chain_id", S),
    ("rcsb_label_asym_id", S), ("rcsb_polymer_entity_id", S), ("uniprot_accession", S),
    ("uniprot_isoform", S), ("mapping_source", S), ("mapping_status", S), ("mapping_coverage", F)],
    ["pocket_instance_id", "pdbbind_auth_chain_id", "uniprot_accession"])
_add("chain_uniprot_mapping_segments", [("segment_id", S), ("pocket_instance_id", S),
    ("uniprot_accession", S), ("pdb_begin", I), ("pdb_end", I), ("uniprot_begin", I), ("uniprot_end", I)], ["segment_id"])
_add("rcsb_ligand_mapping_candidates", [("candidate_id", S), ("ligand_instance_id", S),
    ("rcsb_nonpolymer_entity_id", S), ("label_asym_id", S), ("auth_asym_id", S), ("auth_seq_id", I),
    ("ccd_id", S), ("bird_id", S), ("match_rmsd", F), ("match_method", S),
    ("composition_evidence", S), ("graph_evidence", S), ("rank", I), ("selected", B),
    ("ambiguity_group", S), ("status", S)], ["candidate_id"])
_add("citations", [("citation_id", S), ("title", S), ("journal", S), ("year", I), ("volume", S),
    ("issue", S), ("first_page", S), ("last_page", S), ("doi", S), ("pmid", S), ("crossref_id", S),
    ("publication_status", S), ("source_priority", I), ("bibliographic_fingerprint", S),
    ("metadata_sources", LS), ("conflict_status", S)], ["citation_id"])
_add("pdb_citation_links", [("citation_id", S), ("pdb_id", S), ("source_citation_id", S), ("role", S),
    ("source", S), ("source_priority", I), ("evidence", S)], ["citation_id", "pdb_id", "source", "role"])
_add("citation_authors", [("citation_id", S), ("ordinal", I), ("author_name", S), ("orcid", S), ("source", S)],
    ["citation_id", "source", "ordinal"])
_add("affinity_reference_links", [("measurement_id", S), ("complex_id", S), ("citation_id", S),
    ("candidate_status", S), ("confidence", F), ("evidence_sources", LS), ("evidence_note", S),
    ("automatic_or_manual", S), ("verified_by", S), ("verified_at_utc", TS)],
    ["measurement_id", "citation_id"])
_add("affinity_reference_adjudications", [("measurement_id", S), ("reference_status", S),
    ("selected_citation_id", S), ("rule_version", S), ("confidence", F), ("evidence_summary", S),
    ("adjudicator", S), ("adjudicated_at_utc", TS)], ["measurement_id"])
_add("nearby_nonprotein_components", [("pocket_instance_id", S), ("component_type", S), ("component_id", S),
    ("auth_chain_id", S), ("auth_seq_id", I), ("residue_name", S), ("insertion_code", S), ("atom_count", I),
    ("element_counts_json", S), ("minimum_ligand_distance", F), ("included_in_drugclip_tensor", B),
    ("exclusion_reason", S)], ["pocket_instance_id", "component_id"])
_add("processing_issues", [("issue_id", S), ("stage", S), ("complex_id", S), ("pocket_instance_id", S),
    ("severity", S), ("issue_code", S), ("message", S), ("exception_type", S), ("source_file_id", S),
    ("created_at_utc", TS), ("details_json", S)], ["issue_id"])
_add("lmdb_records", [("library_profile", S), ("lmdb_path", S), ("record_index", I), ("lmdb_key", S),
    ("pocket_instance_id", S), ("pocket_geometry_content_hash", S), ("pocket_derivation_hash", S),
    ("atom_count", I), ("serialized_record_sha256", S), ("logical_record_sha256", S)],
    ["library_profile", "record_index"])


def schema_for(name: str) -> pa.Schema:
    return TABLES[name].schema


ENUMS_BY_TABLE = {
    "source_files": {"validation_status": {"valid", "missing", "empty", "checksum_mismatch", "unreadable"}},
    "binding_measurements": {"parse_status": {"parsed_exact", "parsed_censored", "parsed_approximate",
        "unsupported_measurement_type", "unsupported_unit", "malformed", "missing"}},
    "complexes": {"processing_status": {"accepted", "accepted_with_warnings", "rejected", "not_processed"},
        "geometry_quality_tier": {"A", "B", "C", "rejected", "not_processed"},
        "pocket_comparison_quality": {"concordant", "moderate_difference", "severe_difference", "unavailable", "not_processed"},
        "structure_mapping_quality": {"exact", "aligned", "ambiguous", "unresolved", "unavailable", "not_processed"},
        "bibliography_quality": {"exact", "probable", "unresolved", "unavailable", "not_attempted"}},
    "pockets": {"processing_status": {"accepted", "accepted_with_warnings", "rejected", "not_processed"},
        "geometry_quality_tier": {"A", "B", "C", "rejected", "not_processed"},
        "pocket_comparison_quality": {"concordant", "moderate_difference", "severe_difference", "unavailable", "not_processed"},
        "structure_mapping_quality": {"exact", "aligned", "ambiguous", "unresolved", "unavailable", "not_processed"},
        "bibliography_quality": {"exact", "probable", "unresolved", "unavailable", "not_attempted"}},
    "processing_issues": {"severity": {"info", "warning", "error", "fatal"}},
    "affinity_reference_adjudications": {"reference_status": {"exact_affinity_reference",
        "probable_affinity_reference", "probable_structural_reference", "structural_reference_only",
        "conflicting_references", "reference_unresolved", "no_reference_available", "not_attempted"}},
    "rcsb_ligand_mapping_candidates": {"status": {"exact", "probable", "ambiguous", "unresolved", "not_attempted"}},
}
VOLATILE_BY_TABLE = {
    "source_files": ("modified_time_utc", "downloaded_at_utc"),
    "processing_issues": ("created_at_utc",),
    "affinity_reference_links": ("verified_at_utc",),
    "affinity_reference_adjudications": ("adjudicated_at_utc",),
}
FOREIGN_KEYS_BY_TABLE = {
    "complexes": (("protein_file_id", "source_files", "source_file_id"),
                  ("ligand_sdf_file_id", "source_files", "source_file_id"),
                  ("ligand_mol2_file_id", "source_files", "source_file_id"),
                  ("pdbbind_pocket_file_id", "source_files", "source_file_id")),
    "ligand_instances": (("selected_source_file_id", "source_files", "source_file_id"),),
    "index_record_occurrences": (("complex_id", "complexes", "complex_id"),),
    "protein_chains": (("pocket_instance_id", "pockets", "pocket_instance_id"),),
    "chain_mapping_candidates": (("pocket_instance_id", "pockets", "pocket_instance_id"),),
    "chain_uniprot_mappings": (("pocket_instance_id", "pockets", "pocket_instance_id"),),
    "chain_uniprot_mapping_segments": (("pocket_instance_id", "pockets", "pocket_instance_id"),),
    "rcsb_ligand_mapping_candidates": (("ligand_instance_id", "ligand_instances", "ligand_instance_id"),),
    "pdb_citation_links": (("citation_id", "citations", "citation_id"),),
    "citation_authors": (("citation_id", "citations", "citation_id"),),
    "affinity_reference_links": (("measurement_id", "binding_measurements", "measurement_id"),
                                 ("citation_id", "citations", "citation_id")),
    "affinity_reference_adjudications": (("measurement_id", "binding_measurements", "measurement_id"),
                                         ("selected_citation_id", "citations", "citation_id")),
    "nearby_nonprotein_components": (("pocket_instance_id", "pockets", "pocket_instance_id"),),
    "lmdb_records": (("pocket_instance_id", "pockets", "pocket_instance_id"),),
}
for _name, _spec in list(TABLES.items()):
    TABLES[_name] = replace(
        _spec,
        allowed_enums={column: frozenset(values) for column, values in ENUMS_BY_TABLE.get(_name, {}).items()},
        volatile_columns=VOLATILE_BY_TABLE.get(_name, ()),
        foreign_keys=tuple(dict.fromkeys(_spec.foreign_keys + FOREIGN_KEYS_BY_TABLE.get(_name, ()))),
    )
