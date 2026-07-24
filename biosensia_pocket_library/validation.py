"""Run-level relational, redaction, checksum, and LMDB validation."""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

from .config import BuildConfig
from .hashing import length_frame, normalized_array_bytes, sha256_bytes
from .lmdb_export import validate_lmdb
from .quality import load_quality_effects
from .schemas import TABLES
from .sidecars import read_sidecar, validate_sidecars


def validate_run(run_dir: Path, config: BuildConfig, *, progress: bool = True) -> list[str]:
    errors = validate_sidecars(run_dir / "sidecars")
    for name in TABLES:
        for row in read_sidecar(run_dir / "sidecars", name):
            for key, value in row.items():
                if "pdf" in key.lower():
                    errors.append(f"Forbidden legacy PDF column: {name}.{key}")
                if isinstance(value, str) and re.search(r"(?i)(?:^|\s)\S+\.pdf(?:\s|$)", value):
                    errors.append(f"Unredacted PDF-looking token: {name}.{key}")
                if isinstance(value, list) and any(isinstance(item, str) and ".pdf" in item.lower() for item in value):
                    errors.append(f"Unredacted PDF-looking token in list: {name}.{key}")
    complexes = read_sidecar(run_dir / "sidecars", "complexes")
    ligands = read_sidecar(run_dir / "sidecars", "ligand_instances")
    pockets = read_sidecar(run_dir / "sidecars", "pockets")
    atoms = read_sidecar(run_dir / "sidecars", "pocket_atoms")
    complex_ids = {row["complex_id"] for row in complexes}
    source_file_ids = {row["source_file_id"] for row in read_sidecar(run_dir / "sidecars", "source_files")}
    for complex_row in complexes:
        for column in ("protein_file_id", "ligand_sdf_file_id", "ligand_mol2_file_id", "pdbbind_pocket_file_id"):
            if complex_row[column] is not None and complex_row[column] not in source_file_ids:
                errors.append(f"Unknown source-file reference: complexes.{column}")
    ligand_ids = {row["ligand_instance_id"] for row in ligands}
    atoms_by_pocket: dict[str, list[dict]] = {}
    for row in atoms:
        atoms_by_pocket.setdefault(row["pocket_instance_id"], []).append(row)
    derivations: dict[str, str] = {}
    contents: dict[str, tuple] = {}
    for pocket in pockets:
        if pocket["complex_id"] not in complex_ids or pocket["ligand_instance_id"] not in ligand_ids:
            errors.append(f"Accepted pocket has unknown complex or ligand: {pocket['pocket_instance_id']}")
        if pocket["processing_status"].startswith("accepted") and not atoms_by_pocket.get(pocket["pocket_instance_id"]):
            errors.append(f"Accepted pocket has no atom rows: {pocket['pocket_instance_id']}")
        prior = derivations.setdefault(pocket["pocket_derivation_hash"], pocket["pocket_geometry_content_hash"])
        if prior != pocket["pocket_geometry_content_hash"]:
            errors.append("A derivation hash maps to multiple content hashes")
        source_representation = pocket.get("drugclip_export_view") == "source_pickle"
        exported = sorted((row for row in atoms_by_pocket.get(pocket["pocket_instance_id"], [])
                           if (row.get("included_in_lmdb_source") if source_representation
                               else row["retained_after_crop"])),
                          key=lambda row: row["source_order"] if source_representation else row["export_order"])
        tokens = [row["element"] for row in exported]
        coordinates = np.asarray([[row["x"], row["y"], row["z"]] for row in exported], dtype=np.float32)
        content_schema = b"pocket-content-v2-source-pickle" if source_representation else b"pocket-content-v1"
        regenerated = sha256_bytes(length_frame((content_schema, "\0".join(tokens).encode(),
                                                  normalized_array_bytes(coordinates, "<f4"))))
        if regenerated != pocket["pocket_geometry_content_hash"]:
            errors.append(f"Pocket content hash mismatch: {pocket['pocket_instance_id']}")
        representation = (tuple(tokens), normalized_array_bytes(coordinates, "<f4"))
        previous = contents.setdefault(pocket["pocket_geometry_content_hash"], representation)
        if previous != representation:
            errors.append("Duplicate pocket content hash has inconsistent atoms or coordinates")
    allowed_status = {"accepted", "accepted_with_warnings", "rejected", "not_processed"}
    if any(row["processing_status"] not in allowed_status for row in complexes + pockets):
        errors.append("Unknown geometry processing status")
    enums = {
        "geometry_quality_tier": {"A", "B", "C", "rejected", "not_processed"},
        "pocket_comparison_quality": {"concordant", "moderate_difference", "severe_difference", "unavailable", "not_processed"},
        "structure_mapping_quality": {"exact", "aligned", "ambiguous", "unresolved", "unavailable", "not_processed"},
        "bibliography_quality": {"exact", "probable", "unresolved", "unavailable", "not_attempted"},
    }
    for column, allowed in enums.items():
        if any(row[column] not in allowed for row in complexes + pockets):
            errors.append(f"Unknown {column} value")
    measurement_statuses = {"parsed_exact", "parsed_censored", "parsed_approximate",
                            "unsupported_measurement_type", "unsupported_unit", "malformed", "missing"}
    if any(row["parse_status"] not in measurement_statuses for row in read_sidecar(run_dir / "sidecars", "binding_measurements")):
        errors.append("Unknown binding-measurement parse status")
    issues_table = read_sidecar(run_dir / "sidecars", "processing_issues")
    if any(row["severity"] not in {"info", "warning", "error", "fatal"} for row in issues_table):
        errors.append("Unknown processing-issue severity")
    pocket_ids = {row["pocket_instance_id"] for row in pockets}
    for table_name in ("protein_chains", "chain_mapping_candidates", "chain_uniprot_mappings",
                       "chain_uniprot_mapping_segments", "nearby_nonprotein_components"):
        if any(row["pocket_instance_id"] not in pocket_ids for row in read_sidecar(run_dir / "sidecars", table_name)):
            errors.append(f"Unknown pocket reference in {table_name}")
    candidates = read_sidecar(run_dir / "sidecars", "chain_mapping_candidates")
    selected_groups = {}
    for row in candidates:
        if row["selected"]:
            group = (row["pocket_instance_id"], row["pdbbind_auth_chain_id"])
            selected_groups[group] = selected_groups.get(group, 0) + 1
    if any(count > 1 for count in selected_groups.values()):
        errors.append("Multiple selected chain mappings in one ambiguity group")
    issue_codes = {row["issue_code"] for row in issues_table}
    unknown_codes = issue_codes - set(load_quality_effects(config.quality.rules_file))
    if unknown_codes:
        errors.append(f"Unknown issue codes: {sorted(unknown_codes)}")
    adjudications = read_sidecar(run_dir / "sidecars", "affinity_reference_adjudications")
    measurements = read_sidecar(run_dir / "sidecars", "binding_measurements")
    if sorted(row["measurement_id"] for row in adjudications) != sorted(row["measurement_id"] for row in measurements):
        errors.append("Each binding measurement must have exactly one reference adjudication")
    allowed_references = {"exact_affinity_reference", "probable_affinity_reference", "probable_structural_reference",
                          "structural_reference_only", "conflicting_references", "reference_unresolved",
                          "no_reference_available", "not_attempted"}
    if any(row["reference_status"] not in allowed_references for row in adjudications):
        errors.append("Unknown final affinity-reference status")
    measurement_ids = {row["measurement_id"] for row in measurements}
    citation_ids = {row["citation_id"] for row in read_sidecar(run_dir / "sidecars", "citations")}
    for row in read_sidecar(run_dir / "sidecars", "affinity_reference_links"):
        if row["measurement_id"] not in measurement_ids or row["citation_id"] not in citation_ids:
            errors.append("Affinity-reference link has unknown measurement or citation")
    lmdb_rows = read_sidecar(run_dir / "sidecars", "lmdb_records")
    for profile in sorted({row["library_profile"] for row in lmdb_rows}):
        try:
            validate_lmdb(run_dir, profile, config, progress=progress)
        except Exception as error:
            errors.append(str(error))
    manifest_path = run_dir / "manifest.json"
    if not manifest_path.is_file():
        errors.append("Missing manifest.json")
    else:
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            if manifest.get("run_id") != run_dir.name:
                errors.append("Run-directory name differs from manifest run_id")
        except ValueError:
            errors.append("Malformed manifest.json")
    return errors
