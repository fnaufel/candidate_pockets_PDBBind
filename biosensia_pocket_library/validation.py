"""Run-level relational, redaction, checksum, and LMDB validation."""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from .config import BuildConfig
from .hashing import length_frame, normalized_array_bytes, sha256_bytes
from .lmdb_export import validate_lmdb
from .progress import progress_message, track
from .quality import load_quality_effects
from .schemas import TABLES
from .sidecars import read_sidecar, validate_sidecars


def validate_run(run_dir: Path, config: BuildConfig, *, progress: bool = True) -> list[str]:
    sidecar_dir = run_dir / "sidecars"
    progress_message(f"Validating run: {run_dir}", enabled=progress)
    progress_message("Checking sidecar schemas, keys, enums, and foreign keys", enabled=progress)
    errors = validate_sidecars(sidecar_dir, progress=progress)

    progress_message("Scanning sidecar values for prohibited legacy content", enabled=progress)
    for name, row in track(
        _iter_sidecar_rows(sidecar_dir), description="Scanning sidecar rows",
        total=_sidecar_row_count(sidecar_dir), enabled=progress, unit="row",
    ):
        for key, value in row.items():
            if "pdf" in key.lower():
                errors.append(f"Forbidden legacy PDF column: {name}.{key}")
            if isinstance(value, str) and re.search(r"(?i)(?:^|\s)\S+\.pdf(?:\s|$)", value):
                errors.append(f"Unredacted PDF-looking token: {name}.{key}")
            if isinstance(value, list) and any(isinstance(item, str) and ".pdf" in item.lower() for item in value):
                errors.append(f"Unredacted PDF-looking token in list: {name}.{key}")

    progress_message("Checking pocket geometry and core relationships", enabled=progress)
    complexes = read_sidecar(sidecar_dir, "complexes")
    ligands = read_sidecar(sidecar_dir, "ligand_instances")
    pockets = read_sidecar(sidecar_dir, "pockets")
    atoms = read_sidecar(sidecar_dir, "pocket_atoms")
    complex_ids = {row["complex_id"] for row in complexes}
    source_file_ids = {row["source_file_id"] for row in read_sidecar(sidecar_dir, "source_files")}
    for complex_row in track(
        complexes, description="Checking complex references", total=len(complexes),
        enabled=progress, unit="complex",
    ):
        for column in ("protein_file_id", "ligand_sdf_file_id", "ligand_mol2_file_id", "pdbbind_pocket_file_id"):
            if complex_row[column] is not None and complex_row[column] not in source_file_ids:
                errors.append(f"Unknown source-file reference: complexes.{column}")
    ligand_ids = {row["ligand_instance_id"] for row in ligands}
    atoms_by_pocket: dict[str, list[dict]] = {}
    for row in track(
        atoms, description="Indexing pocket atoms", total=len(atoms), enabled=progress, unit="atom",
    ):
        atoms_by_pocket.setdefault(row["pocket_instance_id"], []).append(row)
    derivations: dict[str, str] = {}
    contents: dict[str, tuple] = {}
    for pocket in track(
        pockets, description="Validating pocket geometry", total=len(pockets),
        enabled=progress, unit="pocket",
    ):
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
    if any(row["parse_status"] not in measurement_statuses for row in read_sidecar(sidecar_dir, "binding_measurements")):
        errors.append("Unknown binding-measurement parse status")
    issues_table = read_sidecar(sidecar_dir, "processing_issues")
    if any(row["severity"] not in {"info", "warning", "error", "fatal"} for row in issues_table):
        errors.append("Unknown processing-issue severity")
    pocket_ids = {row["pocket_instance_id"] for row in pockets}
    relationship_tables = (
        "protein_chains", "chain_mapping_candidates", "chain_uniprot_mappings",
        "chain_uniprot_mapping_segments", "nearby_nonprotein_components",
    )
    for table_name in track(
        relationship_tables, description="Checking pocket relationships",
        total=len(relationship_tables), enabled=progress, unit="table",
    ):
        if any(row["pocket_instance_id"] not in pocket_ids for row in read_sidecar(sidecar_dir, table_name)):
            errors.append(f"Unknown pocket reference in {table_name}")
    candidates = read_sidecar(sidecar_dir, "chain_mapping_candidates")
    selected_groups = {}
    for row in track(
        candidates, description="Checking selected chain mappings", total=len(candidates),
        enabled=progress, unit="mapping",
    ):
        if row["selected"]:
            group = (row["pocket_instance_id"], row["pdbbind_auth_chain_id"])
            selected_groups[group] = selected_groups.get(group, 0) + 1
    if any(count > 1 for count in selected_groups.values()):
        errors.append("Multiple selected chain mappings in one ambiguity group")
    issue_codes = {row["issue_code"] for row in issues_table}
    unknown_codes = issue_codes - set(load_quality_effects(config.quality.rules_file))
    if unknown_codes:
        errors.append(f"Unknown issue codes: {sorted(unknown_codes)}")
    progress_message("Checking affinity and citation relationships", enabled=progress)
    adjudications = read_sidecar(sidecar_dir, "affinity_reference_adjudications")
    measurements = read_sidecar(sidecar_dir, "binding_measurements")
    if sorted(row["measurement_id"] for row in adjudications) != sorted(row["measurement_id"] for row in measurements):
        errors.append("Each binding measurement must have exactly one reference adjudication")
    allowed_references = {"exact_affinity_reference", "probable_affinity_reference", "probable_structural_reference",
                          "structural_reference_only", "conflicting_references", "reference_unresolved",
                          "no_reference_available", "not_attempted"}
    if any(row["reference_status"] not in allowed_references for row in adjudications):
        errors.append("Unknown final affinity-reference status")
    measurement_ids = {row["measurement_id"] for row in measurements}
    citation_ids = {row["citation_id"] for row in read_sidecar(sidecar_dir, "citations")}
    reference_links = read_sidecar(sidecar_dir, "affinity_reference_links")
    for row in track(
        reference_links, description="Checking affinity reference links", total=len(reference_links),
        enabled=progress, unit="link",
    ):
        if row["measurement_id"] not in measurement_ids or row["citation_id"] not in citation_ids:
            errors.append("Affinity-reference link has unknown measurement or citation")
    lmdb_rows = read_sidecar(sidecar_dir, "lmdb_records")
    progress_message("Validating LMDB profiles", enabled=progress)
    for profile in sorted({row["library_profile"] for row in lmdb_rows}):
        progress_message(f"Validating LMDB profile: {profile}", enabled=progress)
        try:
            validate_lmdb(run_dir, profile, config, progress=progress)
        except Exception as error:
            errors.append(str(error))
    progress_message("Checking manifest identity", enabled=progress)
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
    progress_message(f"Validation checks complete: {len(errors)} error(s)", enabled=progress)
    return errors


def _sidecar_row_count(directory: Path) -> int:
    return sum(
        pq.read_metadata(directory / f"{name}.parquet").num_rows
        for name in TABLES
    )


def _iter_sidecar_rows(directory: Path):
    for name in TABLES:
        for row in read_sidecar(directory, name):
            yield name, row
