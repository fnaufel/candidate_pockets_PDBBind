"""Source adapter for DrugCLIP's trusted ``pdb/combine_set`` bundles."""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from .combine_set_mapping import AtomMapping, map_pickle_atoms
from .config import BuildConfig
from .hashing import canonical_json_hash, length_frame, normalized_array_bytes, sha256_bytes, sha256_file, stable_id
from .ligand_parser import ligand_component_rows, parse_ligand
from .models import ProcessingIssue
from .progress import file_progress, track
from .trusted_pickle import CombineSetPickle, drugclip_pocket_token, load_trusted_combine_set_pickle


@dataclass(frozen=True, slots=True)
class CombineSetRecord:
    pdb_id: str
    directory: Path
    pickle_path: Path


def discover_combine_set(root: Path) -> list[CombineSetRecord]:
    if not root.is_dir():
        raise FileNotFoundError(f"DrugCLIP combine_set directory not found: {root}")
    records = [
        CombineSetRecord(path.parent.name.lower(), path.parent.resolve(), path.resolve())
        for path in root.glob("*/data.pkl") if path.is_file()
    ]
    records.sort(key=lambda item: (item.pdb_id, item.directory.as_posix()))
    identifier_counts = Counter(item.pdb_id for item in records)
    duplicates = {pdb_id for pdb_id, count in identifier_counts.items() if count > 1}
    if duplicates:
        raise ValueError(f"Duplicate combine_set directory identifiers: {sorted(duplicates)}")
    return records


def select_combine_set(
    records: Iterable[CombineSetRecord], pdb_ids: Iterable[str] | None = None, limit: int | None = None,
) -> list[CombineSetRecord]:
    wanted = {value.lower() for value in pdb_ids} if pdb_ids else None
    selected = [item for item in records if wanted is None or item.pdb_id in wanted]
    if wanted:
        missing = wanted - {item.pdb_id for item in selected}
        if missing:
            raise ValueError(f"Requested combine_set identifiers not found: {sorted(missing)}")
    return selected[:limit] if limit is not None else selected


def inventory_combine_set(
    records: Iterable[CombineSetRecord], config: BuildConfig, *, progress: bool = True,
) -> tuple[list[dict], dict[str, dict], list[ProcessingIssue]]:
    rows: list[dict] = []
    bundles: dict[str, dict] = {}
    issues: list[ProcessingIssue] = []
    root = config.paths.combine_set_root
    assert root is not None
    for record in track(list(records), description="Inventorying combine_set", enabled=progress):
        actual = sorted(path for path in record.directory.iterdir() if path.is_file())
        files: dict[str, str | Path] = {}
        for path in actual:
            role = _source_role(record.pdb_id, path.name)
            with file_progress(path, description=f"SHA-256 {path.name}",
                               enabled=progress and path.stat().st_size > 50_000_000) as bar:
                digest = sha256_file(path, progress=bar)
            relative = _relative_path(path, config)
            source_id = stable_id("file", f"drugclip_combine_set_{role}", relative, digest)
            rows.append({
                "source_file_id": source_id, "source_kind": f"drugclip_combine_set_{role}",
                "pdb_id": record.pdb_id, "path": relative, "size_bytes": path.stat().st_size,
                "sha256": digest, "modified_time_utc": None, "download_url": None,
                "downloaded_at_utc": None, "http_etag": None, "http_last_modified": None,
                "validation_status": "valid" if path.stat().st_size else "empty", "warning_codes": [],
            })
            files[role] = source_id
            files[f"{role}_path"] = path
            if not path.stat().st_size:
                issues.append(ProcessingIssue("inventory", "EMPTY_SOURCE_FILE", "error",
                                              f"Empty source file {path.name}", source_file_id=source_id))
        if "pickle" not in files:
            raise ValueError(f"Inventory lost required pickle for {record.pdb_id}")
        for role in ("ligand_sdf", "ligand_mol2", "protein_pdb", "pocket_pdb", "pocket6a_pdb"):
            if role not in files:
                issues.append(ProcessingIssue("inventory", "MISSING_COMBINE_SET_NEIGHBOR", "warning",
                                              f"Missing {role} for {record.pdb_id}"))
        bundles[record.pdb_id] = {"record": record, "files": files}
    return rows, bundles, issues


def process_combine_set_record(
    record: CombineSetRecord, bundle: dict, config: BuildConfig, dictionary: set[str],
) -> tuple[dict, dict[str, list[dict]], list[ProcessingIssue], dict]:
    rows: dict[str, list[dict]] = defaultdict(list)
    issues: list[ProcessingIssue] = []
    files = bundle["files"]
    pickle_id = str(files["pickle"])
    pickle_path = Path(files["pickle_path"])
    pickle_sha = sha256_file(pickle_path)
    complex_id = stable_id("complex", config.combine_set.distribution_id, record.pdb_id, pickle_sha)
    complex_row = _complex_row(record, complex_id, files, config)
    pocket_id = None
    try:
        source = load_trusted_combine_set_pickle(
            pickle_path, trusted=config.combine_set.trusted_pickles, dictionary=dictionary
        )
        if source.pocket_name.lower() != record.pdb_id:
            raise ValueError(
                f"Pickle pocket {source.pocket_name!r} differs from directory {record.pdb_id!r}"
            )
        ligand = parse_ligand(
            complex_id, _optional_path(files, "ligand_sdf_path"),
            _optional_path(files, "ligand_mol2_path"), config,
        )
        candidates = [
            (role.removesuffix("_pdb"), Path(files[f"{role}_path"]), str(files[role]))
            for role in ("pocket_pdb", "pocket6a_pdb") if f"{role}_path" in files
        ]
        mapping = map_pickle_atoms(
            source.pocket_atoms, source.pocket_coordinates, candidates, config,
            pdb_id=record.pdb_id, complex_id=complex_id,
        )
        content_hash = _pocket_content_hash(source)
        derivation = canonical_json_hash({
            "schema": "drugclip-combine-set-derivation-v1", "source_pickle_sha256": pickle_sha,
            "pocket_name": source.pocket_name, "preserve_source_geometry": True,
            "content_hash": content_hash,
        })
        pocket_id = stable_id("pocket", complex_id, derivation)
        distances = _ligand_distances(source.pocket_coordinates, ligand)
        mapping_quality = _mapping_quality(mapping, len(source.pocket_atoms))
        comparison_quality = _comparison_quality(mapping, len(source.pocket_atoms), config)
        rows["ligand_instances"].append(_ligand_row(record.pdb_id, complex_id, ligand, files))
        rows["ligand_components"].extend(ligand_component_rows(ligand))
        rows["pockets"].append(_pocket_row(
            record.pdb_id, complex_id, pocket_id, ligand.ligand_instance_id, source, mapping,
            content_hash, derivation, distances, comparison_quality, mapping_quality, config, pickle_id,
        ))
        rows["pocket_atoms"].extend(_atom_rows(
            pocket_id, source, mapping, distances, dictionary, pickle_id,
            config.pocket.distance_cutoff_angstrom,
        ))
        rows["pocket_residues"].extend(_residue_rows(pocket_id, record.pdb_id, mapping, distances))
        rows["protein_chains"].extend(_chain_rows(pocket_id, record.pdb_id, mapping))
        rows["pocket_comparisons"].append(_comparison_row(
            pocket_id, mapping, len(source.pocket_atoms), comparison_quality
        ))
        measurement, adjudication = _measurement_rows(source, complex_id, record.pdb_id)
        rows["binding_measurements"].append(measurement)
        rows["affinity_reference_adjudications"].append(adjudication)
        if mapping.matched_count != len(source.pocket_atoms):
            issues.append(ProcessingIssue(
                "map-source-structures", "PICKLE_ATOM_MAPPING_INCOMPLETE", "warning",
                f"Mapped {mapping.matched_count}/{len(source.pocket_atoms)} pickle atoms to a neighboring pocket PDB",
                complex_id, pocket_id, pickle_id,
                details={"mapping_source_role": mapping.source_role},
            ))
        complex_row.update(
            processing_status="accepted_with_warnings" if issues or ligand.warnings else "accepted",
            geometry_quality_tier="A", pocket_comparison_quality=comparison_quality,
            structure_mapping_quality=mapping_quality,
        )
        event = {"level": "info", "stage": "load-pickles", "code": "POCKET_ACCEPTED",
                 "message": "Trusted pickle pocket accepted", "complex_id": complex_id,
                 "pocket_instance_id": pocket_id,
                 "details": {"atom_count": len(source.pocket_atoms), "mapping_quality": mapping_quality}}
    except Exception as error:
        if config.pipeline.fail_fast:
            raise
        complex_row.update(processing_status="rejected", geometry_quality_tier="rejected",
                           pocket_comparison_quality="unavailable", structure_mapping_quality="unresolved")
        issues.append(ProcessingIssue("load-pickles", "COMBINE_SET_PICKLE_FAILED", "error",
                                      str(error), complex_id, pocket_id, pickle_id,
                                      exception_type=type(error).__name__))
        event = {"level": "error", "stage": "load-pickles", "code": "COMBINE_SET_PICKLE_FAILED",
                 "message": str(error), "complex_id": complex_id}
    return complex_row, rows, issues, event


def _complex_row(record, complex_id, files, config):
    directory = _relative_path(record.directory, config)
    return {
        "complex_id": complex_id, "pdb_id": record.pdb_id,
        "distribution_id": config.combine_set.distribution_id,
        "nominal_complex_set_version": config.combine_set.nominal_release,
        "geometry_origin": "drugclip_combine_set_pickle", "geometry_source_file_id": files["pickle"],
        "structure_processing_version": None, "index_revision_date": None,
        "primary_index_line_number": None, "index_line_redacted": None, "source_line_sha256": None,
        "release_year": None, "resolution_raw": None, "resolution_angstrom": None,
        "experimental_method_hint": None, "ligand_label": None, "index_comment": None,
        "complex_directory": directory, "protein_file_id": files.get("protein_pdb"),
        "ligand_sdf_file_id": files.get("ligand_sdf"), "ligand_mol2_file_id": files.get("ligand_mol2"),
        "pdbbind_pocket_file_id": files.get("pocket_pdb"), "rcsb_entry_status": "not_attempted",
        "processing_status": "not_processed", "geometry_quality_tier": "not_processed",
        "pocket_comparison_quality": "not_processed", "structure_mapping_quality": "not_processed",
        "bibliography_quality": "not_attempted", "warning_count": 0, "error_count": 0,
    }


def _ligand_row(pdb_id, complex_id, ligand, files):
    chemistry = {key: value for key, value in ligand.chemistry.items() if not key.startswith("_")}
    return {
        "ligand_instance_id": ligand.ligand_instance_id, "complex_id": complex_id, "pdb_id": pdb_id,
        "selected_source_format": ligand.source_format,
        "selected_source_file_id": files.get(f"ligand_{ligand.source_format}"),
        "ligand_geometry_content_hash": ligand.content_hash, "ligand_derivation_hash": ligand.derivation_hash,
        "rdkit_parse_status": "parsed", "rdkit_sanitization_status": ligand.chemistry["rdkit_sanitization_status"],
        **chemistry, "sdf_mol2_comparison_status": ligand.comparison["status"],
        "sdf_mol2_coordinate_rmsd": ligand.comparison["coordinate_rmsd"],
        "rcsb_ligand_match_overall_status": "not_processed", "warnings": ligand.warnings,
    }


def _pocket_row(pdb_id, complex_id, pocket_id, ligand_id, source, mapping, content_hash,
                derivation, distances, comparison_quality, mapping_quality, config, pickle_id):
    mapped = [atom for atom in mapping.atoms if atom is not None]
    chains = sorted({atom.auth_chain_id for atom in mapped})
    heavy_count = sum(drugclip_pocket_token(atom) != "H" for atom in source.pocket_atoms)
    contact_mask = distances <= config.pocket.distance_cutoff_angstrom
    return {
        "pocket_instance_id": pocket_id, "complex_id": complex_id, "ligand_instance_id": ligand_id,
        "pdb_id": pdb_id, "pocket_geometry_content_hash": content_hash,
        "pocket_derivation_hash": derivation, "extraction_schema_version": "combine-set-1",
        "geometry_origin": "drugclip_combine_set_pickle", "geometry_source_file_id": pickle_id,
        "derivation_method": "trusted_pickle_source_geometry", "source_geometry_atom_count": len(source.pocket_atoms),
        "source_geometry_heavy_atom_count": heavy_count,
        "distance_cutoff_angstrom": config.pocket.distance_cutoff_angstrom,
        "selected_model_id": mapped[0].model_id if mapped else None, "model_count": None,
        "altloc_policy": "source_pickle_order", "hydrogen_policy": "preserve_source",
        "contact_atom_count": int(np.sum(contact_mask)), "residue_expanded_atom_count": len(source.pocket_atoms),
        "exported_atom_count": len(source.pocket_atoms),
        "contact_residue_count": len({atom.residue_key for index, atom in enumerate(mapping.atoms)
                                      if atom is not None and contact_mask[index]}),
        "drugclip_export_view": "source_pickle", "contributing_chain_count": len(chains),
        "contributing_auth_chain_ids": chains, "minimum_ligand_distance_min": float(distances.min()),
        "minimum_ligand_distance_mean": float(distances.mean()),
        "minimum_ligand_distance_median": float(np.median(distances)),
        "minimum_ligand_distance_max": float(distances.max()), "crop_applied": False,
        "crop_max_atoms": config.pocket.max_pocket_atoms, "maximum_retained_ligand_distance": None,
        "minimum_discarded_ligand_distance": None, "all_elements_supported": True,
        "processing_status": "accepted", "geometry_quality_tier": "A",
        "pocket_comparison_quality": comparison_quality, "structure_mapping_quality": mapping_quality,
        "bibliography_quality": "not_attempted", "warning_codes": [], "error_codes": [],
        "lmdb_profile_memberships": ["default", "tier-a", "tiers-ab", "all-usable"],
    }


def _atom_rows(pocket_id, source, mapping, distances, dictionary, pickle_id, distance_cutoff):
    rows = []
    for index, (raw_token, coordinate, mapped, status, distance) in enumerate(zip(
        source.pocket_atoms, source.pocket_coordinates, mapping.atoms, mapping.statuses, distances, strict=True
    )):
        token = drugclip_pocket_token(raw_token)
        atom_key = stable_id("combine-set-atom", pocket_id, index, raw_token, coordinate.tolist())
        rows.append({
            "pocket_instance_id": pocket_id, "pdbbind_atom_key": atom_key, "source_atom_key": atom_key,
            "geometry_source_file_id": pickle_id, "source_order": index, "model_id": mapped.model_id if mapped else None,
            "record_type": mapped.record_type if mapped else None,
            "auth_chain_id": mapped.auth_chain_id if mapped else None,
            "auth_residue_number": mapped.auth_residue_number if mapped else None,
            "insertion_code": mapped.insertion_code if mapped else None,
            "residue_name": mapped.residue_name if mapped else None, "atom_name": mapped.atom_name if mapped else None,
            "altloc": mapped.altloc if mapped else None, "element": raw_token,
            "occupancy": mapped.occupancy if mapped else None, "b_factor": mapped.b_factor if mapped else None,
            "x": float(coordinate[0]), "y": float(coordinate[1]), "z": float(coordinate[2]),
            "minimum_ligand_distance": float(distance),
            "in_contact_atom_view": bool(distance <= distance_cutoff),
            "in_residue_expanded_atom_view": True, "retained_after_crop": True, "export_order": index,
            "element_supported_by_drugclip": token in dictionary, "rcsb_atom_mapping_status": "not_processed",
            "rcsb_label_asym_id": None, "rcsb_label_seq_id": None, "rcsb_atom_id": None,
            "rcsb_polymer_entity_id": None, "source_mapping_status": status, "included_in_lmdb_source": True,
        })
    return rows


def _residue_rows(pocket_id, pdb_id, mapping, distances):
    grouped: dict[tuple, list[int]] = defaultdict(list)
    for index, atom in enumerate(mapping.atoms):
        if atom is not None:
            grouped[(atom.model_id, *atom.residue_key)].append(index)
    return [{
        "pocket_instance_id": pocket_id, "pdb_id": pdb_id, "model_id": key[0], "auth_chain_id": key[1],
        "auth_residue_number": key[2], "insertion_code": key[3], "residue_name": key[4],
        "minimum_ligand_distance": float(min(distances[index] for index in indices)),
        "selected_atom_count": len(indices),
        "total_heavy_atom_count": sum(mapping.atoms[index].element != "H" for index in indices),
        "rcsb_mapping_status": "not_processed", "rcsb_label_asym_id": None, "rcsb_label_seq_id": None,
        "rcsb_polymer_entity_id": None,
    } for key, indices in sorted(grouped.items())]


def _chain_rows(pocket_id, pdb_id, mapping):
    grouped: dict[str, list] = defaultdict(list)
    for atom in mapping.atoms:
        if atom is not None:
            grouped[atom.auth_chain_id].append(atom)
    return [{
        "pocket_instance_id": pocket_id, "pdb_id": pdb_id, "pdbbind_auth_chain_id": chain,
        "selected_atom_count": len(atoms), "selected_residue_count": len({atom.residue_key for atom in atoms}),
        "rcsb_mapping_status": "not_processed", "warnings": [],
    } for chain, atoms in sorted(grouped.items())]


def _comparison_row(pocket_id, mapping, atom_count, quality):
    matched = mapping.matched_count
    row = {
        "pocket_instance_id": pocket_id, "comparison_view": mapping.source_role or "neighbor_pocket",
        "pdbbind_pocket_file_id": mapping.source_file_id,
        "comparison_status": "compared" if mapping.source_file_id else "unavailable",
        "left_geometry_role": "drugclip_combine_set_pickle",
        "right_geometry_role": mapping.source_role or "neighbor_pocket_unavailable",
        "reextracted_atom_count": atom_count, "pdbbind_atom_count": matched,
        "reextracted_heavy_atom_count": None, "pdbbind_heavy_atom_count": None,
        "common_atom_exact_count": sum(status == "exact_ordered" for status in mapping.statuses),
        "common_atom_fallback_count": sum(status == "spatial_unique" for status in mapping.statuses),
        "only_reextracted_atom_count": atom_count - matched, "only_pdbbind_atom_count": 0,
        "atom_jaccard": matched / atom_count, "reextracted_residue_count": None,
        "pdbbind_residue_count": None, "common_residue_count": None,
        "only_reextracted_residue_count": None, "only_pdbbind_residue_count": None,
        "residue_jaccard": None, "reextracted_chain_ids": [],
        "pdbbind_chain_ids": sorted({atom.auth_chain_id for atom in mapping.atoms if atom is not None}),
        "chain_sets_equal": None, "reextracted_subset_of_pdbbind": matched == atom_count,
        "pdbbind_subset_of_reextracted": True, "common_atom_coordinate_rmsd": None,
        "common_atom_max_coordinate_difference": mapping.maximum_distance,
        "reextracted_maximum_ligand_distance": None, "pdbbind_maximum_ligand_distance": None,
        "reextracted_mean_ligand_distance": None, "pdbbind_mean_ligand_distance": None,
        "reextracted_median_ligand_distance": None, "pdbbind_median_ligand_distance": None,
        "reextracted_p95_ligand_distance": None, "pdbbind_p95_ligand_distance": None,
        "warning_codes": [] if quality == "concordant" else ["PICKLE_ATOM_MAPPING_INCOMPLETE"],
    }
    return row


def _measurement_rows(source, complex_id, pdb_id):
    label = source.label
    kind = label[0] if isinstance(label, tuple) and len(label) == 2 else None
    try:
        value = float(label[1]) if isinstance(label, tuple) and len(label) == 2 else None
    except (TypeError, ValueError):
        value = None
    measurement_id = stable_id("measurement", complex_id, kind, value)
    row = {
        "measurement_id": measurement_id, "complex_id": complex_id, "pdb_id": pdb_id,
        "measurement_type_raw": str(kind) if kind is not None else "missing",
        "measurement_type_normalized": str(kind) if kind is not None else None,
        "relation_raw": "=" if value is not None else None, "relation_normalized": "=" if value is not None else None,
        "value_raw": str(value) if value is not None else None, "value_numeric": value,
        "unit_raw": None, "unit_normalized": None, "value_molar": None, "value_inverse_molar": None,
        "p_measurement_name": str(kind) if value is not None else None, "p_relation": "=" if value is not None else None,
        "p_value": value, "normalization_kind": "drugclip_source_label" if value is not None else None,
        "measurement_raw": json.dumps(label, default=str), "parse_status": "parsed_exact" if value is not None else "missing",
        "parse_warning_codes": [], "source_index_line_number": None,
    }
    adjudication = {
        "measurement_id": measurement_id, "reference_status": "not_attempted", "selected_citation_id": None,
        "rule_version": "1", "confidence": None,
        "evidence_summary": "DrugCLIP source label retained without bibliographic inference",
        "adjudicator": "automatic-v1", "adjudicated_at_utc": None,
    }
    return row, adjudication


def _pocket_content_hash(source: CombineSetPickle) -> str:
    return sha256_bytes(length_frame((
        b"pocket-content-v2-source-pickle", "\0".join(source.pocket_atoms).encode("utf-8"),
        normalized_array_bytes(source.pocket_coordinates, "<f4"),
    )))


def _ligand_distances(coordinates, ligand):
    heavy = ligand.coordinates[ligand.atomic_numbers > 1]
    return np.min(np.linalg.norm(coordinates[:, None, :] - heavy[None, :, :], axis=2), axis=1)


def _mapping_quality(mapping, total):
    if mapping.matched_count == total and all(status == "exact_ordered" for status in mapping.statuses):
        return "exact"
    if mapping.matched_count == total:
        return "aligned"
    if mapping.matched_count:
        return "ambiguous"
    return "unresolved"


def _comparison_quality(mapping, total, config):
    if not mapping.source_file_id:
        return "unavailable"
    ratio = mapping.matched_count / total
    if ratio < config.comparison.atom_jaccard_severe_minimum:
        return "severe_difference"
    if ratio < config.comparison.atom_jaccard_moderate_minimum:
        return "moderate_difference"
    return "concordant"


def _source_role(pdb_id: str, name: str) -> str:
    lower = name.lower()
    exact = {
        "data.pkl": "pickle", f"{pdb_id}_ligand.sdf": "ligand_sdf",
        f"{pdb_id}_ligand.mol2": "ligand_mol2", f"{pdb_id}_protein.pdb": "protein_pdb",
        f"{pdb_id}_pocket.pdb": "pocket_pdb", f"{pdb_id}_pocket6a.pdb": "pocket6a_pdb",
    }
    return exact.get(lower, "extra")


def _relative_path(path: Path, config: BuildConfig) -> str:
    try:
        return path.resolve().relative_to(config.project_root.resolve()).as_posix()
    except ValueError:
        root = config.paths.combine_set_root
        if root is not None:
            try:
                return "external/DrugCLIP/data/pdb/combine_set/" + path.resolve().relative_to(root.resolve()).as_posix()
            except ValueError:
                pass
        return path.resolve().as_posix()


def _optional_path(files: dict, key: str) -> Path | None:
    value = files.get(key)
    return Path(value) if value is not None else None
