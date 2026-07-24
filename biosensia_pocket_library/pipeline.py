"""End-to-end deterministic candidate-pocket build orchestration."""

from __future__ import annotations

import dataclasses
import csv
import json
import os
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np

from .config import BuildConfig
from .drugclip_contract import read_dictionary, verify_drugclip_contract
from .event_logging import EventLogger
from .hashing import canonical_json_hash, sha256_file, stable_id
from .index_parser import parse_index
from .identity import verify_dataset_identity
from .inventory import inventory_sources
from .ligand_parser import ligand_component_rows, parse_ligand
from .lmdb_export import export_lmdb
from .manifest import complete_stage, create_manifest, utc_now, write_manifest, write_resolved_config
from .models import ProcessingIssue
from .pocket_comparison import compare_pdbbind_pocket
from .pocket_extractor import distance_statistics, extract_pocket
from .progress import file_progress, track
from .protein_parser import parse_protein
from .quality import classify_geometry
from .reference_links import merge_affinity_reference_links
from .rcsb import download_mmcif_files, enrich_from_mmcif
from .reporting import generate_reports
from .scrub import scrub
from .sidecars import write_sidecars
from .validation import validate_run

DISTRIBUTION_ID = "pdbbind-2020-v2024p-20250804"


def build_library(config: BuildConfig, *, pdb_ids: Iterable[str] | None = None, limit: int | None = None,
                  year_from: int | None = None, year_to: int | None = None, resume: bool = False,
                  overwrite_run: bool = False, export: bool = True, progress: bool = True) -> Path:
    """Build a run. Bootstrap is read-only until every identity fingerprint is known."""
    index_path = config.paths.index_dir / "INDEX_general_PL.2020R1.lst"
    all_records, occurrences, index_summary = parse_index(index_path, DISTRIBUTION_ID)
    if index_summary["declared_count"] != index_summary["physical_data_line_count"]:
        from .exceptions import SourceIntegrityError
        raise SourceIntegrityError(
            f"Declared PL count {index_summary['declared_count']} differs from parsed data-line count "
            f"{index_summary['physical_data_line_count']}"
        )
    dataset_identity = verify_dataset_identity(config.paths.index_dir, index_summary["declared_count"])
    selected = _select(all_records, pdb_ids, limit, year_from, year_to)
    if not selected:
        raise ValueError("Selection contains no complexes")
    selection_spec = {"pdb_ids": sorted(record.pdb_id for record in selected), "limit": limit,
                      "year_from": year_from, "year_to": year_to}
    selection_fingerprint = canonical_json_hash(selection_spec)
    contract = verify_drugclip_contract(config, progress=progress)
    source_rows_raw, discovery, bootstrap_issues = inventory_sources(
        config.paths.index_dir, config.paths.complex_root, selected, progress=progress)
    source_rows, source_id_map = _source_rows(source_rows_raw, selected, config)
    for item in discovery.values():
        for key, value in list(item["files"].items()):
            if value in source_id_map:
                item["files"][key] = source_id_map[value]
    source_fingerprint = canonical_json_hash({"distribution_id": DISTRIBUTION_ID,
        "files": sorted((row["source_file_id"], row["sha256"], row["validation_status"]) for row in source_rows)})
    run_id = (f"pb20-v24p-20250804-v1-{config.semantic_hash[:8]}-{source_fingerprint[:8]}-"
              f"{selection_fingerprint[:8]}-{contract['drugclip_contract_fingerprint'][:8]}")
    run_dir = config.paths.output_root / run_id
    if run_dir.exists():
        existing = _read_manifest(run_dir)
        compatible = existing and all(existing.get(key) == value for key, value in {
            "semantic_config_hash": config.semantic_hash, "source_fingerprint": source_fingerprint,
            "selection_fingerprint": selection_fingerprint,
            "drugclip_contract_fingerprint": contract["drugclip_contract_fingerprint"]}.items())
        if resume and compatible and existing.get("status") == "complete":
            return run_dir
        if resume and compatible:
            pass
        elif overwrite_run:
            backup = run_dir.with_name(f"{run_dir.name}.backup-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}")
            os.replace(run_dir, backup)
        else:
            raise FileExistsError(f"Run already exists: {run_dir}; use --resume or --overwrite-run")
    run_dir.mkdir(parents=True, exist_ok=True)
    manifest = create_manifest(config, run_id, source_fingerprint, selection_fingerprint, contract,
                               index_summary, dataset_identity)
    manifest["dataset"]["discovered_complex_directory_count"] = sum(bool(item["complex_directory"]) for item in discovery.values())
    write_resolved_config(run_dir / "config.resolved.toml", config)
    write_manifest(run_dir, manifest)
    logger = EventLogger(run_dir)
    logger.emit("info", "bootstrap-identity", "RUN_STARTED", f"Starting {run_id}")
    complete_stage(run_dir, manifest, "bootstrap-identity", "1", {"selection": selection_fingerprint},
                   [run_dir / "config.resolved.toml"], progress=progress)
    rows = _empty_rows()
    rows["source_files"] = source_rows
    selected_ids = {record.complex_id for record in selected}
    rows["index_record_occurrences"] = [row for row in occurrences if row["complex_id"] in selected_ids]
    rows["binding_measurements"] = [dataclasses.asdict(record.measurement) for record in selected]
    issues = list(bootstrap_issues)
    dictionary = read_dictionary(config.paths.drugclip_dictionary)
    chains_by_pdb: dict[str, list[dict]] = defaultdict(list)
    source_hash_by_id = {row["source_file_id"]: row["sha256"] for row in source_rows}
    arguments = [(record, discovery[record.complex_id], config, dictionary, source_hash_by_id)
                 for record in selected]
    results = (_bounded_thread_map(_process_record_geometry, arguments, config.pipeline.workers)
               if config.pipeline.workers > 1 else
               map(lambda values: _process_record_geometry(*values), arguments))
    for complex_row, local_rows, local_issues, event in track(
        results, description="Extracting candidate pockets", total=len(selected), enabled=progress
    ):
        rows["complexes"].append(complex_row)
        for name, values in local_rows.items():
            rows[name].extend(values)
        issues.extend(local_issues)
        chains_by_pdb[complex_row["pdb_id"]].extend(local_rows.get("protein_chains", []))
        logger.emit(**event)
    # Optional RCSB enrichment remains downstream of and cannot modify geometry decisions/hashes.
    if config.rcsb.download_mmcif:
        try:
            cached, external_inventory, download_failures = download_mmcif_files(
                [record.pdb_id for record in selected], config, progress=progress)
            rows["source_files"].extend(_normalize_external_rows(external_inventory, config))
            for failure in download_failures:
                issues.append(ProcessingIssue("download-rcsb", "RCSB_ENRICHMENT_FAILED", "warning",
                                              failure["error"], complex_id=None,
                                              exception_type=failure["exception_type"], details=failure))
            cached_ids = set(cached)
            for complex_row in rows["complexes"]:
                if complex_row["pdb_id"] not in cached_ids:
                    complex_row["rcsb_entry_status"] = "unavailable"
                    complex_row["structure_mapping_quality"] = "unavailable"
            for pocket_row in rows["pockets"]:
                if pocket_row["pdb_id"] not in cached_ids:
                    pocket_row["structure_mapping_quality"] = "unavailable"
            for pdb_id, path in track(sorted(cached.items()), description="Enriching from RCSB", total=len(cached), enabled=progress):
                ligand_inputs = []
                labels_by_complex = {record.complex_id: record.ligand_label for record in selected if record.pdb_id == pdb_id}
                for ligand_row in rows["ligand_instances"]:
                    if ligand_row["pdb_id"] == pdb_id:
                        ligand_inputs.append((ligand_row, labels_by_complex.get(ligand_row["complex_id"])))
                enriched = enrich_from_mmcif(pdb_id, path, config.rcsb.download_compressed,
                                             chains_by_pdb[pdb_id], ligand_inputs)
                for name, values in enriched.items():
                    rows[name].extend(values)
                for ligand_row in rows["ligand_instances"]:
                    if ligand_row["pdb_id"] != pdb_id:
                        continue
                    candidates_for_ligand = [item for item in enriched["rcsb_ligand_mapping_candidates"]
                                             if item["ligand_instance_id"] == ligand_row["ligand_instance_id"]]
                    ligand_row["rcsb_ligand_match_overall_status"] = (
                        "probable" if any(item["selected"] for item in candidates_for_ligand)
                        else "ambiguous" if any(item["status"] == "ambiguous" for item in candidates_for_ligand)
                        else "unresolved"
                    )
                selected_chain_keys = {(item["pocket_instance_id"], item["pdbbind_auth_chain_id"])
                                       for item in enriched["chain_mapping_candidates"] if item["selected"]}
                candidate_chain_keys = {(item["pocket_instance_id"], item["pdbbind_auth_chain_id"])
                                        for item in enriched["chain_mapping_candidates"]}
                def mapping_quality_for(pocket_id):
                    required = {(item["pocket_instance_id"], item["pdbbind_auth_chain_id"])
                                for item in chains_by_pdb[pdb_id] if item["pocket_instance_id"] == pocket_id}
                    return ("exact" if required and required <= selected_chain_keys else
                            "ambiguous" if required & candidate_chain_keys else "unresolved")
                for row in rows["complexes"]:
                    if row["pdb_id"] == pdb_id:
                        row["rcsb_entry_status"] = "current"
                        related = [item for item in rows["pockets"] if item["complex_id"] == row["complex_id"]]
                        row["structure_mapping_quality"] = mapping_quality_for(related[0]["pocket_instance_id"]) if related else "unresolved"
                for chain_row in rows["protein_chains"]:
                    if chain_row["pdb_id"] != pdb_id:
                        continue
                    key = (chain_row["pocket_instance_id"], chain_row["pdbbind_auth_chain_id"])
                    chain_row["rcsb_mapping_status"] = ("exact_identifier_match" if key in selected_chain_keys
                                                        else "ambiguous" if key in candidate_chain_keys
                                                        else "unresolved")
                for pocket_row in rows["pockets"]:
                    if pocket_row["pdb_id"] == pdb_id:
                        pocket_row["structure_mapping_quality"] = mapping_quality_for(pocket_row["pocket_instance_id"])
        except Exception as error:
            issues.append(ProcessingIssue("download-rcsb", "RCSB_ENRICHMENT_FAILED", "warning", str(error),
                                          exception_type=type(error).__name__))
            logger.emit("warning", "download-rcsb", "RCSB_ENRICHMENT_FAILED", str(error))
    _normalize_enrichment_rows(rows)
    citations_by_pdb = defaultdict(list)
    for link in rows["pdb_citation_links"]:
        citations_by_pdb[link["pdb_id"]].append(link)
    for record in selected:
        citation_links = citations_by_pdb[record.pdb_id]
        primary = sorted({link["citation_id"] for link in citation_links if link["role"] == "primary"})
        for link in citation_links:
            rows["affinity_reference_links"].append({
                "measurement_id": record.measurement.measurement_id, "complex_id": record.complex_id,
                "citation_id": link["citation_id"],
                "candidate_status": "probable_structural_reference" if link["role"] == "primary" else "structural_reference_only",
                "confidence": 0.60 if link["role"] == "primary" else 0.30,
                "evidence_sources": [link["source"]],
                "evidence_note": "Depositor citation; not asserted as the affinity-measurement source",
                "automatic_or_manual": "automatic", "verified_by": None, "verified_at_utc": None,
            })
        if len(primary) == 1:
            status, selected_citation, bibliography_quality = "probable_structural_reference", primary[0], "probable"
        elif len(primary) > 1:
            status, selected_citation, bibliography_quality = "conflicting_references", None, "unresolved"
        elif citation_links:
            status, selected_citation, bibliography_quality = "reference_unresolved", None, "unresolved"
        elif config.bibliography.external_enrichment_enabled:
            status, selected_citation, bibliography_quality = "no_reference_available", None, "unavailable"
        else:
            status, selected_citation, bibliography_quality = "not_attempted", None, "not_attempted"
        rows["affinity_reference_adjudications"].append({"measurement_id": record.measurement.measurement_id,
            "reference_status": status, "selected_citation_id": selected_citation, "rule_version": "1",
            "confidence": 0.60 if selected_citation else None,
            "evidence_summary": "Structural-citation evidence only; no automatic affinity-source assertion",
            "adjudicator": "automatic-v1", "adjudicated_at_utc": None})
        for complex_row in rows["complexes"]:
            if complex_row["complex_id"] == record.complex_id:
                complex_row["bibliography_quality"] = bibliography_quality
        for pocket_row in rows["pockets"]:
            if pocket_row["complex_id"] == record.complex_id:
                pocket_row["bibliography_quality"] = bibliography_quality
    _apply_reference_overrides(rows, config, issues)
    rows["affinity_reference_links"] = merge_affinity_reference_links(rows["affinity_reference_links"])
    rows["processing_issues"] = [_issue_row(issue) for issue in issues]
    warning_counts = Counter(issue.complex_id for issue in issues if issue.complex_id and issue.severity == "warning")
    error_counts = Counter(issue.complex_id for issue in issues if issue.complex_id and issue.severity in {"error", "fatal"})
    for row in rows["complexes"]:
        row["warning_count"] = warning_counts[row["complex_id"]]
        row["error_count"] = error_counts[row["complex_id"]]
        if row["processing_status"] == "accepted" and row["warning_count"]:
            row["processing_status"] = "accepted_with_warnings"
    sidecar_results = write_sidecars(run_dir / "sidecars", rows, progress=progress)
    stage_outputs = {
        "inventory": ["source_files"],
        "parse-index": ["complexes", "binding_measurements", "index_record_occurrences"],
        "parse-structures": ["ligand_instances", "ligand_components"],
        "extract-pockets": ["pockets", "pocket_atoms", "pocket_residues", "protein_chains", "nearby_nonprotein_components"],
        "compare-pockets": ["pocket_comparisons", "pocket_atom_differences"],
        "download-rcsb": ["source_files"],
        "map-structures": ["chain_mapping_candidates", "chain_uniprot_mappings",
                           "chain_uniprot_mapping_segments", "rcsb_ligand_mapping_candidates"],
        "enrich-citations": ["citations", "citation_authors", "pdb_citation_links",
                             "affinity_reference_links", "affinity_reference_adjudications"],
        "quality-control": ["complexes", "pockets", "processing_issues"],
    }
    complete_stage(run_dir, manifest, "check-drugclip-contract", "1",
                   {"contract": contract["drugclip_contract_fingerprint"]}, [], progress=progress)
    for stage, names in stage_outputs.items():
        complete_stage(run_dir, manifest, stage, "1", {"source": source_fingerprint},
                       [Path(sidecar_results[name]["path"]) for name in names], progress=progress)
    complete_stage(run_dir, manifest, "write-sidecars", "1", {"source": source_fingerprint},
                   [Path(item["path"]) for item in sidecar_results.values()], progress=progress)
    lmdb_metadata = None
    if export:
        lmdb_metadata, lmdb_rows = export_lmdb(run_dir, config, contract, "default",
                                               overwrite=(config.lmdb.overwrite or resume), progress=progress)
        rows["lmdb_records"] = lmdb_rows
        manifest.setdefault("lmdb_profiles", {})["default"] = {
            **lmdb_metadata, "path": Path(lmdb_metadata["path"]).relative_to(run_dir).as_posix()
        }
        sidecar_results = write_sidecars(run_dir / "sidecars", rows, progress=progress)
        complete_stage(run_dir, manifest, "export-lmdb", "1", {"pockets": sidecar_results["pockets"]["logical_sha256"]},
                       [Path(lmdb_metadata["path"]), Path(lmdb_metadata["path"] + ".profile.json")], progress=progress)
    errors = validate_run(run_dir, config, progress=progress)
    if errors:
        manifest["status"] = "validation_failed"
        manifest["counts"]["validation_errors"] = len(errors)
        write_manifest(run_dir, manifest)
        raise ValueError("Run validation failed: " + "; ".join(errors))
    complete_stage(run_dir, manifest, "validate", "1", {}, [], progress=progress)
    report_results = generate_reports(run_dir, manifest)
    complete_stage(run_dir, manifest, "report", "1", {}, [Path(item["path"]) for item in report_results.values()],
                   progress=progress)
    manifest["counts"] = {"selected_complexes": len(selected), "accepted_pockets": len(rows["pockets"]),
                          "rejected_complexes": sum(row["processing_status"] == "rejected" for row in rows["complexes"]),
                          "lmdb_records": len(rows["lmdb_records"])}
    manifest["sidecar_artifacts"] = {
        name: {**item, "path": Path(item["path"]).relative_to(run_dir).as_posix()}
        for name, item in sidecar_results.items()
    }
    manifest["output_files"] = _artifact_inventory(run_dir, progress=progress)
    manifest["completed_at_utc"] = utc_now()
    manifest["status"] = "complete" if export else "sidecars_complete"
    write_manifest(run_dir, manifest)
    logger.emit("info", "report", "RUN_COMPLETED", f"Completed {run_id}")
    return run_dir


def _select(records, pdb_ids, limit, year_from, year_to):
    wanted = {value.lower() for value in pdb_ids} if pdb_ids else None
    values = [record for record in records if (wanted is None or record.pdb_id in wanted)
              and (year_from is None or record.release_year >= year_from)
              and (year_to is None or record.release_year <= year_to)]
    values.sort(key=lambda record: record.pdb_id)
    return values[:limit] if limit is not None else values


def _bounded_thread_map(operation, arguments, workers):
    """Run at most ``workers`` calls and yield them as they complete.

    ``Executor.map`` eagerly submits the whole input and yields in input order. A
    slow early complex can therefore retain every later completed result in RAM.
    This scheduler bounds both running and completed futures and has no
    head-of-line blocking.
    """
    iterator = iter(arguments)
    executor = ThreadPoolExecutor(max_workers=workers, thread_name_prefix="pocket")
    pending: set[Future] = set()

    def submit_one() -> bool:
        try:
            values = next(iterator)
        except StopIteration:
            return False
        pending.add(executor.submit(operation, *values))
        return True

    try:
        for _ in range(workers):
            if not submit_one():
                break
        while pending:
            completed, pending = wait(pending, return_when=FIRST_COMPLETED)
            # Future completion order is immaterial: sidecars are canonically sorted.
            for future in completed:
                yield future.result()
                submit_one()
    finally:
        for future in pending:
            future.cancel()
        executor.shutdown(wait=True, cancel_futures=True)


def _process_record_geometry(record, info, config, dictionary, source_hash_by_id):
    local_rows: dict[str, list[dict]] = defaultdict(list)
    local_issues: list[ProcessingIssue] = []
    complex_row = _complex_base(record, info, config)
    processing_stage = "ligand"
    try:
        files = info["files"]
        protein_path = _path(files.get("protein_pdb_path"))
        sdf_path = _path(files.get("ligand_sdf_path"))
        mol2_path = _path(files.get("ligand_mol2_path"))
        if protein_path is None:
            processing_stage = "protein"
            raise FileNotFoundError("Missing protein PDB")
        ligand = parse_ligand(record.complex_id, sdf_path, mol2_path, config)
        processing_stage = "protein"
        protein = parse_protein(protein_path, config, pdb_id=record.pdb_id, complex_id=record.complex_id)
        processing_stage = "extraction"
        pocket = extract_pocket(record.complex_id, record.pdb_id, ligand, protein, dictionary, config,
                                protein_source_sha256=source_hash_by_id.get(files.get("protein_pdb")))
        classify_geometry(pocket, ligand, protein, config)
        comparison_rows, difference_rows, comparison_quality = compare_pdbbind_pocket(
            pocket, ligand, _path(files.get("pocket_pdb_path")), files.get("pocket_pdb"), config)
        if comparison_quality == "unavailable":
            local_issues.append(ProcessingIssue("compare-pockets", "PDBBIND_POCKET_COMPARISON_UNAVAILABLE",
                                                "warning", "PDBbind-provided pocket comparison is unavailable",
                                                record.complex_id, pocket.pocket_instance_id))
        local_rows["pocket_comparisons"].extend(comparison_rows)
        local_rows["pocket_atom_differences"].extend(difference_rows)
        local_rows["ligand_instances"].append(_ligand_row(record, ligand, files))
        local_rows["ligand_components"].extend(ligand_component_rows(ligand))
        local_rows["pockets"].append(_pocket_row(
            pocket, protein, comparison_quality, config, files.get("protein_pdb")
        ))
        local_rows["pocket_atoms"].extend(_atom_rows(pocket, dictionary))
        local_rows["pocket_residues"].extend(_residue_rows(pocket))
        local_rows["protein_chains"].extend(_chain_rows(pocket))
        local_rows["nearby_nonprotein_components"].extend(_nearby_rows(pocket, ligand, protein, config))
        complex_row.update(processing_status="accepted_with_warnings" if pocket.warning_codes else "accepted",
                           geometry_quality_tier=pocket.geometry_quality_tier,
                           pocket_comparison_quality=comparison_quality)
        for code in pocket.warning_codes:
            local_issues.append(ProcessingIssue("quality-control", code, "warning",
                                                code.replace("_", " ").title(), record.complex_id,
                                                pocket.pocket_instance_id))
        if protein.discarded_altloc_count:
            local_issues.append(ProcessingIssue("parse-structures", "ALTERNATE_LOCATIONS_DISCARDED", "info",
                                                "Alternate-location atoms were discarded deterministically",
                                                record.complex_id, pocket.pocket_instance_id,
                                                details={"count": protein.discarded_altloc_count}))
        if protein.inferred_element_count:
            local_issues.append(ProcessingIssue("parse-structures", "PROTEIN_ELEMENTS_INFERRED", "info",
                                                "Protein elements were inferred from atom names",
                                                record.complex_id, pocket.pocket_instance_id,
                                                details={"count": protein.inferred_element_count}))
        event = {"level": "info", "stage": "extract-pockets", "code": "POCKET_ACCEPTED",
                 "message": "Pocket geometry accepted", "complex_id": record.complex_id,
                 "pocket_instance_id": pocket.pocket_instance_id,
                 "details": {"tier": pocket.geometry_quality_tier, "atom_count": len(pocket.exported_atoms)}}
    except Exception as error:
        if config.pipeline.fail_fast:
            raise
        complex_row.update(processing_status="rejected", geometry_quality_tier="rejected")
        issue_code = {"ligand": "LIGAND_BOTH_FORMATS_FAILED", "protein": "PROTEIN_PARSE_FAILED",
                      "extraction": "POCKET_EXTRACTION_FAILED"}[processing_stage]
        local_issues.append(ProcessingIssue("parse-structures", issue_code, "error",
                                            str(error), record.complex_id, exception_type=type(error).__name__))
        event = {"level": "error", "stage": "parse-structures", "code": issue_code,
                 "message": str(error), "complex_id": record.complex_id}
    return complex_row, local_rows, local_issues, event


def _empty_rows():
    from .schemas import TABLES
    return {name: [] for name in TABLES}


def _source_rows(raw, records, config):
    pdb_by_complex = {record.complex_id: record.pdb_id for record in records}
    root = config.project_root.resolve()
    result = []
    identifier_map = {}
    for row in raw:
        path = Path(row["path"])
        try:
            normalized = path.relative_to(root).as_posix()
        except ValueError:
            try:
                normalized = "data/raw/P-L/" + path.relative_to(config.paths.complex_root.resolve()).as_posix()
            except ValueError:
                normalized = path.as_posix()
        name = path.name.lower()
        kind = "pdbbind_index" if row["role"] == "index" else (
            "pdbbind_protein" if name.endswith("_protein.pdb") else
            "pdbbind_ligand_sdf" if name.endswith("_ligand.sdf") else
            "pdbbind_ligand_mol2" if name.endswith("_ligand.mol2") else
            "pdbbind_pocket" if name.endswith("_pocket.pdb") else "pdbbind_extra")
        source_id = stable_id("file", kind, normalized, row["sha256"])
        identifier_map[row["source_file_id"]] = source_id
        result.append({"source_file_id": source_id, "source_kind": kind,
            "pdb_id": pdb_by_complex.get(row["complex_id"]), "path": normalized,
            "size_bytes": row["size_bytes"], "sha256": row["sha256"], "modified_time_utc": None,
            "download_url": None, "downloaded_at_utc": None, "http_etag": None, "http_last_modified": None,
            "validation_status": "missing" if row["role"] == "missing" else ("valid" if row["size_bytes"] else "empty"),
            "warning_codes": []})
    return result, identifier_map


def _complex_base(record, info, config):
    files = info["files"]
    directory = info["complex_directory"]
    if directory:
        try:
            suffix = Path(directory).relative_to(config.paths.complex_root.resolve()).as_posix()
            directory = f"data/raw/P-L/{suffix}"
        except ValueError:
            directory = Path(directory).as_posix()
    return {"complex_id": record.complex_id, "pdb_id": record.pdb_id, "distribution_id": record.distribution_id,
        "geometry_origin": "pdbbind_reextracted", "geometry_source_file_id": files.get("protein_pdb"),
        "nominal_complex_set_version": "2020", "structure_processing_version": "2024",
        "index_revision_date": "2025-08-04", "primary_index_line_number": record.primary_index_line_number,
        "index_line_redacted": record.index_line_redacted, "source_line_sha256": record.source_line_sha256,
        "release_year": record.release_year, "resolution_raw": record.resolution_raw,
        "resolution_angstrom": record.resolution_angstrom, "experimental_method_hint": record.experimental_method_hint,
        "ligand_label": record.ligand_label, "index_comment": record.index_comment,
        "complex_directory": directory, "protein_file_id": files.get("protein_pdb"),
        "ligand_sdf_file_id": files.get("ligand_sdf"), "ligand_mol2_file_id": files.get("ligand_mol2"),
        "pdbbind_pocket_file_id": files.get("pocket_pdb"), "rcsb_entry_status": "not_attempted",
        "processing_status": "not_processed", "geometry_quality_tier": "not_processed",
        "pocket_comparison_quality": "not_processed", "structure_mapping_quality": "not_processed",
        "bibliography_quality": "not_attempted", "warning_count": 0, "error_count": 0}


def _normalize_external_rows(rows, config):
    normalized = []
    for row in rows:
        value = dict(row)
        path = Path(value["path"])
        try:
            value["path"] = path.relative_to(config.project_root).as_posix()
        except ValueError:
            value["path"] = path.as_posix()
        normalized.append(value)
    return normalized


def _normalize_enrichment_rows(rows):
    citations = {}
    for row in rows["citations"]:
        prior = citations.get(row["citation_id"])
        if prior is None:
            citations[row["citation_id"]] = row
            continue
        comparable = {key: value for key, value in row.items() if key not in {"metadata_sources", "conflict_status"}}
        prior_comparable = {key: value for key, value in prior.items() if key not in {"metadata_sources", "conflict_status"}}
        if comparable != prior_comparable:
            prior["conflict_status"] = "metadata_conflict"
        prior["metadata_sources"] = sorted(set(prior["metadata_sources"]) | set(row["metadata_sources"]))
    rows["citations"] = list(citations.values())
    for name, keys in (
        ("citation_authors", ("citation_id", "source", "ordinal")),
        ("pdb_citation_links", ("citation_id", "pdb_id", "source", "role")),
        ("chain_uniprot_mappings", ("pocket_instance_id", "pdbbind_auth_chain_id", "uniprot_accession")),
    ):
        unique = {}
        for row in rows[name]:
            unique.setdefault(tuple(row[key] for key in keys), row)
        rows[name] = list(unique.values())


def _apply_reference_overrides(rows, config, issues):
    path = config.bibliography.manual_overrides_path
    if not path.is_file():
        return
    if path.suffix.lower() == ".parquet":
        import pyarrow.parquet as pq
        overrides = pq.read_table(path).to_pylist()
    else:
        with path.open(newline="", encoding="utf-8") as handle:
            overrides = list(csv.DictReader(handle))
    allowed = {"exact_affinity_reference", "probable_affinity_reference", "probable_structural_reference",
               "structural_reference_only", "conflicting_references", "reference_unresolved",
               "no_reference_available", "not_attempted"}
    by_measurement = defaultdict(list)
    for override in overrides:
        by_measurement[override.get("measurement_id")].append(override)
    adjudications = {row["measurement_id"]: row for row in rows["affinity_reference_adjudications"]}
    citation_ids = {row["citation_id"] for row in rows["citations"]}
    for measurement_id, values in by_measurement.items():
        if measurement_id not in adjudications:
            issues.append(ProcessingIssue("enrich-citations", "INVALID_REFERENCE_OVERRIDE", "fatal",
                                          "Override names an unknown measurement", details={"measurement_id": measurement_id}))
            continue
        normalized = {(item.get("citation_id") or None, item.get("reference_status")) for item in values}
        if len(normalized) != 1:
            adjudications[measurement_id].update(reference_status="conflicting_references",
                                                  selected_citation_id=None, confidence=None,
                                                  evidence_summary="Conflicting manual overrides",
                                                  adjudicator="manual-override-v1")
            issues.append(ProcessingIssue("enrich-citations", "CONFLICTING_REFERENCE_OVERRIDES", "fatal",
                                          "Conflicting manual overrides", details={"measurement_id": measurement_id}))
            continue
        citation_id, status = next(iter(normalized))
        if status not in allowed or (citation_id and citation_id not in citation_ids):
            issues.append(ProcessingIssue("enrich-citations", "INVALID_REFERENCE_OVERRIDE", "fatal",
                                          "Override has invalid status or citation", details={"measurement_id": measurement_id}))
            continue
        chosen = values[0]
        measurement_row = next(row for row in rows["binding_measurements"] if row["measurement_id"] == measurement_id)
        if chosen.get("complex_id") and chosen.get("complex_id") != measurement_row["complex_id"]:
            issues.append(ProcessingIssue("enrich-citations", "INVALID_REFERENCE_OVERRIDE", "fatal",
                                          "Override complex does not own measurement", details={"measurement_id": measurement_id}))
            continue
        adjudications[measurement_id].update(reference_status=status, selected_citation_id=citation_id,
                                              confidence=1.0, evidence_summary=chosen.get("evidence_note"),
                                              adjudicator=chosen.get("verified_by") or "manual",
                                              adjudicated_at_utc=_parse_datetime(chosen.get("verified_at_utc")))
        if citation_id:
            complex_id = measurement_row["complex_id"]
            rows["affinity_reference_links"] = [
                row for row in rows["affinity_reference_links"]
                if not (row["measurement_id"] == measurement_id and row["citation_id"] == citation_id)
            ]
            rows["affinity_reference_links"].append({"measurement_id": measurement_id,
                "complex_id": complex_id, "citation_id": citation_id, "candidate_status": status,
                "confidence": 1.0, "evidence_sources": ["manual_override"],
                "evidence_note": chosen.get("evidence_note"), "automatic_or_manual": "manual",
                "verified_by": chosen.get("verified_by"),
                "verified_at_utc": _parse_datetime(chosen.get("verified_at_utc"))})
        bibliography_quality = ("exact" if status == "exact_affinity_reference" else
                                "probable" if status in {"probable_affinity_reference", "probable_structural_reference"}
                                else "unavailable" if status == "no_reference_available" else
                                "not_attempted" if status == "not_attempted" else "unresolved")
        for table_name in ("complexes", "pockets"):
            for row in rows[table_name]:
                if row["complex_id"] == measurement_row["complex_id"]:
                    row["bibliography_quality"] = bibliography_quality


def _parse_datetime(value):
    if not value:
        return None
    try:
        parsed = datetime.fromisoformat(str(value).replace("Z", "+00:00"))
        return parsed if parsed.tzinfo is not None else parsed.replace(tzinfo=timezone.utc)
    except ValueError:
        return None


def _ligand_row(record, ligand, files):
    chemistry = {key: value for key, value in ligand.chemistry.items() if not key.startswith("_")}
    return {"ligand_instance_id": ligand.ligand_instance_id, "complex_id": record.complex_id, "pdb_id": record.pdb_id,
        "selected_source_format": ligand.source_format,
        "selected_source_file_id": files.get(f"ligand_{ligand.source_format}"),
        "ligand_geometry_content_hash": ligand.content_hash, "ligand_derivation_hash": ligand.derivation_hash,
        "rdkit_parse_status": "parsed", "rdkit_sanitization_status": "sanitized",
        **chemistry, "sdf_mol2_comparison_status": ligand.comparison["status"],
        "sdf_mol2_coordinate_rmsd": ligand.comparison["coordinate_rmsd"],
        "rcsb_ligand_match_overall_status": "not_processed", "warnings": ligand.warnings}


def _pocket_row(pocket, protein, comparison_quality, config, geometry_source_file_id=None):
    chains = sorted({atom.auth_chain_id for atom in pocket.exported_atoms})
    return {"pocket_instance_id": pocket.pocket_instance_id, "complex_id": pocket.complex_id,
        "ligand_instance_id": pocket.ligand_instance_id, "pdb_id": pocket.pdb_id,
        "pocket_geometry_content_hash": pocket.content_hash, "pocket_derivation_hash": pocket.derivation_hash,
        "extraction_schema_version": config.pipeline.extraction_version,
        "geometry_origin": "pdbbind_reextracted", "geometry_source_file_id": geometry_source_file_id,
        "derivation_method": "ligand_distance_contact_atom",
        "source_geometry_atom_count": len(pocket.residue_expanded_atoms),
        "source_geometry_heavy_atom_count": sum(atom.element != "H" for atom in pocket.residue_expanded_atoms),
        "distance_cutoff_angstrom": config.pocket.distance_cutoff_angstrom,
        "selected_model_id": protein.selected_model_id, "model_count": protein.model_count,
        "altloc_policy": config.structure.altloc_policy,
        "hydrogen_policy": "include" if config.pocket.include_protein_hydrogens else "exclude",
        "contact_atom_count": len(pocket.contact_atoms), "residue_expanded_atom_count": len(pocket.residue_expanded_atoms),
        "exported_atom_count": len(pocket.exported_atoms), "contact_residue_count": len(pocket.contact_residues),
        "drugclip_export_view": "contact_atom",
        "contributing_chain_count": len(chains), "contributing_auth_chain_ids": chains,
        **distance_statistics(pocket), "crop_applied": pocket.crop_applied,
        "crop_max_atoms": config.pocket.max_pocket_atoms,
        "maximum_retained_ligand_distance": pocket.max_retained_distance,
        "minimum_discarded_ligand_distance": pocket.min_discarded_distance,
        "all_elements_supported": "UNSUPPORTED_ATOM_EXCLUDED" not in pocket.warning_codes,
        "processing_status": "accepted_with_warnings" if pocket.warning_codes else "accepted",
        "geometry_quality_tier": pocket.geometry_quality_tier, "pocket_comparison_quality": comparison_quality,
        "structure_mapping_quality": "not_processed", "bibliography_quality": "not_attempted",
        "warning_codes": pocket.warning_codes, "error_codes": pocket.error_codes,
        "lmdb_profile_memberships": ["default"] if pocket.geometry_quality_tier in set(config.lmdb.include_geometry_quality_tiers) else []}


def _atom_rows(pocket, dictionary):
    contact_keys = {atom.pdbbind_atom_key for atom in pocket.contact_atoms}
    retained = {atom.pdbbind_atom_key: index for index, atom in enumerate(pocket.exported_atoms)}
    exported_by_key = {atom.pdbbind_atom_key: atom for atom in pocket.exported_atoms}
    rows = []
    for atom in pocket.residue_expanded_atoms:
        rows.append({"pocket_instance_id": pocket.pocket_instance_id, "pdbbind_atom_key": atom.pdbbind_atom_key,
            "source_atom_key": atom.pdbbind_atom_key, "geometry_source_file_id": None,
            "source_order": atom.source_order, "model_id": atom.model_id, "record_type": atom.record_type,
            "auth_chain_id": atom.auth_chain_id, "auth_residue_number": atom.auth_residue_number,
            "insertion_code": atom.insertion_code, "residue_name": atom.residue_name, "atom_name": atom.atom_name,
            "altloc": atom.altloc,
            "element": exported_by_key.get(atom.pdbbind_atom_key, atom).element,
            "occupancy": atom.occupancy, "b_factor": atom.b_factor,
            "x": atom.x, "y": atom.y, "z": atom.z,
            "minimum_ligand_distance": pocket.minimum_distances[atom.pdbbind_atom_key],
            "in_contact_atom_view": atom.pdbbind_atom_key in contact_keys,
            "in_residue_expanded_atom_view": True, "retained_after_crop": atom.pdbbind_atom_key in retained,
            "export_order": retained.get(atom.pdbbind_atom_key),
            "element_supported_by_drugclip": exported_by_key.get(atom.pdbbind_atom_key, atom).element in dictionary,
            "rcsb_atom_mapping_status": "not_processed", "rcsb_label_asym_id": None,
            "rcsb_label_seq_id": None, "rcsb_atom_id": None, "rcsb_polymer_entity_id": None,
            "source_mapping_status": "native", "included_in_lmdb_source": atom.pdbbind_atom_key in retained})
    return rows


def _residue_rows(pocket):
    rows = []
    for residue in pocket.contact_residues:
        atoms = [atom for atom in pocket.residue_expanded_atoms if atom.residue_key == residue]
        chain, number, insertion, name = residue
        rows.append({"pocket_instance_id": pocket.pocket_instance_id, "pdb_id": pocket.pdb_id,
            "model_id": atoms[0].model_id, "auth_chain_id": chain, "auth_residue_number": number,
            "insertion_code": insertion, "residue_name": name,
            "minimum_ligand_distance": min(pocket.minimum_distances[atom.pdbbind_atom_key] for atom in atoms),
            "selected_atom_count": len(atoms), "total_heavy_atom_count": len(atoms),
            "rcsb_mapping_status": "not_processed", "rcsb_label_asym_id": None,
            "rcsb_label_seq_id": None, "rcsb_polymer_entity_id": None})
    return rows


def _chain_rows(pocket):
    rows = []
    for chain in sorted({atom.auth_chain_id for atom in pocket.exported_atoms}):
        atoms = [atom for atom in pocket.exported_atoms if atom.auth_chain_id == chain]
        rows.append({"pocket_instance_id": pocket.pocket_instance_id, "pdb_id": pocket.pdb_id,
            "pdbbind_auth_chain_id": chain, "selected_atom_count": len(atoms),
            "selected_residue_count": len({atom.residue_key for atom in atoms}),
            "rcsb_mapping_status": "not_processed", "warnings": []})
    return rows


def _nearby_rows(pocket, ligand, protein, config):
    heavy = ligand.coordinates[ligand.atomic_numbers > 1]
    grouped = defaultdict(list)
    for atom in protein.excluded_atoms:
        grouped[(atom.auth_chain_id, atom.auth_residue_number, atom.insertion_code, atom.residue_name)].append(atom)
    rows = []
    for (chain, number, insertion, name), atoms in grouped.items():
        distances = [float(np.min(np.linalg.norm(heavy - atom.coordinate, axis=1))) for atom in atoms if atom.element != "H"]
        if not distances or min(distances) > config.quality.nearby_nonprotein_cutoff_angstrom:
            continue
        counts = Counter(atom.element for atom in atoms)
        component_id = stable_id("component", pocket.pocket_instance_id, chain, number, insertion, name)
        rows.append({"pocket_instance_id": pocket.pocket_instance_id,
            "component_type": "water" if name in {"HOH", "WAT", "DOD"} else "nonprotein",
            "component_id": component_id, "auth_chain_id": chain, "auth_seq_id": number,
            "residue_name": name, "insertion_code": insertion, "atom_count": len(atoms),
            "element_counts_json": json.dumps(counts, sort_keys=True), "minimum_ligand_distance": min(distances),
            "included_in_drugclip_tensor": False, "exclusion_reason": "NONPROTEIN_COMPONENT"})
    return rows


def _issue_row(issue):
    details = scrub(issue.details)
    return {"issue_id": stable_id("issue", issue.stage, issue.complex_id, issue.pocket_instance_id,
                                   issue.issue_code, canonical_json_hash(details)),
        "stage": issue.stage, "complex_id": issue.complex_id, "pocket_instance_id": issue.pocket_instance_id,
        "severity": issue.severity, "issue_code": issue.issue_code, "message": scrub(issue.message),
        "exception_type": issue.exception_type, "source_file_id": issue.source_file_id,
        "created_at_utc": datetime.now(timezone.utc), "details_json": json.dumps(details, sort_keys=True)}


def _path(value):
    return Path(value) if value else None


def _read_manifest(run_dir):
    try:
        return json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _artifact_inventory(run_dir, *, progress=True):
    excluded = {run_dir / "manifest.json", run_dir / "sidecars/source_files.parquet"}
    values = []
    for path in sorted(item for item in run_dir.rglob("*") if item.is_file()):
        if (path in excluded or "checkpoints" in path.parts or "logs" in path.parts
                or path.name.endswith("-lock") or path.suffix == ".tmp"):
            continue
        with file_progress(path, description=f"Final checksum {path.name}",
                           enabled=progress and path.stat().st_size > 50_000_000) as bar:
            digest = sha256_file(path, progress=bar)
        values.append({"path": path.relative_to(run_dir).as_posix(), "size_bytes": path.stat().st_size,
                       "sha256": digest})
    return values
