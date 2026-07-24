"""End-to-end builder for DrugCLIP ``pdb/combine_set`` source bundles."""

from __future__ import annotations

import json
import os
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from .combine_set_source import (
    discover_combine_set,
    inventory_combine_set,
    process_combine_set_record,
    select_combine_set,
)
from .config import BuildConfig
from .drugclip_contract import read_dictionary, verify_drugclip_contract
from .event_logging import EventLogger
from .hashing import canonical_json_hash
from .lmdb_export import export_lmdb
from .manifest import complete_stage, create_manifest, utc_now, write_manifest, write_resolved_config
from .models import ProcessingIssue
from .pipeline import (_artifact_inventory, _bounded_thread_map, _empty_rows, _issue_row,
                       _normalize_external_rows)
from .progress import track
from .reporting import generate_reports
from .rcsb import download_mmcif_files, enrich_from_mmcif
from .sidecars import write_sidecars
from .validation import validate_run


def build_combine_set_library(
    config: BuildConfig, *, pdb_ids: Iterable[str] | None = None, limit: int | None = None,
    resume: bool = False, overwrite_run: bool = False, export: bool = True, progress: bool = True,
) -> Path:
    if not config.combine_set.trusted_pickles:
        from .exceptions import ConfigurationError
        raise ConfigurationError(
            "combine_set builds require combine_set.trusted_pickles=true because pickle can execute code"
        )
    root = config.paths.combine_set_root
    assert root is not None
    discovered = discover_combine_set(root)
    selected = select_combine_set(discovered, pdb_ids, limit)
    if not selected:
        raise ValueError("Selection contains no combine_set complexes")
    selection_spec = {"pdb_ids": [item.pdb_id for item in selected], "limit": limit}
    selection_fingerprint = canonical_json_hash(selection_spec)
    contract = verify_drugclip_contract(config, progress=progress)
    source_rows, bundles, bootstrap_issues = inventory_combine_set(selected, config, progress=progress)
    source_fingerprint = canonical_json_hash({
        "distribution_id": config.combine_set.distribution_id,
        "files": sorted((row["source_file_id"], row["sha256"], row["validation_status"])
                        for row in source_rows),
    })
    run_id = (f"dc-combine-v1-{config.semantic_hash[:8]}-{source_fingerprint[:8]}-"
              f"{selection_fingerprint[:8]}-{contract['drugclip_contract_fingerprint'][:8]}")
    run_dir = config.paths.combine_set_output_root / run_id
    if run_dir.exists():
        existing = _read_manifest(run_dir)
        compatible = existing and all(existing.get(key) == value for key, value in {
            "semantic_config_hash": config.semantic_hash, "source_fingerprint": source_fingerprint,
            "selection_fingerprint": selection_fingerprint,
            "drugclip_contract_fingerprint": contract["drugclip_contract_fingerprint"],
        }.items())
        if resume and compatible and existing.get("status") == "complete":
            return run_dir
        if not (resume and compatible):
            if overwrite_run:
                backup = run_dir.with_name(
                    f"{run_dir.name}.backup-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
                )
                os.replace(run_dir, backup)
            else:
                raise FileExistsError(f"Run already exists: {run_dir}; use --resume or --overwrite-run")
    run_dir.mkdir(parents=True, exist_ok=True)
    dataset_identity = {
        "distribution_id": config.combine_set.distribution_id,
        "nominal_release": config.combine_set.nominal_release,
        "source_kind": "drugclip_combine_set",
        "trusted_pickle_execution_authorized": True,
        "discovered_record_count": len(discovered),
    }
    index_summary = {"declared_count": len(discovered), "unique_complex_count": len(discovered)}
    manifest = create_manifest(
        config, run_id, source_fingerprint, selection_fingerprint, contract,
        index_summary, dataset_identity,
    )
    manifest["dataset"].update({
        "index_file": None, "index_declared_complex_count": len(discovered),
        "index_parsed_complex_count": len(discovered), "discovered_complex_directory_count": len(discovered),
        "selected_record_count": len(selected), "geometry_origin": "drugclip_combine_set_pickle",
    })
    manifest["pipeline_source"] = "drugclip_combine_set"
    write_resolved_config(run_dir / "config.resolved.toml", config)
    write_manifest(run_dir, manifest)
    logger = EventLogger(run_dir)
    logger.emit("info", "bootstrap-identity", "RUN_STARTED", f"Starting {run_id}")
    complete_stage(run_dir, manifest, "bootstrap-identity", "2", {"selection": selection_fingerprint},
                   [run_dir / "config.resolved.toml"], progress=progress)

    rows = _empty_rows()
    rows["source_files"] = source_rows
    issues = list(bootstrap_issues)
    dictionary = read_dictionary(config.paths.drugclip_dictionary)
    arguments = [(record, bundles[record.pdb_id], config, dictionary) for record in selected]
    results = (_bounded_thread_map(process_combine_set_record, arguments, config.pipeline.workers)
               if config.pipeline.workers > 1 else
               map(lambda values: process_combine_set_record(*values), arguments))
    for complex_row, local_rows, local_issues, event in track(
        results, description="Loading trusted combine_set pickles", total=len(selected), enabled=progress
    ):
        rows["complexes"].append(complex_row)
        for name, values in local_rows.items():
            rows[name].extend(values)
        issues.extend(local_issues)
        logger.emit(**event)
    if config.rcsb.download_mmcif:
        _enrich_from_rcsb(rows, config, issues, progress=progress)
    rows["processing_issues"] = [_issue_row(issue) for issue in issues]
    warning_counts = Counter(issue.complex_id for issue in issues
                             if issue.complex_id and issue.severity == "warning")
    error_counts = Counter(issue.complex_id for issue in issues
                           if issue.complex_id and issue.severity in {"error", "fatal"})
    for row in rows["complexes"]:
        row["warning_count"] = warning_counts[row["complex_id"]]
        row["error_count"] = error_counts[row["complex_id"]]
        if row["processing_status"] == "accepted" and row["warning_count"]:
            row["processing_status"] = "accepted_with_warnings"

    sidecar_results = write_sidecars(run_dir / "sidecars", rows, progress=progress)
    complete_stage(run_dir, manifest, "check-drugclip-contract", "2",
                   {"contract": contract["drugclip_contract_fingerprint"]}, [], progress=progress)
    stage_outputs = {
        "inventory": ["source_files"], "parse-index": ["complexes", "binding_measurements"],
        "parse-structures": ["ligand_instances", "ligand_components"],
        "extract-pockets": ["pockets", "pocket_atoms", "pocket_residues", "protein_chains"],
        "compare-pockets": ["pocket_comparisons"], "download-rcsb": [],
        "map-structures": ["protein_chains"], "enrich-citations": ["affinity_reference_adjudications"],
        "quality-control": ["complexes", "pockets", "processing_issues"],
    }
    for stage, names in stage_outputs.items():
        complete_stage(run_dir, manifest, stage, "2", {"source": source_fingerprint},
                       [Path(sidecar_results[name]["path"]) for name in names], progress=progress)
    complete_stage(run_dir, manifest, "write-sidecars", "2", {"source": source_fingerprint},
                   [Path(item["path"]) for item in sidecar_results.values()], progress=progress)
    lmdb_rows = []
    if export:
        lmdb_metadata, lmdb_rows = export_lmdb(
            run_dir, config, contract, "default", overwrite=(config.lmdb.overwrite or resume), progress=progress
        )
        rows["lmdb_records"] = lmdb_rows
        manifest.setdefault("lmdb_profiles", {})["default"] = {
            **lmdb_metadata, "path": Path(lmdb_metadata["path"]).relative_to(run_dir).as_posix(),
        }
        sidecar_results = write_sidecars(run_dir / "sidecars", rows, progress=progress)
        complete_stage(run_dir, manifest, "export-lmdb", "2",
                       {"pockets": sidecar_results["pockets"]["logical_sha256"]},
                       [Path(lmdb_metadata["path"]), Path(lmdb_metadata["path"] + ".profile.json")],
                       progress=progress)
    errors = validate_run(run_dir, config, progress=progress)
    if errors:
        manifest["status"] = "validation_failed"
        manifest["counts"]["validation_errors"] = len(errors)
        write_manifest(run_dir, manifest)
        raise ValueError("Run validation failed: " + "; ".join(errors))
    complete_stage(run_dir, manifest, "validate", "2", {}, [], progress=progress)
    report_results = generate_reports(run_dir, manifest)
    complete_stage(run_dir, manifest, "report", "2", {},
                   [Path(item["path"]) for item in report_results.values()], progress=progress)
    manifest["counts"] = {
        "selected_complexes": len(selected), "accepted_pockets": len(rows["pockets"]),
        "rejected_complexes": sum(row["processing_status"] == "rejected" for row in rows["complexes"]),
        "lmdb_records": len(lmdb_rows),
    }
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


def _read_manifest(run_dir: Path) -> dict | None:
    try:
        return json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _enrich_from_rcsb(rows: dict[str, list[dict]], config: BuildConfig, issues: list[ProcessingIssue],
                      *, progress: bool) -> None:
    """Reuse the established mmCIF mapper after pickle-to-chain mapping."""
    pdb_ids = sorted({row["pdb_id"] for row in rows["complexes"] if row["processing_status"].startswith("accepted")})
    if not pdb_ids:
        return
    try:
        cached, external_inventory, failures = download_mmcif_files(pdb_ids, config, progress=progress)
        rows["source_files"].extend(_normalize_external_rows(external_inventory, config))
        for failure in failures:
            issues.append(ProcessingIssue(
                "download-rcsb", "RCSB_ENRICHMENT_FAILED", "warning", failure["error"],
                exception_type=failure["exception_type"], details=failure,
            ))
        chains_by_pdb: dict[str, list[dict]] = defaultdict(list)
        for chain in rows["protein_chains"]:
            chains_by_pdb[chain["pdb_id"]].append(chain)
        for pdb_id in pdb_ids:
            if pdb_id not in cached:
                for table_name in ("complexes", "pockets"):
                    for row in rows[table_name]:
                        if row["pdb_id"] == pdb_id:
                            row["structure_mapping_quality"] = "unavailable"
                continue
            ligand_inputs = [(row, None) for row in rows["ligand_instances"] if row["pdb_id"] == pdb_id]
            enriched = enrich_from_mmcif(
                pdb_id, cached[pdb_id], config.rcsb.download_compressed,
                chains_by_pdb[pdb_id], ligand_inputs,
            )
            for name, values in enriched.items():
                rows[name].extend(values)
            selected = {(item["pocket_instance_id"], item["pdbbind_auth_chain_id"])
                        for item in enriched["chain_mapping_candidates"] if item["selected"]}
            candidates = {(item["pocket_instance_id"], item["pdbbind_auth_chain_id"])
                          for item in enriched["chain_mapping_candidates"]}
            for chain in chains_by_pdb[pdb_id]:
                key = (chain["pocket_instance_id"], chain["pdbbind_auth_chain_id"])
                chain["rcsb_mapping_status"] = (
                    "exact_identifier_match" if key in selected else "ambiguous" if key in candidates else "unresolved"
                )
            for pocket in (row for row in rows["pockets"] if row["pdb_id"] == pdb_id):
                required = {(chain["pocket_instance_id"], chain["pdbbind_auth_chain_id"])
                            for chain in chains_by_pdb[pdb_id]
                            if chain["pocket_instance_id"] == pocket["pocket_instance_id"]}
                quality = ("exact" if required and required <= selected else
                           "ambiguous" if required & candidates else "unresolved")
                pocket["structure_mapping_quality"] = quality
                for complex_row in rows["complexes"]:
                    if complex_row["complex_id"] == pocket["complex_id"]:
                        complex_row["rcsb_entry_status"] = "current"
                        complex_row["structure_mapping_quality"] = quality
            for ligand in (row for row in rows["ligand_instances"] if row["pdb_id"] == pdb_id):
                related = [item for item in enriched["rcsb_ligand_mapping_candidates"]
                           if item["ligand_instance_id"] == ligand["ligand_instance_id"]]
                ligand["rcsb_ligand_match_overall_status"] = (
                    "probable" if any(item["selected"] for item in related) else
                    "ambiguous" if any(item["status"] == "ambiguous" for item in related) else "unresolved"
                )
    except Exception as error:
        issues.append(ProcessingIssue("download-rcsb", "RCSB_ENRICHMENT_FAILED", "warning", str(error),
                                      exception_type=type(error).__name__))
