"""Finalize a sidecar-complete or interrupted run without rebuilding geometry."""

from __future__ import annotations

import json
from pathlib import Path

from .config import BuildConfig
from .drugclip_contract import verify_drugclip_contract
from .lmdb_export import export_lmdb
from .manifest import complete_stage, git_identity, project_code_fingerprint, utc_now, write_manifest
from .pipeline import _artifact_inventory
from .reporting import generate_reports
from .sidecars import read_sidecar, write_sidecars
from .validation import validate_run


def finalize_run(run_dir: Path, config: BuildConfig, *, progress: bool = True) -> dict:
    """Publish the required default LMDB and atomically finalize run metadata."""
    manifest_path = run_dir / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("run_id") != run_dir.name:
        raise ValueError("Run-directory name differs from manifest run_id")
    if not (run_dir / "sidecars").is_dir():
        raise FileNotFoundError(f"Missing sidecars directory: {run_dir / 'sidecars'}")

    finalization_started = utc_now()
    finalizer_commit, finalizer_dirty = git_identity(config.project_root)
    finalizer_provenance = {
        "schema_version": "1", "finalized_by_command": "finalize",
        "started_at_utc": finalization_started, "git_commit": finalizer_commit,
        "git_dirty": finalizer_dirty,
        "code_dirty_state_fingerprint": project_code_fingerprint(config.project_root),
        "original_status": manifest.get("status"),
    }
    original_identity = {
        key: manifest.get(key) for key in (
            "run_id", "semantic_config_hash", "source_fingerprint", "selection_fingerprint",
            "drugclip_contract_fingerprint", "git_commit", "code_dirty_state_fingerprint",
        )
    }
    _cleanup_lmdb_temporaries(run_dir)
    contract = verify_drugclip_contract(config, progress=progress)
    _verify_library_contract_compatible(manifest, contract)
    metadata, lmdb_rows = export_lmdb(
        run_dir, config, contract, "default", overwrite=True, progress=progress
    )
    rows = {name: read_sidecar(run_dir / "sidecars", name) for name in _table_names()}
    rows["lmdb_records"] = [
        row for row in rows["lmdb_records"] if row["library_profile"] != "default"
    ] + lmdb_rows
    sidecar_results = write_sidecars(run_dir / "sidecars", rows, progress=progress)
    manifest.setdefault("lmdb_profiles", {})["default"] = {
        **metadata, "path": Path(metadata["path"]).relative_to(run_dir).as_posix()
    }
    manifest["finalization"] = finalizer_provenance
    complete_stage(
        run_dir, manifest, "export-lmdb", "2",
        {"pockets": sidecar_results["pockets"]["logical_sha256"],
         "library_contract": contract["drugclip_contract_fingerprint"]},
        [Path(metadata["path"]), Path(metadata["path"] + ".profile.json")], progress=progress,
        execution_provenance=finalizer_provenance,
    )

    errors = validate_run(run_dir, config, progress=progress)
    if errors:
        manifest["status"] = "validation_failed"
        manifest.setdefault("counts", {})["validation_errors"] = len(errors)
        manifest["sidecar_artifacts"] = _relative_sidecar_results(sidecar_results, run_dir)
        manifest["output_files"] = _artifact_inventory(run_dir, progress=progress)
        write_manifest(run_dir, manifest)
        raise ValueError("Run validation failed: " + "; ".join(errors))
    complete_stage(run_dir, manifest, "validate", "2", {"recovery": "finalize"}, [], progress=progress,
                   execution_provenance=finalizer_provenance)

    report_results = generate_reports(run_dir, manifest)
    complete_stage(
        run_dir, manifest, "report", "2", {"recovery": "finalize"},
        [Path(item["path"]) for item in report_results.values()], progress=progress,
        execution_provenance=finalizer_provenance,
    )
    complexes = rows["complexes"]
    pockets = rows["pockets"]
    manifest["counts"] = {
        "selected_complexes": len(complexes),
        "accepted_pockets": sum(row["processing_status"].startswith("accepted") for row in pockets),
        "rejected_complexes": sum(row["processing_status"] == "rejected" for row in complexes),
        "lmdb_records": len(rows["lmdb_records"]),
        "default_lmdb_records": len(lmdb_rows),
    }
    manifest["sidecar_artifacts"] = _relative_sidecar_results(sidecar_results, run_dir)
    _cleanup_lmdb_temporaries(run_dir)
    manifest["output_files"] = _artifact_inventory(run_dir, progress=progress)
    manifest["completed_at_utc"] = utc_now()
    manifest["status"] = "complete"
    manifest["finalization"]["completed_at_utc"] = manifest["completed_at_utc"]
    if any(manifest.get(key) != value for key, value in original_identity.items()):
        raise RuntimeError("Finalization attempted to alter immutable run identity")
    write_manifest(run_dir, manifest)
    return manifest


def _cleanup_lmdb_temporaries(run_dir: Path) -> None:
    lmdb_dir = run_dir / "lmdb"
    if not lmdb_dir.is_dir():
        return
    for path in lmdb_dir.glob(".*.tmp"):
        if path.is_file():
            path.unlink()


def _relative_sidecar_results(results: dict, run_dir: Path) -> dict:
    return {
        name: {**item, "path": Path(item["path"]).relative_to(run_dir).as_posix()}
        for name, item in results.items()
    }


def _verify_library_contract_compatible(manifest: dict, current: dict) -> None:
    recorded_fingerprint = manifest.get("drugclip_library_contract_fingerprint")
    if recorded_fingerprint and recorded_fingerprint != current["drugclip_library_contract_fingerprint"]:
        raise ValueError("Current DrugCLIP library contract differs from the run manifest")
    recorded = manifest.get("drugclip_contract", {})
    mismatches = [
        key for key in ("dictionary_sha256", "task_sha256", "loader_sha256", "helper_sha256")
        if recorded.get(key) is not None and recorded.get(key) != current.get(key)
    ]
    if mismatches:
        raise ValueError(f"Current DrugCLIP library contract differs in: {', '.join(mismatches)}")


def _table_names():
    from .schemas import TABLES
    return TABLES
