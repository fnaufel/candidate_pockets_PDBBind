"""Plan, prefetch, and apply RCSB enrichment as a post-build workflow."""

from __future__ import annotations

import copy
import json
import os
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from .config import BuildConfig
from .hashing import (
    atomic_write_bytes,
    canonical_json_bytes,
    canonical_json_hash,
    sha256_file,
)
from .manifest import (
    complete_stage,
    git_identity,
    project_code_fingerprint,
    utc_now,
    write_manifest,
    write_resolved_config,
)
from .pipeline import _artifact_inventory, _issue_row, _normalize_external_rows
from .rcsb import (
    cached_inventory,
    download_mmcif_files,
    resolve_cached_mmcif_files,
)
from .rcsb_enrichment import ENRICHMENT_TABLES, apply_rcsb_enrichment
from .reporting import generate_reports
from .schemas import TABLES
from .sidecars import read_sidecar, write_sidecars
from .validation import validate_run


REQUEST_SCHEMA_VERSION = "rcsb-request-v1"
CACHE_MANIFEST_SCHEMA_VERSION = "rcsb-prefetch-manifest-v1"
ENRICHMENT_SCHEMA_VERSION = "rcsb-postbuild-enrichment-v1"

CHANGED_TABLES = {
    "source_files",
    "complexes",
    "ligand_instances",
    "pockets",
    "protein_chains",
    "affinity_reference_links",
    "affinity_reference_adjudications",
    "processing_issues",
    *ENRICHMENT_TABLES,
}
ENRICHMENT_INPUT_TABLES = CHANGED_TABLES | {"binding_measurements"}


def plan_rcsb_request(run_dir: Path, output: Path, *, overwrite: bool = False) -> dict:
    """Write a deterministic request file from a completed library's accepted rows."""
    run_dir = run_dir.resolve()
    manifest = _load_completed_manifest(run_dir)
    if output.exists() and not overwrite:
        raise FileExistsError(f"Request already exists: {output}; use --overwrite")
    pdb_ids = _accepted_pdb_ids(run_dir)
    if not pdb_ids:
        raise ValueError("The run contains no accepted PDB entries to enrich")
    metadata = {
        "schema_version": REQUEST_SCHEMA_VERSION,
        "source_run_id": manifest["run_id"],
        "source_manifest_sha256": sha256_file(run_dir / "manifest.json"),
        "source_identity": _parent_identity(manifest),
        "source_sidecar_logical_sha256": _relevant_sidecar_hashes(manifest),
        "requested_pdb_id_count": len(pdb_ids),
    }
    content = b"# " + canonical_json_bytes(metadata) + b"\npdb_id\n"
    content += "".join(f"{pdb_id}\n" for pdb_id in pdb_ids).encode("ascii")
    atomic_write_bytes(output, content)
    return {
        **metadata,
        "request_path": output.resolve().as_posix(),
        "request_sha256": sha256_file(output),
        "requested_pdb_ids": pdb_ids,
    }


def read_rcsb_request(path: Path) -> dict:
    """Parse and validate a deterministic RCSB request TSV."""
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise FileNotFoundError(f"Cannot read RCSB request: {path}") from error
    if len(lines) < 2 or not lines[0].startswith("# ") or lines[1] != "pdb_id":
        raise ValueError(f"Malformed RCSB request: {path}")
    try:
        metadata = json.loads(lines[0][2:])
    except ValueError as error:
        raise ValueError(f"Malformed RCSB request metadata: {path}") from error
    if metadata.get("schema_version") != REQUEST_SCHEMA_VERSION:
        raise ValueError(f"Unsupported RCSB request schema: {metadata.get('schema_version')}")
    required_metadata = {
        "source_run_id", "source_manifest_sha256", "source_identity",
        "source_sidecar_logical_sha256", "requested_pdb_id_count",
    }
    if not required_metadata <= set(metadata):
        raise ValueError("RCSB request metadata is incomplete")
    if (not isinstance(metadata["source_run_id"], str)
            or not isinstance(metadata["source_identity"], dict)
            or not isinstance(metadata["source_sidecar_logical_sha256"], dict)
            or not _is_sha256(metadata["source_manifest_sha256"])):
        raise ValueError("RCSB request provenance metadata is invalid")
    pdb_ids = [value.strip().lower() for value in lines[2:] if value.strip()]
    if not pdb_ids:
        raise ValueError("RCSB request contains no PDB identifiers")
    if any(len(value) != 4 or not value.isalnum() for value in pdb_ids):
        raise ValueError("RCSB request contains an invalid PDB identifier")
    if pdb_ids != sorted(set(pdb_ids)):
        raise ValueError("RCSB request PDB identifiers must be unique and sorted")
    if metadata.get("requested_pdb_id_count") != len(pdb_ids):
        raise ValueError("RCSB request count differs from its metadata")
    return {
        **metadata,
        "request_path": path.resolve().as_posix(),
        "request_sha256": sha256_file(path),
        "requested_pdb_ids": pdb_ids,
    }


def prefetch_rcsb_request(
    request_path: Path,
    cache_dir: Path,
    config: BuildConfig,
    *,
    manifest_path: Path | None = None,
    refresh: bool = False,
    allow_partial: bool = False,
    progress: bool = True,
) -> dict:
    """Download requested mmCIF files, then independently verify the cache."""
    cache_dir = cache_dir.resolve()
    if config.pipeline.offline:
        raise ValueError("RCSB prefetch requires online mode")
    if config.paths.external_cache_dir.resolve() != cache_dir:
        raise ValueError("Prefetch cache directory differs from the resolved configuration")
    request = read_rcsb_request(request_path)
    requested = request["requested_pdb_ids"]
    _, _, download_failures = download_mmcif_files(
        requested, config, refresh=refresh, progress=progress
    )
    _, records, cache_failures, cache_fingerprint = resolve_cached_mmcif_files(
        requested,
        cache_dir,
        compressed=config.rcsb.download_compressed,
        progress=progress,
    )
    download_by_id = {item["pdb_id"].lower(): item for item in download_failures}
    failures = [download_by_id.get(item["pdb_id"], item) for item in cache_failures]
    result = {
        "schema_version": CACHE_MANIFEST_SCHEMA_VERSION,
        "request": request,
        "cache_fingerprint": cache_fingerprint,
        "download_compressed": config.rcsb.download_compressed,
        "record_count": len(records),
        "records": records,
        "failure_count": len(failures),
        "failures": failures,
        "completed_at_utc": utc_now(),
    }
    destination = manifest_path or (
        cache_dir / "manifests/rcsb_mmcif" / f"{request['request_sha256']}.json"
    )
    atomic_write_bytes(destination, canonical_json_bytes(result) + b"\n")
    result["manifest_path"] = destination.resolve().as_posix()
    if failures and not allow_partial:
        missing = ", ".join(item["pdb_id"] for item in failures[:10])
        suffix = " ..." if len(failures) > 10 else ""
        raise ValueError(
            f"RCSB prefetch is incomplete ({len(failures)} missing: {missing}{suffix}); "
            f"manifest written to {destination}"
        )
    return result


def enrich_library_from_cache(
    parent_run_dir: Path,
    cache_dir: Path,
    config: BuildConfig,
    *,
    output_root: Path | None = None,
    allow_partial: bool = False,
    resume: bool = False,
    overwrite_run: bool = False,
    progress: bool = True,
) -> Path:
    """Create an immutable derived run enriched strictly from a validated local cache."""
    parent_run_dir = parent_run_dir.resolve()
    cache_dir = cache_dir.resolve()
    if not config.pipeline.offline:
        raise ValueError("Cache-only RCSB enrichment requires offline mode")
    if config.paths.external_cache_dir.resolve() != cache_dir:
        raise ValueError("Enrichment cache directory differs from the resolved configuration")
    parent_manifest = _load_completed_manifest(parent_run_dir)
    parent_errors = validate_run(parent_run_dir, config, progress=progress)
    if parent_errors:
        raise ValueError("Parent run validation failed: " + "; ".join(parent_errors))
    _require_unenriched_parent(parent_run_dir)
    pdb_ids = _accepted_pdb_ids(parent_run_dir)
    cached, cache_records, cache_failures, cache_fingerprint = resolve_cached_mmcif_files(
        pdb_ids,
        cache_dir,
        compressed=config.rcsb.download_compressed,
        progress=progress,
    )
    if cache_failures and not allow_partial:
        missing = ", ".join(item["pdb_id"] for item in cache_failures[:10])
        suffix = " ..." if len(cache_failures) > 10 else ""
        raise ValueError(f"RCSB cache is incomplete ({len(cache_failures)} missing: {missing}{suffix})")
    if not cached:
        raise ValueError("RCSB cache contains no valid requested mmCIF files")

    parent_identity = _parent_identity(parent_manifest)
    enrichment_spec = {
        "schema_version": ENRICHMENT_SCHEMA_VERSION,
        "parent_identity": parent_identity,
        "source_sidecar_logical_sha256": _relevant_sidecar_hashes(parent_manifest),
        "cache_fingerprint": cache_fingerprint,
        "allow_partial": allow_partial,
        "download_compressed": config.rcsb.download_compressed,
    }
    enrichment_identity = canonical_json_hash(enrichment_spec)
    run_id = f"{parent_run_dir.name}-rcsb-v1-{enrichment_identity[:12]}"
    destination_root = (output_root or parent_run_dir.parent).resolve()
    run_dir = destination_root / run_id
    if run_dir.exists():
        existing = _read_manifest(run_dir)
        compatible = existing and existing.get("enrichment_identity") == enrichment_identity
        if resume and compatible and existing.get("status") == "complete":
            return run_dir
        if overwrite_run:
            backup = run_dir.with_name(
                f"{run_dir.name}.backup-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
            )
            os.replace(run_dir, backup)
        else:
            raise FileExistsError(f"Derived run already exists: {run_dir}; use --resume or --overwrite-run")

    run_dir.mkdir(parents=True)
    _clone_artifacts(parent_run_dir / "sidecars", run_dir / "sidecars")
    _clone_artifacts(parent_run_dir / "lmdb", run_dir / "lmdb")
    manifest = _derived_manifest(
        parent_manifest,
        parent_run_dir,
        run_id,
        config,
        enrichment_spec,
        enrichment_identity,
        cache_records,
        cache_failures,
        cache_fingerprint,
    )
    write_resolved_config(run_dir / "config.resolved.toml", config)
    write_manifest(run_dir, manifest)
    complete_stage(
        run_dir,
        manifest,
        "bootstrap-identity",
        "rcsb-postbuild-1",
        {"parent": canonical_json_hash(parent_identity), "enrichment": enrichment_identity},
        [run_dir / "config.resolved.toml"],
        progress=progress,
    )
    complete_stage(
        run_dir,
        manifest,
        "download-rcsb",
        "cache-only-1",
        {"cache_fingerprint": cache_fingerprint},
        [],
        progress=progress,
    )

    rows = {name: read_sidecar(run_dir / "sidecars", name) for name in ENRICHMENT_INPUT_TABLES}
    rows["source_files"].extend(_normalize_external_rows(cached_inventory(cache_records, cache_dir), config))
    new_issues = apply_rcsb_enrichment(
        rows,
        cached,
        config,
        failures=cache_failures,
        progress=progress,
        strict=not allow_partial,
    )
    rows["processing_issues"].extend(_issue_row(issue) for issue in new_issues)
    warning_counts = Counter(
        row["complex_id"] for row in rows["processing_issues"]
        if row["complex_id"] and row["severity"] == "warning"
    )
    error_counts = Counter(
        row["complex_id"] for row in rows["processing_issues"]
        if row["complex_id"] and row["severity"] in {"error", "fatal"}
    )
    for row in rows["complexes"]:
        row["warning_count"] = warning_counts[row["complex_id"]]
        row["error_count"] = error_counts[row["complex_id"]]
        if row["processing_status"] == "accepted" and row["warning_count"]:
            row["processing_status"] = "accepted_with_warnings"
    changed_results = write_sidecars(
        run_dir / "sidecars", rows, progress=progress, table_names=CHANGED_TABLES
    )
    sidecar_results = _sidecar_results_from_parent(parent_manifest, run_dir)
    sidecar_results.update(changed_results)
    manifest["sidecar_artifacts"] = _relative_sidecar_results(sidecar_results, run_dir)
    _update_lmdb_profile_metadata(run_dir, manifest, sidecar_results)

    complete_stage(
        run_dir,
        manifest,
        "map-structures",
        "rcsb-postbuild-1",
        {"cache_fingerprint": cache_fingerprint},
        [Path(changed_results[name]["path"]) for name in sorted(ENRICHMENT_TABLES - {
            "citations", "citation_authors", "pdb_citation_links"
        })],
        progress=progress,
    )
    complete_stage(
        run_dir,
        manifest,
        "enrich-citations",
        "rcsb-postbuild-1",
        {"cache_fingerprint": cache_fingerprint},
        [Path(changed_results[name]["path"]) for name in (
            "citations", "citation_authors", "pdb_citation_links",
            "affinity_reference_links", "affinity_reference_adjudications",
        )],
        progress=progress,
    )
    complete_stage(
        run_dir,
        manifest,
        "quality-control",
        "rcsb-postbuild-1",
        {"allow_partial": allow_partial},
        [Path(changed_results[name]["path"]) for name in ("complexes", "pockets", "processing_issues")],
        progress=progress,
    )
    complete_stage(
        run_dir,
        manifest,
        "write-sidecars",
        "rcsb-postbuild-1",
        {"parent_run_id": parent_manifest["run_id"]},
        [Path(item["path"]) for item in sidecar_results.values()],
        progress=progress,
    )
    complete_stage(
        run_dir,
        manifest,
        "export-lmdb",
        "reused-lmdb-1",
        {"parent_run_id": parent_manifest["run_id"]},
        _lmdb_profile_outputs(run_dir, manifest),
        progress=progress,
    )

    errors = validate_run(run_dir, config, progress=progress)
    if errors:
        manifest["status"] = "validation_failed"
        manifest.setdefault("counts", {})["validation_errors"] = len(errors)
        write_manifest(run_dir, manifest)
        raise ValueError("Derived run validation failed: " + "; ".join(errors))
    complete_stage(run_dir, manifest, "validate", "rcsb-postbuild-1", {}, [], progress=progress)
    manifest["counts"] = {
        "selected_complexes": len(rows["complexes"]),
        "accepted_pockets": sum(
            row["processing_status"].startswith("accepted") for row in rows["pockets"]
        ),
        "rejected_complexes": sum(
            row["processing_status"] == "rejected" for row in rows["complexes"]
        ),
        "lmdb_records": len(read_sidecar(run_dir / "sidecars", "lmdb_records")),
        "rcsb_cached_records": len(cache_records),
        "rcsb_cache_failures": len(cache_failures),
    }
    report_results = generate_reports(run_dir, manifest)
    complete_stage(
        run_dir,
        manifest,
        "report",
        "rcsb-postbuild-1",
        {},
        [Path(item["path"]) for item in report_results.values()],
        progress=progress,
    )
    manifest["output_files"] = _artifact_inventory(run_dir, progress=progress)
    manifest["completed_at_utc"] = utc_now()
    manifest["status"] = "complete"
    write_manifest(run_dir, manifest)
    return run_dir


def _load_completed_manifest(run_dir: Path) -> dict:
    manifest = _read_manifest(run_dir)
    if manifest is None:
        raise FileNotFoundError(f"Missing or malformed manifest: {run_dir / 'manifest.json'}")
    if manifest.get("run_id") != run_dir.name:
        raise ValueError("Run-directory name differs from manifest run_id")
    if manifest.get("status") != "complete":
        raise ValueError("RCSB post-build enrichment requires a complete run with an LMDB")
    return manifest


def _read_manifest(run_dir: Path) -> dict | None:
    try:
        return json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _is_sha256(value) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(
        character in "0123456789abcdef" for character in value.lower()
    )


def _accepted_pdb_ids(run_dir: Path) -> list[str]:
    return sorted({
        row["pdb_id"] for row in read_sidecar(run_dir / "sidecars", "complexes")
        if row["processing_status"].startswith("accepted")
    })


def _parent_identity(manifest: dict) -> dict:
    return {
        key: manifest.get(key) for key in (
            "run_id",
            "pipeline_source",
            "semantic_config_hash",
            "source_fingerprint",
            "selection_fingerprint",
            "drugclip_library_contract_fingerprint",
        )
    }


def _relevant_sidecar_hashes(manifest: dict) -> dict:
    names = ("complexes", "pockets", "protein_chains", "ligand_instances", "binding_measurements")
    artifacts = manifest.get("sidecar_artifacts", {})
    missing = [name for name in names if not artifacts.get(name, {}).get("logical_sha256")]
    if missing:
        raise ValueError(f"Parent manifest lacks logical sidecar hashes: {', '.join(missing)}")
    return {name: artifacts[name]["logical_sha256"] for name in names}


def _require_unenriched_parent(run_dir: Path) -> None:
    populated = [name for name in sorted(ENRICHMENT_TABLES) if read_sidecar(run_dir / "sidecars", name)]
    rcsb_sources = [
        row for row in read_sidecar(run_dir / "sidecars", "source_files")
        if row["source_kind"] in {"rcsb_mmcif", "rcsb_api_cache"}
    ]
    if populated or rcsb_sources:
        detail = ", ".join(populated) if populated else "RCSB source files"
        raise ValueError(f"Parent run is already RCSB-enriched ({detail})")


def _clone_artifacts(source: Path, destination: Path) -> None:
    if not source.is_dir():
        raise FileNotFoundError(f"Missing artifact directory: {source}")
    for path in sorted(item for item in source.rglob("*") if item.is_file()):
        target = destination / path.relative_to(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        try:
            os.link(path, target)
        except OSError:
            shutil.copy2(path, target)


def _derived_manifest(
    parent: dict,
    parent_run_dir: Path,
    run_id: str,
    config: BuildConfig,
    enrichment_spec: dict,
    enrichment_identity: str,
    cache_records: list[dict],
    cache_failures: list[dict],
    cache_fingerprint: str,
) -> dict:
    manifest = copy.deepcopy(parent)
    commit, dirty = git_identity(config.project_root)
    parent_manifest_sha256 = sha256_file(parent_run_dir / "manifest.json")
    manifest.update({
        "run_id": run_id,
        "git_commit": commit,
        "git_dirty": dirty,
        "code_dirty_state_fingerprint": project_code_fingerprint(config.project_root),
        "semantic_config_hash": canonical_json_hash({
            "parent": parent["semantic_config_hash"], "enrichment": enrichment_spec
        }),
        "operational_config_hash": config.operational_hash,
        "source_fingerprint": canonical_json_hash({
            "parent": parent["source_fingerprint"], "rcsb_cache": cache_fingerprint
        }),
        "started_at_utc": utc_now(),
        "completed_at_utc": None,
        "status": "running",
        "configuration": config.as_dict(),
        "stage_statuses": {},
        "counts": {},
        "output_files": [],
        "derivation_kind": "post_build_rcsb_enrichment",
        "enrichment_identity": enrichment_identity,
        "parent_run": {
            "run_id": parent["run_id"],
            "manifest_sha256": parent_manifest_sha256,
            "identity": _parent_identity(parent),
            "source_sidecar_logical_sha256": _relevant_sidecar_hashes(parent),
        },
        "rcsb_cache": {
            "schema_version": "rcsb-cache-provenance-v1",
            "cache_fingerprint": cache_fingerprint,
            "download_compressed": config.rcsb.download_compressed,
            "record_count": len(cache_records),
            "failure_count": len(cache_failures),
            "records": cache_records,
            "failures": cache_failures,
            "network_access_used": False,
        },
    })
    return manifest


def _sidecar_results_from_parent(parent: dict, run_dir: Path) -> dict:
    artifacts = parent.get("sidecar_artifacts", {})
    if set(artifacts) != set(TABLES):
        raise ValueError("Parent manifest does not describe every required sidecar")
    return {
        name: {**item, "path": (run_dir / "sidecars" / f"{name}.parquet").as_posix()}
        for name, item in artifacts.items()
    }


def _relative_sidecar_results(results: dict, run_dir: Path) -> dict:
    return {
        name: {**item, "path": Path(item["path"]).relative_to(run_dir).as_posix()}
        for name, item in results.items()
    }


def _update_lmdb_profile_metadata(run_dir: Path, manifest: dict, sidecars: dict) -> None:
    for profile, recorded in manifest.get("lmdb_profiles", {}).items():
        lmdb_path = run_dir / recorded["path"]
        profile_path = lmdb_path.with_suffix(lmdb_path.suffix + ".profile.json")
        metadata = json.loads(profile_path.read_text(encoding="utf-8"))
        metadata["source_sidecar_hashes"] = {
            "pockets.parquet": sidecars["pockets"]["sha256"],
            "pocket_atoms.parquet": sidecars["pocket_atoms"]["sha256"],
        }
        metadata["post_build_enrichment"] = {
            "schema_version": ENRICHMENT_SCHEMA_VERSION,
            "parent_run_id": manifest["parent_run"]["run_id"],
            "rcsb_cache_fingerprint": manifest["rcsb_cache"]["cache_fingerprint"],
        }
        atomic_write_bytes(profile_path, canonical_json_bytes(metadata) + b"\n")
        recorded["profile_metadata_sha256"] = sha256_file(profile_path)
        recorded["source_sidecar_hashes"] = metadata["source_sidecar_hashes"]
        recorded["lmdb_physical_sha256"] = sha256_file(lmdb_path)


def _lmdb_profile_outputs(run_dir: Path, manifest: dict) -> list[Path]:
    outputs = []
    for recorded in manifest.get("lmdb_profiles", {}).values():
        lmdb_path = run_dir / recorded["path"]
        outputs.extend((lmdb_path, lmdb_path.with_suffix(lmdb_path.suffix + ".profile.json")))
    return outputs
