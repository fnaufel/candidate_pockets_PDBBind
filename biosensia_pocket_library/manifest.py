"""Run identity, manifest, checkpoint, and resolved-configuration helpers."""

from __future__ import annotations

import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .config import BuildConfig
from .constants import MANIFEST_SCHEMA_VERSION, PIPELINE_NAME, PIPELINE_VERSION
from .hashing import atomic_write_bytes, canonical_json_bytes, canonical_json_hash, sha256_file
from .quality import load_quality_effects
from .progress import file_progress


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def create_manifest(config: BuildConfig, run_id: str, source_fingerprint: str,
                    selection_fingerprint: str, contract: dict, index_summary: dict,
                    dataset_identity: dict) -> dict:
    commit, dirty = git_identity(Path.cwd())
    code_fingerprint = project_code_fingerprint(config.project_root)
    dependencies = {}
    for name in ("numpy", "rdkit", "lmdb", "pyarrow", "gemmi", "httpx", "tenacity", "tqdm", "scipy"):
        try:
            dependencies[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            dependencies[name] = None
    return {
        "manifest_schema_version": MANIFEST_SCHEMA_VERSION, "pipeline_name": PIPELINE_NAME,
        "pipeline_version": PIPELINE_VERSION, "git_commit": commit, "git_dirty": dirty,
        "code_dirty_state_fingerprint": code_fingerprint,
        "run_id": run_id, "semantic_config_hash": config.semantic_hash,
        "operational_config_hash": config.operational_hash, "source_fingerprint": source_fingerprint,
        "selection_fingerprint": selection_fingerprint,
        "drugclip_contract_fingerprint": contract["drugclip_contract_fingerprint"],
        "drugclip_library_contract_fingerprint": contract["drugclip_contract_fingerprint"],
        "started_at_utc": utc_now(), "completed_at_utc": None, "status": "running",
        "python_version": sys.version, "platform": platform.platform(), "dependency_versions": dependencies,
        "drugclip_contract": contract,
        "dataset": {**dataset_identity,
                    "index_file": f"{config.as_dict()['paths']['index_dir']}/INDEX_general_PL.2020R1.lst",
                    "index_declared_complex_count": index_summary.get("declared_count"),
                    "index_parsed_complex_count": index_summary.get("unique_complex_count"),
                    "discovered_complex_directory_count": None},
        "configuration": config.as_dict(),
        "secret_availability": {
            config.bibliography.contact_email_env: bool(os.environ.get(config.bibliography.contact_email_env)),
            config.bibliography.pubmed_api_key_env: bool(os.environ.get(config.bibliography.pubmed_api_key_env)),
        },
        "quality_rules": {"resolved": load_quality_effects(config.quality.rules_file),
                          "sha256": sha256_file(config.quality.rules_file) if config.quality.rules_file.is_file() else None},
        "stage_statuses": {}, "counts": {}, "output_files": [],
    }


def write_manifest(run_dir: Path, manifest: dict) -> None:
    atomic_write_bytes(run_dir / "manifest.json", canonical_json_bytes(manifest) + b"\n")


def write_resolved_config(path: Path, config: BuildConfig) -> None:
    lines: list[str] = []
    for section, values in config.as_dict().items():
        lines.append(f"[{section}]")
        for key, value in values.items():
            lines.append(f"{key} = {_toml_value(value)}")
        lines.append("")
    atomic_write_bytes(path, "\n".join(lines).encode("utf-8"))


def complete_stage(run_dir: Path, manifest: dict, stage: str, version: str,
                   inputs: dict[str, Any], outputs: list[Path], *, progress: bool = True) -> None:
    output_hashes = {}
    for path in outputs:
        if not path.is_file():
            continue
        with file_progress(path, description=f"Checkpointing {path.name}",
                           enabled=progress and path.stat().st_size > 50_000_000) as bar:
            output_hashes[path.relative_to(run_dir).as_posix()] = sha256_file(path, progress=bar)
    marker = {"stage": stage, "stage_version": version,
              "semantic_config_hash": manifest["semantic_config_hash"],
              "source_fingerprint": manifest["source_fingerprint"], "git_commit": manifest["git_commit"],
              "git_dirty": manifest["git_dirty"],
              "code_dirty_state_fingerprint": manifest["code_dirty_state_fingerprint"],
              "input_hashes": inputs, "output_hashes": output_hashes}
    marker_path = run_dir / "checkpoints" / stage / "complete.json"
    atomic_write_bytes(marker_path, canonical_json_bytes(marker) + b"\n")
    manifest["stage_statuses"][stage] = {"status": "complete", "marker": marker_path.relative_to(run_dir).as_posix()}
    write_manifest(run_dir, manifest)


def git_identity(root: Path) -> tuple[str | None, bool | None]:
    try:
        commit = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"], check=True,
                                capture_output=True, text=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "-C", str(root), "status", "--porcelain"], check=True,
                                    capture_output=True, text=True).stdout.strip())
        return commit, dirty
    except (OSError, subprocess.CalledProcessError):
        return None, None


def project_code_fingerprint(root: Path) -> str:
    included: list[tuple[str, str]] = []
    excluded_roots = {".git", ".venv", "data", "BioSensIA-DC", "__pycache__"}
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if not path.is_file() or any(part in excluded_roots for part in relative.parts):
            continue
        if path.suffix in {".py", ".toml", ".md", ".qmd"} or path.name in {"uv.lock"}:
            included.append((relative.as_posix(), sha256_file(path)))
    return canonical_json_hash(included)


def _toml_value(value: Any) -> str:
    if value is None:
        return '"auto"'
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return repr(value)
    if isinstance(value, str):
        return json.dumps(value, ensure_ascii=False)
    if isinstance(value, list):
        return "[" + ", ".join(_toml_value(item) for item in value) + "]"
    if isinstance(value, dict):
        return "{" + ", ".join(f"{key} = {_toml_value(item)}" for key, item in sorted(value.items())) + "}"
    return json.dumps(str(value))
