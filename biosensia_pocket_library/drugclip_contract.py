"""Verification and fingerprinting of the linked DrugCLIP contract."""

from __future__ import annotations

import subprocess
from pathlib import Path

from .config import BuildConfig
from .exceptions import SourceIntegrityError
from .hashing import canonical_json_hash, sha256_file
from .progress import file_progress


def verify_drugclip_contract(config: BuildConfig, *, progress: bool = True) -> dict:
    drugclip = config.paths.drugclip_dir
    biosensia = config.paths.biosensia_root
    if drugclip is None or not drugclip.is_dir():
        raise SourceIntegrityError(f"DrugCLIP directory not found: {drugclip}")
    candidates = {
        "dictionary": config.paths.drugclip_dictionary,
        "checkpoint": config.paths.drugclip_checkpoint,
        "task": drugclip / "unimol/tasks/drugclip.py",
        "helper": (biosensia / "lmdb_helpers.py") if biosensia else None,
    }
    loader_paths = (
        drugclip / "unimol/data/lmdb_dataset.py",
        drugclip / "unimol/data/affinity_dataset.py",
        drugclip / "unimol/data/remove_hydrogen_dataset.py",
        drugclip / "unimol/data/cropping_dataset.py",
        drugclip / "unimol/data/normalize_dataset.py",
    )
    required = ("dictionary", "task", "helper")
    for name in required:
        if candidates[name] is None or not candidates[name].is_file():
            raise SourceIntegrityError(f"DrugCLIP contract file missing ({name}): {candidates[name]}")
    for path in loader_paths:
        if not path.is_file():
            raise SourceIntegrityError(f"DrugCLIP loader contract file missing: {path}")
    hashes: dict[str, str | None] = {}
    for name, path in candidates.items():
        if path is None or not path.is_file():
            hashes[f"{name}_sha256"] = None
            continue
        with file_progress(path, description=f"Hashing DrugCLIP {name}",
                           enabled=progress and path.stat().st_size > 50_000_000) as bar:
            hashes[f"{name}_sha256"] = sha256_file(path, progress=bar)
    loader_hashes = {path.relative_to(drugclip).as_posix(): sha256_file(path) for path in loader_paths}
    hashes["loader_file_sha256"] = loader_hashes
    hashes["loader_sha256"] = canonical_json_hash(loader_hashes)
    tokens = read_dictionary(candidates["dictionary"])
    revision, dirty = _git_identity(biosensia)
    try:
        link_path = drugclip.relative_to(config.project_root).as_posix()
    except ValueError:
        link_path = drugclip.as_posix()
    result = {"link_path": link_path, "biosensia_commit": revision,
              "biosensia_dirty": dirty, "dictionary_tokens": sorted(tokens), **hashes}
    fingerprint_fields = {key: result.get(key) for key in (
        "biosensia_commit", "biosensia_dirty", "task_sha256", "loader_sha256",
        "helper_sha256", "dictionary_sha256", "checkpoint_sha256")}
    result["drugclip_contract_fingerprint"] = canonical_json_hash(fingerprint_fields)
    return result


def read_dictionary(path: Path) -> set[str]:
    tokens: set[str] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped:
            tokens.add(stripped.split()[0])
    return tokens


def _git_identity(root: Path | None) -> tuple[str | None, bool | None]:
    if root is None:
        return None, None
    try:
        revision = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                                  check=True, capture_output=True, text=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "-C", str(root), "status", "--porcelain"],
                                    check=True, capture_output=True, text=True).stdout.strip())
        return revision, dirty
    except (OSError, subprocess.CalledProcessError):
        return None, None
