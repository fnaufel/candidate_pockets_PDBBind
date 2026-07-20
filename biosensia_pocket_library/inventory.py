"""Deterministic complex discovery and checksum inventory."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable

from .constants import EXPECTED_COMPLEX_SUFFIXES, INDEX_FILES, YEAR_RANGES
from .hashing import sha256_file, stable_id
from .models import IndexRecord, ProcessingIssue
from .progress import file_progress, track


def discover_complex_directories(root: Path) -> dict[str, list[Path]]:
    found: dict[str, list[Path]] = {}
    for year_range in YEAR_RANGES:
        parent = root / year_range
        if not parent.is_dir():
            continue
        with os.scandir(parent) as entries:
            for entry in entries:
                if entry.is_dir(follow_symlinks=True) and len(entry.name) == 4:
                    found.setdefault(entry.name.lower(), []).append(Path(entry.path).resolve())
    return found


def inventory_sources(
    index_dir: Path, complex_root: Path, records: Iterable[IndexRecord], *, progress: bool = True
) -> tuple[list[dict], dict[str, dict], list[ProcessingIssue]]:
    records = list(records)
    directories = discover_complex_directories(complex_root)
    rows: list[dict] = []
    complexes: dict[str, dict] = {}
    issues: list[ProcessingIssue] = []
    for name in track(INDEX_FILES, description="Hashing index files", total=len(INDEX_FILES), enabled=progress):
        path = index_dir / name
        if path.is_file():
            rows.append(_file_row(path, "index", None, progress))
        else:
            issues.append(ProcessingIssue("inventory", "MISSING_INDEX_FILE", "fatal", f"Missing {name}"))
    for record in track(records, description="Inventorying complexes", total=len(records), enabled=progress):
        locations = directories.get(record.pdb_id, [])
        info = {"complex_directory": None, "files": {}}
        if not locations:
            issues.append(ProcessingIssue("inventory", "MISSING_COMPLEX_DIRECTORY", "fatal",
                                          "No complex directory", complex_id=record.complex_id))
            for suffix in EXPECTED_COMPLEX_SUFFIXES:
                expected = f"{record.pdb_id}{suffix}"
                virtual_path = complex_root / "_missing" / record.pdb_id / expected
                missing_id = stable_id("file", "missing", virtual_path.as_posix())
                key = expected.removeprefix(record.pdb_id).lstrip("_").replace(".", "_")
                rows.append({"source_file_id": missing_id, "complex_id": record.complex_id,
                             "role": "missing", "path": virtual_path.as_posix(),
                             "size_bytes": 0, "sha256": None})
                info["files"][key] = missing_id
            complexes[record.complex_id] = info
            continue
        if len(locations) > 1:
            issues.append(ProcessingIssue("inventory", "DUPLICATE_COMPLEX_DIRECTORY", "fatal",
                                          "Multiple complex directories", complex_id=record.complex_id,
                                          details={"paths": [path.as_posix() for path in locations]}))
            for directory in sorted(locations):
                for path in sorted(item for item in directory.iterdir() if item.is_file()):
                    rows.append(_file_row(path, "duplicate_complex", record.complex_id, progress))
            complexes[record.complex_id] = info
            continue
        directory = locations[0]
        info["complex_directory"] = directory.as_posix()
        expected_names = {f"{record.pdb_id}{suffix}" for suffix in EXPECTED_COMPLEX_SUFFIXES}
        actual = {path.name: path for path in directory.iterdir() if path.is_file()}
        lower_actual = {name.lower(): name for name in actual}
        for expected in sorted(expected_names):
            path = actual.get(expected)
            if path is None and expected.lower() in lower_actual:
                path = actual[lower_actual[expected.lower()]]
                issues.append(ProcessingIssue("inventory", "SOURCE_FILENAME_CASE_MISMATCH", "warning",
                                              f"Expected {expected}, got {path.name}", record.complex_id))
            key = expected.removeprefix(record.pdb_id).lstrip("_").replace(".", "_")
            if path is None:
                missing_path = directory / expected
                missing_id = stable_id("file", "missing", missing_path.as_posix())
                rows.append({"source_file_id": missing_id, "complex_id": record.complex_id,
                             "role": "missing", "path": missing_path.as_posix(),
                             "size_bytes": 0, "sha256": None})
                severity = "fatal" if expected.endswith("_protein.pdb") else "warning"
                issues.append(ProcessingIssue("inventory", "MISSING_EXPECTED_FILE", severity,
                                              f"Missing {expected}", record.complex_id,
                                              details={"expected_name": expected}))
                info["files"][key] = missing_id
                continue
            row = _file_row(path, "complex", record.complex_id, progress)
            rows.append(row)
            info["files"][key] = row["source_file_id"]
            info["files"][f"{key}_path"] = path.as_posix()
            if row["size_bytes"] == 0:
                issues.append(ProcessingIssue("inventory", "EMPTY_SOURCE_FILE", "error",
                                              f"Empty {path.name}", record.complex_id,
                                              source_file_id=row["source_file_id"]))
        for name, path in sorted(actual.items()):
            if name not in expected_names and name.lower() not in {x.lower() for x in expected_names}:
                row = _file_row(path, "extra", record.complex_id, progress)
                rows.append(row)
                issues.append(ProcessingIssue("inventory", "EXTRA_COMPLEX_FILE", "info",
                                              f"Extra file {name}", record.complex_id,
                                              source_file_id=row["source_file_id"]))
        complexes[record.complex_id] = info
    return sorted(rows, key=lambda row: (row["path"], row["source_file_id"])), complexes, issues


def _file_row(path: Path, role: str, complex_id: str | None, show_progress: bool) -> dict:
    with file_progress(path, description=f"SHA-256 {path.name}", enabled=show_progress and path.stat().st_size > 50_000_000) as bar:
        digest = sha256_file(path, progress=bar)
    return {"source_file_id": stable_id("file", digest, path.name), "complex_id": complex_id,
            "role": role, "path": path.resolve().as_posix(), "size_bytes": path.stat().st_size,
            "sha256": digest}
