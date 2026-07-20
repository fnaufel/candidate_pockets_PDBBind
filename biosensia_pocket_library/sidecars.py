"""Atomic explicit-schema Parquet I/O and logical validation."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from .hashing import canonical_json_hash, sha256_file
from .progress import track
from .schemas import TABLES


def write_sidecars(directory: Path, rows_by_table: dict[str, list[dict[str, Any]]], *, progress: bool = True) -> dict[str, dict]:
    directory.mkdir(parents=True, exist_ok=True)
    results: dict[str, dict] = {}
    for name in track(sorted(TABLES), description="Writing sidecars", total=len(TABLES), enabled=progress):
        rows = rows_by_table.get(name, [])
        spec = TABLES[name]
        canonical = sorted(rows, key=lambda row: tuple(_sort_value(row.get(key)) for key in spec.sort_by))
        field_names = set(spec.schema.names)
        extras = set().union(*(set(row) for row in canonical)) - field_names if canonical else set()
        if extras:
            raise ValueError(f"Unknown columns for {name}: {sorted(extras)}")
        normalized = [{field.name: row.get(field.name) for field in spec.schema} for row in canonical]
        table = pa.Table.from_pylist(normalized, schema=spec.schema)
        destination = directory / f"{name}.parquet"
        temporary = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
        try:
            pq.write_table(table, temporary, compression="zstd", version="2.6")
            os.replace(temporary, destination)
        finally:
            temporary.unlink(missing_ok=True)
        logical_rows = [{field.name: row.get(field.name) for field in spec.schema
                         if field.name not in spec.volatile_columns} for row in canonical]
        results[name] = {"path": destination.as_posix(), "row_count": len(canonical),
                         "sha256": sha256_file(destination), "logical_sha256": canonical_json_hash(logical_rows)}
    return results


def read_sidecar(directory: Path, name: str) -> list[dict]:
    return pq.read_table(directory / f"{name}.parquet", schema=TABLES[name].schema).to_pylist()


def validate_sidecars(directory: Path) -> list[str]:
    errors: list[str] = []
    loaded: dict[str, list[dict]] = {}
    for name, spec in TABLES.items():
        path = directory / f"{name}.parquet"
        if not path.is_file():
            errors.append(f"Missing sidecar {name}")
            continue
        actual = pq.read_schema(path)
        if actual.remove_metadata() != spec.schema.remove_metadata():
            errors.append(f"Schema mismatch: {name}")
        rows = pq.read_table(path).to_pylist()
        loaded[name] = rows
        for column, allowed in (spec.allowed_enums or {}).items():
            invalid_values = {row[column] for row in rows if row[column] is not None and row[column] not in allowed}
            if invalid_values:
                errors.append(f"Unknown enum values in {name}.{column}: {sorted(invalid_values)}")
        keys = [tuple(row.get(key) for key in spec.primary_key) for row in rows]
        if len(keys) != len(set(keys)):
            errors.append(f"Duplicate primary key: {name}")
        sorted_keys = [tuple(_sort_value(row.get(key)) for key in spec.sort_by) for row in rows]
        if sorted_keys != sorted(sorted_keys):
            errors.append(f"Noncanonical row order: {name}")
    for name, spec in TABLES.items():
        for local, target_table, target in spec.foreign_keys:
            if name not in loaded or target_table not in loaded:
                continue
            valid = {row[target] for row in loaded[target_table]}
            invalid = [row[local] for row in loaded[name] if row.get(local) is not None and row[local] not in valid]
            if invalid:
                errors.append(f"Foreign-key violation: {name}.{local} -> {target_table}.{target}")
    return errors


def _sort_value(value: Any):
    return (value is None, value if not isinstance(value, list) else tuple(value))
