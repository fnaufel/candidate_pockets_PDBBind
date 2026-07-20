"""Canonical hashing and atomic-file helpers."""

from __future__ import annotations

import hashlib
import json
import os
import struct
from pathlib import Path
from typing import Any, Iterable

import numpy as np


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(
    path: str | Path,
    chunk_size: int = 1024 * 1024,
    progress: Any | None = None,
) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
            if progress is not None:
                progress.update(len(chunk))
    return digest.hexdigest()


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
        default=_json_default,
    ).encode("utf-8")


def canonical_json_hash(value: Any) -> str:
    return sha256_bytes(canonical_json_bytes(value))


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def length_frame(parts: Iterable[bytes]) -> bytes:
    output = bytearray()
    for part in parts:
        output.extend(struct.pack("<Q", len(part)))
        output.extend(part)
    return bytes(output)


def normalized_array_bytes(array: np.ndarray, dtype: str) -> bytes:
    normalized = np.asarray(array, dtype=np.dtype(dtype), order="C").copy()
    if np.issubdtype(normalized.dtype, np.floating):
        normalized[normalized == 0] = 0
        if not np.isfinite(normalized).all():
            raise ValueError("Cannot hash nonfinite array")
    return normalized.tobytes(order="C")


def atomic_write_bytes(path: str | Path, data: bytes) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def stable_id(prefix: str, *parts: Any, length: int = 24) -> str:
    digest = canonical_json_hash(parts)
    return f"{prefix}:{digest[:length]}"
