"""Consistent, disableable progress reporting for long operations."""

from __future__ import annotations

import sys
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TypeVar

from tqdm.auto import tqdm

T = TypeVar("T")


def track(
    values: Iterable[T], *, description: str, total: int | None = None, enabled: bool = True,
    unit: str = "item",
) -> Iterator[T]:
    yield from tqdm(values, desc=description, total=total, disable=not enabled, unit=unit)


def progress_message(value: str, *, enabled: bool = True) -> None:
    if enabled:
        tqdm.write(value, file=sys.stderr)


@contextmanager
def file_progress(path: Path, *, description: str, enabled: bool = True):
    bar = tqdm(
        total=path.stat().st_size,
        desc=description,
        unit="B",
        unit_scale=True,
        disable=not enabled,
    )
    try:
        yield bar
    finally:
        bar.close()
