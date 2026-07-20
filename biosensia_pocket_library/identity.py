"""Verification of the mixed-version PDBbind special-distribution identity."""

from __future__ import annotations

import re
from pathlib import Path

from .exceptions import SourceIntegrityError
from .hashing import sha256_file


def verify_dataset_identity(index_dir: Path, declared_pl_count: int | None) -> dict:
    readme = index_dir / "README"
    if not readme.is_file():
        raise SourceIntegrityError("PDBbind distribution README is missing")
    text = readme.read_text(encoding="utf-8", errors="replace")
    count_match = re.search(r"(\d+)\s+protein-ligand complexes", text, re.I)
    readme_count = int(count_match.group(1)) if count_match else None
    if "v2020" not in text or "v2024" not in text or not re.search(r"re[- ]?process", text, re.I):
        raise SourceIntegrityError("README does not identify the v2020/v2024-reprocessed special distribution")
    if declared_pl_count is None or readme_count is None or declared_pl_count != readme_count:
        raise SourceIntegrityError(f"PL count identity mismatch (README={readme_count}, index={declared_pl_count})")
    exact_date = re.search(r"Latest update:\s*Aug(?:ust)?\s+(?:4th|4),?\s+2025", text, re.I)
    month_date = re.search(r"Latest update:\s*Aug\s+2025", text, re.I)
    if not exact_date and not month_date:
        raise SourceIntegrityError("README does not declare the August 2025 index revision")
    return {"name": "PDBbind", "distribution_id": "pdbbind-2020-v2024p-20250804",
            "distribution_label": "PDBbind v2020/v2024-reprocessed special distribution",
            "nominal_complex_set_version": "2020", "structure_processing_version": "2024",
            "index_revision_date": "2025-08-04", "source_readme_sha256": sha256_file(readme)}
