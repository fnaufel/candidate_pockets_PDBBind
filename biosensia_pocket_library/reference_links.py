"""Canonical aggregation for affinity-reference candidate evidence."""

from __future__ import annotations

from collections import defaultdict


_STATUS_PRIORITY = {
    "not_attempted": 0,
    "no_reference_available": 1,
    "reference_unresolved": 2,
    "structural_reference_only": 3,
    "probable_structural_reference": 4,
    "probable_affinity_reference": 5,
    "exact_affinity_reference": 6,
    "conflicting_references": 7,
}


def merge_affinity_reference_links(rows: list[dict]) -> list[dict]:
    """Return one deterministic row per declared primary key."""
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in rows:
        grouped[(row["measurement_id"], row["citation_id"])].append(row)
    merged = []
    for key in sorted(grouped):
        values = grouped[key]
        complex_ids = {row["complex_id"] for row in values}
        if len(complex_ids) != 1:
            raise ValueError(f"Reference candidate {key} maps to multiple complexes")
        chosen = max(values, key=lambda row: (
            row["automatic_or_manual"] == "manual",
            _STATUS_PRIORITY.get(row["candidate_status"], -1),
            row["confidence"] if row["confidence"] is not None else -1.0,
            row["candidate_status"],
        ))
        notes = sorted({row["evidence_note"] for row in values if row.get("evidence_note")})
        verified_by = sorted({row["verified_by"] for row in values if row.get("verified_by")})
        verified_at = [row["verified_at_utc"] for row in values if row.get("verified_at_utc") is not None]
        merged.append({
            **chosen,
            "evidence_sources": sorted({source for row in values for source in row["evidence_sources"]}),
            "evidence_note": " | ".join(notes) if notes else None,
            "automatic_or_manual": "manual" if any(
                row["automatic_or_manual"] == "manual" for row in values
            ) else "automatic",
            "verified_by": " | ".join(verified_by) if verified_by else None,
            "verified_at_utc": max(verified_at) if verified_at else None,
        })
    return merged
