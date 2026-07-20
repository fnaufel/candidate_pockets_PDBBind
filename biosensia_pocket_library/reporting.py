"""Deterministic machine-readable and Markdown build reports."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from .hashing import atomic_write_bytes, canonical_json_bytes, sha256_file
from .sidecars import read_sidecar


def generate_reports(run_dir: Path, manifest: dict) -> dict[str, dict]:
    sidecars = run_dir / "sidecars"
    report_dir = run_dir / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    complexes = read_sidecar(sidecars, "complexes")
    pockets = read_sidecar(sidecars, "pockets")
    issues = read_sidecar(sidecars, "processing_issues")
    comparisons = read_sidecar(sidecars, "pocket_comparisons")
    adjudications = read_sidecar(sidecars, "affinity_reference_adjudications")
    lmdb_rows = read_sidecar(sidecars, "lmdb_records")
    summary = {
        "run_id": manifest["run_id"], "declared_complex_count": manifest["dataset"]["index_declared_complex_count"],
        "selected_complex_count": len(complexes),
        "successfully_processed_count": sum(row["processing_status"].startswith("accepted") for row in complexes),
        "rejected_count": sum(row["processing_status"] == "rejected" for row in complexes),
        "geometry_quality_counts": dict(sorted(Counter(row["geometry_quality_tier"] for row in complexes).items())),
        "pocket_comparison_quality_counts": dict(sorted(Counter(row["pocket_comparison_quality"] for row in complexes).items())),
        "structure_mapping_quality_counts": dict(sorted(Counter(row["structure_mapping_quality"] for row in complexes).items())),
        "bibliography_quality_counts": dict(sorted(Counter(row["bibliography_quality"] for row in complexes).items())),
        "ligand_sdf_success_count": sum(row["selected_source_format"] == "sdf" for row in read_sidecar(sidecars, "ligand_instances")),
        "ligand_mol2_fallback_count": sum(row["selected_source_format"] == "mol2" for row in read_sidecar(sidecars, "ligand_instances")),
        "both_ligand_formats_failed_count": sum(row["issue_code"] == "LIGAND_BOTH_FORMATS_FAILED" for row in issues),
        "cropped_pocket_count": sum(row["crop_applied"] for row in pockets),
        "multi_chain_pocket_count": sum(row["contributing_chain_count"] > 1 for row in pockets),
        "reference_status_counts": dict(sorted(Counter(row["reference_status"] for row in adjudications).items())),
        "lmdb_profile_counts": dict(sorted(Counter(row["library_profile"] for row in lmdb_rows).items())),
        "source_fingerprint": manifest["source_fingerprint"],
        "drugclip_contract_fingerprint": manifest["drugclip_contract_fingerprint"],
        "logical_digest_excluded_columns": {
            "source_files": ["modified_time_utc", "downloaded_at_utc"],
            "processing_issues": ["created_at_utc"],
            "affinity_reference_links": ["verified_at_utc"],
            "affinity_reference_adjudications": ["adjudicated_at_utc"],
        },
    }
    summary_json = report_dir / "build_summary.json"
    atomic_write_bytes(summary_json, canonical_json_bytes(summary) + b"\n")
    markdown = _markdown(summary)
    summary_md = report_dir / "build_summary.md"
    atomic_write_bytes(summary_md, markdown.encode("utf-8"))
    _count_table(report_dir / "quality_counts.parquet", "dimension_value",
                 Counter(row["geometry_quality_tier"] for row in complexes))
    _count_table(report_dir / "failure_counts.parquet", "issue_code",
                 Counter(row["issue_code"] for row in issues if row["severity"] in {"error", "fatal"}))
    _count_table(report_dir / "pocket_size_distribution.parquet", "exported_atom_count",
                 Counter(str(row["exported_atom_count"]) for row in pockets))
    _count_table(report_dir / "pocket_comparison_distribution.parquet", "comparison_status",
                 Counter(row["comparison_status"] for row in comparisons))
    for filename, column, table_name, field in (
        ("chain_mapping_status_counts.parquet", "mapping_status", "chain_mapping_candidates", "mapping_status"),
        ("uniprot_mapping_status_counts.parquet", "mapping_status", "chain_uniprot_mappings", "mapping_status"),
        ("reference_status_counts.parquet", "reference_status", "affinity_reference_adjudications", "reference_status"),
    ):
        _count_table(report_dir / filename, column, Counter(row[field] for row in read_sidecar(sidecars, table_name)))
    return {path.name: {"path": path.as_posix(), "sha256": sha256_file(path)} for path in sorted(report_dir.iterdir()) if path.is_file()}


def _count_table(path: Path, column: str, counts: Counter) -> None:
    schema = pa.schema([(column, pa.string()), ("count", pa.int64())],
                       metadata={b"schema_name": path.stem.encode(), b"semantic_version": b"1.0.0"})
    table = pa.Table.from_pylist([{column: str(key), "count": value} for key, value in sorted(counts.items())], schema=schema)
    pq.write_table(table, path, compression="zstd")


def _markdown(summary: dict) -> str:
    lines = ["# Candidate pocket library build summary", "", f"Run ID: `{summary['run_id']}`", ""]
    for key, value in summary.items():
        if key == "run_id":
            continue
        label = key.replace("_", " ").capitalize()
        lines.append(f"- {label}: `{json.dumps(value, sort_keys=True) if isinstance(value, dict) else value}`")
    lines.extend(["", "BioSensIA-DC embedding caches must be invalidated or namespaced by the LMDB logical checksum.", ""])
    return "\n".join(lines)
