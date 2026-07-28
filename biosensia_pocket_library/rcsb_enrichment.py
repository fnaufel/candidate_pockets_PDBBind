"""Shared application of verified RCSB mmCIF enrichment to library rows."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from .config import BuildConfig
from .models import ProcessingIssue
from .progress import track
from .rcsb import enrich_from_mmcif


ENRICHMENT_TABLES = {
    "chain_mapping_candidates",
    "chain_uniprot_mappings",
    "chain_uniprot_mapping_segments",
    "rcsb_ligand_mapping_candidates",
    "citations",
    "citation_authors",
    "pdb_citation_links",
}


def apply_rcsb_enrichment(
    rows: dict[str, list[dict]],
    cached: dict[str, Path],
    config: BuildConfig,
    *,
    failures: list[dict] | None = None,
    progress: bool = True,
    strict: bool = False,
) -> list[ProcessingIssue]:
    """Apply already-resolved mmCIF files without performing network I/O."""
    issues: list[ProcessingIssue] = []
    accepted_pdb_ids = sorted({
        row["pdb_id"] for row in rows["complexes"]
        if row["processing_status"].startswith("accepted")
    })
    failures_by_id = {item["pdb_id"].lower(): item for item in failures or []}
    chains_by_pdb: dict[str, list[dict]] = defaultdict(list)
    complexes_by_pdb: dict[str, list[dict]] = defaultdict(list)
    pockets_by_pdb: dict[str, list[dict]] = defaultdict(list)
    ligands_by_pdb: dict[str, list[dict]] = defaultdict(list)
    for row in rows["protein_chains"]:
        chains_by_pdb[row["pdb_id"]].append(row)
    for row in rows["complexes"]:
        complexes_by_pdb[row["pdb_id"]].append(row)
    for row in rows["pockets"]:
        pockets_by_pdb[row["pdb_id"]].append(row)
    for row in rows["ligand_instances"]:
        ligands_by_pdb[row["pdb_id"]].append(row)

    for pdb_id in track(
        accepted_pdb_ids,
        description="Applying cached RCSB enrichment",
        total=len(accepted_pdb_ids),
        enabled=progress,
    ):
        if pdb_id not in cached:
            failure = failures_by_id.get(pdb_id, {
                "pdb_id": pdb_id,
                "error": "No valid cached RCSB mmCIF file",
                "exception_type": "RcsbCacheMiss",
            })
            if strict:
                raise ValueError(f"RCSB cache is incomplete for {pdb_id}: {failure['error']}")
            _mark_unavailable(pdb_id, complexes_by_pdb, pockets_by_pdb)
            issues.extend(_failure_issues(pdb_id, failure, complexes_by_pdb, pockets_by_pdb))
            continue
        labels = {row["complex_id"]: row.get("ligand_label") for row in complexes_by_pdb[pdb_id]}
        ligand_inputs = [(row, labels.get(row["complex_id"])) for row in ligands_by_pdb[pdb_id]]
        try:
            enriched = enrich_from_mmcif(
                pdb_id,
                cached[pdb_id],
                config.rcsb.download_compressed,
                chains_by_pdb[pdb_id],
                ligand_inputs,
            )
        except Exception as error:
            if strict:
                raise ValueError(f"Failed to parse cached RCSB mmCIF for {pdb_id}: {error}") from error
            failure = {
                "pdb_id": pdb_id,
                "error": str(error),
                "exception_type": type(error).__name__,
            }
            _mark_unavailable(pdb_id, complexes_by_pdb, pockets_by_pdb)
            issues.extend(_failure_issues(pdb_id, failure, complexes_by_pdb, pockets_by_pdb))
            continue
        for name, values in enriched.items():
            rows[name].extend(values)
        _apply_mapping_statuses(
            pdb_id, enriched, chains_by_pdb, complexes_by_pdb, pockets_by_pdb, ligands_by_pdb
        )

    _normalize_enrichment_rows(rows)
    _apply_structural_citation_evidence(rows)
    return issues


def _apply_mapping_statuses(
    pdb_id: str,
    enriched: dict[str, list[dict]],
    chains_by_pdb: dict[str, list[dict]],
    complexes_by_pdb: dict[str, list[dict]],
    pockets_by_pdb: dict[str, list[dict]],
    ligands_by_pdb: dict[str, list[dict]],
) -> None:
    selected = {
        (item["pocket_instance_id"], item["pdbbind_auth_chain_id"])
        for item in enriched["chain_mapping_candidates"] if item["selected"]
    }
    candidates = {
        (item["pocket_instance_id"], item["pdbbind_auth_chain_id"])
        for item in enriched["chain_mapping_candidates"]
    }
    for chain in chains_by_pdb[pdb_id]:
        key = (chain["pocket_instance_id"], chain["pdbbind_auth_chain_id"])
        chain["rcsb_mapping_status"] = (
            "exact_identifier_match" if key in selected else
            "ambiguous" if key in candidates else "unresolved"
        )
    quality_by_complex: dict[str, str] = {}
    for pocket in pockets_by_pdb[pdb_id]:
        required = {
            (chain["pocket_instance_id"], chain["pdbbind_auth_chain_id"])
            for chain in chains_by_pdb[pdb_id]
            if chain["pocket_instance_id"] == pocket["pocket_instance_id"]
        }
        quality = (
            "exact" if required and required <= selected else
            "ambiguous" if required & candidates else "unresolved"
        )
        pocket["structure_mapping_quality"] = quality
        quality_by_complex[pocket["complex_id"]] = quality
    for complex_row in complexes_by_pdb[pdb_id]:
        complex_row["rcsb_entry_status"] = "current"
        complex_row["structure_mapping_quality"] = quality_by_complex.get(
            complex_row["complex_id"], "unresolved"
        )
    for ligand in ligands_by_pdb[pdb_id]:
        related = [
            item for item in enriched["rcsb_ligand_mapping_candidates"]
            if item["ligand_instance_id"] == ligand["ligand_instance_id"]
        ]
        ligand["rcsb_ligand_match_overall_status"] = (
            "probable" if any(item["selected"] for item in related) else
            "ambiguous" if any(item["status"] == "ambiguous" for item in related) else
            "unresolved"
        )


def _mark_unavailable(pdb_id, complexes_by_pdb, pockets_by_pdb) -> None:
    for row in complexes_by_pdb[pdb_id]:
        row["rcsb_entry_status"] = "unavailable"
        row["structure_mapping_quality"] = "unavailable"
    for row in pockets_by_pdb[pdb_id]:
        row["structure_mapping_quality"] = "unavailable"


def _failure_issues(pdb_id, failure, complexes_by_pdb, pockets_by_pdb):
    pocket_by_complex = {row["complex_id"]: row["pocket_instance_id"] for row in pockets_by_pdb[pdb_id]}
    return [
        ProcessingIssue(
            "download-rcsb",
            "RCSB_ENRICHMENT_FAILED",
            "warning",
            failure["error"],
            complex_id=row["complex_id"],
            pocket_instance_id=pocket_by_complex.get(row["complex_id"]),
            exception_type=failure.get("exception_type"),
            details=failure,
        )
        for row in complexes_by_pdb[pdb_id]
        if row["processing_status"].startswith("accepted")
    ]


def _normalize_enrichment_rows(rows):
    citations = {}
    for row in rows["citations"]:
        prior = citations.get(row["citation_id"])
        if prior is None:
            citations[row["citation_id"]] = row
            continue
        excluded = {"metadata_sources", "conflict_status"}
        if ({key: value for key, value in row.items() if key not in excluded}
                != {key: value for key, value in prior.items() if key not in excluded}):
            prior["conflict_status"] = "metadata_conflict"
        prior["metadata_sources"] = sorted(set(prior["metadata_sources"]) | set(row["metadata_sources"]))
    rows["citations"] = list(citations.values())
    for name, keys in (
        ("citation_authors", ("citation_id", "source", "ordinal")),
        ("pdb_citation_links", ("citation_id", "pdb_id", "source", "role")),
        ("chain_uniprot_mappings", ("pocket_instance_id", "pdbbind_auth_chain_id", "uniprot_accession")),
    ):
        unique = {}
        for row in rows[name]:
            unique.setdefault(tuple(row[key] for key in keys), row)
        rows[name] = list(unique.values())


def _apply_structural_citation_evidence(rows) -> None:
    links_by_pdb: dict[str, list[dict]] = defaultdict(list)
    for link in rows["pdb_citation_links"]:
        links_by_pdb[link["pdb_id"]].append(link)
    existing_links = {
        (row["measurement_id"], row["citation_id"]) for row in rows["affinity_reference_links"]
    }
    adjudications = {row["measurement_id"]: row for row in rows["affinity_reference_adjudications"]}
    complexes = {row["complex_id"]: row for row in rows["complexes"]}
    pockets_by_complex: dict[str, list[dict]] = defaultdict(list)
    for pocket in rows["pockets"]:
        pockets_by_complex[pocket["complex_id"]].append(pocket)
    for measurement in rows["binding_measurements"]:
        citation_links = links_by_pdb[measurement["pdb_id"]]
        primary = sorted({link["citation_id"] for link in citation_links if link["role"] == "primary"})
        for link in citation_links:
            key = (measurement["measurement_id"], link["citation_id"])
            if key in existing_links:
                continue
            rows["affinity_reference_links"].append({
                "measurement_id": measurement["measurement_id"],
                "complex_id": measurement["complex_id"],
                "citation_id": link["citation_id"],
                "candidate_status": (
                    "probable_structural_reference" if link["role"] == "primary"
                    else "structural_reference_only"
                ),
                "confidence": 0.60 if link["role"] == "primary" else 0.30,
                "evidence_sources": [link["source"]],
                "evidence_note": "Depositor citation; not asserted as the affinity-measurement source",
                "automatic_or_manual": "automatic",
                "verified_by": None,
                "verified_at_utc": None,
            })
            existing_links.add(key)
        if len(primary) == 1:
            status, selected, quality = "probable_structural_reference", primary[0], "probable"
        elif len(primary) > 1:
            status, selected, quality = "conflicting_references", None, "unresolved"
        elif citation_links:
            status, selected, quality = "reference_unresolved", None, "unresolved"
        else:
            status, selected, quality = "no_reference_available", None, "unavailable"
        adjudication = adjudications.get(measurement["measurement_id"])
        if adjudication is not None and adjudication["reference_status"] == "not_attempted":
            adjudication.update({
                "reference_status": status,
                "selected_citation_id": selected,
                "rule_version": "rcsb-postbuild-v1",
                "confidence": 0.60 if selected else None,
                "evidence_summary": "Structural-citation evidence only; no automatic affinity-source assertion",
                "adjudicator": "automatic-rcsb-postbuild-v1",
                "adjudicated_at_utc": None,
            })
        complex_row = complexes[measurement["complex_id"]]
        if complex_row["bibliography_quality"] == "not_attempted":
            complex_row["bibliography_quality"] = quality
            for pocket in pockets_by_complex[measurement["complex_id"]]:
                pocket["bibliography_quality"] = quality
