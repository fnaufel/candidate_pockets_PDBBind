"""Content-addressed RCSB mmCIF cache and geometry-independent enrichment."""

from __future__ import annotations

import gzip
import json
import os
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import gemmi
import httpx
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

from .config import BuildConfig
from .hashing import atomic_write_bytes, canonical_json_bytes, sha256_bytes, sha256_file, stable_id
from .progress import track


def download_mmcif_files(pdb_ids: list[str], config: BuildConfig, *, refresh: bool = False,
                         progress: bool = True) -> tuple[dict[str, Path], list[dict], list[dict]]:
    cache = config.paths.external_cache_dir
    refs = cache / "request_index/rcsb_mmcif"
    objects = cache / "objects/sha256"
    refs.mkdir(parents=True, exist_ok=True)
    objects.mkdir(parents=True, exist_ok=True)
    resolved: dict[str, Path] = {}
    inventory: list[dict] = []
    failures: list[dict] = []
    with httpx.Client(timeout=config.rcsb.timeout_seconds, follow_redirects=True) as client:
        for pdb_id in track(sorted(set(pdb_ids)), description="Caching RCSB mmCIF", total=len(set(pdb_ids)), enabled=progress):
            reference = refs / f"{pdb_id.lower()}.json"
            cached = _resolve_reference(reference, objects, pdb_id, config.rcsb.download_compressed)
            if cached is not None and not refresh:
                resolved[pdb_id] = cached
                inventory.append(_cache_inventory(pdb_id, cached, None, None, None))
                inventory.append(_cache_metadata_inventory(pdb_id, reference))
                continue
            if not refresh:
                negative = _active_negative(reference)
                if negative is not None:
                    failures.append(negative)
                    continue
            if config.pipeline.offline:
                continue
            compressed = config.rcsb.download_compressed
            suffix = ".cif.gz" if compressed else ".cif"
            url = f"https://files.rcsb.org/download/{pdb_id.upper()}{suffix}"
            try:
                response = _fetch(client, url, config.rcsb.maximum_retries)
                content = response.content
                plain = gzip.decompress(content) if compressed else content
                _validate_mmcif(plain, pdb_id)
            except Exception as error:
                failure = {"pdb_id": pdb_id, "error": str(error),
                           "exception_type": type(error).__name__, "url": url}
                failures.append(failure)
                envelope = {"pdb_id": pdb_id, "request_method": "GET", "url": url,
                            "normalized_parameters": {}, "request_body_sha256": None,
                            "selected_request_headers": {}, "response_status": None,
                            "response_headers": {}, "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                            "payload_sha256": None, "parser_schema_version": "rcsb-mmcif-cache-v1",
                            "error_classification": type(error).__name__, "error": str(error),
                            "expires_at_utc": (datetime.now(timezone.utc) + timedelta(days=1)).isoformat()}
                atomic_write_bytes(reference, canonical_json_bytes(envelope) + b"\n")
                continue
            digest = sha256_bytes(content)
            object_path = objects / digest[:2] / digest
            if not object_path.exists():
                atomic_write_bytes(object_path, content)
            ref_data = {"pdb_id": pdb_id, "sha256": digest, "payload_sha256": digest,
                        "compressed": compressed, "request_method": "GET", "url": url,
                        "normalized_parameters": {}, "request_body_sha256": None,
                        "selected_request_headers": {}, "response_status": response.status_code,
                        "response_headers": {key: response.headers.get(key) for key in ("etag", "last-modified", "content-type")},
                        "etag": response.headers.get("etag"), "last_modified": response.headers.get("last-modified"),
                        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                        "parser_schema_version": "rcsb-mmcif-cache-v1", "error_classification": None}
            atomic_write_bytes(reference, canonical_json_bytes(ref_data) + b"\n")
            resolved[pdb_id] = object_path
            inventory.append(_cache_inventory(pdb_id, object_path, url, response.headers.get("etag"),
                                              response.headers.get("last-modified")))
            inventory.append(_cache_metadata_inventory(pdb_id, reference))
            if config.rcsb.requests_per_second > 0:
                time.sleep(1.0 / config.rcsb.requests_per_second)
    return resolved, inventory, failures


def enrich_from_mmcif(pdb_id: str, path: Path, compressed: bool, chain_rows: list[dict],
                      ligand_inputs: list[tuple[dict, str | None]] | None = None) -> dict[str, list[dict]]:
    content = path.read_bytes()
    if compressed:
        content = gzip.decompress(content)
    document = gemmi.cif.read_string(content.decode("utf-8", errors="replace"))
    block = document.sole_block()
    result = {"chain_mapping_candidates": [], "chain_uniprot_mappings": [],
              "chain_uniprot_mapping_segments": [], "citations": [], "pdb_citation_links": [],
              "citation_authors": [], "rcsb_ligand_mapping_candidates": []}
    _map_chains(block, pdb_id, chain_rows, result)
    _extract_uniprot(block, pdb_id, chain_rows, result)
    _map_ligands(block, ligand_inputs or [], result)
    _extract_citations(block, pdb_id, result)
    return result


def _fetch(client: httpx.Client, url: str, attempts: int) -> httpx.Response:
    @retry(stop=stop_after_attempt(attempts), wait=wait_exponential(min=1, max=20),
           retry=retry_if_exception_type((httpx.HTTPError, OSError)), reraise=True)
    def request():
        response = client.get(url)
        response.raise_for_status()
        return response
    return request()


def _validate_mmcif(content: bytes, pdb_id: str) -> None:
    try:
        document = gemmi.cif.read_string(content.decode("utf-8", errors="strict"))
        block = document.sole_block()
    except Exception as error:
        raise ValueError(f"Invalid RCSB mmCIF for {pdb_id}: {error}") from error
    block_name = block.name.lower().removeprefix("data_")
    if block_name != pdb_id.lower():
        raise ValueError(f"mmCIF identity {block.name!r} does not match {pdb_id}")


def _resolve_reference(reference: Path, objects: Path, pdb_id: str, expected_compressed: bool) -> Path | None:
    if not reference.is_file():
        return None
    try:
        data = json.loads(reference.read_text(encoding="utf-8"))
        if bool(data.get("compressed")) != expected_compressed:
            return None
        path = objects / data["sha256"][:2] / data["sha256"]
        if not path.is_file() or sha256_file(path) != data["sha256"]:
            return None
        content = path.read_bytes()
        plain = gzip.decompress(content) if data.get("compressed") else content
        _validate_mmcif(plain, pdb_id)
        return path
    except (ValueError, KeyError, OSError, gzip.BadGzipFile):
        return None


def _active_negative(reference: Path) -> dict | None:
    if not reference.is_file():
        return None
    try:
        data = json.loads(reference.read_text(encoding="utf-8"))
        if data.get("payload_sha256") is not None or not data.get("expires_at_utc"):
            return None
        expires = datetime.fromisoformat(data["expires_at_utc"].replace("Z", "+00:00"))
        if expires <= datetime.now(timezone.utc):
            return None
        return {"pdb_id": data.get("pdb_id"), "error": data.get("error", "cached negative response"),
                "exception_type": data.get("error_classification", "CachedNegativeResponse"),
                "url": data.get("url")}
    except (OSError, ValueError, TypeError):
        return None


def _cache_inventory(pdb_id, path, url, etag, modified):
    digest = path.name
    return {"source_file_id": stable_id("file", "rcsb_mmcif", pdb_id, digest), "source_kind": "rcsb_mmcif",
            "pdb_id": pdb_id, "path": path.as_posix(), "size_bytes": path.stat().st_size, "sha256": digest,
            "modified_time_utc": None, "download_url": url, "downloaded_at_utc": None,
            "http_etag": etag, "http_last_modified": modified, "validation_status": "valid", "warning_codes": []}


def _cache_metadata_inventory(pdb_id, path):
    digest = sha256_file(path)
    return {"source_file_id": stable_id("file", "rcsb_api_cache", pdb_id, digest),
            "source_kind": "rcsb_api_cache", "pdb_id": pdb_id, "path": path.as_posix(),
            "size_bytes": path.stat().st_size, "sha256": digest, "modified_time_utc": None,
            "download_url": None, "downloaded_at_utc": None, "http_etag": None,
            "http_last_modified": None, "validation_status": "valid", "warning_codes": []}


def _values(block, prefix, names):
    try:
        table = block.find(prefix, names)
        return [[gemmi.cif.as_string(str(value)) if str(value) not in {".", "?"} else None
                 for value in row] for row in table]
    except (RuntimeError, ValueError):
        return []


def _map_chains(block, pdb_id, chain_rows, result):
    entity_types = {entity_id: entity_type for entity_id, entity_type in
                    _values(block, "_entity_poly.", ["entity_id", "type"])}
    descriptions = {entity_id: description for entity_id, description in
                    _values(block, "_entity.", ["id", "pdbx_description"])}
    atom_rows = _values(block, "_atom_site.", ["group_PDB", "auth_asym_id", "label_asym_id", "label_entity_id"])
    by_auth: dict[str, dict[tuple[str, str], int]] = {}
    for group, auth, label, entity in atom_rows:
        if group == "ATOM" and auth and label:
            counts = by_auth.setdefault(auth, {})
            counts[(label, entity or "")] = counts.get((label, entity or ""), 0) + 1
    for chain in chain_rows:
        candidates = sorted(by_auth.get(chain["pdbbind_auth_chain_id"], {}).items(), key=lambda item: (-item[1], item[0]))
        for rank, ((label, entity), count) in enumerate(candidates, 1):
            entity_type = entity_types.get(entity)
            is_protein = bool(entity_type and "polypeptide" in entity_type.lower())
            selected = len(candidates) == 1 and is_protein
            result["chain_mapping_candidates"].append({
                "chain_mapping_candidate_id": stable_id("chainmap", chain["pocket_instance_id"], label, entity),
                "pocket_instance_id": chain["pocket_instance_id"], "pdb_id": pdb_id,
                "pdbbind_auth_chain_id": chain["pdbbind_auth_chain_id"], "rcsb_label_asym_id": label,
                "rcsb_auth_asym_id": chain["pdbbind_auth_chain_id"], "rcsb_polymer_entity_id": entity,
                "entity_type": entity_type, "entity_description": descriptions.get(entity),
                "organism_name": None, "taxonomy_id": None,
                "sequence_length": None, "mapping_method": "auth_chain_identifier",
                "mapping_status": "exact_identifier_match" if selected else ("not_a_protein" if not is_protein else "ambiguous"),
                "candidate_rank": rank, "selected": selected, "identifier_match_count": count,
                "atom_match_count": 0,
                "residue_coverage": min(1.0, count / max(1, chain["selected_atom_count"])),
                "alignment_atom_count": None, "alignment_rmsd": None,
                "transform_json": None, "evidence_codes": ["AUTH_ASYM_ID_MATCH"], "warning_codes": [],
            })


def _extract_citations(block, pdb_id, result):
    rows = _values(block, "_citation.", ["id", "title", "journal_abbrev", "year", "journal_volume",
                                         "page_first", "page_last", "pdbx_database_id_DOI", "pdbx_database_id_PubMed"])
    for source_id, title, journal, year, volume, first, last, doi, pmid in rows:
        fingerprint = stable_id("bib", title, journal, year, volume, first)
        citation_id = stable_id("citation", doi.lower() if doi else None, pmid, fingerprint)
        result["citations"].append({"citation_id": citation_id, "title": title, "journal": journal,
            "year": int(year) if year and year.isdigit() else None, "volume": volume, "issue": None,
            "first_page": first, "last_page": last, "doi": doi, "pmid": pmid, "crossref_id": None,
            "publication_status": "published", "source_priority": 1, "bibliographic_fingerprint": fingerprint,
            "metadata_sources": ["rcsb_mmcif"], "conflict_status": "none"})
        result["pdb_citation_links"].append({"citation_id": citation_id, "pdb_id": pdb_id,
            "source_citation_id": source_id, "role": "primary" if source_id == "primary" else "reference",
            "source": "rcsb_mmcif", "source_priority": 1, "evidence": "_citation"})
    authors = _values(block, "_citation_author.", ["citation_id", "name", "ordinal"])
    citation_lookup = {row["source_citation_id"]: row["citation_id"] for row in result["pdb_citation_links"]}
    for source_id, name, ordinal in authors:
        if source_id in citation_lookup:
            result["citation_authors"].append({"citation_id": citation_lookup[source_id],
                "ordinal": int(ordinal) if ordinal and ordinal.isdigit() else 0, "author_name": name,
                "orcid": None, "source": "rcsb_mmcif"})


def _extract_uniprot(block, pdb_id, chain_rows, result):
    references = _values(block, "_struct_ref.", ["id", "db_name", "pdbx_db_accession"])
    accessions = {ref_id: accession for ref_id, database, accession in references
                  if database and database.upper() in {"UNP", "UNIPROT"} and accession}
    alignments = _values(block, "_struct_ref_seq.", [
        "align_id", "ref_id", "pdbx_strand_id", "seq_align_beg", "seq_align_end",
        "db_align_beg", "db_align_end",
    ])
    summaries: dict[tuple[str, str, str], dict] = {}
    for chain in chain_rows:
        chain_id = chain["pdbbind_auth_chain_id"]
        for align_id, ref_id, strand_ids, pdb_begin, pdb_end, db_begin, db_end in alignments:
            if ref_id not in accessions or chain_id not in {item.strip() for item in (strand_ids or "").split(",")}:
                continue
            accession = accessions[ref_id]
            segment_id = stable_id("uniprot-segment", chain["pocket_instance_id"], accession, align_id, chain_id)
            result["chain_uniprot_mapping_segments"].append({
                "segment_id": segment_id, "pocket_instance_id": chain["pocket_instance_id"],
                "uniprot_accession": accession, "pdb_begin": _integer(pdb_begin), "pdb_end": _integer(pdb_end),
                "uniprot_begin": _integer(db_begin), "uniprot_end": _integer(db_end),
            })
            selected_candidate = next((item for item in result["chain_mapping_candidates"]
                                       if item["pocket_instance_id"] == chain["pocket_instance_id"]
                                       and item["pdbbind_auth_chain_id"] == chain_id and item["selected"]), None)
            summary_key = (chain["pocket_instance_id"], chain_id, accession)
            summaries.setdefault(summary_key, {
                "pocket_instance_id": chain["pocket_instance_id"], "pdb_id": pdb_id,
                "pdbbind_auth_chain_id": chain_id,
                "rcsb_label_asym_id": selected_candidate["rcsb_label_asym_id"] if selected_candidate else None,
                "rcsb_polymer_entity_id": selected_candidate["rcsb_polymer_entity_id"] if selected_candidate else None,
                "uniprot_accession": accession,
                "uniprot_isoform": None, "mapping_source": "struct_ref_seq",
                "mapping_status": "mapped", "mapping_coverage": None,
            })
    result["chain_uniprot_mappings"].extend(summaries.values())


def _map_ligands(block, ligand_inputs, result):
    instances = []
    for label_asym, entity_id, monomer, auth_chain, auth_seq in _values(
        block, "_pdbx_nonpoly_scheme.", ["asym_id", "entity_id", "mon_id", "pdb_strand_id", "pdb_seq_num"]
    ):
        instances.append({"label_asym_id": label_asym, "entity_id": entity_id, "ccd_id": monomer,
                          "auth_asym_id": auth_chain, "auth_seq_id": _integer(auth_seq), "bird_id": None,
                          "kind": "nonpolymer"})
    for label_asym, entity_id, monomer, auth_chain, auth_seq in _values(
        block, "_pdbx_branch_scheme.", ["asym_id", "entity_id", "mon_id", "pdb_asym_id", "pdb_seq_num"]
    ):
        instances.append({"label_asym_id": label_asym, "entity_id": entity_id, "ccd_id": monomer,
                          "auth_asym_id": auth_chain, "auth_seq_id": _integer(auth_seq), "bird_id": None,
                          "kind": "branched"})
    for ligand, ligand_label in ligand_inputs:
        label_tokens = {token.strip().upper() for token in (ligand_label or "").split("/") if token.strip()}
        ranked = sorted(instances, key=lambda item: (
            0 if item["ccd_id"] in label_tokens else 1,
            item["ccd_id"] or "", item["label_asym_id"] or "", item["auth_seq_id"] or -1,
        ))
        exact_count = sum(item["ccd_id"] in label_tokens for item in ranked)
        ambiguity_group = stable_id("ligand-map-group", ligand["ligand_instance_id"])
        for rank, item in enumerate(ranked, 1):
            label_match = item["ccd_id"] in label_tokens
            selected = label_match and exact_count == 1
            result["rcsb_ligand_mapping_candidates"].append({
                "candidate_id": stable_id("ligand-map", ligand["ligand_instance_id"], item["entity_id"],
                                          item["label_asym_id"], item["auth_seq_id"], item["ccd_id"]),
                "ligand_instance_id": ligand["ligand_instance_id"],
                "rcsb_nonpolymer_entity_id": item["entity_id"], "label_asym_id": item["label_asym_id"],
                "auth_asym_id": item["auth_asym_id"], "auth_seq_id": item["auth_seq_id"],
                "ccd_id": item["ccd_id"], "bird_id": item["bird_id"], "match_rmsd": None,
                "match_method": "pdbbind_ligand_label" if label_match else "enumerated_instance",
                "composition_evidence": "not_computed", "graph_evidence": "not_computed", "rank": rank,
                "selected": selected, "ambiguity_group": ambiguity_group,
                "status": "probable" if selected else ("ambiguous" if label_match else "unresolved"),
            })


def _integer(value):
    try:
        return int(value) if value is not None else None
    except ValueError:
        return None
