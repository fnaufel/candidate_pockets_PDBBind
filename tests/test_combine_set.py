from pathlib import Path

import gzip
import json
import lmdb
import numpy as np
import pickle
import pyarrow.parquet as pq
import pytest

from biosensia_pocket_library.combine_set_pipeline import build_combine_set_library
from biosensia_pocket_library.cli import main
from biosensia_pocket_library.config import load_config
from biosensia_pocket_library.exceptions import ConfigurationError, ParseError
from biosensia_pocket_library.hashing import atomic_write_bytes, canonical_json_bytes, sha256_bytes, sha256_file
from biosensia_pocket_library.rcsb_workflow import (
    enrich_library_from_cache,
    plan_rcsb_request,
    prefetch_rcsb_request,
    read_rcsb_request,
)
from biosensia_pocket_library.trusted_pickle import load_trusted_combine_set_pickle
from biosensia_pocket_library.validation import validate_run


def _pdb_atom(serial, name, element, x):
    return (f"ATOM  {serial:5d} {name:>4s} ALA A   1    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2s}\n")


def _fixture(root: Path):
    drugclip = root / "BioSensIA-DC/external/DrugCLIP"
    for relative in (
        "unimol/tasks/drugclip.py", "unimol/data/lmdb_dataset.py",
        "unimol/data/affinity_dataset.py", "unimol/data/remove_hydrogen_dataset.py",
        "unimol/data/cropping_dataset.py", "unimol/data/normalize_dataset.py",
    ):
        path = drugclip / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"# {relative}\n")
    (drugclip / "data").mkdir(exist_ok=True)
    (drugclip / "data/dict_pkt.txt").write_text("[PAD]\nC\nN\nO\nS\nH\n")
    helper = root / "BioSensIA-DC/lmdb_helpers.py"
    helper.write_text("# helper\n")

    bundle = drugclip / "data/pdb/combine_set/1abc"
    bundle.mkdir(parents=True)
    pocket_coordinates = [[2.0, 0.0, 0.0], [2.2, 0.0, 0.0], [3.0, 0.0, 0.0]]
    with (bundle / "data.pkl").open("wb") as handle:
        pickle.dump({
            "atoms": ["C"], "coordinates": [[[0.0, 0.0, 0.0]]],
            "pocket_atoms": ["C", "H", "N"], "pocket_coordinates": pocket_coordinates,
            "pocket": "1abc", "label": ("-logKd/Ki", 7.5),
        }, handle, protocol=4)
    sdf = """ligand
test

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""
    (bundle / "1abc_ligand.sdf").write_text(sdf)
    pocket = (_pdb_atom(1, "CA", "C", 2.0) + _pdb_atom(2, "HA", "H", 2.2)
              + _pdb_atom(3, "N", "N", 3.0) + "END\n")
    (bundle / "1abc_pocket.pdb").write_text(pocket)
    (bundle / "1abc_pocket6A.pdb").write_text(pocket)
    (bundle / "1abc_protein.pdb").write_text(pocket)
    return bundle


def _config(root: Path):
    return load_config(project_root=root, overrides={
        "pipeline.offline": True, "pipeline.progress": False,
        "combine_set.trusted_pickles": True, "rcsb.download_mmcif": False,
        "pocket.max_pocket_atoms": 2, "pocket.minimum_pocket_atoms_warning": 1,
    })


def _seed_rcsb_cache(cache: Path):
    content = b"""data_1abc
loop_
_entity_poly.entity_id
_entity_poly.type
1 'polypeptide(L)'
loop_
_entity.id
_entity.pdbx_description
1 'Example protein'
loop_
_atom_site.group_PDB
_atom_site.auth_asym_id
_atom_site.label_asym_id
_atom_site.label_entity_id
ATOM A X 1
ATOM A X 1
loop_
_struct_ref.id
_struct_ref.db_name
_struct_ref.pdbx_db_accession
1 UNP P12345
loop_
_struct_ref_seq.align_id
_struct_ref_seq.ref_id
_struct_ref_seq.pdbx_strand_id
_struct_ref_seq.seq_align_beg
_struct_ref_seq.seq_align_end
_struct_ref_seq.db_align_beg
_struct_ref_seq.db_align_end
1 1 A 1 10 5 14
loop_
_pdbx_nonpoly_scheme.asym_id
_pdbx_nonpoly_scheme.entity_id
_pdbx_nonpoly_scheme.mon_id
_pdbx_nonpoly_scheme.pdb_strand_id
_pdbx_nonpoly_scheme.pdb_seq_num
L 2 LIG B 101
loop_
_citation.id
_citation.title
_citation.journal_abbrev
_citation.year
_citation.journal_volume
_citation.page_first
_citation.page_last
_citation.pdbx_database_id_DOI
_citation.pdbx_database_id_PubMed
primary 'Example structure' JTEST 2020 1 1 9 10.1/example 1234
loop_
_citation_author.citation_id
_citation_author.name
_citation_author.ordinal
primary 'Doe, J.' 1
"""
    payload = gzip.compress(content, mtime=0)
    digest = sha256_bytes(payload)
    object_path = cache / "objects/sha256" / digest[:2] / digest
    atomic_write_bytes(object_path, payload)
    reference = {
        "pdb_id": "1abc", "sha256": digest, "payload_sha256": digest,
        "compressed": True, "request_method": "GET",
        "url": "https://files.rcsb.org/download/1ABC.cif.gz",
        "normalized_parameters": {}, "request_body_sha256": None,
        "selected_request_headers": {}, "response_status": 200,
        "response_headers": {"content-type": "application/gzip"},
        "etag": "fixture", "last_modified": None,
        "retrieved_at_utc": "2026-01-01T00:00:00+00:00",
        "parser_schema_version": "rcsb-mmcif-cache-v1", "error_classification": None,
    }
    atomic_write_bytes(
        cache / "request_index/rcsb_mmcif/1abc.json",
        canonical_json_bytes(reference) + b"\n",
    )


def test_trusted_pickle_loader_requires_explicit_authorization(tmp_path: Path):
    bundle = _fixture(tmp_path)
    with pytest.raises(ConfigurationError, match="trusted_pickles=true"):
        load_trusted_combine_set_pickle(
            bundle / "data.pkl", trusted=False, dictionary={"C", "H", "N"}
        )


def test_trusted_pickle_loader_rejects_nonfinite_coordinates(tmp_path: Path):
    path = tmp_path / "bad.pkl"
    with path.open("wb") as handle:
        pickle.dump({"pocket": "bad", "pocket_atoms": ["C"],
                     "pocket_coordinates": [[np.nan, 0.0, 0.0]]}, handle)
    with pytest.raises(ParseError, match="non-finite"):
        load_trusted_combine_set_pickle(path, trusted=True, dictionary={"C"})


def test_combine_set_end_to_end_preserves_raw_loader_input(tmp_path: Path):
    _fixture(tmp_path)
    config = _config(tmp_path)
    run_dir = build_combine_set_library(config, pdb_ids=["1abc"], progress=False)

    assert run_dir.name.startswith("dc-combine-v1-")
    assert validate_run(run_dir, config, progress=False) == []
    pocket = pq.read_table(run_dir / "sidecars/pockets.parquet").to_pylist()[0]
    assert pocket["geometry_origin"] == "drugclip_combine_set_pickle"
    assert pocket["drugclip_export_view"] == "source_pickle"
    assert pocket["source_geometry_atom_count"] == 3
    assert pocket["exported_atom_count"] == 3
    assert pocket["structure_mapping_quality"] == "exact"

    environment = lmdb.open(str(run_dir / "lmdb/candidate_pockets.lmdb"), subdir=False,
                            readonly=True, lock=False)
    try:
        with environment.begin() as transaction:
            record = pickle.loads(transaction.get(b"0"))
    finally:
        environment.close()
    assert record["pocket_atoms"] == ["C", "H", "N"]
    assert record["pocket_coordinates"].dtype == np.dtype("<f4")
    assert record["pocket_coordinates"].shape == (3, 3)


def test_combine_set_builder_refuses_implicit_pickle_trust(tmp_path: Path):
    _fixture(tmp_path)
    config = load_config(project_root=tmp_path, overrides={"pipeline.progress": False})
    with pytest.raises(ConfigurationError, match="trusted_pickles=true"):
        build_combine_set_library(config, limit=1, progress=False)


def test_offline_run_can_be_enriched_immutably_from_verified_cache(tmp_path: Path):
    _fixture(tmp_path)
    config = _config(tmp_path)
    parent = build_combine_set_library(config, pdb_ids=["1abc"], progress=False)
    parent_manifest_path = parent / "manifest.json"
    parent_sidecar_hashes = {
        path.name: sha256_file(path) for path in (parent / "sidecars").glob("*.parquet")
    }
    request_path = tmp_path / "rcsb-request.tsv"
    planned = plan_rcsb_request(parent, request_path)
    parsed = read_rcsb_request(request_path)
    assert planned["requested_pdb_ids"] == ["1abc"]
    assert parsed["source_run_id"] == parent.name

    cache = tmp_path / "shared-rcsb-cache"
    _seed_rcsb_cache(cache)
    online_config = load_config(project_root=tmp_path, overrides={
        "pipeline.offline": False, "pipeline.progress": False,
        "paths.external_cache_dir": cache, "rcsb.download_mmcif": True,
        "combine_set.trusted_pickles": True, "pocket.max_pocket_atoms": 2,
        "pocket.minimum_pocket_atoms_warning": 1,
    })
    prefetched = prefetch_rcsb_request(request_path, cache, online_config, progress=False)
    assert prefetched["record_count"] == 1
    assert prefetched["failure_count"] == 0
    assert Path(prefetched["manifest_path"]).is_file()

    offline_enrichment_config = load_config(project_root=tmp_path, overrides={
        "pipeline.offline": True, "pipeline.progress": False,
        "paths.external_cache_dir": cache, "rcsb.download_mmcif": True,
        "combine_set.trusted_pickles": True, "pocket.max_pocket_atoms": 2,
        "pocket.minimum_pocket_atoms_warning": 1,
    })
    derived = enrich_library_from_cache(
        parent, cache, offline_enrichment_config, output_root=tmp_path / "derived", progress=False
    )

    assert derived != parent
    assert validate_run(derived, offline_enrichment_config, progress=False) == []
    assert sha256_file(parent_manifest_path) == planned["source_manifest_sha256"]
    assert parent_sidecar_hashes == {
        path.name: sha256_file(path) for path in (parent / "sidecars").glob("*.parquet")
    }
    parent_manifest = json.loads(parent_manifest_path.read_text())
    derived_manifest = json.loads((derived / "manifest.json").read_text())
    assert derived_manifest["parent_run"]["run_id"] == parent.name
    assert derived_manifest["rcsb_cache"]["network_access_used"] is False
    assert (
        derived_manifest["lmdb_profiles"]["default"]["lmdb_physical_sha256"]
        == parent_manifest["lmdb_profiles"]["default"]["lmdb_physical_sha256"]
    )
    mappings = pq.read_table(derived / "sidecars/chain_uniprot_mappings.parquet").to_pylist()
    citations = pq.read_table(derived / "sidecars/citations.parquet").to_pylist()
    assert mappings[0]["uniprot_accession"] == "P12345"
    assert citations[0]["doi"] == "10.1/example"

    cache_object = cache / prefetched["records"][0]["object_path"]
    cache_object.write_bytes(b"tampered")
    with pytest.raises(ValueError, match="cache is incomplete"):
        enrich_library_from_cache(
            parent,
            cache,
            offline_enrichment_config,
            output_root=tmp_path / "tampered-derived",
            progress=False,
        )


def test_offline_enrichment_rejects_an_incomplete_cache_before_creating_a_run(tmp_path: Path):
    _fixture(tmp_path)
    config = _config(tmp_path)
    parent = build_combine_set_library(config, pdb_ids=["1abc"], progress=False)
    empty_cache = tmp_path / "empty-cache"
    enrichment_config = load_config(project_root=tmp_path, overrides={
        "pipeline.offline": True, "pipeline.progress": False,
        "paths.external_cache_dir": empty_cache, "rcsb.download_mmcif": True,
        "combine_set.trusted_pickles": True, "pocket.max_pocket_atoms": 2,
        "pocket.minimum_pocket_atoms_warning": 1,
    })
    with pytest.raises(ValueError, match="cache is incomplete"):
        enrich_library_from_cache(
            parent, empty_cache, enrichment_config, output_root=tmp_path / "derived", progress=False
        )
    assert not (tmp_path / "derived").exists()


def test_post_build_rcsb_cli_commands_run_end_to_end_offline_after_prefetch(
    tmp_path: Path, monkeypatch
):
    _fixture(tmp_path)
    parent = build_combine_set_library(_config(tmp_path), pdb_ids=["1abc"], progress=False)
    monkeypatch.chdir(tmp_path)
    request = tmp_path / "cli-rcsb-request.tsv"
    cache = tmp_path / "cli-rcsb-cache"
    output_root = tmp_path / "cli-derived"

    assert main([
        "plan-rcsb", "--run-dir", str(parent), "--output", str(request),
    ]) == 0
    _seed_rcsb_cache(cache)
    assert main([
        "prefetch-rcsb", "--config", str(parent / "config.resolved.toml"),
        "--request", str(request), "--cache-dir", str(cache), "--no-progress",
    ]) == 0
    assert main([
        "enrich-from-cache", "--run-dir", str(parent), "--cache-dir", str(cache),
        "--output-root", str(output_root), "--no-progress",
    ]) == 0

    derived_runs = list(output_root.iterdir())
    assert len(derived_runs) == 1
    manifest = json.loads((derived_runs[0] / "manifest.json").read_text())
    assert manifest["status"] == "complete"
    assert manifest["parent_run"]["run_id"] == parent.name
