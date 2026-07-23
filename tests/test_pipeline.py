from pathlib import Path

import json
import threading
import time

import pyarrow.parquet as pq

from biosensia_pocket_library.config import load_config
from biosensia_pocket_library.drugclip_contract import verify_drugclip_contract
from biosensia_pocket_library.finalization import finalize_run
from biosensia_pocket_library.manifest import PIPELINE_STAGE_ORDER
from biosensia_pocket_library.pipeline import _bounded_thread_map, build_library
from biosensia_pocket_library.reporting import generate_reports
from biosensia_pocket_library.validation import validate_run


def _pdb_atom(serial, name, residue, number, x, element="C"):
    return (f"ATOM  {serial:5d} {name:>4s} {residue:>3s} A{number:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2s}\n")


def _fixture(root: Path):
    index = root / "data/raw/index"
    complex_dir = root / "data/raw/P-L/1981-2000/1abc"
    index.mkdir(parents=True)
    complex_dir.mkdir(parents=True)
    for name in ("INDEX_general_NL.2020R1.lst", "INDEX_general_PN.2020R1.lst", "INDEX_general_PP.2020R1.lst"):
        (index / name).write_text("# empty\n")
    (index / "README").write_text(
        "This special data package contains 1 protein-ligand complexes in PDBbind v2020; "
        "structures were reprocessed in v2024. Latest update: Aug 2025\n"
    )
    (index / "INDEX_general_PL.2020R1.lst").write_text(
        "# 1 protein-ligand complexes in total\n# Latest update: Aug 2025\n"
        "1abc 2.00 2001 Kd<10uM // hidden.PDF (LIG) synthetic\n")
    sdf = """ligand
test

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""
    (complex_dir / "1abc_ligand.sdf").write_text(sdf)
    protein = "HEADER TEST\n" + _pdb_atom(1, "CA", "ALA", 1, 2.2) + _pdb_atom(2, "CB", "ALA", 1, 3.2) + "END\n"
    (complex_dir / "1abc_protein.pdb").write_text(protein)
    (complex_dir / "1abc_pocket.pdb").write_text(protein)
    biosensia = root / "BioSensIA-DC"
    drugclip = biosensia / "external/DrugCLIP"
    for relative in ("unimol/tasks/drugclip.py", "unimol/data/lmdb_dataset.py",
                     "unimol/data/affinity_dataset.py", "unimol/data/remove_hydrogen_dataset.py",
                     "unimol/data/cropping_dataset.py", "unimol/data/normalize_dataset.py"):
        path = drugclip / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"# {relative}\n")
    (drugclip / "data").mkdir(parents=True)
    (drugclip / "data/dict_pkt.txt").write_text("[PAD]\nC\nN\nO\nS\nH\n")
    (biosensia / "lmdb_helpers.py").write_text("# helper\n")


def test_synthetic_end_to_end_build(tmp_path: Path):
    _fixture(tmp_path)
    config = load_config(project_root=tmp_path, overrides={
        "pipeline.offline": True, "pipeline.progress": False, "rcsb.download_mmcif": False,
        "pocket.minimum_pocket_atoms_warning": 1,
    })
    run_dir = build_library(config, pdb_ids=["1abc"], progress=False)
    assert run_dir.name.startswith("pb20-v24p-20250804-v1-")
    assert validate_run(run_dir, config, progress=False) == []
    complex_row = pq.read_table(run_dir / "sidecars/complexes.parquet").to_pylist()[0]
    assert complex_row["processing_status"].startswith("accepted")
    assert ".pdf" not in complex_row["index_line_redacted"].lower()
    pocket = pq.read_table(run_dir / "sidecars/pockets.parquet").to_pylist()[0]
    assert pocket["drugclip_export_view"] == "contact_atom"
    assert pocket["pocket_instance_id"].startswith(complex_row["complex_id"] + ":")
    assert (run_dir / "lmdb/candidate_pockets.lmdb").is_file()
    checkpoint_directories = sorted(path.name for path in (run_dir / "checkpoints").iterdir())
    assert checkpoint_directories == [
        f"{ordinal:03d}_{stage}" for ordinal, stage in enumerate(PIPELINE_STAGE_ORDER)
    ]


def test_worker_count_does_not_change_logical_outputs(tmp_path: Path):
    roots = [tmp_path / "one", tmp_path / "two"]
    manifests = []
    for workers, root in zip((1, 2), roots, strict=True):
        _fixture(root)
        config = load_config(project_root=root, overrides={
            "pipeline.offline": True, "pipeline.progress": False, "pipeline.workers": workers,
            "rcsb.download_mmcif": False, "pocket.minimum_pocket_atoms_warning": 1,
        })
        run_dir = build_library(config, pdb_ids=["1abc"], progress=False)
        manifests.append(json.loads((run_dir / "manifest.json").read_text()))
    assert manifests[0]["semantic_config_hash"] == manifests[1]["semantic_config_hash"]
    assert manifests[0]["operational_config_hash"] != manifests[1]["operational_config_hash"]
    assert {key: value["logical_sha256"] for key, value in manifests[0]["sidecar_artifacts"].items()} == {
        key: value["logical_sha256"] for key, value in manifests[1]["sidecar_artifacts"].items()
    }


def test_bounded_thread_map_avoids_order_blocking_and_eager_submission():
    started = []
    lock = threading.Lock()

    def operation(value):
        with lock:
            started.append(value)
        time.sleep(0.15 if value == 0 else 0.01)
        return value

    results = _bounded_thread_map(operation, [(value,) for value in range(10)], workers=3)
    first = next(results)
    assert first in {1, 2}
    assert len(started) <= 3
    assert sorted([first, *results]) == list(range(10))


def test_library_contract_does_not_hash_or_depend_on_encoder_checkpoint(tmp_path: Path):
    _fixture(tmp_path)
    checkpoint = tmp_path / "BioSensIA-DC/external/DrugCLIP/checkpoint_best.pt"
    checkpoint.write_bytes(b"first encoder")
    config = load_config(project_root=tmp_path)
    alternate_config = load_config(project_root=tmp_path, overrides={
        "paths.drugclip_checkpoint": tmp_path / "another-checkpoint.pt"
    })

    first = verify_drugclip_contract(config, progress=False)
    checkpoint.write_bytes(b"different encoder")
    second = verify_drugclip_contract(config, progress=False)

    assert first["drugclip_library_contract_fingerprint"] == second["drugclip_library_contract_fingerprint"]
    assert config.semantic_hash == alternate_config.semantic_hash
    assert config.operational_hash == alternate_config.operational_hash
    assert "checkpoint_sha256" not in first
    assert first["encoder_checkpoint"]["verification_status"] == "not_evaluated"


def test_legacy_manifest_remains_valid(tmp_path: Path):
    _fixture(tmp_path)
    config = load_config(project_root=tmp_path, overrides={
        "pipeline.offline": True, "pipeline.progress": False, "rcsb.download_mmcif": False,
        "pocket.minimum_pocket_atoms_warning": 1,
    })
    run_dir = build_library(config, pdb_ids=["1abc"], progress=False)
    manifest_path = run_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest.pop("drugclip_library_contract_fingerprint")
    manifest["drugclip_contract"].pop("drugclip_library_contract_fingerprint")
    manifest["drugclip_contract"].pop("contract_schema_version")
    manifest["drugclip_contract"].pop("encoder_checkpoint")
    manifest["drugclip_contract"]["checkpoint_sha256"] = "legacy-checkpoint-hash"
    manifest_path.write_text(json.dumps(manifest))

    assert validate_run(run_dir, config, progress=False) == []
    generate_reports(run_dir, manifest)
    summary = json.loads((run_dir / "reports/build_summary.json").read_text())
    assert summary["drugclip_library_contract_fingerprint"] == manifest["drugclip_contract_fingerprint"]


def test_finalize_recovers_interrupted_legacy_run_without_changing_identity(tmp_path: Path):
    _fixture(tmp_path)
    config = load_config(project_root=tmp_path, overrides={
        "pipeline.offline": True, "pipeline.progress": False, "rcsb.download_mmcif": False,
        "pocket.minimum_pocket_atoms_warning": 1,
    })
    run_dir = build_library(config, pdb_ids=["1abc"], export=False, progress=False)
    manifest_path = run_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    identity = {key: manifest[key] for key in (
        "run_id", "semantic_config_hash", "source_fingerprint", "selection_fingerprint",
        "drugclip_contract_fingerprint", "git_commit", "code_dirty_state_fingerprint",
    )}
    manifest.pop("drugclip_library_contract_fingerprint")
    manifest["drugclip_contract"]["checkpoint_sha256"] = "legacy-checkpoint-hash"
    manifest["completed_at_utc"] = None
    manifest["counts"] = {}
    manifest["status"] = "running"
    manifest_path.write_text(json.dumps(manifest))
    temporary = run_dir / "lmdb/.candidate_pockets.lmdb.123.tmp"
    temporary.parent.mkdir(exist_ok=True)
    temporary.write_bytes(b"abandoned")

    finalized = finalize_run(run_dir, config, progress=False)

    assert finalized["status"] == "complete"
    assert finalized["completed_at_utc"]
    assert finalized["counts"]["default_lmdb_records"] == 1
    assert finalized["lmdb_profiles"]["default"]["profile_name"] == "default"
    assert (run_dir / "lmdb/candidate_pockets.lmdb").is_file()
    assert not temporary.exists()
    assert not any(item["path"].endswith(".tmp") for item in finalized["output_files"])
    assert finalized["drugclip_contract"]["checkpoint_sha256"] == "legacy-checkpoint-hash"
    assert {key: finalized[key] for key in identity} == identity
    assert finalized["stage_statuses"]["export-lmdb"]["status"] == "complete"
    assert finalized["stage_statuses"]["validate"]["status"] == "complete"
    assert finalized["stage_statuses"]["report"]["status"] == "complete"
