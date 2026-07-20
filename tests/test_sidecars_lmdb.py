from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from biosensia_pocket_library.config import load_config
from biosensia_pocket_library.lmdb_export import export_lmdb, validate_lmdb
from biosensia_pocket_library.schemas import TABLES
from biosensia_pocket_library.sidecars import write_sidecars


def test_empty_tables_use_explicit_schemas(tmp_path: Path):
    write_sidecars(tmp_path, {}, progress=False)
    for name, spec in TABLES.items():
        path = tmp_path / f"{name}.parquet"
        assert path.is_file()
        assert pq.read_schema(path).remove_metadata() == spec.schema.remove_metadata()


def test_lmdb_round_trip_from_sidecars(tmp_path: Path):
    config = load_config(project_root=tmp_path, overrides={"paths.drugclip_dictionary": tmp_path / "dict.txt"})
    (tmp_path / "dict.txt").write_text("[PAD]\nC\nN\nO\nS\nH\n")
    rows = {name: [] for name in TABLES}
    pocket_id = "pb20v24p-12345678:1abc:0123456789abcdef"
    rows["pockets"] = [{"pocket_instance_id": pocket_id, "complex_id": "c", "ligand_instance_id": "l",
                        "processing_status": "accepted", "geometry_quality_tier": "A",
                        "pocket_geometry_content_hash": "g", "pocket_derivation_hash": "d"}]
    rows["pocket_atoms"] = [{"pocket_instance_id": pocket_id, "pdbbind_atom_key": "a", "source_order": 1,
                             "element": "C", "x": 1.0, "y": 2.0, "z": 3.0,
                             "retained_after_crop": True, "export_order": 0}]
    write_sidecars(tmp_path / "sidecars", rows, progress=False)
    contract = {key: "hash" for key in ("dictionary_sha256", "task_sha256", "loader_sha256",
                                         "helper_sha256", "checkpoint_sha256")}
    metadata, lmdb_rows = export_lmdb(tmp_path, config, contract, progress=False)
    assert metadata["record_count"] == 1
    assert lmdb_rows[0]["pocket_instance_id"] == pocket_id
    rows["lmdb_records"] = lmdb_rows
    write_sidecars(tmp_path / "sidecars", rows, progress=False)
    assert validate_lmdb(tmp_path, "default", config, progress=False) == []
