from pathlib import Path

import lmdb
import numpy as np
import pickle
import pyarrow.parquet as pq
import pytest

from biosensia_pocket_library.combine_set_pipeline import build_combine_set_library
from biosensia_pocket_library.config import load_config
from biosensia_pocket_library.exceptions import ConfigurationError, ParseError
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
