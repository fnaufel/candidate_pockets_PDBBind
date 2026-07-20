from pathlib import Path

import numpy as np

from biosensia_pocket_library.config import load_config
from biosensia_pocket_library.models import LigandGeometry
from biosensia_pocket_library.pocket_extractor import extract_pocket
from biosensia_pocket_library.protein_parser import parse_protein


def _atom(serial, name, alt, residue, chain, number, x, occupancy=1.0, element="C"):
    return (f"ATOM  {serial:5d} {name:>4s}{alt:1s}{residue:>3s} {chain:1s}{number:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{occupancy:6.2f}{20.0:6.2f}          {element:>2s}\n")


def _ligand(tmp_path):
    coords = np.asarray([[0.0, 0.0, 0.0]], dtype=float)
    return LigandGeometry("ligand:x", "sdf", tmp_path / "x.sdf", "a" * 64,
                          np.asarray([6], dtype=np.int16), ["C"], coords,
                          np.asarray([0], dtype=np.int16), "b" * 64, "c" * 64, {}, {})


def test_altloc_selection_boundary_and_contact_export(tmp_path: Path):
    pdb = tmp_path / "protein.pdb"
    pdb.write_text("HEADER TEST\n" +
        _atom(1, "CA", "A", "ALA", "A", 1, 5.0, 0.4) +
        _atom(2, "CA", "B", "ALA", "A", 1, 4.0, 0.6) +
        _atom(3, "CB", "", "ALA", "A", 1, 7.0) +
        _atom(4, "CA", "", "GLY", "A", 2, 6.0) + "END\n")
    config = load_config(project_root=tmp_path, overrides={"pocket.minimum_pocket_atoms_warning": 1})
    protein = parse_protein(pdb, config)
    assert protein.discarded_altloc_count == 1
    pocket = extract_pocket("complex", "1abc", _ligand(tmp_path), protein, {"C", "N", "O", "S", "H"}, config)
    assert len(pocket.contact_atoms) == 2  # 4.0 and inclusive 6.0
    assert len(pocket.residue_expanded_atoms) == 3  # whole ALA residue includes the 7.0 atom
    assert len(pocket.exported_atoms) == 2  # schema v1 exports contact_atom
    assert pocket.exported_atoms[0].x == 4.0


def test_deterministic_crop_and_hash(tmp_path: Path):
    pdb = tmp_path / "protein.pdb"
    pdb.write_text("HEADER TEST\n" + "".join(_atom(i, f"C{i}", "", "ALA", "A", i, float(i)) for i in range(1, 7)) + "END\n")
    config = load_config(project_root=tmp_path, overrides={"pocket.max_pocket_atoms": 3,
                                                           "pocket.minimum_pocket_atoms_warning": 1})
    protein = parse_protein(pdb, config)
    first = extract_pocket("complex", "1abc", _ligand(tmp_path), protein, {"C", "N", "O", "S", "H"}, config)
    second = extract_pocket("complex", "1abc", _ligand(tmp_path), protein, {"C", "N", "O", "S", "H"}, config)
    assert first.crop_applied and [atom.x for atom in first.exported_atoms] == [1.0, 2.0, 3.0]
    assert first.content_hash == second.content_hash
    assert first.derivation_hash == second.derivation_hash
