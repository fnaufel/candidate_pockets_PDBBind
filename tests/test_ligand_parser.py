from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Geometry import Point3D

from biosensia_pocket_library.config import load_config
from biosensia_pocket_library.exceptions import ParseError
from biosensia_pocket_library.ligand_parser import (
    _canonical_content_hash,
    _connectivity_matches,
    _coordinate_element_mapping,
    parse_ligand,
)


def _sdf(x: float = 0.0) -> str:
    return f"""ligand
test

  1  0  0  0  0  0            999 V2000
    {x:0.4f}    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""


def test_ligand_content_hash_is_source_order_independent():
    numbers = np.asarray([8, 6, 7], dtype=np.int16)
    coordinates = np.asarray([[2.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    original = _canonical_content_hash(numbers, coordinates, ((0, 1, 2),))
    permutation = np.asarray([2, 0, 1])
    permuted = _canonical_content_hash(numbers[permutation], coordinates[permutation], ((0, 1, 2),))
    assert original == permuted


def test_multiple_sdf_records_are_rejected(tmp_path: Path):
    path = tmp_path / "ligand.sdf"
    path.write_text(_sdf() + _sdf(1.0))
    config = load_config(project_root=tmp_path)
    with pytest.raises(ParseError, match="multiple SDF records"):
        parse_ligand("complex", path, None, config)


def test_coordinate_mapping_handles_permuted_symmetric_atoms():
    molecule = Chem.AddHs(Chem.MolFromSmiles("CC"))
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for index in range(molecule.GetNumAtoms()):
        conformer.SetAtomPosition(index, Point3D(float(index), float(index % 3), 0.0))
    molecule.AddConformer(conformer)
    order = list(reversed(range(molecule.GetNumAtoms())))
    permuted = Chem.RenumberAtoms(molecule, order)

    mapping, ambiguous = _coordinate_element_mapping(molecule, permuted)

    assert mapping is not None
    assert not ambiguous
    assert _connectivity_matches(molecule, permuted, mapping)
    left = np.asarray(molecule.GetConformer().GetPositions())
    right = np.asarray(permuted.GetConformer().GetPositions())
    assert np.allclose(left, right[mapping])
