from pathlib import Path

import pytest

from biosensia_pocket_library.config import load_config
from biosensia_pocket_library.drugclip_contract import read_dictionary
from biosensia_pocket_library.ligand_parser import parse_ligand
from biosensia_pocket_library.pocket_extractor import extract_pocket
from biosensia_pocket_library.protein_parser import parse_protein


@pytest.mark.skipif(not Path("data/raw/P-L/1981-2000/2tpi").is_dir(), reason="local PDBbind data unavailable")
def test_local_2tpi_geometry_contract():
    config = load_config(overrides={"pipeline.offline": True})
    root = Path("data/raw/P-L/1981-2000/2tpi")
    ligand = parse_ligand("local:2tpi", root / "2tpi_ligand.sdf", root / "2tpi_ligand.mol2", config)
    protein = parse_protein(root / "2tpi_protein.pdb", config)
    pocket = extract_pocket("local:2tpi", "2tpi", ligand, protein,
                            read_dictionary(config.paths.drugclip_dictionary), config)
    assert 1 <= len(pocket.exported_atoms) <= 256
    assert all(atom.element != "H" for atom in pocket.exported_atoms)
