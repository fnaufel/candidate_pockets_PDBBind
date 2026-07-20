"""Deterministic RDKit ligand loading, validation, and geometry hashing."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from .config import BuildConfig
from .exceptions import ParseError
from .hashing import canonical_json_hash, length_frame, normalized_array_bytes, sha256_bytes, sha256_file, stable_id
from .models import LigandGeometry


@dataclass(slots=True)
class _Attempt:
    source_format: str
    path: Path
    molecule: Chem.Mol | None
    error: str | None
    sanitization_status: str


def parse_ligand(
    complex_id: str, sdf_path: Path | None, mol2_path: Path | None, config: BuildConfig
) -> LigandGeometry:
    paths = {"sdf": sdf_path, "mol2": mol2_path}
    order = tuple(dict.fromkeys((config.ligand.primary_format, config.ligand.fallback_format)))
    attempts: list[_Attempt] = []
    for source_format in order:
        path = paths.get(source_format)
        if path and path.is_file() and path.stat().st_size:
            mol, error, sanitization = (_load_sdf(path, config) if source_format == "sdf"
                                        else _load_mol2(path, config))
            attempts.append(_Attempt(source_format, path, mol, error, sanitization))
    selected = next((attempt for attempt in attempts if attempt.molecule is not None), None)
    if selected is None:
        reasons = "; ".join(f"{item.source_format}: {item.error}" for item in attempts) or "no usable files"
        raise ParseError(f"Both ligand formats are unusable ({reasons})")
    source_format, source_path, mol = selected.source_format, selected.path, selected.molecule
    assert mol is not None
    if mol.GetNumAtoms() < 1 or mol.GetNumConformers() < 1:
        raise ParseError("Ligand has no atoms or explicit conformer")
    coordinates = np.asarray(mol.GetConformer().GetPositions(), dtype=np.float64)
    if coordinates.shape != (mol.GetNumAtoms(), 3) or not np.isfinite(coordinates).all():
        raise ParseError("Ligand coordinates are absent, malformed, or nonfinite")
    atomic_numbers = np.asarray([atom.GetAtomicNum() for atom in mol.GetAtoms()], dtype=np.int16)
    if np.any(atomic_numbers <= 0):
        raise ParseError("Ligand contains unresolved elements")
    if not np.any(atomic_numbers > 1):
        raise ParseError("Ligand contains no heavy atom")
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    components = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
    if len(components) > 1 and not config.ligand.allow_multiple_components:
        raise ParseError("Selected ligand contains multiple disconnected components")
    component_ids = np.zeros(mol.GetNumAtoms(), dtype=np.int16)
    for index, atom_indices in enumerate(components):
        component_ids[list(atom_indices)] = index
    source_hash = sha256_file(source_path)
    content_hash = _canonical_content_hash(atomic_numbers, coordinates, components)
    policy = {
        "schema": "ligand-derivation-v1", "complex_id": complex_id,
        "selected_source_sha256": source_hash, "source_format": source_format,
        "record_ordinal": 0, "conformer_ordinal": 0,
        "selection_policy": list(order), "sanitize": config.ligand.sanitize,
        "sanitization_status": selected.sanitization_status,
        "component_policy": config.ligand.pocket_defining_component_policy,
        "content_hash": content_hash,
        "original_atomic_numbers": atomic_numbers.tolist(),
        "original_coordinates_f64": coordinates.tolist(),
        "original_component_ids": component_ids.tolist(),
    }
    comparison = _compare_formats(attempts)
    warnings: list[str] = []
    if source_format == "mol2" and config.ligand.primary_format != "mol2":
        warnings.append("LIGAND_MOL2_FALLBACK")
    if comparison["status"] == "disagreement":
        warnings.append("SDF_MOL2_DISAGREEMENT")
    if selected.sanitization_status == "failed_geometry_usable":
        warnings.append("LIGAND_CHEMISTRY_UNUSUAL")
    ligand_id = stable_id("ligand", complex_id, content_hash, canonical_json_hash(policy))
    chemistry = _chemistry(mol, elements, selected.sanitization_status)
    return LigandGeometry(
        ligand_instance_id=ligand_id, source_format=source_format, source_path=source_path,
        source_sha256=source_hash, atomic_numbers=atomic_numbers, elements=elements,
        coordinates=coordinates, component_ids=component_ids, content_hash=content_hash,
        derivation_hash=canonical_json_hash(policy), chemistry=chemistry,
        comparison=comparison, warnings=warnings,
    )


def ligand_component_rows(ligand: LigandGeometry) -> list[dict]:
    rows: list[dict] = []
    component_count = int(ligand.component_ids.max()) + 1 if len(ligand.component_ids) else 0
    heavy = ligand.atomic_numbers > 1
    for component_index in range(component_count):
        indices = np.flatnonzero(ligand.component_ids == component_index)
        coords = ligand.coordinates[indices]
        other = np.flatnonzero(ligand.component_ids != component_index)
        separation = None
        if len(other):
            separation = float(np.min(np.linalg.norm(coords[:, None, :] - ligand.coordinates[other][None, :, :], axis=2)))
        counts: dict[str, int] = {}
        for element in np.asarray(ligand.elements)[indices]:
            counts[str(element)] = counts.get(str(element), 0) + 1
        rows.append({
            "ligand_instance_id": ligand.ligand_instance_id, "component_index": component_index,
            "atom_indices": indices.astype(int).tolist(), "atom_count": len(indices),
            "heavy_atom_count": int(np.sum(heavy[indices])), "element_counts_json": json.dumps(counts, sort_keys=True),
            "formal_charge": int(sum(ligand.chemistry["_atom_formal_charges"][index] for index in indices)),
            "centroid_x": float(coords[:, 0].mean()),
            "centroid_y": float(coords[:, 1].mean()), "centroid_z": float(coords[:, 2].mean()),
            "minimum_other_component_separation": separation, "is_pocket_defining": True,
        })
    return rows


def _load_sdf(path: Path, config: BuildConfig) -> tuple[Chem.Mol | None, str | None, str]:
    try:
        nonempty_records = [block for block in path.read_text(encoding="utf-8", errors="replace").split("$$$$")
                            if block.strip()]
        if len(nonempty_records) > 1 and config.ligand.multiple_sdf_record_policy == "reject":
            return None, "multiple SDF records", "not_attempted"
        supplier = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
        molecules = [mol for mol in supplier if mol is not None]
        if not molecules:
            return None, "RDKit returned no molecule", "not_attempted"
        return _attempt_sanitization(molecules[0], config)
    except Exception as error:  # RDKit surfaces several extension exceptions
        return None, f"{type(error).__name__}: {error}", "not_attempted"


def _load_mol2(path: Path, config: BuildConfig) -> tuple[Chem.Mol | None, str | None, str]:
    try:
        mol = Chem.MolFromMol2File(str(path), sanitize=False, removeHs=False)
        if mol is None:
            return None, "RDKit returned no molecule", "not_attempted"
        return _attempt_sanitization(mol, config)
    except Exception as error:
        return None, f"{type(error).__name__}: {error}", "not_attempted"


def _attempt_sanitization(mol: Chem.Mol, config: BuildConfig) -> tuple[Chem.Mol, str | None, str]:
    if not config.ligand.sanitize:
        return mol, None, "not_requested"
    sanitized = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(sanitized)
        return sanitized, None, "sanitized"
    except Exception as error:
        # Coordinates and atomic identities remain usable; chemistry calls become best-effort.
        return mol, f"{type(error).__name__}: {error}", "failed_geometry_usable"


def _compare_formats(attempts: list[_Attempt]) -> dict:
    loaded = {item.source_format: item.molecule for item in attempts if item.molecule is not None}
    if not {"sdf", "mol2"}.issubset(loaded):
        return {"status": "unavailable", "coordinate_rmsd": None}
    sdf, mol2 = loaded["sdf"], loaded["mol2"]
    assert sdf is not None and mol2 is not None
    metrics = {
        "atom_count_equal": sdf.GetNumAtoms() == mol2.GetNumAtoms(),
        "heavy_atom_count_equal": sdf.GetNumHeavyAtoms() == mol2.GetNumHeavyAtoms(),
        "element_multiset_equal": sorted(atom.GetAtomicNum() for atom in sdf.GetAtoms()) ==
                                  sorted(atom.GetAtomicNum() for atom in mol2.GetAtoms()),
        "formal_charge_equal": Chem.GetFormalCharge(sdf) == Chem.GetFormalCharge(mol2),
        "component_count_equal": len(Chem.GetMolFrags(sdf)) == len(Chem.GetMolFrags(mol2)),
        "bond_count_equal": sdf.GetNumBonds() == mol2.GetNumBonds(),
    }
    rmsd = None
    mapping_ambiguity = False
    mapping_cap_exhausted = False
    if all(metrics.values()):
        try:
            matches = mol2.GetSubstructMatches(sdf, uniquify=False, maxMatches=1025)
        except Exception:
            matches = ()
        mapping_cap_exhausted = len(matches) > 1024
        matches = matches[:1024]
        if matches:
            left = np.asarray(sdf.GetConformer().GetPositions(), dtype=np.float64)
            right = np.asarray(mol2.GetConformer().GetPositions(), dtype=np.float64)
            candidates = []
            for mapping in matches:
                value = float(np.sqrt(np.mean(np.sum((left - right[list(mapping)]) ** 2, axis=1))))
                candidates.append((value, tuple(mapping)))
            candidates.sort()
            rmsd = candidates[0][0]
            mapping_ambiguity = len(candidates) > 1 and np.isclose(candidates[0][0], candidates[1][0], atol=1e-12)
    equivalent = all(metrics.values()) and rmsd is not None and rmsd <= 0.1
    return {"status": "equivalent" if equivalent else "disagreement", "coordinate_rmsd": rmsd,
            "mapping_ambiguity": mapping_ambiguity, "mapping_cap_exhausted": mapping_cap_exhausted,
            **metrics}


def _chemistry(mol: Chem.Mol, elements: list[str], sanitization_status: str) -> dict:
    counts: dict[str, int] = {}
    for element in elements:
        counts[element] = counts.get(element, 0) + 1
    try:
        smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
        isomeric = Chem.MolToSmiles(mol, isomericSmiles=True)
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.InchiToInchiKey(inchi) if inchi else None
    except Exception:
        smiles = isomeric = inchi = inchikey = None
    return {
        "canonical_smiles": smiles, "isomeric_smiles": isomeric, "inchi": inchi,
        "inchikey": inchikey, "molecular_formula": _safe(lambda: rdMolDescriptors.CalcMolFormula(mol)),
        "formal_charge": _safe(lambda: int(Chem.GetFormalCharge(mol))),
        "molecular_weight": _safe(lambda: float(Descriptors.MolWt(mol))),
        "atom_count": mol.GetNumAtoms(), "heavy_atom_count": mol.GetNumHeavyAtoms(),
        "component_count": len(Chem.GetMolFrags(mol)), "element_counts_json": json.dumps(counts, sort_keys=True),
        "stereochemistry_status": "perceived" if sanitization_status == "sanitized" else "unresolved",
        "rdkit_sanitization_status": sanitization_status,
        "_atom_formal_charges": [int(atom.GetFormalCharge()) for atom in mol.GetAtoms()],
    }


def _canonical_content_hash(atomic_numbers, coordinates, components) -> str:
    canonical_components = []
    for component in components:
        ordered = sorted(component, key=lambda index: (
            int(atomic_numbers[index]), normalized_array_bytes(coordinates[index], "<f8"), int(index)
        ))
        canonical_components.append(tuple(ordered))
    canonical_components.sort(key=lambda indices: tuple(
        (int(atomic_numbers[index]), normalized_array_bytes(coordinates[index], "<f8")) for index in indices
    ))
    canonical_indices = [index for component in canonical_components for index in component]
    boundaries = np.cumsum([len(component) for component in canonical_components], dtype=np.int64)
    return sha256_bytes(length_frame((
        b"ligand-content-v1", normalized_array_bytes(atomic_numbers[canonical_indices], "<i2"),
        normalized_array_bytes(coordinates[canonical_indices], "<f8"), normalized_array_bytes(boundaries, "<i8"),
    )))


def _safe(operation):
    try:
        return operation()
    except Exception:
        return None
