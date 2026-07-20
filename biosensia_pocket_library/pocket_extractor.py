"""Local PDBbind-only pocket extraction and deterministic DrugCLIP export view."""

from __future__ import annotations

import statistics

import numpy as np

from .config import BuildConfig
from .hashing import canonical_json_hash, length_frame, normalized_array_bytes, sha256_bytes, stable_id
from .models import ExtractedPocket, LigandGeometry, ParsedProtein, ProteinAtom


def extract_pocket(
    complex_id: str, pdb_id: str, ligand: LigandGeometry, protein: ParsedProtein,
    supported_tokens: set[str], config: BuildConfig, protein_source_sha256: str | None = None,
) -> ExtractedPocket:
    ligand_mask = ligand.atomic_numbers > 1 if config.pocket.distance_uses_ligand_heavy_atoms else np.ones(len(ligand.atomic_numbers), bool)
    ligand_coordinates = ligand.coordinates[ligand_mask]
    atoms = [atom for atom in protein.atoms if config.pocket.include_protein_hydrogens or atom.element != "H"]
    if not atoms or not len(ligand_coordinates):
        raise ValueError("No atoms available for pocket distance calculation")
    coordinates = np.asarray([atom.coordinate for atom in atoms])
    minimum = _minimum_distances(coordinates, ligand_coordinates)
    min_by_key = {atom.pdbbind_atom_key: float(distance) for atom, distance in zip(atoms, minimum, strict=True)}
    contact = [atom for atom, distance in zip(atoms, minimum, strict=True) if distance <= config.pocket.distance_cutoff_angstrom]
    if len(contact) < config.pocket.minimum_pocket_atoms_hard:
        raise ValueError("No pocket protein atoms selected")
    residue_keys = {atom.residue_key for atom in contact}
    expanded = [atom for atom in atoms if atom.residue_key in residue_keys]
    supported: list[ProteinAtom] = []
    warnings = list(ligand.warnings)
    unsupported = [atom for atom in contact if config.elements.explicit_mappings.get(atom.element, atom.element) not in supported_tokens]
    if unsupported and config.elements.unsupported_policy == "reject":
        raise ValueError(f"Unsupported pocket elements: {sorted({atom.element for atom in unsupported})}")
    # Schema v1 exports contact_atom, while residue_expanded_atom remains in sidecars.
    for atom in contact:
        mapped = config.elements.explicit_mappings.get(atom.element, atom.element)
        if mapped in supported_tokens:
            if mapped != atom.element:
                atom = _replace_element(atom, mapped)
                warnings.append("EXPLICIT_ELEMENT_MAPPING_APPLIED")
            supported.append(atom)
        else:
            warnings.append("UNSUPPORTED_ATOM_EXCLUDED")
    if not supported:
        raise ValueError("No supported protein heavy atoms remain for DrugCLIP export")
    ordered = sorted(supported, key=lambda atom: (
        min_by_key[atom.pdbbind_atom_key], atom.auth_chain_id, atom.auth_residue_number,
        atom.insertion_code, atom.residue_name, atom.atom_name, atom.altloc, atom.source_order,
    ))
    crop = len(ordered) > config.pocket.max_pocket_atoms
    exported = ordered[: config.pocket.max_pocket_atoms]
    if crop:
        warnings.append("DETERMINISTIC_CROP_APPLIED")
    if len(exported) < config.pocket.minimum_pocket_atoms_warning:
        warnings.append("VERY_SMALL_POCKET")
    if any(atom.residue_name in config.pocket.modified_residue_allowlist for atom in exported):
        warnings.append("MODIFIED_RESIDUE_INCLUDED")
    if residue_keys & protein.missing_atom_residues:
        warnings.append("LOCAL_MISSING_ATOM_RECORD")
    export_elements = [atom.element for atom in exported]
    export_coords = np.asarray([atom.coordinate for atom in exported], dtype=np.float32)
    content_hash = sha256_bytes(length_frame((b"pocket-content-v1", "\0".join(export_elements).encode(),
                                              normalized_array_bytes(export_coords, "<f4"))))
    derivation = {
        "schema": "pocket-derivation-v1", "distribution_id": "pdbbind-2020-v2024p-20250804",
        "complex_id": complex_id, "protein_source_sha256": protein_source_sha256,
        "ligand_content_hash": ligand.content_hash,
        "protein_atom_keys": [atom.pdbbind_atom_key for atom in protein.atoms],
        "cutoff": config.pocket.distance_cutoff_angstrom,
        "heavy_ligand_only": config.pocket.distance_uses_ligand_heavy_atoms,
        "model_policy": config.structure.model_policy, "altloc_policy": config.structure.altloc_policy,
        "hydrogen_policy": config.pocket.include_protein_hydrogens,
        "classification_policy": config.pocket.polymer_classification_policy,
        "export_view": "contact_atom", "max_atoms": config.pocket.max_pocket_atoms,
        "ordered_atom_identities": [atom.pdbbind_atom_key for atom in exported],
        "ordered_source_coordinates_f64": [[atom.x, atom.y, atom.z] for atom in exported],
        "ordered_minimum_distances_f64": [min_by_key[atom.pdbbind_atom_key] for atom in exported],
        "content_hash": content_hash,
    }
    derivation_hash = canonical_json_hash(derivation)
    pocket_id = f"{complex_id}:{derivation_hash[:16]}"
    return ExtractedPocket(
        pocket_instance_id=pocket_id, complex_id=complex_id, ligand_instance_id=ligand.ligand_instance_id,
        pdb_id=pdb_id, content_hash=content_hash, derivation_hash=derivation_hash,
        contact_atoms=contact, residue_expanded_atoms=expanded, exported_atoms=exported,
        minimum_distances=min_by_key, contact_residues=sorted(residue_keys), crop_applied=crop,
        max_retained_distance=max((min_by_key[a.pdbbind_atom_key] for a in exported), default=None),
        min_discarded_distance=min((min_by_key[a.pdbbind_atom_key] for a in ordered[len(exported):]), default=None),
        geometry_quality_tier="A", warning_codes=sorted(set(warnings)), error_codes=[],
    )


def distance_statistics(pocket: ExtractedPocket) -> dict[str, float | None]:
    values = [pocket.minimum_distances[atom.pdbbind_atom_key] for atom in pocket.exported_atoms]
    return {"minimum_ligand_distance_min": min(values, default=None),
            "minimum_ligand_distance_mean": statistics.mean(values) if values else None,
            "minimum_ligand_distance_median": statistics.median(values) if values else None,
            "minimum_ligand_distance_max": max(values, default=None)}


def _minimum_distances(left: np.ndarray, right: np.ndarray, chunk_size: int = 4096) -> np.ndarray:
    result = np.empty(len(left), dtype=np.float64)
    for start in range(0, len(left), chunk_size):
        block = left[start:start + chunk_size]
        result[start:start + len(block)] = np.min(np.linalg.norm(block[:, None, :] - right[None, :, :], axis=2), axis=1)
    return result


def _replace_element(atom: ProteinAtom, element: str) -> ProteinAtom:
    values = {name: getattr(atom, name) for name in atom.__dataclass_fields__}
    values["element"] = element
    return ProteinAtom(**values)
