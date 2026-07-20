"""PDBbind PDB parser implementing model, alternate-location, and element policy."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import gemmi

from .config import BuildConfig
from .constants import CANONICAL_AMINO_ACIDS
from .exceptions import ParseError
from .hashing import stable_id
from .models import ParsedProtein, ProteinAtom


def parse_protein(path: Path, config: BuildConfig, *, pdb_id: str = "", complex_id: str = "",
                  distribution_id: str = "pdbbind-2020-v2024p-20250804") -> ParsedProtein:
    try:
        structure = gemmi.read_structure(str(path))
        model_count = len(structure)
    except Exception as error:
        raise ParseError(f"Gemmi rejected protein PDB: {type(error).__name__}: {error}") from error
    if model_count < 1:
        raise ParseError("Protein has no coordinate model")
    selected_model = _first_model_number(path)
    parsed: list[ProteinAtom] = []
    current_model = 1
    source_order = 0
    inferred = 0
    missing_atom_residues: set[tuple[str, int, str, str]] = set()
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for line in lines:
        record = line[:6].strip().upper()
        if line.startswith("REMARK 470"):
            missing = _parse_missing_atom_residue(line)
            if missing is not None:
                missing_atom_residues.add(missing)
        if record == "MODEL":
            try:
                current_model = int(line[10:14].strip())
            except ValueError:
                current_model += 1
            continue
        if record == "ENDMDL":
            continue
        if record not in {"ATOM", "HETATM"} or current_model != selected_model:
            continue
        source_order += 1
        try:
            residue_name = line[17:20].strip().upper()
            atom_name = line[12:16].strip()
            element = line[76:78].strip().title() if len(line) >= 78 else ""
            if not element:
                element = _infer_element(atom_name)
                inferred += 1
            chain = line[21:22].strip()
            residue_number = int(line[22:26].strip())
            insertion = line[26:27].strip()
            serial = int(line[6:11].strip()) if line[6:11].strip() else None
            occupancy = float(line[54:60]) if line[54:60].strip() else None
            b_factor = float(line[60:66]) if line[60:66].strip() else None
            locally_protein = record == "ATOM" and residue_name in CANONICAL_AMINO_ACIDS
            if record == "HETATM" and residue_name in config.pocket.modified_residue_allowlist:
                locally_protein = config.pocket.include_allowlisted_polymer_hetatm
            key = stable_id("pdbbind-atom", distribution_id, complex_id, pdb_id, selected_model,
                            record, chain, residue_number, insertion, residue_name, atom_name,
                            line[16:17].strip(), element)
            parsed.append(ProteinAtom(
                pdbbind_atom_key=key, source_order=source_order, model_id=selected_model,
                record_type=record, auth_chain_id=chain, auth_residue_number=residue_number,
                insertion_code=insertion, residue_name=residue_name, atom_name=atom_name,
                altloc=line[16:17].strip(), element=element, occupancy=occupancy, b_factor=b_factor,
                x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54]),
                serial_number=serial, locally_classified_protein=locally_protein,
                exclusion_reason=None if locally_protein else "NONPROTEIN_COMPONENT",
            ))
        except (ValueError, IndexError) as error:
            raise ParseError(f"Malformed atom record at source order {source_order}") from error
    if not parsed:
        raise ParseError("Selected protein model has no atom records")
    chosen, discarded = _select_altlocs(parsed)
    protein = [atom for atom in chosen if atom.locally_classified_protein]
    excluded = [atom for atom in chosen if not atom.locally_classified_protein]
    return ParsedProtein(protein, excluded, model_count, selected_model, discarded, inferred,
                         missing_atom_residues=missing_atom_residues)


def _first_model_number(path: Path) -> int:
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            try:
                return int(line[10:14].strip())
            except ValueError:
                return 1
        if line.startswith(("ATOM  ", "HETATM")):
            return 1
    return 1


def _select_altlocs(atoms: list[ProteinAtom]) -> tuple[list[ProteinAtom], int]:
    groups: dict[str, list[ProteinAtom]] = defaultdict(list)
    for atom in atoms:
        group_key = (atom.model_id, atom.record_type, atom.auth_chain_id, atom.auth_residue_number,
                     atom.insertion_code, atom.residue_name, atom.atom_name, atom.element)
        groups[str(group_key)].append(atom)
    selected: list[ProteinAtom] = []
    for candidates in groups.values():
        selected.append(sorted(candidates, key=lambda atom: (
            -(atom.occupancy if atom.occupancy is not None else -1.0),
            0 if not atom.altloc else 1, 0 if atom.altloc == "A" else 1,
            atom.altloc, atom.source_order,
        ))[0])
    selected.sort(key=lambda atom: atom.source_order)
    return selected, len(atoms) - len(selected)


def _infer_element(atom_name: str) -> str:
    letters = "".join(char for char in atom_name if char.isalpha())
    if not letters:
        return ""
    # PDB atom names are right-justified for one-letter elements; protein CA is carbon.
    return letters[0].upper()


def _parse_missing_atom_residue(line: str) -> tuple[str, int, str, str] | None:
    # PDB REMARK 470 data rows: residue, chain, sequence number/insertion, missing atoms.
    import re
    match = re.match(r"^REMARK 470\s+([A-Za-z0-9]{3})\s+(\S?)\s*(-?\d+)([A-Za-z]?)\s+\S+", line)
    if not match:
        return None
    residue, chain, number, insertion = match.groups()
    return chain, int(number), insertion, residue.upper()
