"""PDBbind PL-index parser with mandatory legacy-token redaction."""

from __future__ import annotations

import math
import re
from pathlib import Path

from .exceptions import ParseError, SourceIntegrityError
from .hashing import canonical_json_hash, sha256_bytes, stable_id
from .models import BindingMeasurement, IndexRecord
from .scrub import remove_pdf_tokens

DATA_RE = re.compile(r"^\s*(\S+)\s+(\S+)\s+(\d{4})\s+(.+?)(?:\s*//\s*(.*))?$")
MEASUREMENT_RE = re.compile(
    r"^\s*([A-Za-z][A-Za-z0-9]*)\s*(<=|>=|=|<|>|~|≈)?\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(\S*)\s*$"
)
CONCENTRATION_FACTORS = {
    "M": 1.0, "mM": 1e-3, "uM": 1e-6, "µM": 1e-6,
    "nM": 1e-9, "pM": 1e-12, "fM": 1e-15,
}
KNOWN_CONCENTRATIONS = {"KD": "Kd", "KI": "Ki", "IC50": "IC50", "EC50": "EC50"}
REVERSE_RELATION = {"<": ">", ">": "<", "<=": ">=", ">=": "<=", "=": "="}


def parse_measurement(raw: str, complex_id: str, pdb_id: str, line_number: int) -> BindingMeasurement:
    lexical = raw.strip()
    textual_approximate = bool(re.match(r"(?i)^(?:approx(?:imately)?|ca\.?)\s+", lexical))
    if textual_approximate:
        lexical = re.sub(r"(?i)^(?:approx(?:imately)?|ca\.?)\s+", "", lexical, count=1)
    match = MEASUREMENT_RE.match(lexical)
    base = dict(
        measurement_id=stable_id("measurement", complex_id, 0, raw), complex_id=complex_id,
        pdb_id=pdb_id, measurement_type_raw="", measurement_type_normalized=None,
        relation_raw=None, relation_normalized=None, value_raw=None, value_numeric=None,
        unit_raw=None, unit_normalized=None, value_molar=None, value_inverse_molar=None,
        p_measurement_name=None, p_relation=None, p_value=None, normalization_kind=None,
        measurement_raw=raw, parse_status="malformed", parse_warning_codes=[],
        source_index_line_number=line_number,
    )
    if not raw.strip():
        base["parse_status"] = "missing"
        return BindingMeasurement(**base)
    if not match:
        return BindingMeasurement(**base)
    kind, relation, value_raw, unit = match.groups()
    base.update(measurement_type_raw=kind, relation_raw=relation or "=", value_raw=value_raw,
                unit_raw=unit or None)
    try:
        value = float(value_raw)
    except (ValueError, OverflowError):
        return BindingMeasurement(**base)
    if not math.isfinite(value) or value <= 0:
        base["parse_warning_codes"] = ["NONPOSITIVE_OR_NONFINITE_MEASUREMENT"]
        return BindingMeasurement(**base)
    base["value_numeric"] = value
    normalized_kind = KNOWN_CONCENTRATIONS.get(kind.upper())
    approximate = textual_approximate or relation in {"~", "≈"}
    if normalized_kind:
        base["measurement_type_normalized"] = normalized_kind
        if unit not in CONCENTRATION_FACTORS:
            base["parse_status"] = "unsupported_unit"
            return BindingMeasurement(**base)
        molar = value * CONCENTRATION_FACTORS[unit]
        normalized_relation = "approximately" if approximate else (relation or "=")
        p_relation = "approximately" if approximate else REVERSE_RELATION.get(relation or "=", "=")
        base.update(unit_normalized="M", value_molar=molar, normalization_kind="concentration",
                    p_measurement_name=f"p{normalized_kind}", p_value=-math.log10(molar),
                    relation_normalized=normalized_relation, p_relation=p_relation)
    elif kind.upper() == "KA":
        base["measurement_type_normalized"] = "Ka"
        inverse = _inverse_molar_factor(unit)
        if inverse is None:
            base["parse_status"] = "unsupported_unit"
            return BindingMeasurement(**base)
        value_inverse = value * inverse
        base.update(unit_normalized="M^-1", value_inverse_molar=value_inverse,
                    normalization_kind="association_constant", p_measurement_name="pKa_association",
                    p_value=math.log10(value_inverse), relation_normalized=("approximately" if approximate else relation or "="),
                    p_relation=("approximately" if approximate else relation or "="))
    else:
        base["parse_status"] = "unsupported_measurement_type"
        return BindingMeasurement(**base)
    base["parse_status"] = "parsed_approximate" if approximate else (
        "parsed_censored" if (relation or "=") != "=" else "parsed_exact"
    )
    return BindingMeasurement(**base)


def _inverse_molar_factor(unit: str) -> float | None:
    compact = unit.replace("−", "-").replace("⁻", "-").replace(" ", "")
    match = re.fullmatch(r"(M|mM|uM|µM|nM|pM|fM)(?:\^-?1|-1)", compact)
    if not match:
        return None
    return 1.0 / CONCENTRATION_FACTORS[match.group(1)]


def parse_index(path: Path, distribution_id: str) -> tuple[list[IndexRecord], list[dict], dict]:
    raw_bytes = path.read_bytes()
    distribution_hash = canonical_json_hash({"distribution_id": distribution_id})[:8]
    records: dict[str, IndexRecord] = {}
    canonical_lines: dict[str, str] = {}
    occurrences: list[dict] = []
    declared_count: int | None = None
    data_lines = 0
    for line_number, line_bytes in enumerate(raw_bytes.splitlines(keepends=False), 1):
        line = line_bytes.decode("utf-8", errors="replace")
        if line.lstrip().startswith("#"):
            count_match = re.search(r"#\s*(\d+)\s+protein-ligand complexes", line, re.I)
            if count_match:
                declared_count = int(count_match.group(1))
            continue
        if not line.strip():
            continue
        data_lines += 1
        match = DATA_RE.match(line)
        if not match:
            raise ParseError(f"Malformed index line {line_number}")
        pdb_id, resolution, year_text, measurement_raw, post = match.groups()
        measurement_raw = measurement_raw.strip()
        pdb_id = pdb_id.lower()
        redacted_post = remove_pdf_tokens(post or "")
        ligand_label, comment, warning = _parse_post(redacted_post)
        redacted = f"{pdb_id} {resolution} {year_text} {measurement_raw} // {redacted_post}".rstrip()
        digest = sha256_bytes(line_bytes)
        complex_id = f"pb20v24p-{distribution_hash}:{pdb_id}"
        measurement = parse_measurement(measurement_raw, complex_id, pdb_id, line_number)
        if warning:
            measurement.parse_warning_codes.append(warning)
        record = IndexRecord(
            complex_id=complex_id, distribution_id=distribution_id, pdb_id=pdb_id,
            resolution_raw=resolution,
            resolution_angstrom=float(resolution) if _is_float(resolution) else None,
            experimental_method_hint=None if _is_float(resolution) else resolution,
            release_year=int(year_text), ligand_label=ligand_label, index_comment=comment,
            index_line_redacted=redacted, source_line_sha256=digest,
            primary_index_line_number=line_number, occurrence_line_numbers=[line_number],
            measurement=measurement,
        )
        occurrences.append({"complex_id": complex_id, "source_index_line_number": line_number,
                            "source_line_sha256": digest})
        canonical = canonical_json_hash({"pdb_id": pdb_id, "resolution": resolution, "year": year_text,
                                         "measurement": measurement_raw, "post": redacted_post})
        if pdb_id in records:
            if canonical_lines[pdb_id] != canonical:
                raise SourceIntegrityError(f"Conflicting duplicate PDB ID {pdb_id}")
            records[pdb_id].occurrence_line_numbers.append(line_number)
        else:
            records[pdb_id] = record
            canonical_lines[pdb_id] = canonical
    summary = {"declared_count": declared_count, "physical_data_line_count": data_lines,
               "unique_complex_count": len(records), "equivalent_duplicate_count": data_lines - len(records)}
    return sorted(records.values(), key=lambda item: item.pdb_id), occurrences, summary


def _is_float(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return False
    return True


def _parse_post(text: str) -> tuple[str | None, str | None, str | None]:
    start = text.find("(")
    if start < 0:
        return None, text or None, None
    depth = 0
    for position in range(start, len(text)):
        if text[position] == "(":
            depth += 1
        elif text[position] == ")":
            depth -= 1
            if depth == 0:
                label = text[start + 1:position].strip()
                comment = (text[:start] + " " + text[position + 1:]).strip()
                return label or None, comment or None, None
    return None, text or None, "UNBALANCED_LIGAND_LABEL"
