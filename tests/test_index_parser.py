from pathlib import Path

import pytest

from biosensia_pocket_library.exceptions import SourceIntegrityError
from biosensia_pocket_library.index_parser import parse_index, parse_measurement
from biosensia_pocket_library.scrub import remove_pdf_tokens, scrub


def test_pdf_tokens_are_removed_everywhere_after_separator():
    text = r"C:\refs\ONE.PDF, (LIG) note /tmp/two.pdf; another.pdf"
    redacted = remove_pdf_tokens(text)
    assert ".pdf" not in redacted.lower()
    assert redacted == "(LIG) note"
    assert ".pdf" not in scrub({"message": text})["message"].lower()
    assert remove_pdf_tokens("one.pdf / (LIG)") == "(LIG)"
    assert remove_pdf_tokens("/ one.pdf (LIG)") == "(LIG)"


@pytest.mark.parametrize(
    ("raw", "status", "molar", "inverse", "p_relation"),
    [
        ("Kd=49uM", "parsed_exact", 49e-6, None, "="),
        ("Ki<=17nM", "parsed_censored", 17e-9, None, ">="),
        ("IC50≈5 nM", "parsed_approximate", 5e-9, None, "approximately"),
        ("approximately Kd=5nM", "parsed_approximate", 5e-9, None, "approximately"),
        ("Ka=1.2e6M-1", "parsed_exact", None, 1.2e6, "="),
    ],
)
def test_binding_measurement_normalization(raw, status, molar, inverse, p_relation):
    value = parse_measurement(raw, "complex", "1abc", 8)
    assert value.parse_status == status
    assert value.value_molar == pytest.approx(molar) if molar is not None else value.value_molar is None
    assert value.value_inverse_molar == pytest.approx(inverse) if inverse is not None else value.value_inverse_molar is None
    assert value.p_relation == p_relation


def test_index_redaction_label_and_equivalent_duplicate(tmp_path: Path):
    path = tmp_path / "index.lst"
    path.write_text("# 2 protein-ligand complexes in total\n1ABC 2.0 2001 Kd=1uM // x.pdf (LIG) note\n"
                    "1abc 2.0 2001 Kd=1uM // Y.PDF (LIG) note\n", encoding="utf-8")
    records, occurrences, summary = parse_index(path, "distribution")
    assert len(records) == 1
    assert len(occurrences) == 2
    assert records[0].ligand_label == "LIG"
    assert records[0].index_comment == "note"
    assert ".pdf" not in records[0].index_line_redacted.lower()
    assert summary["equivalent_duplicate_count"] == 1


def test_conflicting_duplicate_is_rejected(tmp_path: Path):
    path = tmp_path / "index.lst"
    path.write_text("1abc 2.0 2001 Kd=1uM // a.pdf (LIG)\n1abc 2.0 2001 Kd=2uM // a.pdf (LIG)\n")
    with pytest.raises(SourceIntegrityError):
        parse_index(path, "distribution")
