from biosensia_pocket_library.reference_links import merge_affinity_reference_links


def test_reference_candidates_merge_sources_and_keep_strongest_role():
    base = {
        "measurement_id": "m", "complex_id": "c", "citation_id": "citation",
        "confidence": 0.3, "evidence_note": "Depositor citation",
        "automatic_or_manual": "automatic", "verified_by": None, "verified_at_utc": None,
    }
    rows = [
        {**base, "candidate_status": "structural_reference_only", "evidence_sources": ["citation"]},
        {**base, "candidate_status": "probable_structural_reference", "confidence": 0.6,
         "evidence_sources": ["primary_citation"]},
    ]

    merged = merge_affinity_reference_links(rows)

    assert len(merged) == 1
    assert merged[0]["candidate_status"] == "probable_structural_reference"
    assert merged[0]["confidence"] == 0.6
    assert merged[0]["evidence_sources"] == ["citation", "primary_citation"]
