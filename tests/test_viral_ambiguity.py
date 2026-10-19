from viral_platform.analysis.viral_gene_detection import (
    apply_virus_assignments,
    derive_detection_state,
    find_ambiguous_features,
)


def _per_virus():
    return {
        "SARS-CoV-2": {"M": {"CoV2_M", "zika_M"}, "ORF8": {"CoV2_ORF8"}},
        "Zika": {"M": {"CoV2_M", "zika_M"}},
    }


def test_shared_features_are_ambiguous():
    ambiguous = find_ambiguous_features(_per_virus())
    assert set(ambiguous) == {"CoV2_M", "zika_M"}


def test_assignment_splits_same_named_genes():
    result = apply_virus_assignments(
        _per_virus(), {"CoV2_M": "SARS-CoV-2", "zika_M": "Zika"}
    )
    assert result["SARS-CoV-2"]["M (SARS-CoV-2)"] == {"CoV2_M"}
    assert result["Zika"]["M (Zika)"] == {"zika_M"}
    assert "M" not in result["Zika"]
    state = derive_detection_state(result)
    assert "M" not in state["matched_features_by_gene"]


def test_single_assignment_keeps_plain_name():
    per_virus = {"A": {"M": {"M"}}, "B": {"M": {"M"}}}
    result = apply_virus_assignments(per_virus, {"M": "A"})
    assert result == {"A": {"M": {"M"}}}
