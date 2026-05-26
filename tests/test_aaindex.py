import pytest
from pydantic import ValidationError

from bio_mcp.core.sequence import MAX_PROTEIN_SEQUENCE_LENGTH, ProteinSequenceError
from bio_mcp.schemas import AaindexLookupInput, AaindexSequenceFeaturesInput
from bio_mcp.tools.aaindex import aaindex_lookup_tool, aaindex_sequence_features_tool


def test_aaindex_lookup_by_id():
    output = aaindex_lookup_tool(AaindexLookupInput(index_id="KYTJ820101"))

    assert output.matches
    record = output.matches[0]
    assert record.accession == "KYTJ820101"
    assert record.values is not None
    assert len(record.values) == 20
    assert record.values["A"] == pytest.approx(1.8)


def test_aaindex_lookup_by_query():
    output = aaindex_lookup_tool(AaindexLookupInput(query="hydrophilicity"))

    assert [record.accession for record in output.matches] == ["HOPT810101"]


def test_aaindex_lookup_requires_id_or_query():
    with pytest.raises(ValidationError, match="Provide index_id or query"):
        AaindexLookupInput()


def test_aaindex_lookup_unknown_returns_warning():
    output = aaindex_lookup_tool(AaindexLookupInput(index_id="NOPE000000"))

    assert output.matches == []
    assert "No packaged AAindex records matched the lookup." in output.warnings


def test_aaindex_sequence_features_valid_sequence():
    output = aaindex_sequence_features_tool(
        AaindexSequenceFeaturesInput(text="ACD", index_id="KYTJ820101")
    )

    assert output.index_id == "KYTJ820101"
    assert output.sequence == "ACD"
    assert [item.value for item in output.per_residue_values] == [1.8, 2.5, -3.5]
    assert output.summary_stats.mean == pytest.approx((1.8 + 2.5 - 3.5) / 3)
    assert output.summary_stats.min == pytest.approx(-3.5)
    assert output.summary_stats.max == pytest.approx(2.5)
    assert output.missing_residues == []


def test_aaindex_sequence_features_rejects_invalid_sequence():
    with pytest.raises(ProteinSequenceError, match="Invalid protein sequence characters"):
        aaindex_sequence_features_tool(
            AaindexSequenceFeaturesInput(text="ACD*", index_id="KYTJ820101")
        )


def test_aaindex_sequence_features_rejects_empty_sequence():
    with pytest.raises(ProteinSequenceError, match="empty"):
        aaindex_sequence_features_tool(
            AaindexSequenceFeaturesInput(text="", index_id="KYTJ820101")
        )


def test_aaindex_sequence_features_rejects_too_long_sequence():
    text = "A" * (MAX_PROTEIN_SEQUENCE_LENGTH + 1)

    with pytest.raises(ProteinSequenceError, match="exceeds maximum"):
        aaindex_sequence_features_tool(
            AaindexSequenceFeaturesInput(text=text, index_id="KYTJ820101")
        )


def test_aaindex_sequence_features_unknown_index():
    with pytest.raises(ValueError, match="was not found"):
        aaindex_sequence_features_tool(
            AaindexSequenceFeaturesInput(text="ACD", index_id="NOPE000000")
        )
