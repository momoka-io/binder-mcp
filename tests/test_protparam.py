import pytest

from bio_mcp.core.sequence import MAX_PROTEIN_SEQUENCE_LENGTH, ProteinSequenceError
from bio_mcp.schemas import ProteinSequenceInput
from bio_mcp.tools.protparam import protparam_analyze_tool


def test_protparam_valid_sequence():
    output = protparam_analyze_tool(ProteinSequenceInput(text="ACDEFWY"))

    assert output.length == 7
    assert output.amino_acid_composition["A"] == 1
    assert output.amino_acid_composition["R"] == 0
    assert output.amino_acid_percent["A"] == pytest.approx(1 / 7)
    assert output.molecular_weight is not None
    assert output.aromaticity == pytest.approx(3 / 7)
    assert output.gravy is not None
    assert output.provenance.execution_mode == "local"


def test_protparam_rejects_invalid_characters():
    with pytest.raises(ProteinSequenceError, match="Invalid protein sequence characters"):
        protparam_analyze_tool(ProteinSequenceInput(text="ACDX"))


def test_protparam_rejects_empty_input():
    with pytest.raises(ProteinSequenceError, match="empty"):
        protparam_analyze_tool(ProteinSequenceInput(text=">empty\n"))


def test_protparam_rejects_too_long_sequence():
    text = "A" * (MAX_PROTEIN_SEQUENCE_LENGTH + 1)

    with pytest.raises(ProteinSequenceError, match="exceeds maximum"):
        protparam_analyze_tool(ProteinSequenceInput(text=text))
