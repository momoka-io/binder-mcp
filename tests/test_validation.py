import pytest

from bio_mcp.core.sequence import MAX_PROTEIN_SEQUENCE_LENGTH, ProteinSequenceError
from bio_mcp.schemas import ProteinSequenceInput
from bio_mcp.tools.validation import validate_protein_sequence_tool


def test_validate_raw_sequence():
    output = validate_protein_sequence_tool(ProteinSequenceInput(text="acdefg"))

    assert output.cleaned_sequence == "ACDEFG"
    assert output.length == 6
    assert output.is_valid is True
    assert output.invalid_characters == []
    assert "converted to uppercase" in " ".join(output.warnings)


def test_validate_fasta_with_invalid_characters():
    output = validate_protein_sequence_tool(ProteinSequenceInput(text=">seq1\nACD-XZ*\n"))

    assert output.cleaned_sequence == "ACD"
    assert output.length == 3
    assert output.is_valid is False
    assert {item.character for item in output.invalid_characters} == {"-", "X", "Z", "*"}
    assert output.provenance.execution_mode == "local"


def test_validate_empty_input():
    output = validate_protein_sequence_tool(ProteinSequenceInput(text=" \n>empty\n"))

    assert output.cleaned_sequence == ""
    assert output.length == 0
    assert output.is_valid is False
    assert "No standard amino acid residues were found." in output.warnings


def test_validate_rejects_too_long_sequence():
    text = "A" * (MAX_PROTEIN_SEQUENCE_LENGTH + 1)

    with pytest.raises(ProteinSequenceError, match="exceeds maximum"):
        validate_protein_sequence_tool(ProteinSequenceInput(text=text))
