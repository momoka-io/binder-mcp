"""Protein sequence validation tool."""

from __future__ import annotations

from bio_mcp import __version__
from bio_mcp.core.sequence import parse_protein_sequence
from bio_mcp.schemas import (
    InvalidCharacter,
    ProteinSequenceInput,
    ToolProvenance,
    ValidateProteinSequenceOutput,
)

TOOL_NAME = "validate_protein_sequence"


def validate_protein_sequence_tool(
    params: ProteinSequenceInput,
) -> ValidateProteinSequenceOutput:
    """Validate raw protein sequence or FASTA text without running remote services."""

    validation = parse_protein_sequence(params.text, max_length=params.max_length)
    warnings = list(validation.warnings)
    invalid_characters = [
        InvalidCharacter(
            character=report.character,
            count=report.count,
            positions=list(report.positions),
        )
        for report in validation.invalid_characters
    ]
    return ValidateProteinSequenceOutput(
        cleaned_sequence=validation.cleaned_sequence,
        length=validation.length,
        invalid_characters=invalid_characters,
        is_valid=validation.is_valid,
        warnings=warnings,
        provenance=ToolProvenance(
            tool_name=TOOL_NAME,
            wrapper_version=__version__,
            backend="bio_mcp.core.sequence",
            execution_mode="local",
            warnings=warnings,
        ),
    )
