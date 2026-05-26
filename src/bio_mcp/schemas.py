"""Pydantic schemas for bio-mcp tool inputs and outputs."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field, model_validator

from bio_mcp.core.sequence import (
    MAX_ALIGNMENT_SEQUENCE_LENGTH,
    MAX_ALIGNMENT_SEQUENCE_COUNT,
    MAX_PROTEIN_SEQUENCE_LENGTH,
    STANDARD_AMINO_ACIDS,
)


class ToolProvenance(BaseModel):
    """Provenance metadata returned by tool calls."""

    tool_name: str
    wrapper_version: str
    backend: str | None = None
    backend_version: str | None = None
    execution_mode: Literal["local"] = "local"
    warnings: list[str] = Field(default_factory=list)


class InvalidCharacter(BaseModel):
    """Invalid sequence character report."""

    character: str
    count: int
    positions: list[int]


class ProteinSequenceInput(BaseModel):
    """Raw protein sequence or FASTA text."""

    text: str = Field(..., description="Raw protein sequence or FASTA text.")
    max_length: int = Field(
        default=MAX_PROTEIN_SEQUENCE_LENGTH,
        ge=1,
        le=MAX_PROTEIN_SEQUENCE_LENGTH,
        description="Maximum parsed sequence length accepted by this local tool.",
    )


class ValidateProteinSequenceOutput(BaseModel):
    """Validation output for a protein sequence."""

    cleaned_sequence: str
    length: int
    invalid_characters: list[InvalidCharacter]
    is_valid: bool
    warnings: list[str] = Field(default_factory=list)
    provenance: ToolProvenance


class ProtParamOutput(BaseModel):
    """Protein property analysis output."""

    length: int
    amino_acid_composition: dict[str, int]
    amino_acid_percent: dict[str, float]
    molecular_weight: float | None
    theoretical_pi: float | None
    aromaticity: float | None
    instability_index: float | None
    gravy: float | None
    unsupported_metrics: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    provenance: ToolProvenance


class AaindexLookupInput(BaseModel):
    """AAindex lookup by accession or text query."""

    index_id: str | None = Field(
        default=None,
        description="AAindex accession/id, for example KYTJ820101.",
    )
    query: str | None = Field(
        default=None,
        description="Case-insensitive text search across packaged AAindex metadata.",
    )
    limit: int = Field(default=10, ge=1, le=25)

    @model_validator(mode="after")
    def require_lookup_value(self) -> "AaindexLookupInput":
        if not (self.index_id and self.index_id.strip()) and not (self.query and self.query.strip()):
            raise ValueError("Provide index_id or query.")
        return self


class AaindexRecord(BaseModel):
    """AAindex metadata and values for one record."""

    accession: str
    title: str
    description: str | None = None
    authors: str | None = None
    journal: str | None = None
    values: dict[str, float] | None = None
    value_order: str = STANDARD_AMINO_ACIDS


class AaindexLookupOutput(BaseModel):
    """AAindex lookup output."""

    query: str
    matches: list[AaindexRecord]
    warnings: list[str] = Field(default_factory=list)
    provenance: ToolProvenance


class AaindexSequenceFeaturesInput(ProteinSequenceInput):
    """Protein sequence and AAindex accession for per-residue features."""

    index_id: str = Field(..., description="AAindex accession/id with numeric amino acid values.")


class ResidueFeatureValue(BaseModel):
    """Per-residue AAindex value."""

    position: int
    residue: str
    value: float | None


class SummaryStats(BaseModel):
    """Summary statistics over available numeric values."""

    mean: float | None
    min: float | None
    max: float | None


class AaindexSequenceFeaturesOutput(BaseModel):
    """AAindex per-sequence feature output."""

    index_id: str
    sequence: str
    length: int
    per_residue_values: list[ResidueFeatureValue]
    summary_stats: SummaryStats
    missing_residues: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    provenance: ToolProvenance


class MafftAlignInput(BaseModel):
    """Input for local MAFFT alignment."""

    fasta_text: str = Field(..., description="Protein FASTA text to align locally.")
    mode: Literal["auto", "localpair", "globalpair", "genafpair"] = Field(
        default="auto",
        description="MAFFT alignment mode.",
    )
    timeout_sec: int = Field(default=60, ge=1, le=600)


class ClustaloAlignInput(BaseModel):
    """Input for local Clustal Omega alignment."""

    fasta_text: str = Field(..., description="Protein FASTA text to align locally.")
    timeout_sec: int = Field(default=60, ge=1, le=600)


class AlignmentOutput(BaseModel):
    """Output from a local multiple-sequence alignment wrapper."""

    aligned_fasta: str | None
    number_of_sequences: int
    command_version: str | None
    stderr_summary: str
    warnings: list[str] = Field(default_factory=list)
    error: str | None = None
    max_sequences: int = MAX_ALIGNMENT_SEQUENCE_COUNT
    max_sequence_length: int = MAX_ALIGNMENT_SEQUENCE_LENGTH
    provenance: ToolProvenance


class OptionalDependencyStatus(BaseModel):
    """Availability and version for an optional dependency."""

    installed: bool
    version: str | None = None


class CommandDependencyStatus(BaseModel):
    """Availability, path, and version for a local CLI binary."""

    available: bool
    path: str | None = None
    version: str | None = None
    error: str | None = None


class HealthOutput(BaseModel):
    """Server health/status output."""

    package_version: str
    python_version: str
    optional_dependencies: dict[str, OptionalDependencyStatus]
    cli_binaries: dict[str, CommandDependencyStatus]
    available_tools: list[str]
    execution_mode: Literal["local"] = "local"
    warnings: list[str] = Field(default_factory=list)
    provenance: ToolProvenance
