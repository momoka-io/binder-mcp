"""Pydantic schemas for bio-mcp tool inputs and outputs."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

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
        description="Case-insensitive text search across AAindex metadata.",
    )
    limit: int = Field(default=20, ge=1, le=500)


class AaindexRecord(BaseModel):
    """AAindex metadata and values for one record."""

    accession: str
    title: str
    description: str | None = None
    pmid: str | None = None
    authors: str | None = None
    journal: str | None = None
    correlations: dict[str, float | None] = Field(default_factory=dict)
    values: dict[str, float | None] | None = None
    value_order: str = STANDARD_AMINO_ACIDS


class AaindexLookupOutput(BaseModel):
    """AAindex lookup output."""

    query: str
    matches: list[AaindexRecord]
    backend_source: Literal["env_file", "packaged_fixture"]
    backend_path: str | None = None
    record_count: int
    warnings: list[str] = Field(default_factory=list)
    error: str | None = None
    provenance: ToolProvenance


class AaindexListInput(BaseModel):
    """AAindex lightweight listing input."""

    query: str | None = Field(
        default=None,
        description="Optional case-insensitive text search across AAindex metadata.",
    )
    limit: int = Field(default=50, ge=1, le=500)


class AaindexMetadataRecord(BaseModel):
    """Lightweight AAindex metadata without amino acid values."""

    accession: str
    title: str
    description: str | None = None
    pmid: str | None = None
    authors: str | None = None
    journal: str | None = None
    correlations: dict[str, float | None] = Field(default_factory=dict)


class AaindexListOutput(BaseModel):
    """AAindex listing output."""

    records: list[AaindexMetadataRecord]
    record_count: int
    backend_source: Literal["env_file", "packaged_fixture"]
    backend_path: str | None = None
    warnings: list[str] = Field(default_factory=list)
    error: str | None = None
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
    backend_source: Literal["env_file", "packaged_fixture"]
    backend_path: str | None = None
    record_count: int
    warnings: list[str] = Field(default_factory=list)
    error: str | None = None
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


class BlastpLocalInput(BaseModel):
    """Input for local BLASTP search."""

    query_fasta: str = Field(..., description="Protein FASTA query text to search locally.")
    db_path: str = Field(..., description="Local BLAST protein database prefix.")
    evalue: float = Field(default=1e-5, gt=0)
    max_target_seqs: int = Field(default=10, ge=1, le=10_000)
    timeout_sec: int = Field(default=120, ge=1, le=3_600)


class PsiblastLocalInput(BlastpLocalInput):
    """Input for local PSI-BLAST search."""

    num_iterations: int = Field(default=3, ge=1, le=10)
    timeout_sec: int = Field(default=300, ge=1, le=3_600)


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


class BlastTabularHit(BaseModel):
    """One parsed BLAST outfmt 6 hit."""

    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float


class BlastpLocalOutput(BaseModel):
    """Output from a local BLASTP wrapper."""

    hits: list[BlastTabularHit]
    raw_tabular: str
    number_of_hits: int
    query_count: int
    db_path: str
    command_version: str | None
    stderr_summary: str
    warnings: list[str] = Field(default_factory=list)
    error: str | None = None
    provenance: ToolProvenance


class PsiblastLocalOutput(BlastpLocalOutput):
    """Output from a local PSI-BLAST wrapper."""

    num_iterations: int


class OptionalDependencyStatus(BaseModel):
    """Availability and version for an optional dependency."""

    installed: bool
    version: str | None = None


class CommandDependencyStatus(BaseModel):
    """Availability, path, and version for a local CLI binary."""

    available: bool
    path: str | None = None
    version: str | None = None
    resolution_source: Literal["env_override", "PATH", "missing"] = "missing"
    error: str | None = None


class AaindexBackendStatus(BaseModel):
    """Availability details for the configured AAindex1 backend."""

    backend_source: Literal["env_file", "packaged_fixture"]
    backend_path: str | None = None
    available: bool
    record_count: int
    error: str | None = None
    warnings: list[str] = Field(default_factory=list)


class HealthOutput(BaseModel):
    """Server health/status output."""

    package_version: str
    python_version: str
    optional_dependencies: dict[str, OptionalDependencyStatus]
    cli_binaries: dict[str, CommandDependencyStatus]
    aaindex_backend: AaindexBackendStatus
    available_tools: list[str]
    execution_mode: Literal["local"] = "local"
    warnings: list[str] = Field(default_factory=list)
    provenance: ToolProvenance
