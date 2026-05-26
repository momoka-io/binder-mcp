"""Small packaged AAindex lookup and sequence feature extraction."""

from __future__ import annotations

import json
from functools import lru_cache
from importlib import resources
from statistics import mean

from bio_mcp import __version__
from bio_mcp.core.sequence import require_valid_protein_sequence
from bio_mcp.schemas import (
    AaindexLookupInput,
    AaindexLookupOutput,
    AaindexRecord,
    AaindexSequenceFeaturesInput,
    AaindexSequenceFeaturesOutput,
    ResidueFeatureValue,
    SummaryStats,
    ToolProvenance,
)

LOOKUP_TOOL_NAME = "aaindex_lookup"
FEATURES_TOOL_NAME = "aaindex_sequence_features"


@lru_cache(maxsize=1)
def _load_records() -> tuple[AaindexRecord, ...]:
    data_path = resources.files("bio_mcp.data").joinpath("aaindex_minimal.json")
    with data_path.open("r", encoding="utf-8") as handle:
        records = json.load(handle)
    return tuple(AaindexRecord.model_validate(record) for record in records)


def aaindex_lookup_tool(params: AaindexLookupInput) -> AaindexLookupOutput:
    """Look up AAindex records from the packaged Phase 1 fixture."""

    query = (params.index_id or params.query or "").strip()
    query_upper = query.upper()
    warnings: list[str] = []
    records = _load_records()

    if params.index_id and params.index_id.strip():
        matches = [record for record in records if record.accession.upper() == query_upper]
    else:
        query_lower = query.lower()
        matches = [
            record
            for record in records
            if query_lower in record.accession.lower()
            or query_lower in record.title.lower()
            or (record.description is not None and query_lower in record.description.lower())
            or (record.authors is not None and query_lower in record.authors.lower())
        ]

    matches = matches[: params.limit]
    if not matches:
        warnings.append("No packaged AAindex records matched the lookup.")

    warnings.append("Using a small packaged AAindex fixture; full AAindex parsing is planned later.")
    return AaindexLookupOutput(
        query=query,
        matches=matches,
        warnings=warnings,
        provenance=ToolProvenance(
            tool_name=LOOKUP_TOOL_NAME,
            wrapper_version=__version__,
            backend="packaged-aaindex-fixture",
            execution_mode="local",
            warnings=warnings,
        ),
    )


def aaindex_sequence_features_tool(
    params: AaindexSequenceFeaturesInput,
) -> AaindexSequenceFeaturesOutput:
    """Extract per-residue AAindex values and summary stats for one protein sequence."""

    validation = require_valid_protein_sequence(params.text, max_length=params.max_length)
    sequence = validation.cleaned_sequence
    lookup = aaindex_lookup_tool(AaindexLookupInput(index_id=params.index_id, limit=1))
    warnings = [*validation.warnings, *lookup.warnings]

    if not lookup.matches:
        raise ValueError(f"AAindex id {params.index_id!r} was not found in the packaged fixture.")

    record = lookup.matches[0]
    if record.values is None:
        raise ValueError(f"AAindex id {params.index_id!r} does not include numeric amino acid values.")

    per_residue_values: list[ResidueFeatureValue] = []
    numeric_values: list[float] = []
    missing_residues: set[str] = set()
    for position, residue in enumerate(sequence, start=1):
        value = record.values.get(residue)
        if value is None:
            missing_residues.add(residue)
        else:
            numeric_values.append(value)
        per_residue_values.append(
            ResidueFeatureValue(position=position, residue=residue, value=value)
        )

    if numeric_values:
        summary_stats = SummaryStats(
            mean=mean(numeric_values),
            min=min(numeric_values),
            max=max(numeric_values),
        )
    else:
        summary_stats = SummaryStats(mean=None, min=None, max=None)
        warnings.append("No AAindex values were available for residues in this sequence.")

    if missing_residues:
        warnings.append("Some residues did not have values in the selected AAindex record.")

    return AaindexSequenceFeaturesOutput(
        index_id=record.accession,
        sequence=sequence,
        length=len(sequence),
        per_residue_values=per_residue_values,
        summary_stats=summary_stats,
        missing_residues=sorted(missing_residues),
        warnings=warnings,
        provenance=ToolProvenance(
            tool_name=FEATURES_TOOL_NAME,
            wrapper_version=__version__,
            backend="packaged-aaindex-fixture",
            execution_mode="local",
            warnings=warnings,
        ),
    )
