"""AAindex lookup, listing, and sequence feature extraction."""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources
from pathlib import Path
from statistics import mean

from bio_mcp import __version__
from bio_mcp.core.aaindex_parser import AaindexParseError, parse_aaindex1_file
from bio_mcp.core.sequence import require_valid_protein_sequence
from bio_mcp.schemas import (
    AaindexBackendStatus,
    AaindexListInput,
    AaindexListOutput,
    AaindexLookupInput,
    AaindexLookupOutput,
    AaindexMetadataRecord,
    AaindexRecord,
    AaindexSequenceFeaturesInput,
    AaindexSequenceFeaturesOutput,
    ResidueFeatureValue,
    SummaryStats,
    ToolProvenance,
)

LOOKUP_TOOL_NAME = "aaindex_lookup"
FEATURES_TOOL_NAME = "aaindex_sequence_features"
LIST_TOOL_NAME = "aaindex_list"
AAINDEX1_ENV_VAR = "BIO_MCP_AAINDEX1_PATH"
PACKAGED_FIXTURE_WARNING = (
    "Using a small packaged AAindex fixture. Set BIO_MCP_AAINDEX1_PATH to a local "
    "AAindex1 flat file for full local AAindex1 support."
)


@dataclass(frozen=True)
class AaindexBackend:
    """Resolved AAindex backend records and status."""

    records: tuple[AaindexRecord, ...]
    backend_source: str
    backend_path: str | None
    warnings: tuple[str, ...] = ()
    error: str | None = None

    @property
    def record_count(self) -> int:
        return len(self.records)

    @property
    def available(self) -> bool:
        return self.error is None and bool(self.records)


@lru_cache(maxsize=1)
def _load_packaged_records() -> tuple[AaindexRecord, ...]:
    data_path = resources.files("bio_mcp.data").joinpath("aaindex_minimal.json")
    with data_path.open("r", encoding="utf-8") as handle:
        records = json.load(handle)
    return tuple(AaindexRecord.model_validate(record) for record in records)


def resolve_aaindex_backend() -> AaindexBackend:
    """Resolve AAindex records from env-configured AAindex1 or packaged fallback."""

    env_path = os.environ.get(AAINDEX1_ENV_VAR)
    if env_path and env_path.strip():
        path = Path(env_path).expanduser()
        try:
            records = parse_aaindex1_file(path)
        except AaindexParseError as exc:
            return AaindexBackend(
                records=(),
                backend_source="env_file",
                backend_path=str(path),
                warnings=(str(exc),),
                error=str(exc),
            )
        return AaindexBackend(
            records=records,
            backend_source="env_file",
            backend_path=str(path),
        )

    return AaindexBackend(
        records=_load_packaged_records(),
        backend_source="packaged_fixture",
        backend_path=None,
        warnings=(PACKAGED_FIXTURE_WARNING,),
    )


def aaindex_backend_status() -> AaindexBackendStatus:
    """Return health/status details for the configured AAindex backend."""

    backend = resolve_aaindex_backend()
    return AaindexBackendStatus(
        backend_source=backend.backend_source,  # type: ignore[arg-type]
        backend_path=backend.backend_path,
        available=backend.available,
        record_count=backend.record_count,
        error=backend.error,
        warnings=list(backend.warnings),
    )


def aaindex_lookup_tool(params: AaindexLookupInput) -> AaindexLookupOutput:
    """Look up AAindex records by accession or text query from the local backend."""

    backend = resolve_aaindex_backend()
    query = (params.index_id or params.query or "").strip()
    warnings = list(backend.warnings)

    if backend.error:
        return _lookup_output(
            query=query,
            matches=[],
            backend=backend,
            warnings=warnings,
            error=backend.error,
        )

    index_id = (params.index_id or "").strip()
    text_query = (params.query or "").strip()
    if not index_id and not text_query:
        error = "Provide index_id or query."
        warnings.append(error)
        return _lookup_output(
            query=query,
            matches=[],
            backend=backend,
            warnings=warnings,
            error=error,
        )

    matches: list[AaindexRecord] = []
    seen: set[str] = set()
    if index_id:
        wanted = index_id.upper()
        for record in backend.records:
            if record.accession.upper() == wanted:
                matches.append(record)
                seen.add(record.accession.upper())
                break

    if text_query:
        lowered = text_query.lower()
        for record in backend.records:
            accession = record.accession.upper()
            if accession in seen:
                continue
            if _record_matches_query(record, lowered):
                matches.append(record)
                seen.add(accession)

    matches = matches[: params.limit]
    if not matches:
        warnings.append("No AAindex records matched the lookup.")

    return _lookup_output(
        query=query,
        matches=matches,
        backend=backend,
        warnings=warnings,
        error=None,
    )


def aaindex_list_tool(params: AaindexListInput) -> AaindexListOutput:
    """List lightweight AAindex metadata from the local backend."""

    backend = resolve_aaindex_backend()
    warnings = list(backend.warnings)
    records: list[AaindexRecord]
    if backend.error:
        records = []
    elif params.query and params.query.strip():
        query = params.query.strip().lower()
        records = [record for record in backend.records if _record_matches_query(record, query)]
    else:
        records = list(backend.records)

    records = records[: params.limit]
    metadata_records = [
        AaindexMetadataRecord(
            accession=record.accession,
            title=record.title,
            description=record.description,
            pmid=record.pmid,
            authors=record.authors,
            journal=record.journal,
            correlations=record.correlations,
        )
        for record in records
    ]
    return AaindexListOutput(
        records=metadata_records,
        record_count=backend.record_count,
        backend_source=backend.backend_source,  # type: ignore[arg-type]
        backend_path=backend.backend_path,
        warnings=warnings,
        error=backend.error,
        provenance=ToolProvenance(
            tool_name=LIST_TOOL_NAME,
            wrapper_version=__version__,
            backend=f"aaindex1:{backend.backend_source}",
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
    backend = resolve_aaindex_backend()
    warnings = [*validation.warnings, *backend.warnings]

    if backend.error:
        return _features_error_output(
            params=params,
            sequence=sequence,
            backend=backend,
            warnings=warnings,
            error=backend.error,
        )

    record = _find_record_by_accession(backend.records, params.index_id)
    if record is None:
        error = f"AAindex id {params.index_id!r} was not found in the configured backend."
        warnings.append(error)
        return _features_error_output(
            params=params,
            sequence=sequence,
            backend=backend,
            warnings=warnings,
            error=error,
        )
    if record.values is None:
        error = f"AAindex id {params.index_id!r} does not include amino acid values."
        warnings.append(error)
        return _features_error_output(
            params=params,
            sequence=sequence,
            backend=backend,
            warnings=warnings,
            error=error,
        )

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
        backend_source=backend.backend_source,  # type: ignore[arg-type]
        backend_path=backend.backend_path,
        record_count=backend.record_count,
        warnings=warnings,
        error=None,
        provenance=ToolProvenance(
            tool_name=FEATURES_TOOL_NAME,
            wrapper_version=__version__,
            backend=f"aaindex1:{backend.backend_source}",
            execution_mode="local",
            warnings=warnings,
        ),
    )


def _lookup_output(
    *,
    query: str,
    matches: list[AaindexRecord],
    backend: AaindexBackend,
    warnings: list[str],
    error: str | None,
) -> AaindexLookupOutput:
    return AaindexLookupOutput(
        query=query,
        matches=matches,
        backend_source=backend.backend_source,  # type: ignore[arg-type]
        backend_path=backend.backend_path,
        record_count=backend.record_count,
        warnings=warnings,
        error=error,
        provenance=ToolProvenance(
            tool_name=LOOKUP_TOOL_NAME,
            wrapper_version=__version__,
            backend=f"aaindex1:{backend.backend_source}",
            execution_mode="local",
            warnings=warnings,
        ),
    )


def _features_error_output(
    *,
    params: AaindexSequenceFeaturesInput,
    sequence: str,
    backend: AaindexBackend,
    warnings: list[str],
    error: str,
) -> AaindexSequenceFeaturesOutput:
    return AaindexSequenceFeaturesOutput(
        index_id=params.index_id,
        sequence=sequence,
        length=len(sequence),
        per_residue_values=[],
        summary_stats=SummaryStats(mean=None, min=None, max=None),
        missing_residues=[],
        backend_source=backend.backend_source,  # type: ignore[arg-type]
        backend_path=backend.backend_path,
        record_count=backend.record_count,
        warnings=warnings,
        error=error,
        provenance=ToolProvenance(
            tool_name=FEATURES_TOOL_NAME,
            wrapper_version=__version__,
            backend=f"aaindex1:{backend.backend_source}",
            execution_mode="local",
            warnings=warnings,
        ),
    )


def _find_record_by_accession(
    records: tuple[AaindexRecord, ...], index_id: str
) -> AaindexRecord | None:
    wanted = index_id.upper()
    for record in records:
        if record.accession.upper() == wanted:
            return record
    return None


def _record_matches_query(record: AaindexRecord, query: str) -> bool:
    searchable = [
        record.accession,
        record.title,
        record.description or "",
        record.pmid or "",
        record.authors or "",
        record.journal or "",
        " ".join(record.correlations),
    ]
    return any(query in value.lower() for value in searchable)
