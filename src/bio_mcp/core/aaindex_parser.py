"""Parser for local AAindex1 flat files."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

from bio_mcp.core.sequence import STANDARD_AMINO_ACID_SET, STANDARD_AMINO_ACIDS
from bio_mcp.schemas import AaindexRecord

AAINDEX1_FIELDS = {"H", "D", "R", "A", "T", "J", "C", "I"}
MISSING_VALUE_TOKENS = {"NA", "N/A", "-", "."}
_NUMBER_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?$")


class AaindexParseError(ValueError):
    """Raised when an AAindex1 flat file or entry cannot be parsed."""


def parse_aaindex1_file(path: str | Path) -> tuple[AaindexRecord, ...]:
    """Parse AAindex1 records from a local flat file."""

    file_path = Path(path)
    try:
        text = file_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise AaindexParseError(f"Could not read AAindex1 file {file_path}: {exc}") from exc
    return parse_aaindex1_text(text)


def parse_aaindex1_text(text: str) -> tuple[AaindexRecord, ...]:
    """Parse AAindex1 records from flat-file text."""

    records: list[AaindexRecord] = []
    for entry_number, raw_entry in enumerate(text.split("//"), start=1):
        entry = raw_entry.strip("\n")
        if not entry.strip():
            continue
        records.append(_parse_entry(entry, entry_number=entry_number))
    if not records:
        raise AaindexParseError("AAindex1 file did not contain any records.")
    return tuple(records)


def _parse_entry(entry: str, *, entry_number: int) -> AaindexRecord:
    fields = _collect_fields(entry)
    accession = _field_text(fields, "H").strip()
    if not accession:
        raise AaindexParseError(f"AAindex1 entry {entry_number} is missing required H accession.")

    values = _parse_i_field(fields.get("I", []), accession=accession)
    if not values:
        raise AaindexParseError(f"AAindex1 entry {accession} is missing required I values.")

    title = _field_text(fields, "D").strip() or accession
    return AaindexRecord(
        accession=accession,
        title=title,
        description=title,
        pmid=_field_text(fields, "R") or None,
        authors=_field_text(fields, "A") or None,
        journal=_metadata_text(fields, ("T", "J")) or None,
        correlations=_parse_correlations(fields.get("C", [])),
        values=values,
        value_order=STANDARD_AMINO_ACIDS,
    )


def _collect_fields(entry: str) -> dict[str, list[str]]:
    fields: dict[str, list[str]] = defaultdict(list)
    current_field: str | None = None
    for raw_line in entry.splitlines():
        if not raw_line.strip():
            continue
        if len(raw_line) >= 2 and raw_line[0] in AAINDEX1_FIELDS and raw_line[1].isspace():
            current_field = raw_line[0]
            fields[current_field].append(raw_line[1:].strip())
        elif current_field is not None:
            fields[current_field].append(raw_line.strip())
    return fields


def _field_text(fields: dict[str, list[str]], field: str) -> str:
    return " ".join(part.strip() for part in fields.get(field, []) if part.strip())


def _metadata_text(fields: dict[str, list[str]], field_names: tuple[str, ...]) -> str:
    parts = [_field_text(fields, field) for field in field_names]
    return " ".join(part for part in parts if part)


def _parse_i_field(lines: list[str], *, accession: str) -> dict[str, float | None]:
    if not lines:
        return {}

    first_order: list[str] = []
    second_order: list[str] = []
    explicit_order: list[str] = []
    value_tokens: list[str] = []

    for line in lines:
        tokens = line.split()
        slash_tokens = [token for token in tokens if _aa_pair_token(token)]
        single_aa_tokens = [token.upper() for token in tokens if token.upper() in STANDARD_AMINO_ACID_SET]
        if slash_tokens:
            for token in slash_tokens:
                left, right = token.upper().split("/", 1)
                first_order.append(left)
                second_order.append(right)
            continue
        if len(single_aa_tokens) >= 10 and not any(_is_value_token(token) for token in tokens):
            explicit_order.extend(single_aa_tokens)
            continue
        value_tokens.extend(token for token in tokens if _is_value_token(token))

    order = explicit_order or [*first_order, *second_order]
    if len(order) != 20:
        raise AaindexParseError(
            f"AAindex1 entry {accession} has malformed I amino acid header; expected 20 residues."
        )
    if set(order) != STANDARD_AMINO_ACID_SET:
        raise AaindexParseError(
            f"AAindex1 entry {accession} has malformed I amino acid header; residues are not standard."
        )
    if len(value_tokens) < 20:
        raise AaindexParseError(
            f"AAindex1 entry {accession} has malformed I values; expected 20 values."
        )

    values: dict[str, float | None] = {}
    for residue, raw_value in zip(order, value_tokens[:20], strict=False):
        values[residue] = _parse_value(raw_value)
    return {residue: values[residue] for residue in STANDARD_AMINO_ACIDS}


def _aa_pair_token(token: str) -> bool:
    pieces = token.upper().split("/")
    return (
        len(pieces) == 2
        and len(pieces[0]) == 1
        and len(pieces[1]) == 1
        and pieces[0] in STANDARD_AMINO_ACID_SET
        and pieces[1] in STANDARD_AMINO_ACID_SET
    )


def _is_value_token(token: str) -> bool:
    normalized = token.upper()
    return normalized in MISSING_VALUE_TOKENS or bool(_NUMBER_RE.match(token))


def _parse_value(token: str) -> float | None:
    if token.upper() in MISSING_VALUE_TOKENS:
        return None
    return float(token)


def _parse_correlations(lines: list[str]) -> dict[str, float | None]:
    correlations: dict[str, float | None] = {}
    tokens = " ".join(lines).split()
    index = 0
    while index < len(tokens):
        accession = tokens[index]
        value: float | None = None
        if index + 1 < len(tokens) and _is_value_token(tokens[index + 1]):
            value = _parse_value(tokens[index + 1])
            index += 2
        else:
            index += 1
        correlations[accession] = value
    return correlations
