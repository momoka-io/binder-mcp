"""Protein sequence parsing and validation helpers."""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from textwrap import wrap

STANDARD_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
STANDARD_AMINO_ACID_SET = set(STANDARD_AMINO_ACIDS)
MAX_PROTEIN_SEQUENCE_LENGTH = 10_000
MAX_ALIGNMENT_SEQUENCE_LENGTH = 5_000
MAX_ALIGNMENT_SEQUENCE_COUNT = 100


@dataclass(frozen=True)
class InvalidCharacterReport:
    """A single invalid character and all 1-based parsed sequence positions."""

    character: str
    count: int
    positions: tuple[int, ...]


@dataclass(frozen=True)
class ProteinSequenceValidation:
    """Parsed protein sequence validation details."""

    cleaned_sequence: str
    candidate_sequence: str
    invalid_characters: tuple[InvalidCharacterReport, ...]
    warnings: tuple[str, ...]

    @property
    def length(self) -> int:
        return len(self.cleaned_sequence)

    @property
    def is_valid(self) -> bool:
        return self.length > 0 and not self.invalid_characters


@dataclass(frozen=True)
class FastaRecord:
    """A validated protein FASTA record."""

    header: str
    sequence: str


class ProteinSequenceError(ValueError):
    """Raised when a protein sequence cannot be used for computation."""


def parse_protein_sequence(
    text: str,
    *,
    max_length: int = MAX_PROTEIN_SEQUENCE_LENGTH,
) -> ProteinSequenceValidation:
    """Parse raw sequence or FASTA text into a cleaned standard protein sequence.

    FASTA headers are removed, sequence lines are concatenated, residues are
    uppercased, and non-standard characters are reported rather than silently
    accepted. The returned cleaned sequence contains only the 20 standard amino
    acid one-letter codes.
    """

    if max_length < 1:
        raise ProteinSequenceError("max_length must be at least 1.")

    warnings: list[str] = []
    sequence_chunks: list[str] = []
    fasta_headers = 0
    whitespace_removed = False
    lower_case_found = False

    for raw_line in text.splitlines() or [text]:
        stripped = raw_line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            fasta_headers += 1
            continue
        if any(char.isspace() for char in stripped):
            whitespace_removed = True
        if any(char.islower() for char in stripped):
            lower_case_found = True
        sequence_chunks.append("".join(stripped.split()))

    if fasta_headers == 1:
        warnings.append("Input contained a FASTA header; sequence lines were parsed.")
    elif fasta_headers > 1:
        warnings.append("Input contained multiple FASTA records; sequence lines were concatenated.")

    if whitespace_removed:
        warnings.append("Whitespace inside sequence lines was removed.")
    if lower_case_found:
        warnings.append("Sequence was converted to uppercase.")

    candidate_sequence = "".join(sequence_chunks).upper()
    if len(candidate_sequence) > max_length:
        raise ProteinSequenceError(
            f"Protein sequence length {len(candidate_sequence)} exceeds maximum {max_length}."
        )

    positions_by_character: OrderedDict[str, list[int]] = OrderedDict()
    cleaned: list[str] = []
    for position, residue in enumerate(candidate_sequence, start=1):
        if residue in STANDARD_AMINO_ACID_SET:
            cleaned.append(residue)
        else:
            positions_by_character.setdefault(residue, []).append(position)

    invalid_reports = tuple(
        InvalidCharacterReport(character=character, count=len(positions), positions=tuple(positions))
        for character, positions in positions_by_character.items()
    )
    if invalid_reports:
        warnings.append("Invalid characters were omitted from cleaned_sequence.")

    cleaned_sequence = "".join(cleaned)
    if not cleaned_sequence:
        warnings.append("No standard amino acid residues were found.")

    return ProteinSequenceValidation(
        cleaned_sequence=cleaned_sequence,
        candidate_sequence=candidate_sequence,
        invalid_characters=invalid_reports,
        warnings=tuple(warnings),
    )


def require_valid_protein_sequence(
    text: str,
    *,
    max_length: int = MAX_PROTEIN_SEQUENCE_LENGTH,
) -> ProteinSequenceValidation:
    """Parse a protein sequence and raise a clean error if it is not computable."""

    validation = parse_protein_sequence(text, max_length=max_length)
    if validation.length == 0:
        raise ProteinSequenceError("Protein sequence is empty after parsing.")
    if validation.invalid_characters:
        invalid_summary = ", ".join(
            f"{report.character!r} at position(s) {', '.join(str(pos) for pos in report.positions)}"
            for report in validation.invalid_characters
        )
        raise ProteinSequenceError(f"Invalid protein sequence characters: {invalid_summary}.")
    return validation


def parse_protein_fasta_records(
    fasta_text: str,
    *,
    max_sequences: int = MAX_ALIGNMENT_SEQUENCE_COUNT,
    max_sequence_length: int = MAX_ALIGNMENT_SEQUENCE_LENGTH,
) -> tuple[tuple[FastaRecord, ...], tuple[str, ...]]:
    """Parse and validate protein FASTA records for local alignment wrappers."""

    if max_sequences < 1:
        raise ProteinSequenceError("max_sequences must be at least 1.")
    if max_sequence_length < 1:
        raise ProteinSequenceError("max_sequence_length must be at least 1.")

    records: list[tuple[str, list[str]]] = []
    current_header: str | None = None
    current_lines: list[str] = []
    warnings: list[str] = []
    saw_content = False

    for raw_line in fasta_text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        saw_content = True
        if line.startswith(">"):
            if current_header is not None:
                records.append((current_header, current_lines))
            current_header = line[1:].strip() or f"sequence_{len(records) + 1}"
            current_lines = []
            if len(records) + 1 > max_sequences:
                raise ProteinSequenceError(
                    f"FASTA sequence count exceeds maximum {max_sequences}."
                )
            continue
        if current_header is None:
            raise ProteinSequenceError("Alignment input must be FASTA text with header lines.")
        current_lines.append(line)

    if not saw_content:
        raise ProteinSequenceError("FASTA input is empty.")
    if current_header is None:
        raise ProteinSequenceError("Alignment input must be FASTA text with header lines.")

    records.append((current_header, current_lines))
    if len(records) > max_sequences:
        raise ProteinSequenceError(f"FASTA sequence count exceeds maximum {max_sequences}.")

    parsed_records: list[FastaRecord] = []
    for index, (header, lines) in enumerate(records, start=1):
        sequence_text = "".join(lines)
        validation = parse_protein_sequence(sequence_text, max_length=max_sequence_length)
        warnings.extend(f"Record {index} ({header}): {warning}" for warning in validation.warnings)
        if validation.length == 0:
            raise ProteinSequenceError(f"FASTA record {index} ({header}) is empty.")
        if validation.invalid_characters:
            invalid_summary = ", ".join(
                f"{report.character!r} at position(s) {', '.join(str(pos) for pos in report.positions)}"
                for report in validation.invalid_characters
            )
            raise ProteinSequenceError(
                f"Invalid protein sequence characters in FASTA record {index} ({header}): "
                f"{invalid_summary}."
            )
        parsed_records.append(FastaRecord(header=header, sequence=validation.cleaned_sequence))

    if len(parsed_records) == 1:
        warnings.append("Only one FASTA record was provided; alignment output may be unchanged.")

    return tuple(parsed_records), tuple(warnings)


def render_fasta(records: tuple[FastaRecord, ...], *, line_width: int = 80) -> str:
    """Render validated FASTA records with stable wrapping."""

    lines: list[str] = []
    for record in records:
        lines.append(f">{record.header}")
        lines.extend(wrap(record.sequence, width=line_width) or [""])
    return "\n".join(lines) + "\n"
