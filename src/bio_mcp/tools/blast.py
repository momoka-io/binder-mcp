"""Local BLASTP and PSI-BLAST wrappers."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

from bio_mcp import __version__
from bio_mcp.core.sequence import ProteinSequenceError, parse_protein_fasta_records, render_fasta
from bio_mcp.core.subprocesses import (
    CommandStatus,
    CommandTimeoutError,
    detect_command,
    run_command,
    summarize_stderr,
)
from bio_mcp.schemas import (
    BlastpLocalInput,
    BlastpLocalOutput,
    BlastTabularHit,
    PsiblastLocalInput,
    PsiblastLocalOutput,
    ToolProvenance,
)

BLASTP_TOOL_NAME = "blastp_local"
PSIBLAST_TOOL_NAME = "psiblast_local"
BLASTP_ENV_VAR = "BIO_MCP_BLASTP_BIN"
PSIBLAST_ENV_VAR = "BIO_MCP_PSIBLAST_BIN"
MAKEBLASTDB_ENV_VAR = "BIO_MCP_MAKEBLASTDB_BIN"
BLAST_OUTFMT_FIELDS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
)
BLAST_OUTFMT = "6 " + " ".join(BLAST_OUTFMT_FIELDS)
PROTEIN_DB_EXACT_SUFFIXES = (".pin", ".phr", ".psq")
PROTEIN_DB_VOLUME_SUFFIXES = (".pin", ".phr", ".psq")
PROTEIN_DB_ALIAS_SUFFIX = ".pal"


def detect_blastp() -> CommandStatus:
    """Detect the local BLASTP binary and version."""

    return detect_command("blastp", ["-version"], env_var=BLASTP_ENV_VAR)


def detect_psiblast() -> CommandStatus:
    """Detect the local PSI-BLAST binary and version."""

    return detect_command("psiblast", ["-version"], env_var=PSIBLAST_ENV_VAR)


def detect_makeblastdb() -> CommandStatus:
    """Detect the local makeblastdb binary and version."""

    return detect_command("makeblastdb", ["-version"], env_var=MAKEBLASTDB_ENV_VAR)


def validate_protein_blast_db_prefix(db_path: str) -> str | None:
    """Return a clean error if db_path is not a local protein BLAST database prefix."""

    if not db_path or not db_path.strip():
        return "Database error: db_path must be a local BLAST protein database prefix."

    prefix = Path(db_path).expanduser()
    if prefix.is_dir():
        return (
            f"Database error: db_path '{db_path}' is a directory; provide the BLAST database "
            "prefix, not a directory."
        )

    exact_markers = [prefix.with_name(prefix.name + suffix) for suffix in PROTEIN_DB_EXACT_SUFFIXES]
    if all(marker.is_file() for marker in exact_markers):
        return None

    if prefix.with_name(prefix.name + PROTEIN_DB_ALIAS_SUFFIX).is_file():
        return None

    parent = prefix.parent if str(prefix.parent) else Path(".")
    if parent.exists():
        volume_markers = {
            suffix: list(parent.glob(f"{prefix.name}.*{suffix}"))
            for suffix in PROTEIN_DB_VOLUME_SUFFIXES
        }
        if all(markers for markers in volume_markers.values()):
            return None

    expected = ", ".join(str(marker) for marker in exact_markers)
    return (
        f"Database error: db_path '{db_path}' does not look like a local protein BLAST "
        f"database prefix. Expected marker files such as: {expected}."
    )


def parse_blast_outfmt6(raw_tabular: str) -> list[BlastTabularHit]:
    """Parse BLAST outfmt 6 fields into structured hits."""

    hits: list[BlastTabularHit] = []
    for line_number, line in enumerate(raw_tabular.splitlines(), start=1):
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != len(BLAST_OUTFMT_FIELDS):
            raise ValueError(
                f"BLAST outfmt 6 line {line_number} has {len(fields)} fields; "
                f"expected {len(BLAST_OUTFMT_FIELDS)}."
            )
        try:
            hits.append(
                BlastTabularHit(
                    qseqid=fields[0],
                    sseqid=fields[1],
                    pident=float(fields[2]),
                    length=int(fields[3]),
                    mismatch=int(fields[4]),
                    gapopen=int(fields[5]),
                    qstart=int(fields[6]),
                    qend=int(fields[7]),
                    sstart=int(fields[8]),
                    send=int(fields[9]),
                    evalue=float(fields[10]),
                    bitscore=float(fields[11]),
                )
            )
        except ValueError as exc:
            raise ValueError(f"Could not parse BLAST outfmt 6 line {line_number}: {exc}") from exc
    return hits


def clean_blast_outfmt6(raw_output: str, warnings: list[str]) -> str:
    """Keep BLAST outfmt 6 hit rows and report non-tabular stdout lines."""

    tabular_lines: list[str] = []
    ignored_lines: list[str] = []
    for line in raw_output.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) == len(BLAST_OUTFMT_FIELDS):
            tabular_lines.append(line)
        elif len(fields) == 1:
            ignored_lines.append(line.strip())
        else:
            tabular_lines.append(line)

    if ignored_lines:
        warnings.append(
            "Ignored non-tabular BLAST stdout line(s): "
            + "; ".join(ignored_lines[:3])
            + ("..." if len(ignored_lines) > 3 else "")
        )
    return "\n".join(tabular_lines) + ("\n" if tabular_lines else "")


def blastp_local_tool(params: BlastpLocalInput) -> BlastpLocalOutput:
    """Search a local protein BLAST database with a local BLASTP binary."""

    try:
        records, validation_warnings = parse_protein_fasta_records(params.query_fasta)
    except ProteinSequenceError as exc:
        return _blastp_error_output(
            tool_name=BLASTP_TOOL_NAME,
            backend="blastp",
            query_count=0,
            db_path=params.db_path,
            command_version=None,
            warnings=[],
            error=str(exc),
        )

    warnings = list(validation_warnings)
    db_error = validate_protein_blast_db_prefix(params.db_path)
    if db_error:
        return _blastp_error_output(
            tool_name=BLASTP_TOOL_NAME,
            backend="blastp",
            query_count=len(records),
            db_path=params.db_path,
            command_version=None,
            warnings=warnings,
            error=db_error,
        )

    status = detect_blastp()
    if not status.available:
        return _blastp_dependency_error_output(
            tool_name=BLASTP_TOOL_NAME,
            backend="blastp",
            query_count=len(records),
            db_path=params.db_path,
            status=status,
            warnings=warnings,
        )
    if status.error:
        warnings.append(status.error)

    with TemporaryDirectory(prefix="bio-mcp-blastp-") as tmpdir:
        tmp_path = Path(tmpdir)
        query_path = tmp_path / "query.fasta"
        query_path.write_text(render_fasta(records), encoding="utf-8")
        command = [
            status.path or "blastp",
            "-query",
            str(query_path),
            "-db",
            params.db_path,
            "-evalue",
            str(params.evalue),
            "-max_target_seqs",
            str(params.max_target_seqs),
            "-outfmt",
            BLAST_OUTFMT,
        ]
        try:
            result = run_command(command, cwd=tmp_path, timeout_sec=params.timeout_sec)
        except CommandTimeoutError as exc:
            return _blastp_error_output(
                tool_name=BLASTP_TOOL_NAME,
                backend="blastp",
                query_count=len(records),
                db_path=params.db_path,
                command_version=status.version,
                warnings=warnings,
                error=str(exc),
            )
        except OSError as exc:
            return _blastp_error_output(
                tool_name=BLASTP_TOOL_NAME,
                backend="blastp",
                query_count=len(records),
                db_path=params.db_path,
                command_version=status.version,
                warnings=warnings,
                error=f"Could not run BLASTP: {exc}",
            )

    stderr_summary = summarize_stderr(result.stderr)
    if result.returncode != 0:
        return _blastp_error_output(
            tool_name=BLASTP_TOOL_NAME,
            backend="blastp",
            query_count=len(records),
            db_path=params.db_path,
            command_version=status.version,
            warnings=warnings,
            error=f"BLASTP exited with status {result.returncode}.",
            stderr_summary=stderr_summary,
            raw_tabular=result.stdout,
        )

    raw_tabular = clean_blast_outfmt6(result.stdout, warnings)
    try:
        hits = parse_blast_outfmt6(raw_tabular)
    except ValueError as exc:
        return _blastp_error_output(
            tool_name=BLASTP_TOOL_NAME,
            backend="blastp",
            query_count=len(records),
            db_path=params.db_path,
            command_version=status.version,
            warnings=warnings,
            error=str(exc),
            stderr_summary=stderr_summary,
            raw_tabular=raw_tabular,
        )

    return BlastpLocalOutput(
        hits=hits,
        raw_tabular=raw_tabular,
        number_of_hits=len(hits),
        query_count=len(records),
        db_path=params.db_path,
        command_version=status.version,
        stderr_summary=stderr_summary,
        warnings=warnings,
        error=None,
        provenance=_provenance(BLASTP_TOOL_NAME, "blastp", status.version, warnings),
    )


def psiblast_local_tool(params: PsiblastLocalInput) -> PsiblastLocalOutput:
    """Search a local protein BLAST database with a local PSI-BLAST binary."""

    try:
        records, validation_warnings = parse_protein_fasta_records(params.query_fasta)
    except ProteinSequenceError as exc:
        return _psiblast_error_output(
            query_count=0,
            db_path=params.db_path,
            num_iterations=params.num_iterations,
            command_version=None,
            warnings=[],
            error=str(exc),
        )

    warnings = list(validation_warnings)
    db_error = validate_protein_blast_db_prefix(params.db_path)
    if db_error:
        return _psiblast_error_output(
            query_count=len(records),
            db_path=params.db_path,
            num_iterations=params.num_iterations,
            command_version=None,
            warnings=warnings,
            error=db_error,
        )

    status = detect_psiblast()
    if not status.available:
        error = status.error or "Dependency error: required binary 'psiblast' was not found on PATH."
        output_warnings = [*warnings, error]
        return PsiblastLocalOutput(
            hits=[],
            raw_tabular="",
            number_of_hits=0,
            query_count=len(records),
            db_path=params.db_path,
            command_version=status.version,
            stderr_summary="",
            warnings=output_warnings,
            error=error,
            provenance=_provenance(PSIBLAST_TOOL_NAME, "psiblast", status.version, output_warnings),
            num_iterations=params.num_iterations,
        )
    if status.error:
        warnings.append(status.error)

    with TemporaryDirectory(prefix="bio-mcp-psiblast-") as tmpdir:
        tmp_path = Path(tmpdir)
        query_path = tmp_path / "query.fasta"
        query_path.write_text(render_fasta(records), encoding="utf-8")
        command = [
            status.path or "psiblast",
            "-query",
            str(query_path),
            "-db",
            params.db_path,
            "-num_iterations",
            str(params.num_iterations),
            "-evalue",
            str(params.evalue),
            "-max_target_seqs",
            str(params.max_target_seqs),
            "-outfmt",
            BLAST_OUTFMT,
        ]
        try:
            result = run_command(command, cwd=tmp_path, timeout_sec=params.timeout_sec)
        except CommandTimeoutError as exc:
            return _psiblast_error_output(
                query_count=len(records),
                db_path=params.db_path,
                num_iterations=params.num_iterations,
                command_version=status.version,
                warnings=warnings,
                error=str(exc),
            )
        except OSError as exc:
            return _psiblast_error_output(
                query_count=len(records),
                db_path=params.db_path,
                num_iterations=params.num_iterations,
                command_version=status.version,
                warnings=warnings,
                error=f"Could not run PSI-BLAST: {exc}",
            )

    stderr_summary = summarize_stderr(result.stderr)
    if result.returncode != 0:
        return _psiblast_error_output(
            query_count=len(records),
            db_path=params.db_path,
            num_iterations=params.num_iterations,
            command_version=status.version,
            warnings=warnings,
            error=f"PSI-BLAST exited with status {result.returncode}.",
            stderr_summary=stderr_summary,
            raw_tabular=result.stdout,
        )

    raw_tabular = clean_blast_outfmt6(result.stdout, warnings)
    try:
        hits = parse_blast_outfmt6(raw_tabular)
    except ValueError as exc:
        return _psiblast_error_output(
            query_count=len(records),
            db_path=params.db_path,
            num_iterations=params.num_iterations,
            command_version=status.version,
            warnings=warnings,
            error=str(exc),
            stderr_summary=stderr_summary,
            raw_tabular=raw_tabular,
        )

    return PsiblastLocalOutput(
        hits=hits,
        raw_tabular=raw_tabular,
        number_of_hits=len(hits),
        query_count=len(records),
        db_path=params.db_path,
        command_version=status.version,
        stderr_summary=stderr_summary,
        warnings=warnings,
        error=None,
        provenance=_provenance(PSIBLAST_TOOL_NAME, "psiblast", status.version, warnings),
        num_iterations=params.num_iterations,
    )


def _blastp_dependency_error_output(
    *,
    tool_name: str,
    backend: str,
    query_count: int,
    db_path: str,
    status: CommandStatus,
    warnings: list[str],
) -> BlastpLocalOutput:
    error = status.error or f"Dependency error: required binary '{backend}' was not found on PATH."
    return _blastp_error_output(
        tool_name=tool_name,
        backend=backend,
        query_count=query_count,
        db_path=db_path,
        command_version=status.version,
        warnings=warnings,
        error=error,
    )


def _blastp_error_output(
    *,
    tool_name: str,
    backend: str,
    query_count: int,
    db_path: str,
    command_version: str | None,
    warnings: list[str],
    error: str,
    stderr_summary: str = "",
    raw_tabular: str = "",
) -> BlastpLocalOutput:
    output_warnings = [*warnings, error]
    return BlastpLocalOutput(
        hits=[],
        raw_tabular=raw_tabular,
        number_of_hits=0,
        query_count=query_count,
        db_path=db_path,
        command_version=command_version,
        stderr_summary=stderr_summary,
        warnings=output_warnings,
        error=error,
        provenance=_provenance(tool_name, backend, command_version, output_warnings),
    )


def _psiblast_error_output(
    *,
    query_count: int,
    db_path: str,
    num_iterations: int,
    command_version: str | None,
    warnings: list[str],
    error: str,
    stderr_summary: str = "",
    raw_tabular: str = "",
) -> PsiblastLocalOutput:
    output_warnings = [*warnings, error]
    return PsiblastLocalOutput(
        hits=[],
        raw_tabular=raw_tabular,
        number_of_hits=0,
        query_count=query_count,
        db_path=db_path,
        command_version=command_version,
        stderr_summary=stderr_summary,
        warnings=output_warnings,
        error=error,
        provenance=_provenance(PSIBLAST_TOOL_NAME, "psiblast", command_version, output_warnings),
        num_iterations=num_iterations,
    )


def _provenance(
    tool_name: str,
    backend: str,
    command_version: str | None,
    warnings: list[str],
) -> ToolProvenance:
    return ToolProvenance(
        tool_name=tool_name,
        wrapper_version=__version__,
        backend=backend,
        backend_version=command_version,
        execution_mode="local",
        warnings=warnings,
    )
