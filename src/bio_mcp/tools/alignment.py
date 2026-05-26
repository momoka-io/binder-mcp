"""Local MAFFT and Clustal Omega alignment wrappers."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

from bio_mcp import __version__
from bio_mcp.core.sequence import parse_protein_fasta_records, render_fasta
from bio_mcp.core.subprocesses import (
    CommandStatus,
    CommandTimeoutError,
    detect_command,
    run_command,
    summarize_stderr,
)
from bio_mcp.schemas import AlignmentOutput, ClustaloAlignInput, MafftAlignInput, ToolProvenance

MAFFT_TOOL_NAME = "mafft_align"
CLUSTALO_TOOL_NAME = "clustalo_align"

MAFFT_MODE_ARGS = {
    "auto": ["--auto"],
    "localpair": ["--localpair", "--maxiterate", "1000"],
    "globalpair": ["--globalpair", "--maxiterate", "1000"],
    "genafpair": ["--genafpair", "--maxiterate", "1000"],
}


def detect_mafft() -> CommandStatus:
    """Detect the local MAFFT binary and version."""

    return detect_command("mafft", ["--version"])


def detect_clustalo() -> CommandStatus:
    """Detect the local Clustal Omega binary and version."""

    return detect_command("clustalo", ["--version"])


def mafft_align_tool(params: MafftAlignInput) -> AlignmentOutput:
    """Align protein FASTA records with a local MAFFT binary."""

    records, validation_warnings = parse_protein_fasta_records(params.fasta_text)
    warnings = list(validation_warnings)
    status = detect_mafft()
    if not status.available:
        return _dependency_error_output(
            tool_name=MAFFT_TOOL_NAME,
            backend="mafft",
            number_of_sequences=len(records),
            status=status,
            warnings=warnings,
        )
    if status.error:
        warnings.append(status.error)

    with TemporaryDirectory(prefix="bio-mcp-mafft-") as tmpdir:
        tmp_path = Path(tmpdir)
        input_path = tmp_path / "input.fasta"
        input_path.write_text(render_fasta(records), encoding="utf-8")
        command = [status.path or "mafft", *MAFFT_MODE_ARGS[params.mode], str(input_path)]
        try:
            result = run_command(command, cwd=tmp_path, timeout_sec=params.timeout_sec)
        except CommandTimeoutError as exc:
            return _runtime_error_output(
                tool_name=MAFFT_TOOL_NAME,
                backend="mafft",
                number_of_sequences=len(records),
                status=status,
                warnings=warnings,
                error=str(exc),
            )
        except OSError as exc:
            return _runtime_error_output(
                tool_name=MAFFT_TOOL_NAME,
                backend="mafft",
                number_of_sequences=len(records),
                status=status,
                warnings=warnings,
                error=f"Could not run MAFFT: {exc}",
            )

    stderr_summary = summarize_stderr(result.stderr)
    if result.returncode != 0:
        error = f"MAFFT exited with status {result.returncode}."
        return _runtime_error_output(
            tool_name=MAFFT_TOOL_NAME,
            backend="mafft",
            number_of_sequences=len(records),
            status=status,
            warnings=warnings,
            error=error,
            stderr_summary=stderr_summary,
        )
    if not result.stdout.strip():
        return _runtime_error_output(
            tool_name=MAFFT_TOOL_NAME,
            backend="mafft",
            number_of_sequences=len(records),
            status=status,
            warnings=warnings,
            error="MAFFT completed but did not produce aligned FASTA output.",
            stderr_summary=stderr_summary,
        )

    return AlignmentOutput(
        aligned_fasta=result.stdout,
        number_of_sequences=len(records),
        command_version=status.version,
        stderr_summary=stderr_summary,
        warnings=warnings,
        error=None,
        provenance=ToolProvenance(
            tool_name=MAFFT_TOOL_NAME,
            wrapper_version=__version__,
            backend="mafft",
            backend_version=status.version,
            execution_mode="local",
            warnings=warnings,
        ),
    )


def clustalo_align_tool(params: ClustaloAlignInput) -> AlignmentOutput:
    """Align protein FASTA records with a local Clustal Omega binary."""

    records, validation_warnings = parse_protein_fasta_records(params.fasta_text)
    warnings = list(validation_warnings)
    status = detect_clustalo()
    if not status.available:
        return _dependency_error_output(
            tool_name=CLUSTALO_TOOL_NAME,
            backend="clustalo",
            number_of_sequences=len(records),
            status=status,
            warnings=warnings,
        )
    if status.error:
        warnings.append(status.error)

    with TemporaryDirectory(prefix="bio-mcp-clustalo-") as tmpdir:
        tmp_path = Path(tmpdir)
        input_path = tmp_path / "input.fasta"
        output_path = tmp_path / "aligned.fasta"
        input_path.write_text(render_fasta(records), encoding="utf-8")
        command = [
            status.path or "clustalo",
            "--infile",
            str(input_path),
            "--outfile",
            str(output_path),
            "--outfmt",
            "fasta",
            "--force",
        ]
        try:
            result = run_command(command, cwd=tmp_path, timeout_sec=params.timeout_sec)
        except CommandTimeoutError as exc:
            return _runtime_error_output(
                tool_name=CLUSTALO_TOOL_NAME,
                backend="clustalo",
                number_of_sequences=len(records),
                status=status,
                warnings=warnings,
                error=str(exc),
            )
        except OSError as exc:
            return _runtime_error_output(
                tool_name=CLUSTALO_TOOL_NAME,
                backend="clustalo",
                number_of_sequences=len(records),
                status=status,
                warnings=warnings,
                error=f"Could not run Clustal Omega: {exc}",
            )
        stderr_summary = summarize_stderr(result.stderr)
        if result.returncode != 0:
            return _runtime_error_output(
                tool_name=CLUSTALO_TOOL_NAME,
                backend="clustalo",
                number_of_sequences=len(records),
                status=status,
                warnings=warnings,
                error=f"Clustal Omega exited with status {result.returncode}.",
                stderr_summary=stderr_summary,
            )
        if not output_path.exists() or not output_path.read_text(encoding="utf-8").strip():
            return _runtime_error_output(
                tool_name=CLUSTALO_TOOL_NAME,
                backend="clustalo",
                number_of_sequences=len(records),
                status=status,
                warnings=warnings,
                error="Clustal Omega completed but did not produce aligned FASTA output.",
                stderr_summary=stderr_summary,
            )
        aligned_fasta = output_path.read_text(encoding="utf-8")

    return AlignmentOutput(
        aligned_fasta=aligned_fasta,
        number_of_sequences=len(records),
        command_version=status.version,
        stderr_summary=stderr_summary,
        warnings=warnings,
        error=None,
        provenance=ToolProvenance(
            tool_name=CLUSTALO_TOOL_NAME,
            wrapper_version=__version__,
            backend="clustalo",
            backend_version=status.version,
            execution_mode="local",
            warnings=warnings,
        ),
    )


def _dependency_error_output(
    *,
    tool_name: str,
    backend: str,
    number_of_sequences: int,
    status: CommandStatus,
    warnings: list[str],
) -> AlignmentOutput:
    error = status.error or f"Dependency error: required binary '{backend}' was not found on PATH."
    output_warnings = [*warnings, error]
    return AlignmentOutput(
        aligned_fasta=None,
        number_of_sequences=number_of_sequences,
        command_version=status.version,
        stderr_summary="",
        warnings=output_warnings,
        error=error,
        provenance=ToolProvenance(
            tool_name=tool_name,
            wrapper_version=__version__,
            backend=backend,
            backend_version=status.version,
            execution_mode="local",
            warnings=output_warnings,
        ),
    )


def _runtime_error_output(
    *,
    tool_name: str,
    backend: str,
    number_of_sequences: int,
    status: CommandStatus,
    warnings: list[str],
    error: str,
    stderr_summary: str = "",
) -> AlignmentOutput:
    output_warnings = [*warnings, error]
    return AlignmentOutput(
        aligned_fasta=None,
        number_of_sequences=number_of_sequences,
        command_version=status.version,
        stderr_summary=stderr_summary,
        warnings=output_warnings,
        error=error,
        provenance=ToolProvenance(
            tool_name=tool_name,
            wrapper_version=__version__,
            backend=backend,
            backend_version=status.version,
            execution_mode="local",
            warnings=output_warnings,
        ),
    )
