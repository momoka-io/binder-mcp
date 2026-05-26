import subprocess
from pathlib import Path

import pytest

from bio_mcp.core.sequence import (
    MAX_ALIGNMENT_SEQUENCE_COUNT,
    MAX_ALIGNMENT_SEQUENCE_LENGTH,
    ProteinSequenceError,
)
from bio_mcp.schemas import ClustaloAlignInput, MafftAlignInput
from bio_mcp.tools.alignment import clustalo_align_tool, mafft_align_tool

FASTA = ">seq1\nACD\n>seq2\nACE\n"


def test_mafft_align_success_mocked_subprocess(monkeypatch):
    calls = []

    monkeypatch.setattr("bio_mcp.core.subprocesses.shutil.which", lambda command: "mafft")

    def fake_run(args, **kwargs):
        calls.append(args)
        assert kwargs["shell"] is False
        assert kwargs["capture_output"] is True
        assert Path(kwargs["cwd"]).exists()
        if "--version" in args:
            return subprocess.CompletedProcess(args, 0, stdout="", stderr="v7.520\n")
        assert "--localpair" in args
        return subprocess.CompletedProcess(
            args,
            0,
            stdout=">seq1\nAC-D\n>seq2\nACE-\n",
            stderr="Progressive alignment complete\n",
        )

    monkeypatch.setattr("bio_mcp.core.subprocesses.subprocess.run", fake_run)

    output = mafft_align_tool(MafftAlignInput(fasta_text=FASTA, mode="localpair"))

    assert len(calls) == 2
    assert output.aligned_fasta == ">seq1\nAC-D\n>seq2\nACE-\n"
    assert output.number_of_sequences == 2
    assert output.command_version == "v7.520"
    assert output.stderr_summary == "Progressive alignment complete"
    assert output.error is None
    assert output.provenance.backend == "mafft"


def test_mafft_align_missing_binary_returns_dependency_error(monkeypatch):
    monkeypatch.setattr("bio_mcp.core.subprocesses.shutil.which", lambda command: None)

    output = mafft_align_tool(MafftAlignInput(fasta_text=FASTA))

    assert output.aligned_fasta is None
    assert output.number_of_sequences == 2
    assert output.command_version is None
    assert output.error == "Dependency error: required binary 'mafft' was not found on PATH."
    assert output.error in output.warnings


def test_mafft_align_timeout_returns_error(monkeypatch):
    monkeypatch.setattr("bio_mcp.core.subprocesses.shutil.which", lambda command: "mafft")

    def fake_run(args, **kwargs):
        if "--version" in args:
            return subprocess.CompletedProcess(args, 0, stdout="", stderr="v7.520\n")
        raise subprocess.TimeoutExpired(args, kwargs["timeout"])

    monkeypatch.setattr("bio_mcp.core.subprocesses.subprocess.run", fake_run)

    output = mafft_align_tool(MafftAlignInput(fasta_text=FASTA, timeout_sec=1))

    assert output.aligned_fasta is None
    assert output.error == "Command timed out after 1 seconds: mafft"
    assert output.error in output.warnings


def test_mafft_align_rejects_invalid_sequence():
    with pytest.raises(ProteinSequenceError, match="Invalid protein sequence characters"):
        mafft_align_tool(MafftAlignInput(fasta_text=">seq1\nACD*\n>seq2\nACE\n"))


def test_mafft_align_rejects_empty_input():
    with pytest.raises(ProteinSequenceError, match="FASTA input is empty"):
        mafft_align_tool(MafftAlignInput(fasta_text=" \n"))


def test_mafft_align_rejects_too_long_sequence():
    fasta = f">seq1\n{'A' * (MAX_ALIGNMENT_SEQUENCE_LENGTH + 1)}\n"

    with pytest.raises(ProteinSequenceError, match="exceeds maximum"):
        mafft_align_tool(MafftAlignInput(fasta_text=fasta))


def test_clustalo_align_success_mocked_subprocess(monkeypatch):
    calls = []

    monkeypatch.setattr("bio_mcp.core.subprocesses.shutil.which", lambda command: "clustalo")

    def fake_run(args, **kwargs):
        calls.append(args)
        assert kwargs["shell"] is False
        assert kwargs["capture_output"] is True
        assert Path(kwargs["cwd"]).exists()
        if "--version" in args:
            return subprocess.CompletedProcess(args, 0, stdout="1.2.4\n", stderr="")
        output_path = Path(args[args.index("--outfile") + 1])
        output_path.write_text(">seq1\nAC-D\n>seq2\nACE-\n", encoding="utf-8")
        return subprocess.CompletedProcess(args, 0, stdout="", stderr="clustalo done\n")

    monkeypatch.setattr("bio_mcp.core.subprocesses.subprocess.run", fake_run)

    output = clustalo_align_tool(ClustaloAlignInput(fasta_text=FASTA))

    assert len(calls) == 2
    assert output.aligned_fasta == ">seq1\nAC-D\n>seq2\nACE-\n"
    assert output.number_of_sequences == 2
    assert output.command_version == "1.2.4"
    assert output.stderr_summary == "clustalo done"
    assert output.error is None
    assert output.provenance.backend == "clustalo"


def test_clustalo_align_missing_binary_returns_dependency_error(monkeypatch):
    monkeypatch.setattr("bio_mcp.core.subprocesses.shutil.which", lambda command: None)

    output = clustalo_align_tool(ClustaloAlignInput(fasta_text=FASTA))

    assert output.aligned_fasta is None
    assert output.number_of_sequences == 2
    assert output.command_version is None
    assert output.error == "Dependency error: required binary 'clustalo' was not found on PATH."
    assert output.error in output.warnings


def test_clustalo_align_timeout_returns_error(monkeypatch):
    monkeypatch.setattr("bio_mcp.core.subprocesses.shutil.which", lambda command: "clustalo")

    def fake_run(args, **kwargs):
        if "--version" in args:
            return subprocess.CompletedProcess(args, 0, stdout="1.2.4\n", stderr="")
        raise subprocess.TimeoutExpired(args, kwargs["timeout"])

    monkeypatch.setattr("bio_mcp.core.subprocesses.subprocess.run", fake_run)

    output = clustalo_align_tool(ClustaloAlignInput(fasta_text=FASTA, timeout_sec=1))

    assert output.aligned_fasta is None
    assert output.error == "Command timed out after 1 seconds: clustalo"
    assert output.error in output.warnings


def test_clustalo_align_rejects_invalid_sequence():
    with pytest.raises(ProteinSequenceError, match="Invalid protein sequence characters"):
        clustalo_align_tool(ClustaloAlignInput(fasta_text=">seq1\nACD*\n>seq2\nACE\n"))


def test_clustalo_align_rejects_empty_record():
    with pytest.raises(ProteinSequenceError, match="is empty"):
        clustalo_align_tool(ClustaloAlignInput(fasta_text=">seq1\n\n>seq2\nACE\n"))


def test_clustalo_align_rejects_too_many_sequences():
    fasta = "".join(f">seq{i}\nA\n" for i in range(MAX_ALIGNMENT_SEQUENCE_COUNT + 1))

    with pytest.raises(ProteinSequenceError, match="sequence count exceeds maximum"):
        clustalo_align_tool(ClustaloAlignInput(fasta_text=fasta))
