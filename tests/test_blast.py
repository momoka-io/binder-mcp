import subprocess
from pathlib import Path

from bio_mcp.schemas import BlastpLocalInput, PsiblastLocalInput
from bio_mcp.tools.blast import (
    blastp_local_tool,
    detect_blastp,
    parse_blast_outfmt6,
    psiblast_local_tool,
)

QUERY_FASTA = ">query1\nMKTAYIAKQRQISFVKSHFSRQDILD\n"
RAW_HIT = "query1\tsubject1\t100.000\t26\t0\t0\t1\t26\t1\t26\t1e-20\t55.5\n"


def _make_db_markers(tmp_path: Path) -> str:
    db_prefix = tmp_path / "toydb"
    for suffix in (".pin", ".phr", ".psq"):
        db_prefix.with_name(db_prefix.name + suffix).write_text("marker", encoding="utf-8")
    return str(db_prefix)


def test_parse_blast_outfmt6_into_structured_hits():
    hits = parse_blast_outfmt6(RAW_HIT)

    assert len(hits) == 1
    assert hits[0].qseqid == "query1"
    assert hits[0].sseqid == "subject1"
    assert hits[0].pident == 100.0
    assert hits[0].length == 26
    assert hits[0].evalue == 1e-20
    assert hits[0].bitscore == 55.5


def test_blastp_local_success_mocked_subprocess(monkeypatch, tmp_path):
    calls = []
    db_path = _make_db_markers(tmp_path)

    monkeypatch.delenv("BIO_MCP_BLASTP_BIN", raising=False)
    monkeypatch.setattr("bio_mcp.core.subprocesses.shutil.which", lambda command: "blastp")

    def fake_run(args, **kwargs):
        calls.append(args)
        assert kwargs["shell"] is False
        assert kwargs["capture_output"] is True
        assert Path(kwargs["cwd"]).exists()
        if "-version" in args:
            return subprocess.CompletedProcess(args, 0, stdout="blastp: 2.17.0+\n", stderr="")
        assert "-query" in args
        assert args[args.index("-db") + 1] == db_path
        assert args[args.index("-outfmt") + 1].startswith("6 qseqid sseqid pident")
        return subprocess.CompletedProcess(args, 0, stdout=RAW_HIT, stderr="blastp done\n")

    monkeypatch.setattr("bio_mcp.core.subprocesses.subprocess.run", fake_run)

    output = blastp_local_tool(BlastpLocalInput(query_fasta=QUERY_FASTA, db_path=db_path))

    assert len(calls) == 2
    assert output.error is None
    assert output.raw_tabular == RAW_HIT
    assert output.number_of_hits == 1
    assert output.query_count == 1
    assert output.hits[0].sseqid == "subject1"
    assert output.command_version == "blastp: 2.17.0+"
    assert output.stderr_summary == "blastp done"
    assert output.provenance.backend == "blastp"


def test_psiblast_local_success_mocked_subprocess(monkeypatch, tmp_path):
    calls = []
    db_path = _make_db_markers(tmp_path)

    monkeypatch.delenv("BIO_MCP_PSIBLAST_BIN", raising=False)
    monkeypatch.setattr("bio_mcp.core.subprocesses.shutil.which", lambda command: "psiblast")

    def fake_run(args, **kwargs):
        calls.append(args)
        assert kwargs["shell"] is False
        assert kwargs["capture_output"] is True
        if "-version" in args:
            return subprocess.CompletedProcess(args, 0, stdout="psiblast: 2.17.0+\n", stderr="")
        assert args[args.index("-num_iterations") + 1] == "2"
        assert args[args.index("-db") + 1] == db_path
        return subprocess.CompletedProcess(args, 0, stdout=RAW_HIT, stderr="psiblast done\n")

    monkeypatch.setattr("bio_mcp.core.subprocesses.subprocess.run", fake_run)

    output = psiblast_local_tool(
        PsiblastLocalInput(query_fasta=QUERY_FASTA, db_path=db_path, num_iterations=2)
    )

    assert len(calls) == 2
    assert output.error is None
    assert output.number_of_hits == 1
    assert output.query_count == 1
    assert output.num_iterations == 2
    assert output.hits[0].qseqid == "query1"
    assert output.command_version == "psiblast: 2.17.0+"
    assert output.stderr_summary == "psiblast done"
    assert output.provenance.backend == "psiblast"


def test_blastp_local_missing_binary_returns_dependency_error(monkeypatch, tmp_path):
    db_path = _make_db_markers(tmp_path)
    monkeypatch.delenv("BIO_MCP_BLASTP_BIN", raising=False)
    monkeypatch.setattr("bio_mcp.core.subprocesses.shutil.which", lambda command: None)

    output = blastp_local_tool(BlastpLocalInput(query_fasta=QUERY_FASTA, db_path=db_path))

    assert output.hits == []
    assert output.query_count == 1
    assert output.command_version is None
    assert output.error == "Dependency error: required binary 'blastp' was not found on PATH."
    assert output.error in output.warnings


def test_psiblast_local_missing_binary_returns_dependency_error(monkeypatch, tmp_path):
    db_path = _make_db_markers(tmp_path)
    monkeypatch.delenv("BIO_MCP_PSIBLAST_BIN", raising=False)
    monkeypatch.setattr("bio_mcp.core.subprocesses.shutil.which", lambda command: None)

    output = psiblast_local_tool(PsiblastLocalInput(query_fasta=QUERY_FASTA, db_path=db_path))

    assert output.hits == []
    assert output.query_count == 1
    assert output.command_version is None
    assert output.error == "Dependency error: required binary 'psiblast' was not found on PATH."
    assert output.error in output.warnings


def test_blastp_env_override_valid_executable(monkeypatch, tmp_path):
    fake_blastp = tmp_path / "blastp"
    fake_blastp.write_text("#!/bin/sh\necho 'blastp: test'\n", encoding="utf-8")
    fake_blastp.chmod(0o755)
    monkeypatch.setenv("BIO_MCP_BLASTP_BIN", str(fake_blastp))
    monkeypatch.setattr("bio_mcp.core.subprocesses.shutil.which", lambda command: None)

    status = detect_blastp()

    assert status.available is True
    assert status.path == str(fake_blastp)
    assert status.version == "blastp: test"
    assert status.resolution_source == "env_override"


def test_psiblast_env_override_missing_path(monkeypatch, tmp_path):
    missing = tmp_path / "missing-psiblast"
    db_path = _make_db_markers(tmp_path)
    monkeypatch.setenv("BIO_MCP_PSIBLAST_BIN", str(missing))
    monkeypatch.setattr("bio_mcp.core.subprocesses.shutil.which", lambda command: None)

    output = psiblast_local_tool(PsiblastLocalInput(query_fasta=QUERY_FASTA, db_path=db_path))

    assert output.error is not None
    assert "BIO_MCP_PSIBLAST_BIN points to" in output.error
    assert str(missing) in output.error


def test_blastp_local_invalid_fasta_returns_clean_error(tmp_path):
    output = blastp_local_tool(
        BlastpLocalInput(query_fasta=">query1\nMKT*\n", db_path=str(tmp_path / "toydb"))
    )

    assert output.hits == []
    assert output.query_count == 0
    assert output.error is not None
    assert "Invalid protein sequence characters" in output.error


def test_psiblast_local_invalid_fasta_returns_clean_error(tmp_path):
    output = psiblast_local_tool(
        PsiblastLocalInput(query_fasta=">query1\nMKT*\n", db_path=str(tmp_path / "toydb"))
    )

    assert output.hits == []
    assert output.query_count == 0
    assert output.error is not None
    assert "Invalid protein sequence characters" in output.error


def test_blastp_local_missing_db_path_returns_database_error():
    output = blastp_local_tool(BlastpLocalInput(query_fasta=QUERY_FASTA, db_path=" "))

    assert output.hits == []
    assert output.query_count == 1
    assert output.error == "Database error: db_path must be a local BLAST protein database prefix."
    assert output.error in output.warnings


def test_psiblast_local_db_path_not_found_returns_database_error(tmp_path):
    missing = tmp_path / "missingdb"

    output = psiblast_local_tool(PsiblastLocalInput(query_fasta=QUERY_FASTA, db_path=str(missing)))

    assert output.hits == []
    assert output.query_count == 1
    assert output.error is not None
    assert "does not look like a local protein BLAST database prefix" in output.error


def test_blastp_local_timeout_returns_error(monkeypatch, tmp_path):
    db_path = _make_db_markers(tmp_path)
    monkeypatch.delenv("BIO_MCP_BLASTP_BIN", raising=False)
    monkeypatch.setattr("bio_mcp.core.subprocesses.shutil.which", lambda command: "blastp")

    def fake_run(args, **kwargs):
        if "-version" in args:
            return subprocess.CompletedProcess(args, 0, stdout="blastp: 2.17.0+\n", stderr="")
        raise subprocess.TimeoutExpired(args, kwargs["timeout"])

    monkeypatch.setattr("bio_mcp.core.subprocesses.subprocess.run", fake_run)

    output = blastp_local_tool(
        BlastpLocalInput(query_fasta=QUERY_FASTA, db_path=db_path, timeout_sec=1)
    )

    assert output.error == "Command timed out after 1 seconds: blastp"
    assert output.error in output.warnings
