import os
import subprocess

import pytest

from bio_mcp.schemas import BlastpLocalInput, PsiblastLocalInput
from bio_mcp.tools.blast import (
    blastp_local_tool,
    detect_blastp,
    detect_makeblastdb,
    detect_psiblast,
    psiblast_local_tool,
)


pytestmark = pytest.mark.skipif(
    os.environ.get("BIO_MCP_RUN_INTEGRATION") != "1",
    reason="Set BIO_MCP_RUN_INTEGRATION=1 to run local CLI integration tests.",
)

QUERY_FASTA = ">query1\nMKTAYIAKQRQISFVKSHFSRQDILD\n"
DB_FASTA = (
    ">subject1\nMKTAYIAKQRQISFVKSHFSRQDILD\n"
    ">subject2\nGAVLIPFWYDERKSTNQCMH\n"
)


def _build_tiny_blast_db(tmp_path):
    makeblastdb = detect_makeblastdb()
    if not makeblastdb.available:
        pytest.skip(makeblastdb.error or "makeblastdb binary is not available.")

    fasta_path = tmp_path / "proteins.fasta"
    db_prefix = tmp_path / "proteins_db"
    fasta_path.write_text(DB_FASTA, encoding="utf-8")
    result = subprocess.run(
        [
            makeblastdb.path or "makeblastdb",
            "-in",
            str(fasta_path),
            "-dbtype",
            "prot",
            "-out",
            str(db_prefix),
        ],
        cwd=str(tmp_path),
        text=True,
        capture_output=True,
        timeout=30,
        shell=False,
    )
    assert result.returncode == 0, result.stderr
    return str(db_prefix)


def test_blastp_local_real_cli_with_tiny_database(tmp_path):
    status = detect_blastp()
    if not status.available:
        pytest.skip(status.error or "blastp binary is not available.")
    db_path = _build_tiny_blast_db(tmp_path)

    output = blastp_local_tool(
        BlastpLocalInput(query_fasta=QUERY_FASTA, db_path=db_path, timeout_sec=30)
    )

    assert output.error is None
    assert output.number_of_hits >= 1
    assert output.query_count == 1
    assert output.hits[0].qseqid == "query1"
    assert output.command_version


def test_psiblast_local_real_cli_with_tiny_database(tmp_path):
    status = detect_psiblast()
    if not status.available:
        pytest.skip(status.error or "psiblast binary is not available.")
    db_path = _build_tiny_blast_db(tmp_path)

    output = psiblast_local_tool(
        PsiblastLocalInput(
            query_fasta=QUERY_FASTA,
            db_path=db_path,
            num_iterations=2,
            timeout_sec=60,
        )
    )

    assert output.error is None
    assert output.number_of_hits >= 1
    assert output.query_count == 1
    assert output.num_iterations == 2
    assert output.hits[0].qseqid == "query1"
    assert output.command_version
