import os

import pytest

from bio_mcp.schemas import ClustaloAlignInput, MafftAlignInput
from bio_mcp.tools.alignment import clustalo_align_tool, detect_clustalo, detect_mafft, mafft_align_tool


pytestmark = pytest.mark.skipif(
    os.environ.get("BIO_MCP_RUN_INTEGRATION") != "1",
    reason="Set BIO_MCP_RUN_INTEGRATION=1 to run local CLI integration tests.",
)

FASTA = ">seq1\nACDEFG\n>seq2\nACDFFG\n"


def test_mafft_align_real_cli():
    status = detect_mafft()
    if not status.available:
        pytest.skip(status.error or "MAFFT binary is not available.")

    output = mafft_align_tool(MafftAlignInput(fasta_text=FASTA, timeout_sec=30))

    assert output.error is None
    assert output.aligned_fasta is not None
    assert ">seq1" in output.aligned_fasta
    assert ">seq2" in output.aligned_fasta
    assert output.number_of_sequences == 2
    assert output.command_version


def test_clustalo_align_real_cli():
    status = detect_clustalo()
    if not status.available:
        pytest.skip(status.error or "Clustal Omega binary is not available.")

    output = clustalo_align_tool(ClustaloAlignInput(fasta_text=FASTA, timeout_sec=30))

    assert output.error is None
    assert output.aligned_fasta is not None
    assert ">seq1" in output.aligned_fasta
    assert ">seq2" in output.aligned_fasta
    assert output.number_of_sequences == 2
    assert output.command_version
