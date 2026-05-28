"""MCP server entrypoint for bio-mcp."""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from bio_mcp import __version__
from bio_mcp.core.sequence import MAX_PROTEIN_SEQUENCE_LENGTH
from bio_mcp.schemas import (
    AaindexListInput,
    AaindexLookupInput,
    AaindexSequenceFeaturesInput,
    BlastpLocalInput,
    ClustaloAlignInput,
    MafftAlignInput,
    PsiblastLocalInput,
    ProteinSequenceInput,
)
from bio_mcp.tools.aaindex import (
    aaindex_list_tool,
    aaindex_lookup_tool,
    aaindex_sequence_features_tool,
)
from bio_mcp.tools.alignment import clustalo_align_tool, mafft_align_tool
from bio_mcp.tools.blast import blastp_local_tool, psiblast_local_tool
from bio_mcp.tools.health import bio_mcp_health_tool
from bio_mcp.tools.protparam import protparam_analyze_tool
from bio_mcp.tools.validation import validate_protein_sequence_tool


def build_server():
    """Build and return the FastMCP server instance."""

    try:
        from mcp.server.fastmcp import FastMCP
    except ImportError as exc:  # pragma: no cover - depends on runtime environment
        raise RuntimeError(
            "The 'mcp' package is required to run the MCP server. "
            "Install with: python -m pip install -e ."
        ) from exc

    mcp = FastMCP("bio_mcp")

    @mcp.tool()
    def validate_protein_sequence(
        text: str,
        max_length: int = MAX_PROTEIN_SEQUENCE_LENGTH,
    ) -> dict:
        """Validate raw protein sequence or FASTA text and report invalid characters."""

        return validate_protein_sequence_tool(
            ProteinSequenceInput(text=text, max_length=max_length)
        ).model_dump()

    @mcp.tool()
    def protparam_analyze(
        text: str,
        max_length: int = MAX_PROTEIN_SEQUENCE_LENGTH,
    ) -> dict:
        """Analyze lightweight ProtParam-like protein properties locally."""

        return protparam_analyze_tool(
            ProteinSequenceInput(text=text, max_length=max_length)
        ).model_dump()

    @mcp.tool()
    def aaindex_lookup(
        index_id: str | None = None,
        query: str | None = None,
        limit: int = 20,
    ) -> dict:
        """Look up local AAindex metadata and amino acid values by id or query."""

        return aaindex_lookup_tool(
            AaindexLookupInput(index_id=index_id, query=query, limit=limit)
        ).model_dump()

    @mcp.tool()
    def aaindex_list(
        query: str | None = None,
        limit: int = 50,
    ) -> dict:
        """List local AAindex metadata without amino acid value tables."""

        return aaindex_list_tool(AaindexListInput(query=query, limit=limit)).model_dump()

    @mcp.tool()
    def aaindex_sequence_features(
        text: str,
        index_id: str,
        max_length: int = MAX_PROTEIN_SEQUENCE_LENGTH,
    ) -> dict:
        """Return per-residue AAindex values and summary stats for a protein sequence."""

        return aaindex_sequence_features_tool(
            AaindexSequenceFeaturesInput(
                text=text,
                index_id=index_id,
                max_length=max_length,
            )
        ).model_dump()

    @mcp.tool()
    def mafft_align(
        fasta_text: str,
        mode: str = "auto",
        timeout_sec: int = 60,
    ) -> dict:
        """Align protein FASTA records with the local MAFFT binary."""

        return mafft_align_tool(
            MafftAlignInput(
                fasta_text=fasta_text,
                mode=mode,
                timeout_sec=timeout_sec,
            )
        ).model_dump()

    @mcp.tool()
    def clustalo_align(
        fasta_text: str,
        timeout_sec: int = 60,
    ) -> dict:
        """Align protein FASTA records with the local Clustal Omega binary."""

        return clustalo_align_tool(
            ClustaloAlignInput(
                fasta_text=fasta_text,
                timeout_sec=timeout_sec,
            )
        ).model_dump()

    @mcp.tool()
    def blastp_local(
        query_fasta: str,
        db_path: str,
        evalue: float = 1e-5,
        max_target_seqs: int = 10,
        timeout_sec: int = 120,
    ) -> dict:
        """Search a local protein BLAST database with the local BLASTP binary."""

        return blastp_local_tool(
            BlastpLocalInput(
                query_fasta=query_fasta,
                db_path=db_path,
                evalue=evalue,
                max_target_seqs=max_target_seqs,
                timeout_sec=timeout_sec,
            )
        ).model_dump()

    @mcp.tool()
    def psiblast_local(
        query_fasta: str,
        db_path: str,
        num_iterations: int = 3,
        evalue: float = 1e-5,
        max_target_seqs: int = 10,
        timeout_sec: int = 300,
    ) -> dict:
        """Search a local protein BLAST database with the local PSI-BLAST binary."""

        return psiblast_local_tool(
            PsiblastLocalInput(
                query_fasta=query_fasta,
                db_path=db_path,
                num_iterations=num_iterations,
                evalue=evalue,
                max_target_seqs=max_target_seqs,
                timeout_sec=timeout_sec,
            )
        ).model_dump()

    @mcp.tool()
    def bio_mcp_health() -> dict:
        """Return package, Python, optional dependency, and tool availability status."""

        return bio_mcp_health_tool().model_dump()

    return mcp


def main(argv: Sequence[str] | None = None) -> int:
    """Run the MCP server over stdio."""

    parser = argparse.ArgumentParser(description="Run the local bio-mcp MCP server.")
    parser.add_argument(
        "--version",
        action="version",
        version=f"bio-mcp {__version__}",
    )
    parser.parse_args(argv)

    server = build_server()
    server.run()
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
