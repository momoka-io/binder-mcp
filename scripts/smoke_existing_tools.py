#!/usr/bin/env python
"""Smoke-test existing Phase 1 and Phase 1B bio-mcp tools."""

from __future__ import annotations

import sys
import os
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

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
    AAINDEX1_ENV_VAR,
    aaindex_list_tool,
    aaindex_lookup_tool,
    aaindex_sequence_features_tool,
)
from bio_mcp.tools.alignment import clustalo_align_tool, detect_clustalo, detect_mafft, mafft_align_tool
from bio_mcp.tools.blast import (
    blastp_local_tool,
    detect_blastp,
    detect_makeblastdb,
    detect_psiblast,
    psiblast_local_tool,
)
from bio_mcp.tools.health import bio_mcp_health_tool
from bio_mcp.tools.protparam import protparam_analyze_tool
from bio_mcp.tools.validation import validate_protein_sequence_tool


FASTA = ">seq1\nACDEFG\n>seq2\nACDFFG\n"
BLAST_QUERY_FASTA = ">query1\nMKTAYIAKQRQISFVKSHFSRQDILD\n"
BLAST_DB_FASTA = (
    ">subject1\nMKTAYIAKQRQISFVKSHFSRQDILD\n"
    ">subject2\nGAVLIPFWYDERKSTNQCMH\n"
)


def main() -> int:
    results: list[tuple[str, str, str]] = []

    def check(name: str, func) -> None:
        try:
            message = func()
        except Exception as exc:  # noqa: BLE001 - smoke script reports all failures.
            results.append((name, "FAIL", str(exc)))
            return
        if isinstance(message, tuple):
            status, detail = message
            results.append((name, status, detail))
        else:
            results.append((name, "PASS", message or ""))

    check("validate_protein_sequence valid", _smoke_validate_valid)
    check("validate_protein_sequence invalid", _smoke_validate_invalid)
    check("protparam_analyze", _smoke_protparam)
    check("aaindex_backend", _smoke_aaindex_backend)
    check("aaindex_lookup", _smoke_aaindex_lookup)
    check("aaindex_sequence_features", _smoke_aaindex_features)
    check("bio_mcp_health", _smoke_health)
    check("mafft_align", _smoke_mafft)
    check("clustalo_align", _smoke_clustalo)
    check("blastp_local", _smoke_blastp)
    check("psiblast_local", _smoke_psiblast)

    for name, status, detail in results:
        suffix = f" - {detail}" if detail else ""
        print(f"{status}: {name}{suffix}")

    passed = sum(1 for _, status, _ in results if status == "PASS")
    skipped = sum(1 for _, status, _ in results if status == "SKIP")
    failed = sum(1 for _, status, _ in results if status == "FAIL")
    print(f"SUMMARY: PASS={passed} SKIP={skipped} FAIL={failed}")
    return 1 if failed else 0


def _smoke_validate_valid() -> str:
    output = validate_protein_sequence_tool(ProteinSequenceInput(text="ACDEFG"))
    assert output.is_valid is True
    assert output.length == 6
    return f"length={output.length}"


def _smoke_validate_invalid() -> str:
    output = validate_protein_sequence_tool(ProteinSequenceInput(text="ACD*"))
    assert output.is_valid is False
    assert output.invalid_characters
    return f"invalid={output.invalid_characters[0].character}"


def _smoke_protparam() -> str:
    output = protparam_analyze_tool(ProteinSequenceInput(text="ACDEFWY"))
    assert output.length == 7
    assert output.molecular_weight is not None
    return f"length={output.length}"


def _smoke_aaindex_lookup() -> str:
    index_id = _smoke_aaindex_index_id()
    output = aaindex_lookup_tool(AaindexLookupInput(index_id=index_id))
    assert output.error is None, output.error
    assert output.matches
    assert output.matches[0].accession == index_id
    return (
        f"backend={output.backend_source}, records={output.record_count}, "
        f"id={index_id}, matches={len(output.matches)}"
    )


def _smoke_aaindex_features() -> str:
    index_id, sequence = _smoke_aaindex_feature_input()
    output = aaindex_sequence_features_tool(
        AaindexSequenceFeaturesInput(text=sequence, index_id=index_id)
    )
    assert output.error is None, output.error
    assert output.summary_stats.mean is not None
    return f"backend={output.backend_source}, id={index_id}, mean={output.summary_stats.mean:.3f}"


def _smoke_aaindex_backend() -> str:
    listing = aaindex_list_tool(AaindexListInput(limit=1))
    assert listing.error is None, listing.error
    if AAINDEX1_ENV_VAR in os.environ:
        assert listing.backend_source == "env_file"
    else:
        assert listing.backend_source == "packaged_fixture"
    return f"source={listing.backend_source}, records={listing.record_count}"


def _smoke_aaindex_index_id() -> str:
    listing = aaindex_list_tool(AaindexListInput(limit=1))
    assert listing.error is None, listing.error
    assert listing.records
    return listing.records[0].accession


def _smoke_aaindex_feature_input() -> tuple[str, str]:
    lookup = aaindex_lookup_tool(AaindexLookupInput(index_id=_smoke_aaindex_index_id()))
    assert lookup.error is None, lookup.error
    assert lookup.matches and lookup.matches[0].values
    record = lookup.matches[0]
    residues = [residue for residue, value in record.values.items() if value is not None]
    assert residues
    sequence = "".join(residues[:3]) if len(residues) >= 3 else residues[0]
    return record.accession, sequence


def _smoke_health() -> str:
    output = bio_mcp_health_tool()
    assert "bio_mcp_health" in output.available_tools
    return f"tools={len(output.available_tools)}"


def _smoke_mafft() -> tuple[str, str] | str:
    status = detect_mafft()
    if not status.available:
        return ("SKIP", status.error or "MAFFT is not available.")
    output = mafft_align_tool(MafftAlignInput(fasta_text=FASTA, timeout_sec=30))
    assert output.error is None, output.error
    assert output.aligned_fasta and ">seq1" in output.aligned_fasta
    return f"version={output.command_version}"


def _smoke_clustalo() -> tuple[str, str] | str:
    status = detect_clustalo()
    if not status.available:
        return ("SKIP", status.error or "Clustal Omega is not available.")
    output = clustalo_align_tool(ClustaloAlignInput(fasta_text=FASTA, timeout_sec=30))
    assert output.error is None, output.error
    assert output.aligned_fasta and ">seq1" in output.aligned_fasta
    return f"version={output.command_version}"


def _build_tiny_blast_db(tmp_path: Path) -> tuple[str, str] | tuple[None, str]:
    makeblastdb = detect_makeblastdb()
    if not makeblastdb.available:
        return None, makeblastdb.error or "makeblastdb is not available."
    fasta_path = tmp_path / "proteins.fasta"
    db_prefix = tmp_path / "proteins_db"
    fasta_path.write_text(BLAST_DB_FASTA, encoding="utf-8")
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
    if result.returncode != 0:
        return None, result.stderr.strip() or result.stdout.strip() or "makeblastdb failed."
    return str(db_prefix), makeblastdb.version or "unknown"


def _smoke_blastp() -> tuple[str, str] | str:
    status = detect_blastp()
    if not status.available:
        return ("SKIP", status.error or "BLASTP is not available.")
    with TemporaryDirectory(prefix="bio-mcp-smoke-blastp-") as tmpdir:
        db_path, db_detail = _build_tiny_blast_db(Path(tmpdir))
        if db_path is None:
            return ("SKIP", db_detail)
        output = blastp_local_tool(
            BlastpLocalInput(query_fasta=BLAST_QUERY_FASTA, db_path=db_path, timeout_sec=30)
        )
    assert output.error is None, output.error
    assert output.number_of_hits >= 1
    return f"version={output.command_version}, hits={output.number_of_hits}, makeblastdb={db_detail}"


def _smoke_psiblast() -> tuple[str, str] | str:
    status = detect_psiblast()
    if not status.available:
        return ("SKIP", status.error or "PSI-BLAST is not available.")
    with TemporaryDirectory(prefix="bio-mcp-smoke-psiblast-") as tmpdir:
        db_path, db_detail = _build_tiny_blast_db(Path(tmpdir))
        if db_path is None:
            return ("SKIP", db_detail)
        output = psiblast_local_tool(
            PsiblastLocalInput(
                query_fasta=BLAST_QUERY_FASTA,
                db_path=db_path,
                num_iterations=2,
                timeout_sec=60,
            )
        )
    assert output.error is None, output.error
    assert output.number_of_hits >= 1
    return f"version={output.command_version}, hits={output.number_of_hits}, makeblastdb={db_detail}"


if __name__ == "__main__":
    raise SystemExit(main())
