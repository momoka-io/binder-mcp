# Tool Catalog

This roadmap keeps Phase 1 local-only and lightweight. Biological inputs are not sent to remote services.

| Tool | Status | Notes |
| --- | --- | --- |
| `validate_protein_sequence` | phase1 | Implemented. Parses raw sequence or FASTA text, cleans to standard amino acids, and reports invalid characters. |
| `protparam_analyze` | phase1 | Implemented. Uses Biopython ProtParam when installed; otherwise returns stable local composition, molecular weight, aromaticity, and GRAVY with unsupported metrics marked. |
| `aaindex_lookup` | phase1 | Implemented with a small packaged fixture and an interface intended for a full AAindex backend later. |
| `aaindex_sequence_features` | phase1 | Implemented. Uses packaged AAindex values for per-residue features and summary statistics. |
| `bio_mcp_health` | phase1 | Implemented. Reports package version, Python version, optional dependencies, and available tools. |
| `mafft_align` | phase1B | Implemented. Local CLI wrapper only; requires a `mafft` binary on `PATH`, validates protein FASTA input, uses temporary directories, and enforces timeouts and input limits. |
| `clustalo_align` | phase1B | Implemented. Local CLI wrapper only; requires a `clustalo` binary on `PATH`, validates protein FASTA input, uses temporary directories, and enforces timeouts and input limits. |
| `blastp_local` | phase1B | Implemented. Local BLASTP wrapper only; requires BLAST+ `blastp` and a separately created local protein BLAST database prefix. No remote BLAST or database download. |
| `psiblast_local` | phase1B | Implemented. Local PSI-BLAST wrapper only; requires BLAST+ `psiblast` and a separately created local protein BLAST database prefix. Iterations are limited to 1-10. No remote BLAST or database download. |
| AlphaFold or other GPU structure models | blocked/needs API/license/GPU | Not implemented in Phase 1. Requires separate safety, dependency, and compute planning. |
| Rosetta | blocked/needs API/license/GPU | Not implemented in Phase 1. Requires license and heavyweight local setup. |
| GROMACS / AMBER / CHARMM molecular dynamics | blocked/needs API/license/GPU | Not implemented in Phase 1. Heavy simulation workflows are out of scope. |
| Docking tools | blocked/needs API/license/GPU | Not implemented in Phase 1. Requires separate tool-specific safety and licensing review. |
| Remote web-service wrappers | later | Not implemented. Any future remote wrapper must require `remote_ok: true` and document transmitted data. |
