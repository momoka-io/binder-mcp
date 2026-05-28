# bio-mcp

`bio-mcp` is a local-only Python MCP server for lightweight bioinformatics sequence analysis. Phase 1 implements FASTA/protein validation, ProtParam-like metrics, AAindex1 lookup and sequence features, and health/status reporting. Phase 1B adds lightweight local MAFFT, Clustal Omega, BLASTP, and PSI-BLAST wrappers.

No Phase 1 or Phase 1B tool makes network calls, scrapes websites, invokes heavy/GPU compute tools, or sends biological data to remote services. System packages are never installed automatically.

The repository checkout directory may still be named `binder-mcp` for historical reasons. The project display name is `bio-mcp`, the Python module is `bio_mcp`, the Codex MCP server key is `bio_mcp`, and the CLI command is `bio-mcp`.

## Install

```bash
python -m pip install -e .[dev]
```

For runtime-only use:

```bash
python -m pip install -e .
```

Optional local BLAST+ binaries can be installed outside `bio-mcp`, for example:

```bash
conda install -c bioconda -c conda-forge blast
```

## Run

```bash
python -m bio_mcp.server
```

Help/version checks:

```bash
python -m bio_mcp.server --help
bio-mcp --version
```

## Codex `config.toml` Example

```toml
[mcp_servers.bio_mcp]
command = "python"
args = ["-m", "bio_mcp.server"]
cwd = "E:/gitrepository/binder-mcp"
```

## Tools

### `validate_protein_sequence`

Input:

```json
{
  "text": ">seq1\nACD-XZ*\n"
}
```

Example output excerpt:

```json
{
  "cleaned_sequence": "ACD",
  "length": 3,
  "invalid_characters": [
    { "character": "-", "count": 1, "positions": [4] },
    { "character": "X", "count": 1, "positions": [5] },
    { "character": "Z", "count": 1, "positions": [6] },
    { "character": "*", "count": 1, "positions": [7] }
  ],
  "is_valid": false,
  "warnings": [
    "Input contained a FASTA header; sequence lines were parsed.",
    "Invalid characters were omitted from cleaned_sequence."
  ]
}
```

### `protparam_analyze`

Input:

```json
{
  "text": "ACDEFWY"
}
```

Example output excerpt with Biopython installed:

```json
{
  "length": 7,
  "amino_acid_composition": { "A": 1, "C": 1, "D": 1, "E": 1, "F": 1, "G": 0 },
  "molecular_weight": 932.9943999999999,
  "theoretical_pi": 4.0500284194946286,
  "aromaticity": 0.42857142857142855,
  "instability_index": 98.85714285714286,
  "gravy": -0.3,
  "unsupported_metrics": []
}
```

When Biopython is not installed, `protparam_analyze` uses local fallback metrics and returns `null` for `theoretical_pi` and `instability_index` while listing them in `unsupported_metrics`.

### `aaindex_lookup`

AAindex support has two local backend modes:

- Default mode uses a small packaged fixture containing `KYTJ820101` and `HOPT810101`.
- Full local mode parses a local AAindex1 flat file selected by `BIO_MCP_AAINDEX1_PATH`.

```bash
export BIO_MCP_AAINDEX1_PATH=/path/to/aaindex1
```

`bio-mcp` does not download AAindex automatically. AAindex parsing is local-only, and no biological data is sent to remote services. AAindex1 is supported; AAindex2 and AAindex3 are not implemented yet.

Input:

```json
{
  "index_id": "KYTJ820101"
}
```

Example output excerpt:

```json
{
  "query": "KYTJ820101",
  "matches": [
    {
      "accession": "KYTJ820101",
      "title": "Hydropathy index",
      "description": "Kyte-Doolittle hydropathy scale packaged as a Phase 1 AAindex fixture.",
      "pmid": null,
      "authors": "Kyte J. and Doolittle R.F.",
      "journal": "J. Mol. Biol. 157, 105-132 (1982)",
      "correlations": {},
      "values": {
        "A": 1.8,
        "C": 2.5,
        "D": -3.5
      }
    }
  ],
  "backend_source": "packaged_fixture",
  "backend_path": null,
  "record_count": 2,
  "warnings": [
    "Using a small packaged AAindex fixture. Set BIO_MCP_AAINDEX1_PATH to a local AAindex1 flat file for full local AAindex1 support."
  ],
  "error": null
}
```

Query search:

```json
{
  "query": "hydrophilicity",
  "limit": 20
}
```

### `aaindex_list`

Input:

```json
{
  "query": "hydro",
  "limit": 50
}
```

Example output excerpt:

```json
{
  "records": [
    {
      "accession": "KYTJ820101",
      "title": "Hydropathy index",
      "description": "Kyte-Doolittle hydropathy scale packaged as a Phase 1 AAindex fixture.",
      "pmid": null,
      "authors": "Kyte J. and Doolittle R.F.",
      "journal": "J. Mol. Biol. 157, 105-132 (1982)",
      "correlations": {}
    }
  ],
  "record_count": 2,
  "backend_source": "packaged_fixture",
  "backend_path": null,
  "warnings": [
    "Using a small packaged AAindex fixture. Set BIO_MCP_AAINDEX1_PATH to a local AAindex1 flat file for full local AAindex1 support."
  ],
  "error": null
}
```

### `aaindex_sequence_features`

Input:

```json
{
  "text": "ACD",
  "index_id": "KYTJ820101"
}
```

Example output excerpt:

```json
{
  "index_id": "KYTJ820101",
  "sequence": "ACD",
  "length": 3,
  "per_residue_values": [
    { "position": 1, "residue": "A", "value": 1.8 },
    { "position": 2, "residue": "C", "value": 2.5 },
    { "position": 3, "residue": "D", "value": -3.5 }
  ],
  "summary_stats": {
    "mean": 0.26666666666666666,
    "min": -3.5,
    "max": 2.5
  },
  "missing_residues": [],
  "backend_source": "packaged_fixture",
  "backend_path": null,
  "record_count": 2,
  "warnings": [
    "Using a small packaged AAindex fixture. Set BIO_MCP_AAINDEX1_PATH to a local AAindex1 flat file for full local AAindex1 support."
  ],
  "error": null
}
```

With a full local AAindex1 backend:

```bash
export BIO_MCP_AAINDEX1_PATH=/path/to/aaindex1
```

```json
{
  "text": "ACD",
  "index_id": "KYTJ820101"
}
```

The response uses `backend_source: "env_file"` and reports `backend_path` and the parsed `record_count`.

### `mafft_align`

Requires a local `mafft` binary. The wrapper resolves it from `BIO_MCP_MAFFT_BIN`
first, then from `PATH`.

Input:

```json
{
  "fasta_text": ">seq1\nACDEFG\n>seq2\nACDFFG\n",
  "mode": "auto",
  "timeout_sec": 60
}
```

Missing binary output excerpt:

```json
{
  "aligned_fasta": null,
  "number_of_sequences": 2,
  "command_version": null,
  "stderr_summary": "",
  "warnings": [
    "Dependency error: required binary 'mafft' was not found on PATH."
  ],
  "error": "Dependency error: required binary 'mafft' was not found on PATH."
}
```

Successful execution output excerpt, from a mocked test or a machine with MAFFT installed:

```json
{
  "aligned_fasta": ">seq1\nACDEFG-\n>seq2\nACD-FFG\n",
  "number_of_sequences": 2,
  "command_version": "v7.520",
  "stderr_summary": "Progressive alignment complete",
  "warnings": [],
  "error": null
}
```

Supported MAFFT modes are `auto`, `localpair`, `globalpair`, and `genafpair`.

### `clustalo_align`

Requires a local `clustalo` binary. The wrapper resolves it from
`BIO_MCP_CLUSTALO_BIN` first, then from `PATH`.

Input:

```json
{
  "fasta_text": ">seq1\nACDEFG\n>seq2\nACDFFG\n",
  "timeout_sec": 60
}
```

Missing binary output excerpt:

```json
{
  "aligned_fasta": null,
  "number_of_sequences": 2,
  "command_version": null,
  "stderr_summary": "",
  "warnings": [
    "Dependency error: required binary 'clustalo' was not found on PATH."
  ],
  "error": "Dependency error: required binary 'clustalo' was not found on PATH."
}
```

Successful execution output excerpt, from a mocked test or a machine with Clustal Omega installed:

```json
{
  "aligned_fasta": ">seq1\nACDEFG-\n>seq2\nACD-FFG\n",
  "number_of_sequences": 2,
  "command_version": "1.2.4",
  "stderr_summary": "clustalo done",
  "warnings": [],
  "error": null
}
```

### `blastp_local`

Requires a local BLAST+ `blastp` binary and a local protein BLAST database prefix. The wrapper resolves `blastp` from `BIO_MCP_BLASTP_BIN` first, then from `PATH`. It is local-only, makes no network calls, does not download databases, and is not NCBI remote BLAST.

Create a local protein database separately:

```bash
makeblastdb -in proteins.fasta -dbtype prot -out proteins_db
```

Input:

```json
{
  "query_fasta": ">query1\nMKTAYIAKQRQISFVKSHFSRQDILD\n",
  "db_path": "/path/to/proteins_db",
  "evalue": 1e-5,
  "max_target_seqs": 10,
  "timeout_sec": 120
}
```

Example output excerpt:

```json
{
  "hits": [
    {
      "qseqid": "query1",
      "sseqid": "subject1",
      "pident": 100.0,
      "length": 26,
      "mismatch": 0,
      "gapopen": 0,
      "qstart": 1,
      "qend": 26,
      "sstart": 1,
      "send": 26,
      "evalue": 1e-20,
      "bitscore": 55.5
    }
  ],
  "raw_tabular": "query1\tsubject1\t100.000\t26\t0\t0\t1\t26\t1\t26\t1e-20\t55.5\n",
  "number_of_hits": 1,
  "query_count": 1,
  "db_path": "/path/to/proteins_db",
  "command_version": "blastp: 2.17.0+",
  "stderr_summary": "",
  "warnings": [],
  "error": null
}
```

### `psiblast_local`

Requires a local BLAST+ `psiblast` binary and a local protein BLAST database prefix. The wrapper resolves `psiblast` from `BIO_MCP_PSIBLAST_BIN` first, then from `PATH`. It is local-only, makes no network calls, does not download databases, and is not NCBI remote BLAST. `num_iterations` is constrained to 1-10.

Input:

```json
{
  "query_fasta": ">query1\nMKTAYIAKQRQISFVKSHFSRQDILD\n",
  "db_path": "/path/to/proteins_db",
  "num_iterations": 3,
  "evalue": 1e-5,
  "max_target_seqs": 10,
  "timeout_sec": 300
}
```

Example output excerpt:

```json
{
  "hits": [],
  "raw_tabular": "",
  "number_of_hits": 0,
  "query_count": 1,
  "db_path": "/path/to/proteins_db",
  "num_iterations": 3,
  "command_version": "psiblast: 2.17.0+",
  "stderr_summary": "",
  "warnings": [],
  "error": null
}
```

### `bio_mcp_health`

Example output excerpt:

```json
{
  "package_version": "0.1.0",
  "python_version": "3.11.9",
  "optional_dependencies": {
    "biopython": { "installed": true, "version": "1.87" },
    "mcp": { "installed": true, "version": "1.27.1" }
  },
  "cli_binaries": {
    "mafft": {
      "available": false,
      "path": null,
      "version": null,
      "resolution_source": "missing",
      "error": "Dependency error: required binary 'mafft' was not found on PATH."
    },
    "clustalo": {
      "available": true,
      "path": "/usr/bin/clustalo",
      "version": "1.2.4",
      "resolution_source": "PATH",
      "error": null
    },
    "blastp": {
      "available": true,
      "path": "/home/ttyang/miniconda3/envs/bio-mcp/bin/blastp",
      "version": "blastp: 2.17.0+",
      "resolution_source": "PATH",
      "error": null
    },
    "psiblast": {
      "available": true,
      "path": "/home/ttyang/miniconda3/envs/bio-mcp/bin/psiblast",
      "version": "psiblast: 2.17.0+",
      "resolution_source": "PATH",
      "error": null
    },
    "makeblastdb": {
      "available": true,
      "path": "/home/ttyang/miniconda3/envs/bio-mcp/bin/makeblastdb",
      "version": "makeblastdb: 2.17.0+",
      "resolution_source": "PATH",
      "error": null
    }
  },
  "aaindex_backend": {
    "backend_source": "packaged_fixture",
    "backend_path": null,
    "available": true,
    "record_count": 2,
    "error": null,
    "warnings": [
      "Using a small packaged AAindex fixture. Set BIO_MCP_AAINDEX1_PATH to a local AAindex1 flat file for full local AAindex1 support."
    ]
  },
  "available_tools": [
    "validate_protein_sequence",
    "protparam_analyze",
    "aaindex_lookup",
    "aaindex_list",
    "aaindex_sequence_features",
    "mafft_align",
    "clustalo_align",
    "blastp_local",
    "psiblast_local",
    "bio_mcp_health"
  ],
  "execution_mode": "local"
}
```

## Development

```bash
pytest
python -m bio_mcp.server --help
```

The default pytest suite is local-only and may use mocks for CLI wrapper behavior.
To run real local smoke and CLI integration checks:

```bash
python scripts/smoke_existing_tools.py
BIO_MCP_RUN_INTEGRATION=1 pytest -q
```

Optional binary override examples:

```bash
BIO_MCP_MAFFT_BIN=/home/ttyang/miniconda3/envs/bio-mcp/bin/mafft
BIO_MCP_CLUSTALO_BIN=/home/ttyang/miniconda3/envs/bio-mcp/bin/clustalo
BIO_MCP_BLASTP_BIN=/home/ttyang/miniconda3/envs/bio-mcp/bin/blastp
BIO_MCP_PSIBLAST_BIN=/home/ttyang/miniconda3/envs/bio-mcp/bin/psiblast
BIO_MCP_MAKEBLASTDB_BIN=/home/ttyang/miniconda3/envs/bio-mcp/bin/makeblastdb
```

MAFFT, Clustal Omega, BLASTP, and PSI-BLAST wrappers are local-only and make no network calls.

## Known Phase 1 Limitations

- Full AAindex1 coverage requires a local AAindex1 flat file via `BIO_MCP_AAINDEX1_PATH`.
- AAindex2 and AAindex3 are not implemented yet.
- AAindex data is not downloaded automatically.
- Without Biopython, `protparam_analyze` does not compute theoretical pI or instability index.
- Only the 20 standard amino acid one-letter codes are accepted for computation.
- MAFFT and Clustal Omega wrappers require local binaries via environment override or `PATH`; they are not installed or downloaded by `bio-mcp`.
- BLAST/PSI-BLAST require local BLAST+ binaries via environment override or `PATH`; they are not installed or downloaded by `bio-mcp`.
- BLAST/PSI-BLAST require a local protein database created separately, for example with `makeblastdb -dbtype prot`.
- No nucleotide BLAST wrappers are implemented yet.
- No remote BLAST wrappers are implemented yet.
- Phase 1B does not include AlphaFold, Rosetta, docking, molecular dynamics, GPU tools, or remote services.
