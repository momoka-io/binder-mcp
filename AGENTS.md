# AGENTS.md

## Project goal

Build a Python MCP server named `binder-mcp` that exposes selected bioinformatics tools as safe, typed, testable MCP tools.

The project will be implemented incrementally. Do not attempt to wrap all tools at once.

## Priority order

Phase 1 must be local-only and lightweight:

1. FASTA and protein sequence validation
2. ProtParam-like protein property analysis
3. AAindex lookup and per-sequence AAindex feature extraction
4. Tool health/status endpoint

Phase 1B may add lightweight local CLI wrappers:

1. MAFFT
2. Clustal Omega
3. Local BLAST/PSI-BLAST only if binaries and databases are explicitly configured

Do not implement AlphaFold, Rosetta, GROMACS, AMBER, CHARMM, docking, GPU models, or remote web-service wrappers until Phase 1 is complete and tested.

## Architecture

Use Python 3.11+.

Preferred layout:

- `src/binder_mcp/server.py`: MCP server entrypoint
- `src/binder_mcp/schemas.py`: Pydantic input/output models
- `src/binder_mcp/tools/`: individual tool modules
- `src/binder_mcp/core/`: shared validation, sequence parsing, subprocess helpers
- `tests/`: pytest tests
- `docs/tool_catalog.md`: tool roadmap and status
- `README.md`: install, run, Codex MCP config, examples

Use the official MCP Python SDK / FastMCP style unless there is a strong reason not to.

## Tool design rules

Every MCP tool must:

- Have a clear name, docstring, typed input, and typed structured output.
- Validate input sequences before running computation.
- Return JSON-serializable output.
- Include provenance fields where relevant:
  - tool name
  - wrapper version
  - backend package or binary version, if detectable
  - local/remote execution mode
  - warnings
- Avoid writing outside the project workspace.
- Never use `shell=True` for subprocess calls.
- Use timeouts for subprocesses.
- Use temporary working directories for tools that produce files.
- Return useful error messages instead of raw stack traces.

## Safety and privacy

By default, do not send user sequences, structures, antibodies, epitopes, or other biological data to remote services.

Any remote service wrapper must require an explicit `remote_ok: true` parameter and must document what data is sent.

Do not scrape websites unless the site explicitly allows it or provides a documented API.

Do not provide wet-lab protocols, organism engineering instructions, or autonomous optimization loops. These wrappers are for computational analysis and require human review.

## Testing requirements

For every implemented tool:

- Add unit tests for valid input.
- Add tests for invalid sequence characters.
- Add tests for empty input.
- Add tests for maximum length or timeout behavior where relevant.
- Add at least one snapshot/example output in `README.md` or `docs/examples.md`.

Before considering a task complete, run:

```bash
pytest
python -m bio_mcp.server --help || true