from bio_mcp.core.subprocesses import CommandStatus
from bio_mcp.server import build_server
from bio_mcp.tools.aaindex import AAINDEX1_ENV_VAR
from bio_mcp.tools.health import AVAILABLE_TOOLS, bio_mcp_health_tool


EXPECTED_TOOLS = {
    "validate_protein_sequence",
    "protparam_analyze",
    "aaindex_lookup",
    "aaindex_list",
    "aaindex_sequence_features",
    "mafft_align",
    "clustalo_align",
    "blastp_local",
    "psiblast_local",
    "bio_mcp_health",
}


def test_server_registers_expected_tools():
    server = build_server()

    registered_tools = set(server._tool_manager._tools)

    assert registered_tools == EXPECTED_TOOLS


def test_health_reports_available_tools():
    output = bio_mcp_health_tool()

    assert output.package_version
    assert output.python_version
    assert set(AVAILABLE_TOOLS) == EXPECTED_TOOLS
    assert set(output.available_tools) == EXPECTED_TOOLS
    assert "biopython" in output.optional_dependencies
    assert "mcp" in output.optional_dependencies
    assert "mafft" in output.cli_binaries
    assert "clustalo" in output.cli_binaries
    assert "blastp" in output.cli_binaries
    assert "psiblast" in output.cli_binaries
    assert "makeblastdb" in output.cli_binaries
    assert output.aaindex_backend.available is True
    assert output.aaindex_backend.record_count >= 2
    assert output.execution_mode == "local"


def test_health_reports_bad_aaindex_backend_without_crashing(monkeypatch, tmp_path):
    missing_path = tmp_path / "missing-aaindex1"
    monkeypatch.setenv(AAINDEX1_ENV_VAR, str(missing_path))

    output = bio_mcp_health_tool()

    assert output.aaindex_backend.backend_source == "env_file"
    assert output.aaindex_backend.available is False
    assert output.aaindex_backend.record_count == 0
    assert output.aaindex_backend.error is not None
    assert "AAindex backend error" in " ".join(output.warnings)


def test_health_reports_cli_binary_versions(monkeypatch):
    monkeypatch.setattr(
        "bio_mcp.tools.health.detect_mafft",
        lambda: CommandStatus(
            available=True,
            path="/usr/bin/mafft",
            version="v7.520",
            resolution_source="PATH",
        ),
    )
    monkeypatch.setattr(
        "bio_mcp.tools.health.detect_clustalo",
        lambda: CommandStatus(
            available=False,
            path=None,
            version=None,
            resolution_source="missing",
            error="missing clustalo",
        ),
    )
    monkeypatch.setattr(
        "bio_mcp.tools.health.detect_blastp",
        lambda: CommandStatus(
            available=True,
            path="/usr/bin/blastp",
            version="blastp: 2.17.0+",
            resolution_source="PATH",
        ),
    )
    monkeypatch.setattr(
        "bio_mcp.tools.health.detect_psiblast",
        lambda: CommandStatus(
            available=False,
            path=None,
            version=None,
            resolution_source="missing",
            error="missing psiblast",
        ),
    )
    monkeypatch.setattr(
        "bio_mcp.tools.health.detect_makeblastdb",
        lambda: CommandStatus(
            available=True,
            path="/usr/bin/makeblastdb",
            version="makeblastdb: 2.17.0+",
            resolution_source="PATH",
        ),
    )

    output = bio_mcp_health_tool()

    assert output.cli_binaries["mafft"].available is True
    assert output.cli_binaries["mafft"].version == "v7.520"
    assert output.cli_binaries["mafft"].resolution_source == "PATH"
    assert output.cli_binaries["clustalo"].available is False
    assert output.cli_binaries["clustalo"].resolution_source == "missing"
    assert output.cli_binaries["clustalo"].error == "missing clustalo"
    assert output.cli_binaries["blastp"].available is True
    assert output.cli_binaries["blastp"].version == "blastp: 2.17.0+"
    assert output.cli_binaries["psiblast"].available is False
    assert output.cli_binaries["psiblast"].error == "missing psiblast"
    assert output.cli_binaries["makeblastdb"].available is True
