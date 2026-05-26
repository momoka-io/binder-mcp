from bio_mcp.tools.health import AVAILABLE_TOOLS, bio_mcp_health_tool
from bio_mcp.core.subprocesses import CommandStatus


def test_health_reports_available_tools():
    output = bio_mcp_health_tool()

    assert output.package_version
    assert output.python_version
    assert set(AVAILABLE_TOOLS).issubset(output.available_tools)
    assert "biopython" in output.optional_dependencies
    assert "mcp" in output.optional_dependencies
    assert "mafft" in output.cli_binaries
    assert "clustalo" in output.cli_binaries
    assert output.execution_mode == "local"


def test_health_reports_cli_binary_versions(monkeypatch):
    monkeypatch.setattr(
        "bio_mcp.tools.health.detect_mafft",
        lambda: CommandStatus(available=True, path="/usr/bin/mafft", version="v7.520"),
    )
    monkeypatch.setattr(
        "bio_mcp.tools.health.detect_clustalo",
        lambda: CommandStatus(available=False, path=None, version=None, error="missing clustalo"),
    )

    output = bio_mcp_health_tool()

    assert output.cli_binaries["mafft"].available is True
    assert output.cli_binaries["mafft"].version == "v7.520"
    assert output.cli_binaries["clustalo"].available is False
    assert output.cli_binaries["clustalo"].error == "missing clustalo"
