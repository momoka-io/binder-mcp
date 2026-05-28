"""Health/status tool for the local bio-mcp server."""

from __future__ import annotations

import platform
from importlib import metadata
from importlib.util import find_spec

from bio_mcp import __version__
from bio_mcp.schemas import (
    CommandDependencyStatus,
    HealthOutput,
    OptionalDependencyStatus,
    ToolProvenance,
)
from bio_mcp.tools.alignment import detect_clustalo, detect_mafft
from bio_mcp.tools.aaindex import aaindex_backend_status
from bio_mcp.tools.blast import detect_blastp, detect_makeblastdb, detect_psiblast

TOOL_NAME = "bio_mcp_health"
AVAILABLE_TOOLS = [
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
]
OPTIONAL_DEPENDENCIES = {
    "biopython": "Bio",
    "mcp": "mcp",
}


def _dependency_status(distribution_name: str, import_name: str) -> OptionalDependencyStatus:
    if find_spec(import_name) is None:
        return OptionalDependencyStatus(installed=False, version=None)
    try:
        version = metadata.version(distribution_name)
    except metadata.PackageNotFoundError:
        version = None
    return OptionalDependencyStatus(installed=True, version=version)


def bio_mcp_health_tool() -> HealthOutput:
    """Return local package, Python, optional dependency, and tool availability status."""

    optional_dependencies = {
        name: _dependency_status(name, import_name)
        for name, import_name in OPTIONAL_DEPENDENCIES.items()
    }
    warnings: list[str] = []
    if not optional_dependencies["biopython"].installed:
        warnings.append("Biopython is not installed; protparam_analyze will use fallback metrics.")
    if not optional_dependencies["mcp"].installed:
        warnings.append("mcp is not installed; importable tool functions work but the MCP server cannot run.")

    command_statuses = {
        "mafft": detect_mafft(),
        "clustalo": detect_clustalo(),
        "blastp": detect_blastp(),
        "psiblast": detect_psiblast(),
        "makeblastdb": detect_makeblastdb(),
    }
    cli_binaries = {
        name: CommandDependencyStatus(
            available=status.available,
            path=status.path,
            version=status.version,
            resolution_source=status.resolution_source,
            error=status.error,
        )
        for name, status in command_statuses.items()
    }
    for name, status in command_statuses.items():
        if not status.available and status.error:
            warnings.append(status.error)
        elif status.available and status.error:
            warnings.append(f"{name} was found but version detection reported: {status.error}")

    aaindex_status = aaindex_backend_status()
    warnings.extend(aaindex_status.warnings)
    if aaindex_status.error:
        warnings.append(f"AAindex backend error: {aaindex_status.error}")

    return HealthOutput(
        package_version=__version__,
        python_version=platform.python_version(),
        optional_dependencies=optional_dependencies,
        cli_binaries=cli_binaries,
        aaindex_backend=aaindex_status,
        available_tools=AVAILABLE_TOOLS,
        execution_mode="local",
        warnings=warnings,
        provenance=ToolProvenance(
            tool_name=TOOL_NAME,
            wrapper_version=__version__,
            backend="bio_mcp",
            execution_mode="local",
            warnings=warnings,
        ),
    )
