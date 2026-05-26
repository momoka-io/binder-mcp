"""ProtParam-like protein property analysis."""

from __future__ import annotations

from importlib import metadata
from importlib.util import find_spec

from bio_mcp import __version__
from bio_mcp.core.sequence import STANDARD_AMINO_ACIDS, require_valid_protein_sequence
from bio_mcp.schemas import ProteinSequenceInput, ProtParamOutput, ToolProvenance

TOOL_NAME = "protparam_analyze"

AVERAGE_RESIDUE_MASS = {
    "A": 71.0788,
    "R": 156.1875,
    "N": 114.1038,
    "D": 115.0886,
    "C": 103.1388,
    "E": 129.1155,
    "Q": 128.1307,
    "G": 57.0519,
    "H": 137.1411,
    "I": 113.1594,
    "L": 113.1594,
    "K": 128.1741,
    "M": 131.1926,
    "F": 147.1766,
    "P": 97.1167,
    "S": 87.0782,
    "T": 101.1051,
    "W": 186.2132,
    "Y": 163.1760,
    "V": 99.1326,
}

KYTE_DOOLITTLE = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}

WATER_MASS = 18.0153


def protparam_analyze_tool(params: ProteinSequenceInput) -> ProtParamOutput:
    """Analyze local protein sequence properties using Biopython when available."""

    validation = require_valid_protein_sequence(params.text, max_length=params.max_length)
    sequence = validation.cleaned_sequence
    warnings = list(validation.warnings)

    if find_spec("Bio") is not None:
        try:
            return _analyze_with_biopython(sequence, warnings)
        except Exception as exc:  # pragma: no cover - defensive fallback
            warnings.append(f"Biopython ProtParam failed; used built-in fallback metrics. Reason: {exc}")

    return _analyze_with_fallback(sequence, warnings)


def _composition(sequence: str) -> dict[str, int]:
    return {residue: sequence.count(residue) for residue in STANDARD_AMINO_ACIDS}


def _percentages(composition: dict[str, int], length: int) -> dict[str, float]:
    return {residue: count / length for residue, count in composition.items()}


def _analyze_with_biopython(sequence: str, warnings: list[str]) -> ProtParamOutput:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    analysis = ProteinAnalysis(sequence)
    composition = {residue: int(analysis.count_amino_acids().get(residue, 0)) for residue in STANDARD_AMINO_ACIDS}
    try:
        backend_version = metadata.version("biopython")
    except metadata.PackageNotFoundError:
        backend_version = None

    return ProtParamOutput(
        length=len(sequence),
        amino_acid_composition=composition,
        amino_acid_percent=_percentages(composition, len(sequence)),
        molecular_weight=analysis.molecular_weight(),
        theoretical_pi=analysis.isoelectric_point(),
        aromaticity=analysis.aromaticity(),
        instability_index=analysis.instability_index(),
        gravy=analysis.gravy(),
        unsupported_metrics=[],
        warnings=warnings,
        provenance=ToolProvenance(
            tool_name=TOOL_NAME,
            wrapper_version=__version__,
            backend="biopython.ProteinAnalysis",
            backend_version=backend_version,
            execution_mode="local",
            warnings=warnings,
        ),
    )


def _analyze_with_fallback(sequence: str, warnings: list[str]) -> ProtParamOutput:
    composition = _composition(sequence)
    length = len(sequence)
    unsupported_metrics = ["theoretical_pi", "instability_index"]
    warnings = [
        *warnings,
        "Biopython is not installed; theoretical_pi and instability_index are unsupported.",
    ]
    molecular_weight = sum(AVERAGE_RESIDUE_MASS[residue] for residue in sequence) + WATER_MASS
    aromaticity = (composition["F"] + composition["W"] + composition["Y"]) / length
    gravy = sum(KYTE_DOOLITTLE[residue] for residue in sequence) / length

    return ProtParamOutput(
        length=length,
        amino_acid_composition=composition,
        amino_acid_percent=_percentages(composition, length),
        molecular_weight=molecular_weight,
        theoretical_pi=None,
        aromaticity=aromaticity,
        instability_index=None,
        gravy=gravy,
        unsupported_metrics=unsupported_metrics,
        warnings=warnings,
        provenance=ToolProvenance(
            tool_name=TOOL_NAME,
            wrapper_version=__version__,
            backend="bio_mcp.fallback_protparam",
            execution_mode="local",
            warnings=warnings,
        ),
    )
