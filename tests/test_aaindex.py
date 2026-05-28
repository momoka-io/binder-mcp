from pathlib import Path

import pytest

from bio_mcp.core.aaindex_parser import AaindexParseError, parse_aaindex1_text
from bio_mcp.core.sequence import MAX_PROTEIN_SEQUENCE_LENGTH, ProteinSequenceError
from bio_mcp.schemas import AaindexListInput, AaindexLookupInput, AaindexSequenceFeaturesInput
from bio_mcp.tools.aaindex import (
    AAINDEX1_ENV_VAR,
    aaindex_list_tool,
    aaindex_lookup_tool,
    aaindex_sequence_features_tool,
)

FIXTURE_PATH = Path(__file__).parent / "fixtures" / "aaindex1_sample.txt"


def test_parse_sample_aaindex1_file_with_multiple_records():
    records = parse_aaindex1_text(FIXTURE_PATH.read_text(encoding="utf-8"))

    assert [record.accession for record in records] == [
        "KYTJ820101",
        "HOPT810101",
        "MISS000101",
    ]


def test_parse_kytj820101_values_and_metadata():
    record = parse_aaindex1_text(FIXTURE_PATH.read_text(encoding="utf-8"))[0]

    assert record.accession == "KYTJ820101"
    assert record.title == "Hydropathy index"
    assert record.pmid == "PMID:7108955"
    assert record.authors == "Kyte J. and Doolittle R.F."
    assert record.journal and "J. Mol. Biol." in record.journal
    assert record.correlations == {"HOPT810101": pytest.approx(-0.740)}
    assert record.values is not None
    assert record.values["A"] == pytest.approx(1.8)
    assert record.values["L"] == pytest.approx(3.8)
    assert record.values["V"] == pytest.approx(4.2)


def test_parse_hopt810101_multiline_metadata():
    records = parse_aaindex1_text(FIXTURE_PATH.read_text(encoding="utf-8"))
    record = records[1]

    assert record.accession == "HOPT810101"
    assert record.title == "Hydrophilicity value with continued metadata"
    assert record.values is not None
    assert record.values["D"] == pytest.approx(3.0)


def test_parse_missing_values_as_none():
    record = parse_aaindex1_text(FIXTURE_PATH.read_text(encoding="utf-8"))[2]

    assert record.values is not None
    assert record.values["A"] is None
    assert record.values["N"] is None
    assert record.values["C"] == pytest.approx(3.0)


def test_parse_rejects_malformed_entry_cleanly():
    malformed = """H BAD000101
D Bad entry
I    A/L    R/K
     1.0    2.0
//
"""

    with pytest.raises(AaindexParseError, match="malformed I amino acid header"):
        parse_aaindex1_text(malformed)


def test_aaindex_lookup_by_id_packaged_fixture(monkeypatch):
    monkeypatch.delenv(AAINDEX1_ENV_VAR, raising=False)

    output = aaindex_lookup_tool(AaindexLookupInput(index_id="KYTJ820101"))

    assert output.error is None
    assert output.backend_source == "packaged_fixture"
    assert output.backend_path is None
    assert output.record_count >= 2
    assert output.matches
    record = output.matches[0]
    assert record.accession == "KYTJ820101"
    assert record.values is not None
    assert len(record.values) == 20
    assert record.values["A"] == pytest.approx(1.8)


def test_aaindex_default_packaged_fixture_fallback(monkeypatch):
    monkeypatch.delenv(AAINDEX1_ENV_VAR, raising=False)

    output = aaindex_lookup_tool(AaindexLookupInput(query="hydrophilicity"))

    assert output.backend_source == "packaged_fixture"
    assert [record.accession for record in output.matches] == ["HOPT810101"]
    assert output.warnings


def test_aaindex_env_file_valid_path_works(monkeypatch):
    monkeypatch.setenv(AAINDEX1_ENV_VAR, str(FIXTURE_PATH))

    output = aaindex_lookup_tool(AaindexLookupInput(index_id="HOPT810101"))

    assert output.error is None
    assert output.backend_source == "env_file"
    assert output.backend_path == str(FIXTURE_PATH)
    assert output.record_count == 3
    assert output.matches[0].title == "Hydrophilicity value with continued metadata"


def test_aaindex_env_file_missing_path_returns_clean_error(monkeypatch, tmp_path):
    missing_path = tmp_path / "missing-aaindex1"
    monkeypatch.setenv(AAINDEX1_ENV_VAR, str(missing_path))

    output = aaindex_lookup_tool(AaindexLookupInput(index_id="KYTJ820101"))

    assert output.matches == []
    assert output.backend_source == "env_file"
    assert output.backend_path == str(missing_path)
    assert output.record_count == 0
    assert output.error is not None
    assert "Could not read AAindex1 file" in output.error


def test_aaindex_lookup_search_by_exact_accession(monkeypatch):
    monkeypatch.setenv(AAINDEX1_ENV_VAR, str(FIXTURE_PATH))

    output = aaindex_lookup_tool(AaindexLookupInput(index_id="kytj820101"))

    assert [record.accession for record in output.matches] == ["KYTJ820101"]


def test_aaindex_lookup_search_by_text_query(monkeypatch):
    monkeypatch.setenv(AAINDEX1_ENV_VAR, str(FIXTURE_PATH))

    output = aaindex_lookup_tool(AaindexLookupInput(query="antigenic"))

    assert [record.accession for record in output.matches] == ["HOPT810101"]


def test_aaindex_lookup_limit_is_respected(monkeypatch):
    monkeypatch.setenv(AAINDEX1_ENV_VAR, str(FIXTURE_PATH))

    output = aaindex_lookup_tool(AaindexLookupInput(query="PMID", limit=2))

    assert len(output.matches) == 2


def test_aaindex_lookup_requires_id_or_query():
    output = aaindex_lookup_tool(AaindexLookupInput())

    assert output.matches == []
    assert output.error == "Provide index_id or query."


def test_aaindex_lookup_unknown_returns_warning(monkeypatch):
    monkeypatch.delenv(AAINDEX1_ENV_VAR, raising=False)

    output = aaindex_lookup_tool(AaindexLookupInput(index_id="NOPE000000"))

    assert output.matches == []
    assert "No AAindex records matched the lookup." in output.warnings


def test_aaindex_sequence_features_valid_sequence_packaged_fixture(monkeypatch):
    monkeypatch.delenv(AAINDEX1_ENV_VAR, raising=False)

    output = aaindex_sequence_features_tool(
        AaindexSequenceFeaturesInput(text="ACD", index_id="KYTJ820101")
    )

    assert output.error is None
    assert output.index_id == "KYTJ820101"
    assert output.backend_source == "packaged_fixture"
    assert output.sequence == "ACD"
    assert [item.value for item in output.per_residue_values] == [1.8, 2.5, -3.5]
    assert output.summary_stats.mean == pytest.approx((1.8 + 2.5 - 3.5) / 3)
    assert output.summary_stats.min == pytest.approx(-3.5)
    assert output.summary_stats.max == pytest.approx(2.5)
    assert output.missing_residues == []


def test_aaindex_sequence_features_with_local_sample_file(monkeypatch):
    monkeypatch.setenv(AAINDEX1_ENV_VAR, str(FIXTURE_PATH))

    output = aaindex_sequence_features_tool(
        AaindexSequenceFeaturesInput(text="ALV", index_id="KYTJ820101")
    )

    assert output.error is None
    assert output.backend_source == "env_file"
    assert output.record_count == 3
    assert [item.value for item in output.per_residue_values] == [1.8, 3.8, 4.2]
    assert output.summary_stats.mean == pytest.approx((1.8 + 3.8 + 4.2) / 3)


def test_aaindex_sequence_features_handles_missing_per_residue_values(monkeypatch):
    monkeypatch.setenv(AAINDEX1_ENV_VAR, str(FIXTURE_PATH))

    output = aaindex_sequence_features_tool(
        AaindexSequenceFeaturesInput(text="ANC", index_id="MISS000101")
    )

    assert output.error is None
    assert [item.value for item in output.per_residue_values] == [None, None, 3.0]
    assert output.missing_residues == ["A", "N"]
    assert output.summary_stats.mean == pytest.approx(3.0)


def test_aaindex_sequence_features_all_missing_values_warns(monkeypatch):
    monkeypatch.setenv(AAINDEX1_ENV_VAR, str(FIXTURE_PATH))

    output = aaindex_sequence_features_tool(
        AaindexSequenceFeaturesInput(text="AN", index_id="MISS000101")
    )

    assert output.error is None
    assert output.summary_stats.mean is None
    assert output.summary_stats.min is None
    assert output.summary_stats.max is None
    assert "No AAindex values were available for residues in this sequence." in output.warnings


def test_aaindex_sequence_features_rejects_invalid_sequence():
    with pytest.raises(ProteinSequenceError, match="Invalid protein sequence characters"):
        aaindex_sequence_features_tool(
            AaindexSequenceFeaturesInput(text="ACD*", index_id="KYTJ820101")
        )


def test_aaindex_sequence_features_rejects_empty_sequence():
    with pytest.raises(ProteinSequenceError, match="empty"):
        aaindex_sequence_features_tool(
            AaindexSequenceFeaturesInput(text="", index_id="KYTJ820101")
        )


def test_aaindex_sequence_features_rejects_too_long_sequence():
    text = "A" * (MAX_PROTEIN_SEQUENCE_LENGTH + 1)

    with pytest.raises(ProteinSequenceError, match="exceeds maximum"):
        aaindex_sequence_features_tool(
            AaindexSequenceFeaturesInput(text=text, index_id="KYTJ820101")
        )


def test_aaindex_sequence_features_unknown_index():
    output = aaindex_sequence_features_tool(
        AaindexSequenceFeaturesInput(text="ACD", index_id="NOPE000000")
    )

    assert output.error == "AAindex id 'NOPE000000' was not found in the configured backend."
    assert output.per_residue_values == []


def test_aaindex_list_basic_search(monkeypatch):
    monkeypatch.setenv(AAINDEX1_ENV_VAR, str(FIXTURE_PATH))

    output = aaindex_list_tool(AaindexListInput(query="hydro", limit=1))

    assert output.error is None
    assert output.backend_source == "env_file"
    assert output.record_count == 3
    assert len(output.records) == 1
    assert output.records[0].accession == "KYTJ820101"
    assert not hasattr(output.records[0], "values")
