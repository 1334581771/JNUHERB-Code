"""Tests for the protein–disease link layer (CTD / Open Targets).

The full 1.2 M-row links table is not required here: a small synthetic
links frame is patched into ``jnuherb.query._load_table`` so the tests run
before and after the real ingest lands.
"""

from __future__ import annotations

import pandas as pd
import pytest

import jnuherb.query as query
from jnuherb.analysis import resolve_disease_to_targets
from jnuherb.chat import translate_disease_name


@pytest.fixture()
def links_table(monkeypatch):
    """Patch a minimal protein–disease links table into the query layer."""
    rows = [
        # EGFR & co for lung cancer (OT)
        ("ENSP00000275493", "EGFR", "MONDO:0008903", "lung cancer", "OT",
         "direct;ev=100", "0.7181", ""),
        ("ENSP00000288602", "BRAF", "MONDO:0008903", "lung cancer", "OT",
         "direct;ev=50", "0.6000", ""),
        ("ENSP00000256078", "KRAS", "MONDO:0008903", "lung cancer", "OT",
         "direct;ev=40", "0.5500", ""),
        ("ENSP00000269571", "ERBB2", "MONDO:0008903", "lung cancer", "OT",
         "direct;ev=10", "0.0500", ""),
        # same disease, CTD curated (no numeric score)
        ("ENSP00000275493", "EGFR", "MONDO:0008903", "lung cancer", "CTD",
         "therapeutic", "", "12345678"),
        # diabetes (OT)
        ("ENSP00000381976", "SLC30A8", "MONDO:0005148",
         "type 2 diabetes mellitus", "OT", "direct;ev=20", "0.8717", ""),
        ("ENSP00000355255", "PPARG", "MONDO:0005148",
         "type 2 diabetes mellitus", "OT", "direct;ev=9", "0.0800", ""),
        # alzheimer (CTD inference)
        ("ENSP00000344456", "APP", "MESH:D000544", "Alzheimer Disease", "CTD",
         "inference:65", "65.0", "87654321"),
    ]
    df = pd.DataFrame(rows, columns=[
        "Ensembl_ID", "gene_name", "disease_id", "disease_name", "source",
        "evidence_type", "score", "pmids",
    ])
    real = query._load_table
    monkeypatch.setattr(query, "_load_table",
                        lambda name: df if name == "jnuherb_protein_disease_links"
                        else real(name))
    return df


class TestGetProteinDiseaseLinks:
    def test_by_disease_name(self, links_table):
        out = query.get_protein_disease_links(by_disease="lung cancer")
        assert set(out["gene_name"]) == {"EGFR", "BRAF", "KRAS", "ERBB2"}
        assert set(out["source"]) == {"OT", "CTD"}

    def test_by_disease_id(self, links_table):
        out = query.get_protein_disease_links(by_disease="MONDO:0005148")
        assert len(out) == 2
        assert "PPARG" in set(out["gene_name"])

    def test_min_score_filters_ot_but_keeps_ctd(self, links_table):
        out = query.get_protein_disease_links(by_disease="lung cancer",
                                              min_score=0.1)
        # ERBB2 (0.05) drops; the CTD curated row (empty score) survives
        genes = set(out["gene_name"])
        assert "ERBB2" not in genes
        assert "EGFR" in genes
        assert "CTD" in set(out["source"])

    def test_by_protein(self, links_table):
        out = query.get_protein_disease_links(
            by_protein="ENSP00000275493")
        assert {r["disease_id"] for _, r in out.iterrows()} == {"MONDO:0008903"}

    def test_source_filter(self, links_table):
        out = query.get_protein_disease_links(by_disease="lung cancer",
                                              sources=["OT"])
        assert set(out["source"]) == {"OT"}


class TestGetDiseases:
    def test_keyword_search(self, links_table):
        out = query.get_diseases("lung")
        assert len(out) == 1
        assert out.iloc[0]["disease_id"] == "MONDO:0008903"
        assert out.iloc[0]["n_links"] == 4          # unique proteins

    def test_empty_keyword(self, links_table):
        assert query.get_diseases("").empty

    def test_no_match(self, links_table):
        assert query.get_diseases("zzzzz").empty


class TestResolveDiseaseToTargets:
    def test_resolves_and_ranks(self, links_table):
        r = resolve_disease_to_targets("lung cancer", min_score=0.1,
                                       max_targets=10)
        assert r["disease_id"] == "MONDO:0008903"
        assert r["targets"][0] == "ENSP00000275493"   # EGFR, top score
        assert r["n_total"] == 3                       # ERBB2 filtered out

    def test_max_targets_cap(self, links_table):
        r = resolve_disease_to_targets("lung cancer", min_score=0.0,
                                       max_targets=2)
        assert len(r["targets"]) == 2

    def test_unknown_disease(self, links_table):
        r = resolve_disease_to_targets("flibbertigibbet disease")
        assert r["targets"] == [] and r["n_total"] == 0

    def test_ctd_inference_row_resolves(self, links_table):
        r = resolve_disease_to_targets("Alzheimer Disease", min_score=0.0)
        assert r["targets"] == ["ENSP00000344456"]

    def test_english_query_never_calls_translator(self, links_table,
                                                  monkeypatch):
        """English misses must NOT hit the LLM — only CJK queries do."""
        def _boom(text):  # pragma: no cover - must not be called
            raise AssertionError("translator called for English query")

        monkeypatch.setattr("jnuherb.chat.translate_disease_name", _boom)
        r = resolve_disease_to_targets("no such disease here")
        assert r["targets"] == []
        assert r["translated_from"] is None

    def test_chinese_query_translates_then_resolves(self, links_table,
                                                    monkeypatch):
        monkeypatch.setattr("jnuherb.chat.translate_disease_name",
                            lambda text: "Alzheimer Disease")
        r = resolve_disease_to_targets("阿尔茨海默病", min_score=0.0)
        assert r["targets"] == ["ENSP00000344456"]
        assert r["disease_name"] == "Alzheimer Disease"
        assert r["translated_from"] == "阿尔茨海默病"

    def test_chinese_query_translation_failure_is_silent(self, links_table,
                                                         monkeypatch):
        monkeypatch.setattr("jnuherb.chat.translate_disease_name",
                            lambda text: None)
        r = resolve_disease_to_targets("某种不存在的病")
        assert r["targets"] == [] and r["n_total"] == 0
        assert r["translated_from"] is None

    def test_chinese_query_translation_without_match(self, links_table,
                                                     monkeypatch):
        """A translation that matches nothing yields empty, not an error."""
        monkeypatch.setattr("jnuherb.chat.translate_disease_name",
                            lambda text: "zzzzz")
        r = resolve_disease_to_targets("某种不存在的病")
        assert r["targets"] == []
        assert r["translated_from"] is None


class TestTranslateDiseaseName:
    def test_no_cjk_returns_none_without_api(self, monkeypatch):
        def _boom(*a, **kw):  # pragma: no cover - must not be called
            raise AssertionError("API called for non-Chinese input")

        monkeypatch.setattr("jnuherb.chat._resolve_api_key", _boom)
        assert translate_disease_name("hypertension") is None

    def test_translates_and_strips_punctuation(self, monkeypatch):
        monkeypatch.setattr("jnuherb.chat._resolve_api_key",
                            lambda key: "fake")
        monkeypatch.setattr("jnuherb.chat._call_deepseek",
                            lambda *a, **kw: "Hypertension。")
        assert translate_disease_name("高血压") == "Hypertension"

    def test_unknown_input_returns_none(self, monkeypatch):
        monkeypatch.setattr("jnuherb.chat._resolve_api_key",
                            lambda key: "fake")
        monkeypatch.setattr("jnuherb.chat._call_deepseek",
                            lambda *a, **kw: "UNKNOWN")
        assert translate_disease_name("今天天气不错") is None

    def test_api_failure_returns_none(self, monkeypatch):
        monkeypatch.setattr("jnuherb.chat._resolve_api_key",
                            lambda key: "fake")

        def _raise(*a, **kw):
            raise RuntimeError("network down")

        monkeypatch.setattr("jnuherb.chat._call_deepseek", _raise)
        assert translate_disease_name("高血压") is None
