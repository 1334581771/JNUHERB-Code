"""Tests for the AI consultation module (jnuherb.chat)."""

from __future__ import annotations

import json
from unittest.mock import MagicMock, patch

import pytest

from jnuherb.chat import (
    build_system_prompt,
    consult,
    retrieve_candidates,
    validate_recommendation,
    _extract_symptom_keywords,
    _load_syndrome_profiles,
)


# ---------------------------------------------------------------------------
# Keyword extraction
# ---------------------------------------------------------------------------

def test_extract_symptom_keywords_basic():
    text = "头痛发热，恶寒无汗，脉浮紧"
    kws = _extract_symptom_keywords(text)
    assert "头痛" in kws
    assert "headache" in kws
    assert "脉浮紧" in kws


def test_extract_symptom_keywords_tcm_terms():
    """Newly added TCM-specific terms should be extracted and expanded."""
    text = "往来寒热，胸胁苦满，默默不欲饮食，心烦喜呕，口苦咽干，舌苔薄白，脉弦"
    kws = _extract_symptom_keywords(text)
    assert "往来寒热" in kws
    assert "alternating chills and fever" in kws
    assert "胸胁苦满" in kws
    assert "脉弦" in kws
    assert "wiry pulse" in kws
    assert "口苦" in kws


def test_extract_symptom_keywords_syndrome_patterns():
    """Syndrome pattern names should be extracted."""
    text = "患者有少阳证，表现为寒热往来"
    kws = _extract_symptom_keywords(text)
    assert "少阳证" in kws
    assert "Shaoyang syndrome" in kws


# ---------------------------------------------------------------------------
# Profile loading
# ---------------------------------------------------------------------------

def test_load_syndrome_profiles_contains_classics():
    profiles = _load_syndrome_profiles()
    expected = {
        "小柴胡汤", "右归丸", "四物汤", "独活寄生汤", "六味地黄丸",
        "桃红四物汤", "当归拈痛汤", "芍药甘草汤", "四君子汤", "理中丸", "安宫牛黄丸",
    }
    assert expected.issubset(set(profiles.keys()))


# ---------------------------------------------------------------------------
# Candidate retrieval
# ---------------------------------------------------------------------------

def test_retrieve_candidates_xiao_chai_hu():
    """Classic 小柴胡汤 symptoms should retrieve 小柴胡汤."""
    symptoms = "往来寒热，胸胁苦满，默默不欲饮食，心烦喜呕，口苦咽干，舌苔薄白，脉弦"
    candidates = retrieve_candidates(symptoms, top_k_formulas=10)
    formula_names = [f["name"] for f in candidates["formulas"]]
    assert any("小柴胡汤" in n for n in formula_names)


def test_retrieve_candidates_liu_wei_di_huang():
    """Classic 六味地黄丸 symptoms should retrieve 六味地黄丸."""
    symptoms = "腰膝酸软，头晕耳鸣，盗汗遗精，舌红少苔，脉细数"
    candidates = retrieve_candidates(symptoms, top_k_formulas=10)
    formula_names = [f["name"] for f in candidates["formulas"]]
    assert any("六味地黄丸" in n for n in formula_names)


def test_retrieve_candidates_you_gui_wan_synthetic():
    """右归丸 should now be a REAL DB record after HERB formula enrichment."""
    symptoms = "腰膝酸软，畏寒肢冷，阳痿遗精，舌淡苔白，脉沉细"
    candidates = retrieve_candidates(symptoms, top_k_formulas=10)
    formula_names = [f["name"] for f in candidates["formulas"]]
    assert any("右归丸" in n for n in formula_names)
    you_gui = next(f for f in candidates["formulas"] if f["name"] == "右归丸")
    assert you_gui.get("_synthetic") is False


def test_retrieve_candidates_si_jun_zi():
    """Classic 四君子汤 symptoms should retrieve 四君子汤."""
    symptoms = "食少便溏，神疲乏力，气短懒言，面色萎白，舌淡苔白，脉细弱"
    candidates = retrieve_candidates(symptoms, top_k_formulas=10)
    formula_names = [f["name"] for f in candidates["formulas"]]
    assert any("四君子汤" in n for n in formula_names)


def test_retrieve_candidates_li_zhong():
    """Classic 理中丸 symptoms should retrieve 理中丸."""
    symptoms = "畏寒肢冷，腹部喜温喜按，呕吐清水，下利清谷，舌淡胖苔白滑，脉沉迟"
    candidates = retrieve_candidates(symptoms, top_k_formulas=10)
    formula_names = [f["name"] for f in candidates["formulas"]]
    assert any("理中丸" in n for n in formula_names)


# ---------------------------------------------------------------------------
# Validation layer
# ---------------------------------------------------------------------------

def test_validate_recommendation_match():
    parsed = {
        "syndrome": {"zh": "少阳证", "en": "Shaoyang Syndrome"},
        "formulas": {"zh": ["小柴胡汤"], "en": ["Xiao Chai Hu Tang"]},
    }
    result = validate_recommendation(parsed, {"formulas": []})
    assert result["syndrome_formula_match"]["小柴胡汤"] is True
    assert result["validation_warnings"] == []


def test_validate_recommendation_mismatch():
    parsed = {
        "syndrome": {"zh": "肾阴虚证", "en": "Kidney Yin Deficiency"},
        "formulas": {"zh": ["小柴胡汤"], "en": ["Xiao Chai Hu Tang"]},
    }
    result = validate_recommendation(parsed, {"formulas": []})
    assert result["syndrome_formula_match"]["小柴胡汤"] is False
    assert any("不匹配" in w for w in result["validation_warnings"])


def test_validate_recommendation_synthetic_warning():
    """Missing-DB formulas should trigger a coverage warning."""
    parsed = {
        "syndrome": {"zh": "下焦虚寒证", "en": "Lower Jiao Deficiency Cold"},
        "formulas": {"zh": ["萆薢分清饮"], "en": ["Bixie Fenqing Yin"]},
    }
    candidates = {"formulas": [{"name": "萆薢分清饮", "_synthetic": True}]}
    result = validate_recommendation(parsed, candidates)
    assert result["db_coverage"]["萆薢分清饮"] is False
    assert any("不在数据库中" in w for w in result["validation_warnings"])


# ---------------------------------------------------------------------------
# System prompt builder
# ---------------------------------------------------------------------------

def test_build_system_prompt_includes_profiles():
    candidates = retrieve_candidates(
        "往来寒热，胸胁苦满，脉弦", top_k_formulas=5
    )
    prompt = build_system_prompt(candidates)
    assert "小柴胡汤" in prompt
    assert "少阳证" in prompt
    assert "辨证步骤" in prompt
    assert "示例1" in prompt


# ---------------------------------------------------------------------------
# Consult end-to-end with mocked API
# ---------------------------------------------------------------------------

def test_consult_success_mocked():
    mock_json = {
        "choices": [{
            "message": {
                "content": json.dumps({
                    "syndrome": {"zh": "风寒束表证", "en": "Wind-Cold Exterior Syndrome"},
                    "analysis": {"zh": "分析", "en": "Analysis"},
                    "herbs": {"zh": ["麻黄", "桂枝"], "en": ["Ephedra", "Cinnamon Twig"]},
                    "formulas": {"zh": ["麻黄汤"], "en": ["Ma Huang Tang"]},
                    "explanation": {"zh": "说明", "en": "Explanation"},
                })
            }
        }]
    }

    with patch("jnuherb.chat.requests.post") as mock_post:
        mock_resp = MagicMock()
        mock_resp.status_code = 200
        mock_resp.raise_for_status = MagicMock()
        mock_resp.json = MagicMock(return_value=mock_json)
        mock_post.return_value = mock_resp

        result = consult("头痛发热，恶寒无汗，脉浮紧", api_key="test-key")

    assert result["syndrome"]["zh"] == "风寒束表证"
    assert "麻黄汤" in result["formulas"]["zh"]
    assert "candidates" in result
    assert "validation_warnings" in result


def test_consult_api_timeout():
    with patch("jnuherb.chat.requests.post") as mock_post:
        import requests
        mock_post.side_effect = requests.exceptions.Timeout("timed out")
        with pytest.raises(RuntimeError, match="timed out"):
            consult("symptoms", api_key="test-key")


def test_consult_api_401():
    import requests
    with patch("jnuherb.chat.requests.post") as mock_post:
        mock_resp = MagicMock()
        mock_resp.status_code = 401
        mock_resp.raise_for_status.side_effect = requests.exceptions.HTTPError("Unauthorized")
        mock_post.return_value = mock_resp
        with pytest.raises(RuntimeError, match="authentication failed"):
            consult("symptoms", api_key="test-key")
