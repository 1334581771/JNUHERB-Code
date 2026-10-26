"""Integration tests: chat.consult enforces the daily quota server-side."""

import pytest

import jnuherb.chat as chat_mod
import jnuherb.quota as quota_mod


@pytest.fixture
def usage_db(tmp_path, monkeypatch):
    path = str(tmp_path / "usage.db")
    monkeypatch.setenv("JNUHERB_USAGE_DB", path)
    # Keep these tests fast regardless of the production limit.
    monkeypatch.setattr(quota_mod, "DAILY_QUESTION_LIMIT", 3)
    return path


def _no_api(*args, **kwargs):
    raise AssertionError("DeepSeek API must not be called in this scenario")


def test_consult_consumes_quota_and_marks_result(usage_db, monkeypatch):
    monkeypatch.setattr(chat_mod, "_call_deepseek", _no_api)
    # Bypass both LLM stages: stage-1 OFFTOPIC ends the pipeline after one
    # stage-1 call would have happened — here we stub the stage entirely.
    monkeypatch.setattr(
        chat_mod, "infer_syndrome",
        lambda *a, **k: {"syndrome": {"zh": "OFFTOPIC", "en": ""}},
    )
    result = chat_mod.consult("写一段快速排序代码", api_key="sk-test-dummy")
    assert result["identity_reply"] == chat_mod.BRAND_REPLY
    assert result["quota"]["used"] == 1
    assert result["quota"]["limit"] == 3
    assert result["quota"]["remaining"] == 2


def test_exhausted_pool_refuses_without_api_call(usage_db, monkeypatch):
    monkeypatch.setattr(chat_mod, "_call_deepseek", _no_api)
    for _ in range(3):
        assert quota_mod.consume_quota()["allowed"]
    result = chat_mod.consult("头痛发热，恶寒无汗")
    assert result["quota_exceeded"] == quota_mod.QUOTA_EXHAUSTED_MESSAGE
    # Still exactly at the cap — the refused request consumed nothing.
    assert quota_mod.quota_status()["used"] == 3


def test_pure_identity_question_does_not_consume(usage_db, monkeypatch):
    monkeypatch.setattr(chat_mod, "_call_deepseek", _no_api)
    result = chat_mod.consult("你是什么模型？")
    assert result["identity_reply"] == chat_mod.BRAND_REPLY
    assert result["quota"] is None
    assert quota_mod.quota_status()["used"] == 0
