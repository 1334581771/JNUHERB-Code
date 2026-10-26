"""Tests for jnuherb.quota — server-side daily question limit."""

import threading

import pytest

import jnuherb.quota as quota_mod
from jnuherb.quota import consume_quota, quota_status


@pytest.fixture
def usage_db(tmp_path, monkeypatch):
    """Isolate the usage database into a temp file per test."""
    path = str(tmp_path / "usage.db")
    monkeypatch.setenv("JNUHERB_USAGE_DB", path)
    return path


def test_status_starts_full(usage_db):
    s = quota_status()
    assert s["used"] == 0
    assert s["remaining"] == quota_mod.DAILY_QUESTION_LIMIT


def test_status_does_not_consume(usage_db):
    quota_status()
    quota_status()
    assert quota_status()["used"] == 0


def test_consume_up_to_limit_then_deny(usage_db):
    limit = quota_mod.DAILY_QUESTION_LIMIT
    for i in range(limit):
        v = consume_quota()
        assert v["allowed"] is True
        assert v["used"] == i + 1
        assert v["remaining"] == limit - i - 1
    # Pool exhausted: further requests are refused and count does not grow.
    v = consume_quota()
    assert v["allowed"] is False
    assert v["remaining"] == 0
    assert quota_status()["used"] == limit


def test_new_day_resets(usage_db, monkeypatch):
    monkeypatch.setattr(quota_mod, "_beijing_today", lambda: "2026-10-24")
    consume_quota()
    consume_quota()
    assert quota_status()["used"] == 2
    # Midnight in Beijing: a fresh row starts from zero.
    monkeypatch.setattr(quota_mod, "_beijing_today", lambda: "2026-10-25")
    s = quota_status()
    assert s["day"] == "2026-10-25"
    assert s["used"] == 0
    assert s["remaining"] == quota_mod.DAILY_QUESTION_LIMIT


def test_concurrent_consumption_never_exceeds_limit(usage_db):
    """Simultaneous consumers must not push the daily count past the cap."""
    limit = quota_mod.DAILY_QUESTION_LIMIT
    results: list[dict] = []
    lock = threading.Lock()

    def worker():
        v = consume_quota()
        with lock:
            results.append(v)

    threads = [threading.Thread(target=worker) for _ in range(limit + 20)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    allowed = [v for v in results if v["allowed"]]
    denied = [v for v in results if not v["allowed"]]
    assert len(allowed) == limit
    assert len(denied) == 20
    assert quota_status()["used"] == limit


def test_custom_limit(usage_db):
    v = consume_quota(limit=2)
    assert v["allowed"] and v["remaining"] == 1
    v = consume_quota(limit=2)
    assert v["allowed"] and v["remaining"] == 0
    v = consume_quota(limit=2)
    assert not v["allowed"]
