"""Server-side daily quota for AI (DeepSeek) questions.

The site has no user accounts, so the quota is a **shared daily pool**:
at most :data:`DAILY_QUESTION_LIMIT` AI consultations per calendar day
(Beijing time) across all visitors.  Enforcement lives in
:func:`jnuherb.chat.consult` *before* any DeepSeek API call, so the limit
cannot be bypassed by editing frontend code.

Storage is a small standalone SQLite database (default
``<tempdir>/jnuherb_usage.db``; override with the ``JNUHERB_USAGE_DB``
environment variable).  It is deliberately separate from the main
jnuherb.db so database rebuilds never touch usage state.  Each calendar
day gets its own row, so resetting at 00:00 Beijing time needs no cron
job — a new day simply starts from zero.  Concurrency safety comes from
the single-statement conditional ``UPDATE`` below: SQLite serializes
writers, and the increment only happens when the count is still below
the limit, so simultaneous requests can never push the total past the
cap.

Usage::

    from jnuherb.quota import consume_quota, quota_status

    verdict = consume_quota()            # atomically consume one question
    if not verdict["allowed"]:
        ...                                # deny before calling DeepSeek
    print(verdict["remaining"])            # e.g. 35
    print(quota_status()["remaining"])     # peek without consuming
"""

from __future__ import annotations

import os
import sqlite3
import tempfile
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Optional

#: Maximum AI consultations per day (shared across all visitors).
DAILY_QUESTION_LIMIT: int = 50

#: Exact message shown (and returned) once the daily pool is exhausted.
#: Deliberately omits the cap number so visitors cannot probe the limit.
QUOTA_EXHAUSTED_MESSAGE: str = (
    "We have reached our daily limit of questions. Please try again tomorrow."
)

# Beijing time.  Computed as UTC+8 directly instead of ZoneInfo so the
# code works on minimal images without the tzdata package (China has had
# no DST since 1991, so a fixed offset is correct).
_BEIJING_TZ = timezone(timedelta(hours=8))

_TABLE = "daily_quota"


def _beijing_today() -> str:
    """Current calendar date in Beijing time, ISO format (YYYY-MM-DD)."""
    return datetime.now(_BEIJING_TZ).date().isoformat()


def _db_path() -> str:
    """Resolve the usage database path.

    ``JNUHERB_USAGE_DB`` wins when set (useful for tests and for pinning
    the file to a Docker volume); otherwise fall back to the system temp
    directory, which is always writable inside the container.
    """
    env = os.environ.get("JNUHERB_USAGE_DB")
    if env:
        return env
    return str(Path(tempfile.gettempdir()) / "jnuherb_usage.db")


def _connect() -> sqlite3.Connection:
    """Open the usage database, creating schema and tuning for concurrency."""
    conn = sqlite3.connect(_db_path(), timeout=10)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA busy_timeout=5000")
    conn.execute(
        f"CREATE TABLE IF NOT EXISTS {_TABLE} ("
        "day TEXT PRIMARY KEY, "
        "count INTEGER NOT NULL DEFAULT 0)"
    )
    return conn


def quota_status(limit: Optional[int] = None) -> dict[str, Any]:
    """Peek at today's usage WITHOUT consuming anything.

    Returns a dict with ``day``, ``used``, ``limit`` and ``remaining``.
    Old days' rows are opportunistically pruned on every read so the
    table stays at one row per active day.  ``limit`` defaults to the
    current :data:`DAILY_QUESTION_LIMIT`, read at call time.
    """
    if limit is None:
        limit = DAILY_QUESTION_LIMIT
    day = _beijing_today()
    try:
        conn = _connect()
        try:
            row = conn.execute(
                f"SELECT count FROM {_TABLE} WHERE day=?", (day,)
            ).fetchone()
            conn.execute(f"DELETE FROM {_TABLE} WHERE day<?", (day,))
            conn.commit()
        finally:
            conn.close()
    except sqlite3.Error:
        # Never let a bookkeeping failure block the chat feature.
        return {"day": day, "used": 0, "limit": limit, "remaining": limit}
    used = int(row[0]) if row else 0
    return {
        "day": day,
        "used": used,
        "limit": limit,
        "remaining": max(0, limit - used),
    }


def consume_quota(limit: Optional[int] = None) -> dict[str, Any]:
    """Atomically consume one question from today's pool.

    The conditional ``UPDATE ... WHERE count < ?`` is the whole trick:
    SQLite applies it under its write lock, so with N simultaneous
    requests at most ``limit`` increments succeed — no check-then-act
    race.  The row for today is created on first use; a fresh day is a
    fresh row, which is what resets the counter at 00:00 Beijing time.
    ``limit`` defaults to the current :data:`DAILY_QUESTION_LIMIT`, read
    at call time.

    Returns a dict with ``allowed`` (bool), ``day``, ``used`` (count
    after this call, unchanged when denied), ``limit`` and ``remaining``.
    On any storage error the request is allowed (fail-open) so a broken
    bookkeeping file never takes the consultation feature down.
    """
    if limit is None:
        limit = DAILY_QUESTION_LIMIT
    day = _beijing_today()
    try:
        conn = _connect()
        try:
            conn.execute(
                f"INSERT INTO {_TABLE}(day, count) VALUES(?, 0) "
                "ON CONFLICT(day) DO NOTHING",
                (day,),
            )
            cur = conn.execute(
                f"UPDATE {_TABLE} SET count=count+1 WHERE day=? AND count<?",
                (day, limit),
            )
            allowed = cur.rowcount == 1
            conn.commit()
            row = conn.execute(
                f"SELECT count FROM {_TABLE} WHERE day=?", (day,)
            ).fetchone()
        finally:
            conn.close()
    except sqlite3.Error:
        return {
            "allowed": True, "day": day, "used": -1,
            "limit": limit, "remaining": limit,
        }
    used = int(row[0]) if row else 0
    return {
        "allowed": allowed,
        "day": day,
        "used": used,
        "limit": limit,
        "remaining": max(0, limit - used),
    }
