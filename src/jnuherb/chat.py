"""DeepSeek-powered TCM consultation and recommendation engine.

Uses the DeepSeek API (OpenAI-compatible) to:

1. Analyze symptoms → TCM syndrome differentiation (证型判断).
2. Recommend herbs and formulas based on the syndrome.
3. Cross-reference recommendations with the jnuherb database.

Usage::

    from jnuherb.chat import consult

    result = consult("头痛发热，恶寒无汗，脉浮紧", api_key="sk-...")
    print(result["syndrome"])
    print(result["herbs"])
"""

from __future__ import annotations

import json
import os
import re
from typing import Any, Optional

import requests

# ---------------------------------------------------------------------------
# DeepSeek API configuration
# ---------------------------------------------------------------------------
DEEPSEEK_API_URL = "https://api.deepseek.com/chat/completions"
DEEPSEEK_MODEL = "deepseek-chat"

# ---------------------------------------------------------------------------
# System prompt for TCM syndrome differentiation
# ---------------------------------------------------------------------------
SYSTEM_PROMPT = """你是一位经验丰富的中医专家。用户会描述疾病症状，请你进行中医辨证分析并给出用药建议。

请严格按以下 JSON 格式回复（不要包含其他文字）：
{
  "syndrome": "中医证型名称",
  "analysis": "简要的病机分析（100字以内）",
  "herbs": ["草药中文名1", "草药中文名2", "草药中文名3", ...],
  "formulas": ["方剂名1", "方剂名2", ...],
  "explanation": "用药思路说明（150字以内）"
}

要求：
1. syndrome: 必须是标准的中医证型名（如：风寒束表证、肝郁气滞证等）
2. herbs: 列出 5-10 味最相关的草药（使用中文正名，如柴胡、黄芩等）
3. formulas: 列出 2-5 个经典方剂（使用标准方剂名，如小柴胡汤、麻黄汤等）
4. herbs 和 formulas 中的名称必须是标准名称，以便在中药数据库中检索

请根据用户描述的症状进行辨证分析。"""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def consult(
    symptoms: str,
    api_key: Optional[str] = None,
    model: str = DEEPSEEK_MODEL,
    extra_context: str = "",
) -> dict[str, Any]:
    """Consult the DeepSeek AI for TCM syndrome analysis and recommendations.

    Args:
        symptoms: Patient symptom description (Chinese or English).
        api_key: DeepSeek API key.  If None, reads from env var
            ``DEEPSEEK_API_KEY``.
        model: DeepSeek model name.
        extra_context: Optional additional context (e.g. lab results).

    Returns:
        Dict with keys ``syndrome``, ``analysis``, ``herbs``,
        ``formulas``, ``explanation``, ``raw_response``.

    Raises:
        ValueError: If no API key is available.
        RuntimeError: If the API call fails.
    """
    api_key = _resolve_api_key(api_key)

    user_message = symptoms
    if extra_context:
        user_message = f"{symptoms}\n\n补充信息：{extra_context}"

    response = _call_deepseek(api_key, model, SYSTEM_PROMPT, user_message)

    # Parse JSON from response
    parsed = _parse_response(response)

    parsed["raw_response"] = response
    return parsed


def check_database(
    herbs: list[str],
    formulas: list[str],
) -> dict[str, Any]:
    """Cross-reference AI recommendations against the jnuherb database.

    Args:
        herbs: List of Chinese herb names.
        formulas: List of formula names.

    Returns:
        Dict with keys ``herbs_found``, ``herbs_not_found``,
        ``formulas_found``, ``formulas_not_found``, each containing
        the corresponding database records or names.
    """
    from jnuherb._data import _load_table
    import pandas as pd

    tcm_df = _load_table("jnuherb_tcm")
    fml_df = _load_table("jnuherb_formula")

    # Check herbs by Chinese name
    herbs_found_list: list[dict] = []
    herbs_not_found: list[str] = []
    for herb in herbs:
        match = tcm_df[tcm_df["cn_name"] == herb]
        if len(match) > 0:
            row = match.iloc[0]
            herbs_found_list.append({
                "name": herb,
                "HVMID": row["HVMID"],
                "pinyin": row.get("pinyin_name", ""),
                "en_name": row.get("en_name", ""),
                "in_db": True,
            })
        else:
            herbs_not_found.append(herb)

    # Check formulas by name
    formulas_found_list: list[dict] = []
    formulas_not_found: list[str] = []
    for fml_name in formulas:
        mask = fml_df["name"].str.contains(
            fml_name, na=False, case=False
        )
        if mask.any():
            match = fml_df[mask].iloc[0]
            formulas_found_list.append({
                "name": fml_name,
                "HVPID": match["HVPID"],
                "full_name": match["name"],
                "source": match.get("Source Document", ""),
                "in_db": True,
            })
        else:
            formulas_not_found.append(fml_name)

    return {
        "herbs_found": herbs_found_list,
        "herbs_not_found": herbs_not_found,
        "formulas_found": formulas_found_list,
        "formulas_not_found": formulas_not_found,
    }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _resolve_api_key(api_key: Optional[str]) -> str:
    """Resolve API key from argument, env var, or Streamlit secrets."""
    if api_key:
        return api_key

    # Environment variable
    env_key = os.environ.get("DEEPSEEK_API_KEY")
    if env_key:
        return env_key

    # Try Streamlit secrets
    try:
        import streamlit as st

        st_key = st.secrets.get("DEEPSEEK_API_KEY")
        if st_key:
            return st_key
    except Exception:
        pass

    raise ValueError(
        "DeepSeek API key not found. Provide it as an argument, set the "
        "DEEPSEEK_API_KEY environment variable, or add it to "
        ".streamlit/secrets.toml as DEEPSEEK_API_KEY."
    )


def _call_deepseek(
    api_key: str,
    model: str,
    system_prompt: str,
    user_message: str,
    timeout: int = 60,
) -> str:
    """Call the DeepSeek chat completions API and return the message text."""
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
    }
    payload = {
        "model": model,
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_message},
        ],
        "temperature": 0.3,
        "max_tokens": 2048,
    }

    try:
        resp = requests.post(
            DEEPSEEK_API_URL,
            headers=headers,
            json=payload,
            timeout=timeout,
        )
        resp.raise_for_status()
        data = resp.json()
        return data["choices"][0]["message"]["content"]
    except requests.exceptions.Timeout:
        raise RuntimeError(
            "DeepSeek API request timed out. Please try again."
        )
    except requests.exceptions.HTTPError as e:
        if resp.status_code == 401:
            raise RuntimeError(
                "DeepSeek API authentication failed. Check your API key."
            )
        raise RuntimeError(
            f"DeepSeek API error (HTTP {resp.status_code}): {e}"
        )
    except requests.exceptions.RequestException as e:
        raise RuntimeError(f"DeepSeek API request failed: {e}")
    except (KeyError, IndexError):
        raise RuntimeError(
            f"Unexpected DeepSeek API response format: {resp.text[:200]}"
        )


def _parse_response(text: str) -> dict[str, Any]:
    """Extract structured JSON from the model response.

    Handles cases where the model wraps JSON in ``` fences or
    includes trailing text.
    """
    # Try direct JSON parse first
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass

    # Try extracting from ```json ... ``` fence
    m = re.search(r"```(?:json)?\s*([\s\S]*?)```", text)
    if m:
        try:
            return json.loads(m.group(1))
        except json.JSONDecodeError:
            pass

    # Try finding the first { ... } block
    m = re.search(r"\{[\s\S]*\}", text)
    if m:
        try:
            return json.loads(m.group(0))
        except json.JSONDecodeError:
            pass

    raise RuntimeError(
        f"Could not parse DeepSeek response as JSON. Raw response:\n{text[:500]}"
    )
