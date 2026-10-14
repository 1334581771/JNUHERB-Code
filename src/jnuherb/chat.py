"""DeepSeek-powered TCM consultation and recommendation engine.

Uses the DeepSeek API (OpenAI-compatible) to:

1. Analyze symptoms → TCM syndrome differentiation (证型判断).
2. Recommend herbs and formulas based on the syndrome.
3. Cross-reference recommendations with the jnuherb database.

Usage::

    from jnuherb.chat import consult

    result = consult("头痛发热，恶寒无汗，脉浮紧", api_key="sk-...")
    print(result["syndrome"]["zh"])
    print(result["syndrome"]["en"])
    print(result["herbs"]["zh"])
    print(result["herbs"]["en"])
"""

from __future__ import annotations

import json
import os
import re
from typing import Any, Optional

import numpy as np
import pandas as pd
import requests

# ---------------------------------------------------------------------------
# DeepSeek API configuration
# ---------------------------------------------------------------------------
DEEPSEEK_API_URL = "https://api.deepseek.com/chat/completions"
DEEPSEEK_MODEL = "deepseek-chat"

# ---------------------------------------------------------------------------
# Bilingual symptom keyword expansion (Chinese → English synonyms)
# ---------------------------------------------------------------------------
_SYMPTOM_TRANSLATIONS = {
    "头痛": ["headache"],
    "头疼": ["headache"],
    "发热": ["fever"],
    "发烧": ["fever"],
    "恶寒": ["aversion to cold", "chilly", "chill"],
    "畏寒": ["aversion to cold", "chilly", "chill"],
    "无汗": ["no sweating", "anhidrosis", "without sweating"],
    "汗出": ["sweating", "perspiration"],
    "出汗": ["sweating", "perspiration"],
    "脉浮": ["floating pulse"],
    "脉紧": ["tight pulse"],
    "咳嗽": ["cough"],
    "咳": ["cough"],
    "咽痛": ["sore throat", "throat pain"],
    "咽喉痛": ["sore throat", "throat pain"],
    "鼻塞": ["nasal congestion", "stuffy nose"],
    "流涕": ["nasal discharge", "runny nose"],
    "失眠": ["insomnia"],
    "多梦": ["dream-disturbed sleep", "dream disturbed"],
    "心悸": ["palpitations"],
    "心慌": ["palpitations"],
    "健忘": ["forgetfulness"],
    "乏力": ["fatigue", "lassitude"],
    "疲倦": ["fatigue", "lassitude"],
    "舌淡": ["pale tongue"],
    "苔薄": ["thin coating"],
    "胃脘痛": ["epigastric pain", "stomachache"],
    "胃痛": ["epigastric pain", "stomachache"],
    "嗳气": ["belching", "eructation"],
    "食欲不振": ["poor appetite", "anorexia", "loss of appetite"],
    "纳差": ["poor appetite", "anorexia"],
    "大便不畅": ["constipation"],
    "便秘": ["constipation"],
    "腹泻": ["diarrhea"],
    "泄泻": ["diarrhea"],
    "呕吐": ["vomiting"],
    "恶心": ["nausea"],
    "腹胀": ["abdominal distension", "bloating"],
    "腹痛": ["abdominal pain"],
    "黄疸": ["jaundice"],
    "水肿": ["edema"],
    "眩晕": ["dizziness"],
    "头晕": ["dizziness"],
    "口干": ["dry mouth"],
    "口渴": ["thirst"],
    "口苦": ["bitter taste"],
    "烦躁": ["irritability", "restlessness"],
    "胸闷": ["chest oppression", "chest stuffiness"],
    "胸痛": ["chest pain"],
    "气短": ["shortness of breath"],
    "气喘": ["asthma", "dyspnea"],
    "腰痛": ["low back pain", "lumbar pain"],
    "关节痛": ["joint pain", "arthralgia"],
    "痛经": ["dysmenorrhea"],
    "月经不调": ["irregular menstruation"],
    "白带": ["leukorrhea"],
    "遗精": ["seminal emission"],
    "阳痿": ["impotence"],
    "尿频": ["frequent urination"],
    "尿急": ["urgent urination"],
    "尿痛": ["painful urination"],
    "皮疹": ["rash"],
    "瘙痒": ["itching", "pruritus"],
    "发热恶寒": ["fever and aversion to cold"],
    "寒热往来": ["alternating chills and fever"],
    "脉浮紧": ["floating and tight pulse"],
    "脉浮数": ["floating and rapid pulse"],
    "脉细数": ["thin and rapid pulse"],
    "舌淡苔白": ["pale tongue with white coating"],
    "舌红苔黄": ["red tongue with yellow coating"],
    "恶寒发热": ["aversion to cold and fever"],
    "头痛无汗": ["headache without sweating"],
    "鼻塞流涕": ["nasal congestion and discharge"],
    "心悸失眠": ["palpitations and insomnia"],
    "神疲乏力": ["fatigue and lassitude"],
    # TCM pattern / syndrome names
    "少阳证": ["Shaoyang syndrome"],
    "肾阴虚": ["kidney yin deficiency"],
    "肾阳虚": ["kidney yang deficiency"],
    "肝郁": ["liver qi stagnation"],
    "脾虚": ["spleen qi deficiency"],
    "脾胃虚寒": ["spleen stomach deficiency cold"],
    "血虚": ["blood deficiency"],
    "血瘀": ["blood stasis"],
    "湿热": ["damp-heat"],
    "风寒湿痹": ["wind cold damp bi syndrome"],
    # Pulse signs
    "脉弦": ["wiry pulse"],
    "脉细": ["thin pulse"],
    "脉数": ["rapid pulse"],
    "脉沉": ["deep pulse"],
    "脉弱": ["weak pulse"],
    "脉细数": ["thin and rapid pulse"],
    "脉沉细": ["deep and thin pulse"],
    "脉滑数": ["slippery and rapid pulse"],
    "脉涩": ["choppy pulse"],
    "脉弦涩": ["wiry and choppy pulse"],
    "脉沉迟": ["deep and slow pulse"],
    "脉濡弱": ["soft and weak pulse"],
    "脉弦滑": ["wiry and slippery pulse"],
    # Tongue signs
    "舌红": ["red tongue"],
    "舌淡": ["pale tongue"],
    "舌质红": ["red tongue body"],
    "舌质淡": ["pale tongue body"],
    "舌质暗红": ["dark red tongue body"],
    "舌质紫暗": ["dark purple tongue body"],
    "舌淡胖": ["pale swollen tongue"],
    "舌体胖大": ["swollen tongue body"],
    "齿痕": ["teeth marks"],
    "苔薄白": ["thin white coating"],
    "苔黄": ["yellow coating"],
    "苔白": ["white coating"],
    "苔少": ["scanty coating"],
    "少苔": ["scanty coating"],
    "无苔": ["no coating"],
    "苔黄腻": ["yellow greasy coating"],
    "苔白滑": ["white slippery coating"],
    "舌下络脉": ["sublingual veins"],
    # Classic formula key symptoms
    "往来寒热": ["alternating chills and fever"],
    "胸胁苦满": ["chest and hypochondrium fullness"],
    "默默不欲饮食": ["silent anorexia"],
    "心烦喜呕": ["restlessness and nausea"],
    "口苦": ["bitter taste"],
    "咽干": ["dry throat"],
    "目眩": ["dizziness"],
    "腰膝酸软": ["weakness and soreness of waist and knees"],
    "头晕耳鸣": ["dizziness and tinnitus"],
    "盗汗": ["night sweats"],
    "骨蒸潮热": ["bone-steaming tidal fever"],
    "手足心热": ["heat in palms and soles"],
    "畏寒": ["aversion to cold"],
    "腰膝酸冷": ["cold and sore waist and knees"],
    "阳痿": ["impotence"],
    "遗精": ["seminal emission"],
    "小便清长": ["clear long urination"],
    "夜尿频多": ["nocturia"],
    "五更泄泻": ["dawn diarrhea"],
    "大便溏薄": ["loose stools"],
    "便溏": ["loose stools"],
    "面色萎黄": ["sallow complexion"],
    "面色晦暗": ["dark complexion"],
    "唇甲色淡": ["pale lips and nails"],
    "头晕眼花": ["dizziness and blurred vision"],
    "手足发麻": ["numbness of hands and feet"],
    "小腹隐痛": ["dull lower abdominal pain"],
    "月经量少": ["scanty menstruation"],
    "经血紫暗": ["dark purple menstrual blood"],
    "痛经": ["dysmenorrhea"],
    "肌肤甲错": ["rough dry skin"],
    "关节红肿": ["red swollen joints"],
    "关节灼热": ["burning joints"],
    "肢节烦痛": ["vexing pain in limbs and joints"],
    "肩背沉重": ["heaviness in shoulder and back"],
    "遍身酸痛": ["generalized soreness"],
    "小便黄赤": ["dark yellow urine"],
    "大便粘滞": ["sticky stools"],
    "腹中挛急疼痛": ["spasmodic abdominal pain"],
    "小腿转筋": ["calf cramps"],
    "抽筋": ["muscle cramps"],
    "三叉神经痛": ["trigeminal neuralgia"],
    "坐骨神经痛": ["sciatica"],
    "食少": ["poor appetite"],
    "腹胀": ["abdominal distension"],
    "肠鸣": ["borborygmus"],
    "气短懒言": ["shortness of breath and reluctance to speak"],
    "腹部喜温喜按": ["abdominal pain relieved by warmth and pressure"],
    "腹痛绵绵": ["dull persistent abdominal pain"],
    "呕吐清水": ["vomiting clear fluid"],
    "下利清谷": ["diarrhea with undigested food"],
    "吐涎沫": ["drooling frothy saliva"],
    "高热": ["high fever"],
    "神昏": ["loss of consciousness"],
    "谵语": ["delirium"],
    "面红目赤": ["red face and eyes"],
    "牙关紧闭": ["clenched jaw"],
    "惊厥抽搐": ["convulsions and seizures"],
    "喉间痰鸣": ["phlegm rattle in throat"],
    "口渴引饮": ["extreme thirst"],
    "麻木不仁": ["numbness"],
    "屈伸不利": ["limited flexion and extension"],
    "刺痛": ["sharp pain"],
    "痛有定处": ["fixed location pain"],
    "拒按": ["tender to palpation"],
}

# ---------------------------------------------------------------------------
# Syndrome tokenization (zang-fu / pathogenesis / treatment-method tokens)
# ---------------------------------------------------------------------------
_ZANGFU_TOKENS = [
    "心包", "膀胱", "三焦", "小肠", "大肠",
    "心", "肝", "脾", "肺", "肾", "胃", "胆",
    "上焦", "中焦", "下焦",
    "太阳", "阳明", "少阳", "太阴", "少阴", "厥阴",
    "卫分", "气分", "营分", "血分", "半表半里",
]
_BINGJI_TOKENS = [
    "气血两虚", "气血亏虚", "气虚", "血虚", "阴虚", "阳虚", "津亏", "气陷",
    "气滞", "气郁", "血瘀", "血热", "血寒", "出血",
    "痰湿", "痰热", "痰饮", "湿热", "寒湿", "风寒", "风热", "风湿",
    "火热", "实热", "虚热", "虚火", "食积", "水饮", "水肿",
    "阳亢", "肝风", "内风", "腑实", "积滞",
]
_DEFICIENCY_TOKENS = ("气虚", "血虚", "阴虚", "阳虚", "气血两虚", "气血亏虚", "津亏")

# 临床同义的 GB/T 标准词对：双方在标准中各自独立成词、标准未给
# 别名关系，但检索时应当互通（肝郁气滞证 ≈ 肝气郁结证）。不参与
# 归一化（两者都是标准词，不应互相覆盖），只在匹配时给近义加分。
_SYNDROME_NEAR_EQUIV: list[set[str]] = [
    {"肝气郁结证", "肝郁气滞证"},
    {"瘀阻胸胁证", "气滞血瘀证"},
]


def _syndrome_stem(term: str) -> str:
    """去掉结尾"证"的词干：肾阳虚水泛 -> 肾阳虚（父子证型包含匹配用）。"""
    return term[:-1] if term.endswith("证") else term


def _syndrome_near_equiv(a: str, b: str) -> bool:
    return any(a in grp and b in grp for grp in _SYNDROME_NEAR_EQUIV)
_TREATMENT_TOKENS = [
    "益气", "补气", "健脾", "补脾", "养血", "补血", "滋阴", "养阴",
    "温阳", "助阳", "壮阳", "回阳", "清热", "泻火", "解毒", "化痰",
    "祛痰", "燥湿", "利湿", "利水", "化湿", "祛湿", "活血", "化瘀",
    "疏肝", "理气", "解郁", "解表", "发汗", "温里", "散寒", "消食",
    "导滞", "安神", "镇惊", "息风", "开窍", "固涩", "收敛", "祛风湿",
    "通络", "补虚", "和胃", "降逆", "润肠",
]


def _extract_syndrome_tokens(syndrome_zh: str) -> tuple[list[str], list[str]]:
    """Split a TCM syndrome name into zang-fu organ and pathogenesis tokens.

    Returns ``(zangfu_tokens, bingji_tokens)``. When the syndrome carries a
    bare "虚" without a specific deficiency qualifier, a weak "虚" token is
    added so compatibility checks can bridge it to specific deficiencies.
    """
    if not syndrome_zh:
        return [], []
    zangfu = [t for t in _ZANGFU_TOKENS if t in syndrome_zh]
    bingji = [t for t in _BINGJI_TOKENS if t in syndrome_zh]
    if not bingji and "虚" in syndrome_zh:
        bingji = ["虚"]
    return zangfu, bingji


def _syndrome_compatibility(syndrome_zh: str, pattern_zh: str) -> tuple[str, list[str]]:
    """Classify how compatible a formula's typical pattern is with a syndrome.

    Returns ``(tier, shared_tokens)`` where tier is one of:

    - ``"full"``     – exact or parent/child pattern match (✅)
    - ``"partial"``  – shared zang-fu / pathogenesis tokens; usable as a
                       base formula with modifications (🟡)
    - ``"mismatch"`` – no meaningful overlap (🔴)
    - ``"unusable"`` – the profile's pattern field is an effect description
                       (e.g. 气血双补) rather than a syndrome name, so it
                       cannot be auto-validated (ℹ️)
    """
    if not pattern_zh:
        return "unusable", []
    # GB/T 16751.2 canonicalization: official synonyms and curated clinical
    # aliases resolve to the same standard term before any comparison.
    syndrome_zh = canonical_syndrome(syndrome_zh)
    pattern_zh = canonical_syndrome(pattern_zh)
    if pattern_zh in syndrome_zh or syndrome_zh in pattern_zh:
        return "full", []
    # 父子证型：去"证"词干互为子串（肾阳虚 ⊂ 肾阳虚水泛）或临床同义
    # 标准词对（肝郁气滞 ≈ 肝气郁结）均按 full 计
    s_stem, p_stem = _syndrome_stem(syndrome_zh), _syndrome_stem(pattern_zh)
    if (s_stem and p_stem and (s_stem in p_stem or p_stem in s_stem)) \
            or _syndrome_near_equiv(syndrome_zh, pattern_zh):
        return "full", []
    s_zf, s_bj = _extract_syndrome_tokens(syndrome_zh)
    p_zf, p_bj = _extract_syndrome_tokens(pattern_zh)
    if not p_zf and not p_bj:
        return "unusable", []
    shared_zf = [t for t in s_zf if t in p_zf]
    shared_bj = [t for t in s_bj if t in p_bj]
    score = 3 * len(shared_zf) + 4 * len(shared_bj)
    # Deficiency bridge: a bare "虚" on one side vs a specific deficiency
    # (气虚/血虚/...) on the other counts as weak partial overlap.  This is
    # what lets 四君子汤 (脾胃气虚证) pass as a 🟡 partial match for
    # 心脾两虚证 instead of a hard mismatch.
    if not shared_bj and "虚" in s_bj:
        if "虚" in p_bj or any(t in p_bj for t in _DEFICIENCY_TOKENS):
            score += 2
            shared_bj = ["虚"]
    if score >= 4:
        return "partial", shared_zf + shared_bj
    return "mismatch", []


# ---------------------------------------------------------------------------
# System prompt for TCM syndrome differentiation
# ---------------------------------------------------------------------------
SYSTEM_PROMPT_TEMPLATE = """你是一位经验丰富的中医专家，同时具备用英文向国际读者解释中医概念的能力。用户会描述疾病症状，请你进行中医辨证分析并给出用药建议。

**重要约束：你只能从下方“数据库候选”列表中选择草药和方剂进行推荐。** 治疗思路（explanation）必须围绕你实际推荐的、来自候选列表的草药/方剂展开，不能引入列表外未经验证的药物。如果候选列表为空，请在 explanation 中说明当前数据库未检索到直接匹配的症状记录，不要编造药物。

请严格按以下 JSON 格式回复（不要包含其他文字）。所有文本字段必须同时提供中文（zh）和英文（en）：
{
  "syndrome": {
    "zh": "中医证型名称",
    "en": "English translation of the TCM syndrome (e.g., Wind-Cold Exterior Syndrome)"
  },
  "analysis": {
    "zh": "简要的病机分析（100字以内）",
    "en": "Brief pathogenesis analysis in English (within 150 words)"
  },
  "herbs": {
    "zh": ["草药中文名1", "草药中文名2", "草药中文名3", ...],
    "en": ["English herb name 1", "English herb name 2", "English herb name 3", ...]
  },
  "formulas": {
    "zh": ["方剂名1", "方剂名2", ...],
    "en": ["English formula name 1 (pinyin/translation)", "English formula name 2", ...]
  },
  "explanation": {
    "zh": "按以下结构输出中医辨证论治完整回答（严格基于上述候选草药/方剂，不得引入候选列表外的药物）：\n\n【症状采集】\n  - 主诉：[用户主要症状概括]\n  - 补充四诊信息（如有）：[舌象/脉象/问诊要点]\n\n【八纲定性】\n  - 表/里：...\n  - 寒/热：...\n  - 虚/实：...\n  - 阴/阳总属：...\n\n【脏腑/气血定位】\n  - 病变脏腑：...\n  - 气血津液状态：...\n\n【证型判断】\n  - 主证：[标准证型名]\n  - 兼证（如有）：...\n\n【病机分析】\n  - 病因：...\n  - 病机（核心机制）：...\n\n【治则治法】\n  - 治疗原则：...\n  - 治法（八法）：...\n\n【方药参考】\n  - 代表方剂：[从候选列表中选取，2-5个]\n  - 方义简析（君臣佐使）：...\n  - 加减建议：...\n\n【预防调摄】\n  - 生活建议：...\n  - 饮食宜忌：...",
    "en": "Provide a complete TCM syndrome differentiation and treatment response in English, organized by the same sections: Symptom Collection, Eight-Principles Analysis, Zang-Fu/Qi-Blood Localization, Syndrome Pattern, Pathogenesis, Treatment Principle, Formula Reference (only from candidate list), and Prevention/Lifestyle. Keep each section concise and grounded in the candidate herbs/formulas."
  }
}

辨证步骤（必须遵循）：
1. 先判断八纲：阴阳、表里、寒热、虚实。
2. 再归纳四诊：望诊（舌象、面色）、闻诊（气味、声音）、问诊（症状、病史）、切诊（脉象）。
3. 根据辨证结果，从候选方剂中选择最匹配的经典方剂。
4. 如果候选列表中没有匹配方剂，请明确说明，不要编造。

要求：
1. syndrome.zh: 必须是标准的中医证型名（如：风寒束表证、肝郁气滞证等）。
2. syndrome.en: 使用国际读者能理解的英文证型名（可保留拼音并加括号解释）。
3. herbs.zh: 列出 5-10 味最相关的草药（必须从候选列表中选取，使用中文正名）。
4. herbs.en: 对应草药的常用英文名或拉丁名。
5. formulas.zh: 列出 2-5 个经典方剂（必须从候选列表中选取，使用标准方剂名）。
6. formulas.en: 对应方剂的英文翻译或拼音（如 Xiao Chai Hu Tang）。
7. 英文翻译应准确、专业，便于 SCI 论文读者理解。
8. explanation.zh 必须按上述【症状采集】→【预防调摄】的章节结构输出，并引用你实际推荐的具体草药和方剂名称，说明它们如何针对用户症状协同作用，避免泛泛而谈。"""


# ---------------------------------------------------------------------------
# Classic formula syndrome-differentiation profiles
# ---------------------------------------------------------------------------
# 档案 source_classic 字段中的数据库来源码（非典籍出处），用于经典方判定
_DB_SOURCE_CODES = {"ITCM", "ETCM", "HERB", "TCMID", "CPMCP", "BATMAN", "SYMMAP"}


def _is_classic_source(source_classic: str) -> bool:
    """True when the profile's source looks like a classical book citation.

    Herbarium-style database provenance codes (ITCM/ETCM/HERB/CPMCP...) are
    NOT classic citations, while anything else — 《伤寒论》, 景岳全书, pinyin
    book titles like 'shang han lun' — counts as one.
    """
    sc = (source_classic or "").strip()
    if not sc:
        return False
    head = re.split(r"[：:，,。\s]", sc)[0].strip()
    return head.upper() not in _DB_SOURCE_CODES


_PROFILES_CACHE: Optional[dict[str, Any]] = None
_GBT_CACHE: Optional[dict[str, Any]] = None


def _load_gbt_syndromes() -> dict[str, Any]:
    """Load the GB/T 16751.2-2021 syndrome terminology bundle.

    Returns a dict with ``standard_terms`` (set), ``alias_to_standard`` and
    ``clinical_alias_to_standard`` (official synonyms and curated mappings
    of common clinical variants respectively).  Empty/missing fields are
    tolerated so the module works with older data bundles.
    """
    global _GBT_CACHE
    if _GBT_CACHE is not None:
        return _GBT_CACHE

    data: dict[str, Any] = {
        "standard_set": set(),
        "alias_to_standard": {},
        "clinical_alias_to_standard": {},
    }
    try:
        from importlib.resources import files

        path = files("jnuherb.data").joinpath("gbt16751_syndromes.json")
        with open(str(path), "r", encoding="utf-8") as f:
            raw = json.load(f)
        data["standard_set"] = set(raw.get("standard_terms", []))
        data["alias_to_standard"] = raw.get("alias_to_standard", {})
        data["clinical_alias_to_standard"] = raw.get("clinical_alias_to_standard", {})
    except Exception:
        pass

    _GBT_CACHE = data
    return data


def canonical_syndrome(name: str) -> str:
    """Map a syndrome name to its GB/T 16751.2 standard term when known.

    Resolution order: curated clinical-alias mapping (explicit merges such
    as 心脾气血两虚证 -> 心脾两虚证, which shadow the standard set) →
    official synonym → unchanged.  Unknown names pass through as-is.
    """
    name = (name or "").strip()
    if not name:
        return name
    gbt = _load_gbt_syndromes()
    mapped = gbt["clinical_alias_to_standard"].get(name)
    if mapped:
        return mapped
    mapped = gbt["alias_to_standard"].get(name)
    if mapped:
        return mapped
    return name


def _load_syndrome_profiles() -> dict[str, Any]:
    """Load the bundled classic-formula syndrome-differentiation profiles.

    Returns a dict mapping Chinese formula name to its profile dict.
    """
    global _PROFILES_CACHE
    if _PROFILES_CACHE is not None:
        return _PROFILES_CACHE

    try:
        from importlib.resources import files

        path = files("jnuherb.data").joinpath("formula_syndrome_profiles.json")
        with open(str(path), "r", encoding="utf-8") as f:
            profiles = json.load(f)
    except Exception:
        profiles = []

    _PROFILES_CACHE = {p["name_zh"]: p for p in profiles}
    return _PROFILES_CACHE


def _build_profile_context(candidates: dict[str, Any], profiles: dict[str, Any]) -> str:
    """Build in-context formula profiles for candidates that have them."""
    formula_names = {f.get("name", "") for f in candidates.get("formulas", [])}
    lines = ["=== 经典方剂辨证档案 / Classic formula syndrome profiles ==="]
    found = False
    for name in formula_names:
        profile = profiles.get(name)
        if not profile:
            continue
        found = True
        bg = profile.get("ba_gang", {})
        sz = profile.get("si_zhen", {})
        key_syms = profile.get("key_symptoms", [])[:8]
        contras = profile.get("contraindications", [])[:4]
        lines.append(
            f"【{name}】对应证型：{profile.get('syndrome_pattern_zh', '')} / "
            f"{profile.get('syndrome_pattern_en', '')}"
        )
        lines.append(
            f"  八纲：阴阳={bg.get('yin_yang','')}, 表里={bg.get('interior_exterior','')}, "
            f"寒热={bg.get('cold_heat','')}, 虚实={bg.get('deficiency_excess','')}"
        )
        lines.append(f"  望诊：{sz.get('wang','')}")
        lines.append(f"  闻诊：{sz.get('wen','')}")
        lines.append(f"  问诊：{sz.get('wen_zhen','')}")
        lines.append(f"  切诊：{sz.get('qie','')}")
        lines.append(f"  关键症状：{', '.join(key_syms)}")
        lines.append(f"  禁忌：{', '.join(contras)}")
    if not found:
        lines.append("（当前候选方剂暂无辨证档案）")
    return "\n".join(lines)


def _build_few_shot_examples() -> str:
    """Return static few-shot examples for syndrome-formula reasoning."""
    return """=== 辨证示例 / Syndrome differentiation examples ===

示例1：
症状：往来寒热，胸胁苦满，默默不欲饮食，心烦喜呕，口苦咽干，舌苔薄白，脉弦。
八纲：阳证，半表半里，热证（寒热往来偏热），虚实夹杂。
四诊：望诊见舌质偏红、苔薄白；闻诊见默默不欲饮食；问诊见口苦咽干、目眩、胸胁苦满、心烦喜呕；切诊脉弦。
辨证：少阳证。
推荐方剂：小柴胡汤。

示例2：
症状：腰膝酸软，头晕耳鸣，盗汗遗精，舌红少苔，脉细数。
八纲：阴证，里证，虚证，热证（虚火上炎）。
四诊：望诊见颧红、舌红少苔；问诊见腰膝酸软、头晕耳鸣、盗汗遗精、骨蒸潮热、手足心热；切诊脉沉细数。
辨证：肾阴虚证。
推荐方剂：六味地黄丸。

示例3：
症状：食少便溏，神疲乏力，气短懒言，面色萎白，舌淡苔白，脉细弱。
八纲：阴证，里证，虚证，无明显寒热偏性（偏温）。
四诊：望诊见面色萎白、舌淡嫩有齿痕；闻诊见气短懒言；问诊见食少便溏、神疲乏力、腹胀肠鸣；切诊脉细弱。
辨证：脾胃气虚证。
推荐方剂：四君子汤。"""


def build_system_prompt(candidates: dict[str, Any]) -> str:
    """Assemble the full system prompt with profiles and few-shot examples."""
    profiles = _load_syndrome_profiles()
    parts = [
        SYSTEM_PROMPT_TEMPLATE,
        _build_profile_context(candidates, profiles),
        _build_few_shot_examples(),
    ]
    return "\n\n".join(parts)


# ---------------------------------------------------------------------------
# Stage 1: syndrome differentiation prompt
# ---------------------------------------------------------------------------
_STAGE1_PROMPT = """你是一位经验丰富的中医专家。用户会描述患者的疾病症状，请你完成第一阶段的辨证分析。本阶段只辨证、不推荐任何草药或方剂。

请严格按以下 JSON 格式回复（不要包含其他文字）：
{
  "syndrome": {"zh": "标准中医证型名（如：心脾两虚证、风寒束表证）", "en": "English syndrome name"},
  "ba_gang": {"zh": "八纲归纳（如：阴证，里证，虚证，偏寒）", "en": ""},
  "zang_fu": {"zh": "病变脏腑与气血状态（如：心脾两虚，气血不足）", "en": ""},
  "pathogenesis": {"zh": "100字以内病机分析", "en": "Brief English pathogenesis"},
  "treatment_principle": {"zh": "治法，使用规范术语（如：益气健脾，养血安神）", "en": "English treatment principle"},
  "key_symptoms": {"zh": ["支持该证型的关键症状1", "关键症状2", "关键症状3"], "en": []}
}

辨证要求：
1. 先辨八纲（阴阳、表里、寒热、虚实），再定脏腑气血定位；
2. 证型名必须规范、完整（如 心脾两虚证、脾胃气虚证、肝郁气滞证、肾阴虚证）；
3. 治法用规范术语，以四字结构为主（益气健脾/养血安神/清热化痰/疏肝解郁 等）；
4. 本阶段只输出辨证结论，不要推荐方药。"""


def _bilingual_field(parsed: dict[str, Any], key: str) -> dict[str, str]:
    """Normalize a possibly-missing field to a {'zh': ..., 'en': ...} dict."""
    val = parsed.get(key, "")
    if isinstance(val, dict):
        return {"zh": str(val.get("zh", "")), "en": str(val.get("en", ""))}
    if isinstance(val, str):
        return {"zh": val, "en": ""}
    return {"zh": "", "en": ""}


def infer_syndrome(
    symptoms: str,
    api_key: Optional[str] = None,
    model: str = DEEPSEEK_MODEL,
    extra_context: str = "",
) -> dict[str, Any]:
    """Stage-1 LLM call: differentiate the syndrome without recommending formulas.

    Returns a dict with bilingual ``syndrome``, ``ba_gang``, ``zang_fu``,
    ``pathogenesis``, ``treatment_principle`` and a Chinese
    ``key_symptoms`` list.
    """
    api_key = _resolve_api_key(api_key)

    user_message = f"患者症状描述：{symptoms}"
    if extra_context:
        user_message += f"\n\n补充信息：{extra_context}"

    raw = _call_deepseek(api_key, model, _STAGE1_PROMPT, user_message)
    parsed = _parse_response(raw)

    key_symptoms = parsed.get("key_symptoms", {})
    if isinstance(key_symptoms, dict):
        key_symptoms_zh = key_symptoms.get("zh", []) or []
    elif isinstance(key_symptoms, list):
        key_symptoms_zh = key_symptoms
    else:
        key_symptoms_zh = []

    return {
        "syndrome": _bilingual_field(parsed, "syndrome"),
        "ba_gang": _bilingual_field(parsed, "ba_gang"),
        "zang_fu": _bilingual_field(parsed, "zang_fu"),
        "pathogenesis": _bilingual_field(parsed, "pathogenesis"),
        "treatment_principle": _bilingual_field(parsed, "treatment_principle"),
        "key_symptoms": [str(s) for s in key_symptoms_zh],
        "raw_response": raw,
    }


# ---------------------------------------------------------------------------
# Stage 2: syndrome-driven full-database formula scoring
# ---------------------------------------------------------------------------
def _score_formulas_by_syndrome(
    syndrome_zh: str,
    treatment_zh: str = "",
    top_k: int = 10,
) -> list[dict[str, Any]]:
    """Traverse ALL formulas in the database and score them against an
    inferred syndrome pattern (以证索方).

    Scoring per formula:
      1. Exact profile-pattern match with the syndrome      → +12
      2. Zang-fu organ token overlap (心/脾/肾/...)          → +3 each
      3. Pathogenesis token overlap (气虚/阴虚/湿热/...)     → +4 each
      4. Treatment-method token overlap with the formula's
         function/instruction text (益气/健脾/养血/...)      → +2 each
      5. Classic-formula bonus: formulas whose syndrome
         profile cites a classical source (出自经典医籍)     → +6

    Unlike keyword retrieval this channel sees the whole table, so a formula
    whose wording differs from the patient's description can still surface
    when its recorded syndrome matches the differentiation result.  The
    classic bonus keeps canonical decoctions (归脾汤, 四君子汤, ...) from
    being crowded out of the candidate list by same-syndrome proprietary
    preparations scoring identical pattern points.
    """
    from jnuherb._data import _load_table

    if not syndrome_zh:
        return []

    syndrome_zh = canonical_syndrome(syndrome_zh)
    zangfu, bingji = _extract_syndrome_tokens(syndrome_zh)
    treat_tokens = [
        t for t in _TREATMENT_TOKENS
        if t in treatment_zh or t in syndrome_zh
    ]
    if not (zangfu or bingji) and not treat_tokens:
        return []

    profiles = _load_syndrome_profiles()
    fml_df = _load_table("jnuherb_formula")
    if fml_df is None or fml_df.empty:
        return []

    syn_cols = [c for c in ("Syndromes_zh", "Syndromes_en") if c in fml_df.columns]
    func_cols = [
        c for c in ("name", "Instructions", "Treatment Symptoms", "Formula_category")
        if c in fml_df.columns
    ]
    syn_text = (
        fml_df[syn_cols].fillna("").astype(str).agg(" ".join, axis=1)
        if syn_cols else pd.Series([""] * len(fml_df), index=fml_df.index)
    )
    func_text = (
        fml_df[func_cols].fillna("").astype(str).agg(" ".join, axis=1)
        if func_cols else pd.Series([""] * len(fml_df), index=fml_df.index)
    )

    scored: list[dict[str, Any]] = []
    for idx in fml_df.index:
        name = str(fml_df.at[idx, "name"]) if "name" in fml_df.columns else ""
        profile = profiles.get(name, {}) if name else {}
        pattern_zh = canonical_syndrome(profile.get("syndrome_pattern_zh", ""))
        pattern_en = profile.get("syndrome_pattern_en", "")
        # 库方 Syndromes_zh 字段的每个证名单独归一化后一并纳入匹配语料：
        # HERB 来源的该字段常是旧称/异名（食积证），归一化后才能与
        # 推断证型（饮食积滞证）全匹配
        syn_raw = syn_text.at[idx]
        syn_canon_list = [
            canonical_syndrome(p)
            for p in re.split(r"[;；、/,，]", syn_raw)
            if p.strip()
        ]
        syn_canon = " ".join(syn_canon_list)
        haystack = f"{syn_raw} {syn_canon} {pattern_zh} {pattern_en}"

        score = 0
        hits: list[str] = []
        s_stem = _syndrome_stem(syndrome_zh)
        p_stem = _syndrome_stem(pattern_zh)
        if pattern_zh and (
            pattern_zh in syndrome_zh or syndrome_zh in pattern_zh
            or (s_stem and p_stem and (s_stem in p_stem or p_stem in s_stem))
        ):
            score += 12
            hits.append(f"证型全匹配:{pattern_zh}")
            # canonical 后精确相等（而非父子证型词干包含）再加 +1：
            # 肾阳虚证查询时让直指该证的经典（右归丸 pattern=肾阳不足证
            # →canonical 肾阳虚证）排在父子证型经典（金匮肾气丸=肾阳虚水泛证、
            # 真人养脏汤=脾肾阳虚证）之前，破解 39 个同分经典方内卷
            if pattern_zh == syndrome_zh:
                score += 1
                hits.append("直指该证")
        elif syndrome_zh and any(
            syndrome_zh == t or syndrome_zh in t for t in syn_canon_list
        ):
            # 库方 Syndromes_zh 字段（旧称/异名经归一化后）与推断证型全匹配
            score += 12
            hits.append(f"库证型全匹配:{syndrome_zh}")
        elif pattern_zh and _syndrome_near_equiv(syndrome_zh, pattern_zh):
            # 临床同义标准词对（肝郁气滞 ≈ 肝气郁结）：与全匹配同分，
            # 保证柴胡疏肝散/血府逐瘀汤等以主证身份进入候选池而非
            # 被同证型中成药用功效字段 token 淹没
            score += 12
            hits.append(f"近义证型:{pattern_zh}")
        else:
            for t in zangfu:
                if t in haystack:
                    score += 3
                    hits.append(f"脏腑:{t}")
            for t in bingji:
                if t != "虚" and t in haystack:
                    score += 4
                    hits.append(f"病机:{t}")
        for t in treat_tokens:
            if t in func_text.at[idx]:
                score += 2
                hits.append(f"治法:{t}")
        # 经典方加权：辨证档案出处为经典医籍者 +6，避免被同分中成药淹没
        # （说明书式中成药功效字段几乎必含治法词，经典方档案文本稀疏，
        # 加权不足会被同证型中成药系统性挤出候选池）
        if _is_classic_source(profile.get("source_classic", "")):
            score += 6
            hits.append("经典方")

        if score == 0:
            continue
        rec = fml_df.loc[idx].to_dict()
        rec["match_score"] = score
        rec["match_source"] = "syndrome"
        rec["_match_detail"] = hits
        scored.append(rec)

    # Classic-formula fallback: if a profiled classic formula's pattern
    # exactly matches the inferred syndrome but the formula itself is not in
    # the database, inject a synthetic record (flagged for the UI) so the
    # LLM can still consider the canonical prescription.
    by_name = {r.get("name", ""): r for r in scored}
    for name, profile in profiles.items():
        pattern_zh = canonical_syndrome(profile.get("syndrome_pattern_zh", ""))
        if not pattern_zh:
            continue
        f_stem = _syndrome_stem(pattern_zh)
        if not (pattern_zh in syndrome_zh or syndrome_zh in pattern_zh
                or (s_stem and f_stem and (s_stem in f_stem or f_stem in s_stem))):
            continue
        if name in by_name:
            continue
        db_record = _find_formula_in_db(name, profile.get("hvpid") or None)
        if db_record is None:
            db_record = {
                "HVPID": profile.get("hvpid") or "UNKNOWN",
                "name": name,
                "Prescription Composition": "",
                "Instructions": "",
                "Treatment Symptoms": "",
                "Source Document": profile.get("source_classic", ""),
                "_synthetic": True,
            }
        db_record["match_score"] = 12 + 6  # exact pattern + classic bonus
        db_record["match_source"] = "syndrome"
        db_record["_match_detail"] = [f"证型全匹配:{pattern_zh}", "经典方"]
        scored.append(db_record)

    # 同名去重（库中存在大量重名行，如人参散x26、当归补血汤x6），保留最高分；
    # 同分时经典方优先，避免被排序稳定性交给库行序的中成药占位
    best: dict[str, dict[str, Any]] = {}
    for rec in scored:
        key = rec.get("name", "")
        if not key:
            continue
        if key not in best or rec.get("match_score", 0) > best[key].get("match_score", 0):
            best[key] = rec
    scored = list(best.values())
    scored.sort(key=lambda r: (
        -r.get("match_score", 0),
        "经典方" not in r.get("_match_detail", []),
    ))
    return scored[:top_k]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def consult(
    symptoms: str,
    api_key: Optional[str] = None,
    model: str = DEEPSEEK_MODEL,
    extra_context: str = "",
    top_k_herbs: int = 15,
    top_k_formulas: int = 20,
) -> dict[str, Any]:
    """Consult the DeepSeek AI for TCM syndrome analysis and recommendations.

    Two-stage pipeline (two LLM API calls per consultation):

    1. **辨证 (differentiation)** – the model infers the TCM syndrome,
       zang-fu localization and treatment principle from the symptoms alone.
    2. **选方 (formula selection)** – formulas are retrieved by traversing
       the whole database against the inferred syndrome (以证索方) and
       merging with symptom-keyword candidates; the model then picks the
       final herbs/formulas only from this syndrome-grounded candidate
       list, anchored to the stage-1 differentiation result.

    Recommendations are therefore grounded in the jnuherb database while
    staying consistent with the inferred syndrome pattern.

    Args:
        symptoms: Patient symptom description (Chinese or English).
        api_key: DeepSeek API key.  If None, reads from env var
            ``DEEPSEEK_API_KEY``.
        model: DeepSeek model name.
        extra_context: Optional additional context (e.g. lab results).
        top_k_herbs: Number of herb candidates to retrieve from the database.
        top_k_formulas: Number of formula candidates to retrieve from the database.

    Returns:
        Dict with bilingual keys ``syndrome``, ``analysis``, ``herbs``,
        ``formulas``, ``explanation``. Each is either a dict with ``zh``
        and ``en`` strings (for text fields) or a dict with ``zh`` and
        ``en`` lists (for herbs/formulas). Also includes ``raw_response``
        and ``candidates``.

    Raises:
        ValueError: If no API key is available.
        RuntimeError: If the API call fails.
    """
    api_key = _resolve_api_key(api_key)

    # ------------------------------------------------------------------
    # Stage 1: syndrome differentiation (no formula constraint yet)
    # ------------------------------------------------------------------
    stage1 = infer_syndrome(
        symptoms, api_key=api_key, model=model, extra_context=extra_context
    )

    # ------------------------------------------------------------------
    # Stage 2a: retrieve candidates. Herbs come from symptom keywords;
    # formulas are re-ranked by traversing the WHOLE database against the
    # inferred syndrome (以证索方), then merged with the keyword channel.
    # ------------------------------------------------------------------
    candidates = retrieve_candidates(
        symptoms,
        top_k_herbs=top_k_herbs,
        top_k_formulas=top_k_formulas,
    )
    syndrome_zh = stage1["syndrome"]["zh"]
    treatment_zh = stage1["treatment_principle"]["zh"]
    syndrome_formulas = _score_formulas_by_syndrome(
        syndrome_zh, treatment_zh, top_k=top_k_formulas
    )
    candidates["formulas"] = _merge_formula_scores(
        candidates["formulas"], syndrome_formulas
    )[:top_k_formulas]

    rag_context = _build_rag_context(candidates)

    # ------------------------------------------------------------------
    # Stage 2b: final LLM call — pick formulas from the syndrome-matched
    # candidates, grounded in the confirmed stage-1 differentiation.
    # ------------------------------------------------------------------
    stage1_block = (
        "=== 第一阶段辨证结论（已经确定，请在此基础上有据选方）/ Confirmed stage-1 differentiation ===\n"
        f"证型 Syndrome：{syndrome_zh} / {stage1['syndrome']['en']}\n"
        f"八纲：{stage1['ba_gang']['zh']}\n"
        f"脏腑气血：{stage1['zang_fu']['zh']}\n"
        f"病机：{stage1['pathogenesis']['zh']}\n"
        f"治法 Treatment principle：{treatment_zh} / {stage1['treatment_principle']['en']}\n\n"
        "选方要求：请优先从下方候选方剂中选择【对应证型与上述证型一致，或脏腑病机相合（可作基础方加减）】者；"
        "每个候选方剂附带辨证档案/证型字段，请据此核对后再选。"
        "若候选中同时存在经典汤剂（如X汤/X丸出自经典医籍）与功效相近的中成药，请优先推荐经典汤剂，"
        "中成药可作为衍生替代方案列出。"
        "若候选中确实没有对证方剂，请明确说明，不要勉强推荐证型相悖的方剂。"
    )

    system_prompt = "\n\n".join([stage1_block, build_system_prompt(candidates)])
    user_message = f"患者症状描述：{symptoms}\n\n{rag_context}"
    if extra_context:
        user_message += f"\n\n补充信息：{extra_context}"

    response = _call_deepseek(api_key, model, system_prompt, user_message)

    # Parse JSON from response and normalize to bilingual format
    parsed = _parse_response(response)
    parsed = _normalize_bilingual(parsed)

    # Carry the stage-1 differentiation through so the UI and the validator
    # stay consistent with the syndrome the formulas were selected against.
    parsed["syndrome"] = stage1["syndrome"]
    parsed["ba_gang"] = stage1["ba_gang"]
    parsed["zang_fu"] = stage1["zang_fu"]
    parsed["pathogenesis"] = stage1["pathogenesis"]
    parsed["treatment_principle"] = stage1["treatment_principle"]
    parsed["stage1"] = stage1
    parsed = validate_recommendation(parsed, candidates)

    parsed["raw_response"] = response
    parsed["candidates"] = candidates
    return parsed


def retrieve_candidates(
    symptoms: str,
    top_k_herbs: int = 15,
    top_k_formulas: int = 10,
) -> dict[str, Any]:
    """Retrieve candidate herbs and formulas from the jnuherb database.

    Uses a multi-channel scoring strategy:

    1. Keyword overlap against database text fields (classic RAG).
    2. Syndrome-profile keyword boosting for classic formulas.
    3. Synthetic records for classic formulas missing from the DB.

    Args:
        symptoms: Patient symptom description.
        top_k_herbs: Number of top herbs to return.
        top_k_formulas: Number of top formula candidates to return.

    Returns:
        Dict with keys ``herbs`` and ``formulas``. Each is a list of dicts
        containing record fields plus a ``match_score`` and optionally
        ``match_source`` / ``_synthetic``.
    """
    from jnuherb._data import _load_table

    keywords = _extract_symptom_keywords(symptoms)
    if not keywords:
        return {"herbs": [], "formulas": []}

    # Herb fields to search
    tcm_df = _load_table("jnuherb_tcm")
    herb_text_fields = [
        "cn_name", "pinyin_name", "en_name", "Properties", "Meridians",
        "Function", "Indication", "Clinical_manifestations",
        "Therapeutic_en_class", "Therapeutic_cn_class",
    ]
    herb_records = _score_records(tcm_df, herb_text_fields, keywords)

    # Formula channel 1: database keyword overlap
    # Enrichment columns (Syndromes_zh/en, Indications_en, Formula_alias,
    # Formula_category, Formula_type) are added by scripts/ingest_herb_formula.py;
    # _score_records skips columns absent from the table, so this is graceful.
    fml_df = _load_table("jnuherb_formula")
    formula_text_fields = [
        "name", "Prescription Composition", "Instructions",
        "Treatment Symptoms", "Source Document",
        "Syndromes_zh", "Syndromes_en", "Indications_en",
        "Formula_alias", "Formula_category", "Formula_type",
    ]
    formula_records = _score_records(fml_df, formula_text_fields, keywords)
    for r in formula_records:
        r["match_source"] = "database"

    # Formula channel 2: syndrome-profile boosting + missing-formula fallback
    profile_boosts = _score_formulas_by_profile(symptoms, keywords)

    # Merge and re-rank formula candidates
    merged_formulas = _merge_formula_scores(formula_records, profile_boosts)

    return {
        "herbs": herb_records[:top_k_herbs],
        "formulas": merged_formulas[:top_k_formulas],
    }


def _find_formula_in_db(name: str, hvpid: Optional[str] = None) -> Optional[dict[str, Any]]:
    """Find a formula record by name or HVPID, or return a synthetic record.

    Synthetic records are used for classic formulas that are not yet present
    in the bundled database so that the LLM can still recommend them when the
    syndrome pattern matches, while the downstream UI can flag them.
    """
    from jnuherb._data import _load_table

    fml_df = _load_table("jnuherb_formula")
    if hvpid:
        match = fml_df[fml_df["HVPID"] == hvpid]
        if len(match) > 0:
            row = match.iloc[0].to_dict()
            row["_synthetic"] = False
            return row

    mask = fml_df["name"].str.contains(name, na=False, case=False)
    if mask.any():
        row = fml_df[mask].iloc[0].to_dict()
        row["_synthetic"] = False
        return row

    return None


def _score_formulas_by_profile(
    symptoms: str,
    keywords: list[str],
) -> list[dict[str, Any]]:
    """Score classic formulas by overlap with their syndrome profiles.

    Returns database records (or synthetic records for missing formulas)
    with a profile-derived match_score.
    """
    profiles = _load_syndrome_profiles()
    symptoms_lower = symptoms.lower()
    scored: list[dict[str, Any]] = []

    for name, profile in profiles.items():
        score = 0
        sources: list[str] = []

        # Match against key_symptoms (Chinese terms have higher weight)
        for kw in profile.get("key_symptoms", []):
            if len(kw) <= 1:
                continue
            if kw.lower() in symptoms_lower:
                # Chinese terms are usually short and precise; weight them higher
                weight = 3 if re.search(r"[一-鿿]", kw) else 1
                score += weight

        # Direct syndrome-pattern mention is a strong signal
        pattern_zh = profile.get("syndrome_pattern_zh", "")
        if pattern_zh and pattern_zh in symptoms:
            score += 8
            sources.append("syndrome_pattern")

        if score == 0:
            continue

        db_record = _find_formula_in_db(name, profile.get("hvpid") or None)
        if db_record is None:
            db_record = {
                "HVPID": profile.get("hvpid") or "UNKNOWN",
                "name": name,
                "Prescription Composition": "",
                "Instructions": "",
                "Treatment Symptoms": "",
                "Source Document": profile.get("source_classic", ""),
                "_synthetic": True,
            }
        db_record["match_score"] = score
        db_record["match_source"] = "profile"
        scored.append(db_record)

    return scored


def _merge_formula_scores(
    db_records: list[dict[str, Any]],
    profile_records: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Merge database keyword scores with profile-boosted scores.

    If a formula appears in both channels, scores are summed.
    """
    merged: dict[str, dict[str, Any]] = {}

    for r in db_records:
        name = r.get("name", "")
        if not name:
            continue
        merged[name] = dict(r)

    for r in profile_records:
        name = r.get("name", "")
        if not name:
            continue
        if name in merged:
            existing = merged[name]
            existing["match_score"] = existing.get("match_score", 0) + r.get("match_score", 0)
            existing["match_source"] = "database+profile"
            existing["_synthetic"] = existing.get("_synthetic", False) and r.get("_synthetic", False)
        else:
            merged[name] = dict(r)

    # Sort by score descending
    return sorted(merged.values(), key=lambda x: x.get("match_score", 0), reverse=True)


def _extract_symptom_keywords(text: str) -> list[str]:
    """Extract symptom keywords from user input.

    Matches known Chinese/English symptom terms in the raw text, expands
    Chinese symptoms to their English synonyms, and removes common
    stopwords. This ensures bilingual database fields are covered.
    """
    # Keep Chinese characters, English letters, digits, and whitespace
    cleaned = re.sub(r"[^一-龥a-zA-Z0-9\s]", " ", text)
    tokens = [t.strip() for t in cleaned.split() if len(t.strip()) >= 2]

    stopwords = {
        "的", "了", "和", "或", "是", "在", "有", "我", "你", "他", "她", "它",
        "们", "这", "那", "等", "及", "与", "并", "但", "而", "个", "为", "以",
        "可", "能", "会", "要", "不", "无", "如", "例如", "感觉", "出现", "伴有",
        "明显", "加重", "缓解", "患者", "病人", "自述", "主诉", "出现", "大约",
        "最近", "目前", "现在", "已经", "有些", "有时", "经常", "持续",
    }
    tokens = [t for t in tokens if t.lower() not in stopwords]

    # Extract known symptom terms from raw text (handles multi-character
    # symptoms that may be embedded in longer phrases such as "头痛发热")
    expanded = []
    text_lower = text.lower()
    for cn_term, en_terms in _SYMPTOM_TRANSLATIONS.items():
        if cn_term in text:
            expanded.extend([cn_term] + en_terms)

    # Also include any English tokens directly present in the input
    for t in tokens:
        if re.search(r"[a-zA-Z]", t):
            expanded.append(t)

    # Deduplicate while preserving order
    seen = set()
    out = []
    for t in expanded:
        key = t.lower()
        if key not in seen:
            seen.add(key)
            out.append(t)
    return out


def _score_records(
    df,
    text_columns: list[str],
    keywords: list[str],
) -> list[dict[str, Any]]:
    """Score dataframe rows by keyword overlap across text columns."""
    if df.empty or not keywords:
        return []

    available_cols = [c for c in text_columns if c in df.columns]
    if not available_cols:
        return []

    # Concatenate relevant text columns per row
    text_series = df[available_cols].fillna("").astype(str).agg(" ".join, axis=1)
    text_series = text_series.str.lower()

    scores = np.zeros(len(df), dtype=int)
    for kw in keywords:
        kw_lower = kw.lower()
        scores += text_series.str.contains(kw_lower, regex=False, na=False).astype(int)

    if scores.max() == 0:
        return []

    df = df.copy()
    df["match_score"] = scores
    df = df[df["match_score"] > 0].sort_values("match_score", ascending=False)
    return df.to_dict(orient="records")


def _build_rag_context(candidates: dict[str, Any]) -> str:
    """Build a context string from retrieved candidates for the LLM prompt."""
    lines = ["=== 数据库候选（请严格从中选择推荐）/ Database candidates (recommend ONLY from this list) ==="]

    herb_list = candidates.get("herbs", [])
    formula_list = candidates.get("formulas", [])

    if not herb_list and not formula_list:
        lines.append(
            "未在数据库中检索到与用户症状直接匹配的草药或方剂。"
            "No herbs or formulas directly matching the user's symptoms were found in the database."
        )
        return "\n".join(lines)

    lines.append("\n--- 候选草药 / Candidate herbs ---")
    for i, h in enumerate(herb_list, 1):
        cn = h.get("cn_name", "")
        en = h.get("en_name", "")
        func = h.get("Function", "") or ""
        ind = h.get("Indication", "") or ""
        lines.append(
            f"{i}. {cn} ({en}) - Function/功效: {func[:80]}; Indication/主治: {ind[:80]}"
        )

    lines.append("\n--- 候选方剂 / Candidate formulas ---")
    for i, f in enumerate(formula_list, 1):
        name = f.get("name", "")
        symptoms = f.get("Treatment Symptoms", "") or ""
        composition = f.get("Prescription Composition", "") or ""
        synd_zh = f.get("Syndromes_zh", "") or ""
        synd_en = f.get("Syndromes_en", "") or ""
        ind_en = f.get("Indications_en", "") or ""
        synthetic_flag = " [数据库暂未收录 / not in DB]" if f.get("_synthetic", False) else ""
        extra_bits = []
        if synd_zh:
            extra_bits.append(f"Syndrome/证型: {synd_zh[:60]}"
                              + (f" ({synd_en[:60]})" if synd_en else ""))
        if ind_en:
            extra_bits.append(f"Indication: {ind_en[:80]}")
        extra = "; " + "; ".join(extra_bits) if extra_bits else ""
        lines.append(
            f"{i}. {name}{synthetic_flag} - Symptoms/主治: {symptoms[:100]}; Composition/组成: {composition[:80]}{extra}"
        )

    return "\n".join(lines)


def validate_recommendation(
    parsed: dict[str, Any],
    candidates: dict[str, Any],
) -> dict[str, Any]:
    """Post-process and validate AI recommendations against syndrome profiles.

    Checks:
    1. Whether recommended formulas match the inferred syndrome pattern,
       graded as 高度匹配 (✅) / 部分匹配 (🟡) / 不匹配 (🔴).  Profile
       entries whose "pattern" is actually an effect description (e.g.
       气血双补) are reported as 档案待规范 (ℹ️) instead of a mismatch.
    2. Whether any contraindicated pattern appears.
    3. Whether recommended formulas are present in the database.

    Adds keys ``validation_warnings`` (hard problems only),
    ``validation_notes`` (informational ✅/🟡/ℹ️ messages),
    ``syndrome_formula_match`` and ``db_coverage``.
    """
    profiles = _load_syndrome_profiles()

    syndrome_zh = parsed.get("syndrome", {}).get("zh", "")
    syndrome_en = parsed.get("syndrome", {}).get("en", "")
    formulas_zh = parsed.get("formulas", {}).get("zh", [])

    warnings: list[str] = []
    notes: list[str] = []
    match_results: dict[str, bool] = {}
    db_coverage: dict[str, bool] = {}

    for fml_name in formulas_zh:
        profile = profiles.get(fml_name)

        # Database coverage: query the actual database regardless of candidates
        db_record = _find_formula_in_db(fml_name, profile.get("hvpid") if profile else None)
        in_db = db_record is not None and not db_record.get("_synthetic", False)
        db_coverage[fml_name] = in_db

        if not profile:
            if not in_db:
                warnings.append(
                    f"提示：推荐方剂【{fml_name}】暂不在数据库中，请人工复核。"
                )
            continue

        # Syndrome-pattern compatibility (graded, tolerates variants)
        pattern_zh = profile.get("syndrome_pattern_zh", "")
        pattern_en = profile.get("syndrome_pattern_en", "")
        tier, shared = _syndrome_compatibility(syndrome_zh, pattern_zh)
        if tier == "mismatch" and pattern_en and syndrome_en:
            if pattern_en.lower() in syndrome_en.lower():
                tier, shared = "full", []
        match_results[fml_name] = tier in ("full", "partial")

        if tier == "full":
            notes.append(
                f"✅ 推荐方剂【{fml_name}】与推断证型【{syndrome_zh}】高度匹配"
                f"（典型证型：{pattern_zh}）。"
            )
        elif tier == "partial":
            notes.append(
                f"🟡 推荐方剂【{fml_name}】与推断证型部分匹配："
                f"该方剂典型证型为【{pattern_zh}】，与【{syndrome_zh}】在"
                f"【{'、'.join(shared)}】上一致，可作为基础方加减使用，请结合临床复核。"
            )
        elif tier == "unusable":
            notes.append(
                f"ℹ️ 推荐方剂【{fml_name}】的辨证档案“典型证型”字段为"
                f"【{pattern_zh}】，疑似功效描述而非规范证型名，无法自动校验，"
                f"建议人工核对该档案。"
            )
        else:
            warnings.append(
                f"注意：推荐方剂【{fml_name}】与推断证型可能不匹配。"
                f"该方剂典型对应证型为【{pattern_zh}】，"
                f"但当前推断证型为【{syndrome_zh}】，请结合临床复核。"
            )

        # Contraindications
        for contra in profile.get("contraindications", []):
            if contra in syndrome_zh:
                warnings.append(
                    f"警告：推荐方剂【{fml_name}】禁忌证型包含【{contra}】，"
                    f"但当前推断证型出现此特征，请避免使用。"
                )

        if not in_db:
            warnings.append(
                f"提示：推荐方剂【{fml_name}】暂不在数据库中，后续网络药理学分析可能受限。"
            )

    parsed["validation_warnings"] = warnings
    parsed["validation_notes"] = notes
    parsed["syndrome_formula_match"] = match_results
    parsed["db_coverage"] = db_coverage
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
    max_tokens: int = 8192,
) -> str:
    """Call the DeepSeek chat completions API and return the message text.

    ``max_tokens`` defaults to 8192 (the deepseek-chat ceiling): the stage-2
    consultation must emit a long structured differentiation report, and a
    lower cap truncates the JSON mid-stream, which then fails parsing.
    """
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
        "max_tokens": max_tokens,
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


_CJK_RE = re.compile(r"[一-鿿]")


_DISEASE_TRANSLATE_PROMPT = (
    "You are a medical terminology translator for a biomedical database "
    "lookup. Translate the user's Chinese disease name into its standard "
    "English medical term as used in the MONDO, EFO, MeSH, or Orphanet "
    "disease vocabularies. Reply with ONLY the English disease term — no "
    "explanation, no alternatives, no punctuation. If the input is not a "
    "disease name, reply with exactly: UNKNOWN"
)


def translate_disease_name(
    text: str,
    api_key: Optional[str] = None,
    model: str = DEEPSEEK_MODEL,
    timeout: int = 20,
) -> Optional[str]:
    """Translate a Chinese disease name to its standard English term.

    Used by the Reverse-mode disease lookup: the protein–disease link
    table only carries MONDO/EFO/MeSH/Orphanet English vocabulary, so a
    Chinese query (e.g. 高血压) is translated first (e.g. Hypertension)
    and the English term is what actually gets looked up.

    Args:
        text: User-entered disease query.
        api_key: DeepSeek API key (falls back to env / Streamlit secrets).
        model: DeepSeek model id.
        timeout: HTTP timeout in seconds (lookup box needs snappy UX).

    Returns:
        The English disease term, or None when the input contains no
        CJK characters, the translation fails, or the model reports the
        input is not a disease name. Never raises on API errors — a
        failed translation just means "no Chinese fallback".
    """
    if not _CJK_RE.search(text):
        return None
    try:
        key = _resolve_api_key(api_key)
        raw = _call_deepseek(
            key, model, _DISEASE_TRANSLATE_PROMPT, text.strip(),
            timeout=timeout, max_tokens=30,
        )
    except Exception:
        return None
    term = raw.strip().splitlines()[0].strip().strip("\"'。.，,;；")
    if not term or term.upper() == "UNKNOWN":
        return None
    return term


def _normalize_bilingual(parsed: dict[str, Any]) -> dict[str, Any]:
    """Ensure legacy/single-language responses are upgraded to bilingual format.

    Returns a dict where ``syndrome``, ``analysis``, ``explanation`` are
    ``{"zh": ..., "en": ...}`` and ``herbs`` / ``formulas`` are
    ``{"zh": [...], "en": [...]}``.
    """
    out: dict[str, Any] = {}

    # Text fields
    for key in ("syndrome", "analysis", "explanation"):
        val = parsed.get(key, "")
        if isinstance(val, dict) and "zh" in val and "en" in val:
            out[key] = val
        elif isinstance(val, str):
            out[key] = {"zh": val, "en": val}
        else:
            out[key] = {"zh": "", "en": ""}

    # List fields
    for key in ("herbs", "formulas"):
        val = parsed.get(key, [])
        if isinstance(val, dict) and "zh" in val and "en" in val:
            out[key] = val
        elif isinstance(val, list):
            out[key] = {"zh": val, "en": val}
        else:
            out[key] = {"zh": [], "en": []}

    # Preserve any extra keys
    for key, val in parsed.items():
        if key not in out:
            out[key] = val

    return out


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

    hint = ""
    if not text.rstrip().endswith(("}", "```")):
        hint = (
            "\n\n(Hint: the response looks truncated — the consultation "
            "report may exceed the model's max_tokens limit. If this "
            "persists, try a shorter symptom description.)"
        )
    raise RuntimeError(
        f"Could not parse DeepSeek response as JSON. Raw response:\n{text[:500]}{hint}"
    )
