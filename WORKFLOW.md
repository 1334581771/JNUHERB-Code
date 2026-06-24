# JNUHERB — 项目开发工作流与对话总结

> **日期**: 2026-06-18 ~ 2026-06-19  
> **目标**: 将 R 包 `jnuherb` 转换为 Python 包 `jnuherb`，并构建可部署的中药网络药理学 Web 分析平台

---

## 一、起点：原始 R 包分析

### 原始包概况

| 属性 | 内容 |
|------|------|
| 包名 | `jnuherb` |
| 语言 | R |
| 定位 | 中药网络药理学分析工具 |
| 数据规模 | 7 个 CSV 文件，~578K 条记录，~20MB |
| 导出函数 | 23 个 |

### 原始 R 包功能模块

| 模块 | 文件 | 核心功能 |
|------|------|----------|
| 数据库查询 | `R/get.R` (171行) | 7 个 `get_*` 函数，查询方剂/草药/化合物/蛋白 |
| 核心分析 | `R/analysis.R` (375行) | `from_tcm_or_formula`（正向）、`from_proteins`（反向）、`dfs_filter` |
| 评分优化 | `R/compute.R` (263行) | `herbiv_score`（概率联合评分）、`component`（背包优化） |
| 网络分析 | `R/network.R` (521行) | 网络过滤、度中心性、箱线图/条形图 |
| 可视化 | `R/output.R` (287行) | Cytoscape CSV 导出、D3.js HTML 网络图 |
| 管道 | `R/herb_analysis.R` (261行) | 一体化集成分析 |
| GUI | `R/launch_app.R` + `inst/shiny/app.R` (683行) | Shiny Dashboard |

### 原始数据库

| 数据集 | 记录数 |
|--------|--------|
| `jnuherb_formula.csv` | ~6,200 |
| `jnuherb_tcm.csv` | ~4,700 |
| `jnuherb_chemicals.csv` | ~13,100 |
| `jnuherb_proteins.csv` | ~14,200 |
| `jnuherb_formula_tcm_links.csv` | ~28,300 |
| `jnuherb_tcm_chemical_links.csv` | ~60,600 |
| `jnuherb_chemical_protein_links.csv` | ~450,800 |

---

## 二、Phase 1-9：R → Python 完整转换

### 关键设计决策

| 决策 | 选择 | 理由 |
|------|------|------|
| 包布局 | `src/` layout | PEP 517 现代标准 |
| 数据访问 | `importlib.resources` | Python 3.9+ 标准库 |
| DataFrame | `pandas` | 直接替代 R `data.frame` |
| 静态绘图 | `matplotlib` + `seaborn` | 替代 base R `boxplot()` 和 `ggplot2` |
| 交互绘图 | `plotly` | 替代 `plotly::ggplotly()` |
| Web GUI | `Streamlit` | 最接近 Shiny 的 Python 方案 |
| 构建 | `setuptools` via `pyproject.toml` | 成熟稳定 |

### R → Python 核心翻译模式

```
R                                    Python
──────────────────────────────────────────────────────────
df[df$col %in% items, ]              df[df["col"].isin(items)]
utils::read.csv(path)                pd.read_csv(path, dtype=str)
utils::write.csv(df, path)           df.to_csv(path, index=False)
sapply(x, fn)                        x.apply(fn) / list comprehension
1 - prod(1 - scores)                1 - np.prod(1 - np.array(scores))
order(df$col, decreasing=TRUE)      df.sort_values("col", ascending=False)
is.null(x)                          x is None
nrow(df) / ncol(df)                 len(df) / df.shape[1]
rowMeans(df[, cols, drop=FALSE])    df[cols].mean(axis=1)
system.file("extdata", ...)         importlib.resources.files(...)
paste0(x, y)                        f"{x}{y}"
```

### 最终包结构

```
jnuherb/
├── pyproject.toml              # 包元数据 + 依赖 + CLI入口
├── README.md                   # 使用文档
├── LICENSE                     # MIT
├── Dockerfile                  # Docker 容器化
├── docker-compose.yml          # Docker Compose 编排
├── .dockerignore
├── .streamlit/config.toml      # Streamlit 服务器配置
├── run.sh                      # 一键启动脚本
├── app_entrypoint.py           # 独立入口（Streamlit Cloud）
├── examples/
│   └── deoxyshikonin_example.py
├── tests/
│   ├── test_query.py           # 查询层测试
│   ├── test_analysis.py        # 分析层测试
│   ├── test_score.py           # 评分层测试
│   └── test_network.py         # 网络层测试
└── src/jnuherb/
    ├── __init__.py             # 公共 API 导出
    ├── __main__.py             # python -m jnuherb CLI
    ├── _data.py                # 懒加载 CSV + 缓存
    ├── query.py                # 7 个 get_* 查询函数
    ├── analysis.py             # from_tcm_or_formula, from_proteins, dfs_filter
    ├── score.py                # herbiv_score, component, knapsack
    ├── network.py              # 7 个网络分析/可视化函数
    ├── output.py               # out_for_cyto, vis_network
    ├── pipeline.py             # herb_analysis 一体化工作流
    ├── app.py                  # launch_app + CLI
    ├── chat.py                 # DeepSeek AI 辨证咨询
    └── data/
        ├── streamlit_app.py    # Streamlit 完整 GUI (9个标签页)
        └── *.csv               # 7 个 CSV 数据文件
```

### 导出函数清单 (23个)

```python
# 查询层 (7)
get_formula, get_formula_tcm_links, get_tcm, get_tcm_chem_links,
get_chemicals, get_chem_protein_links, get_proteins

# 分析层 (3)
from_tcm_or_formula, from_proteins, dfs_filter

# 评分层 (2)
herbiv_score, component

# 网络层 (7)
filter_formula_network, filter_target_network, calculate_degree,
degree_summary, plot_degree_boxplot, plot_degree_barchart,
network_analysis_pipeline

# 输出层 (2)
out_for_cyto, vis_network

# 管道 (1)
herb_analysis

# GUI (1)
launch_app
```

---

## 三、Bug 修复记录

### Bug 1: 显示名 → JNUHERB 大写
- **问题**: 用户要求包显示名为 `JNUHERB`（非 `jnuherb`）
- **解决**: pip 包名保持小写（PEP 8 兼容），所有用户可见文本改为大写

### Bug 2: `score.py` PerformanceWarning
- **问题**: 逐列 `df[col] = ...` 导致 DataFrame 碎片化警告
- **解决**: 重构为 `groupby` 预建索引 + `pd.concat` 批量赋值

### Bug 3: `query.py` 单字符串输入崩溃
- **问题**: `.isin("string")` 报错 `TypeError: only list-like objects are allowed`
- **解决**: 添加 `_ensure_list()` 统一将输入转为 list

### Bug 4: 蛋白 ID 类型不匹配无提示
- **问题**: 用户输入 gene name 但选 `Ensembl_ID` 模式，报"没找到链接"但不说明原因
- **解决**: 在 `from_proteins` 中添加蛋白查找为空的提前检查，给出针对性错误提示

### Bug 5: 示例 Ensembl ID 不在数据库中
- **问题**: Streamlit 界面展示的示例 ID 来自旧 R 包文档，不在当前数据库
- **解决**: 替换为已验证的 ID（`ENSP00000381588` 等），并添加数据库限制说明

### Bug 6: `from_tcm_or_formula` 空列表崩溃
- **问题**: 方剂名搜索返回空时，`tcm_or_formula_id[0]` 触发 `IndexError`
- **解决**: 添加空列表检查，抛出明确错误

### Bug 7: D3.js 网络可视化渲染失败
- **问题**: 大型网络的 D3.js HTML 在 Streamlit iframe 中无法渲染
- **解决**: 换成 Plotly 轻量级 Top-100 核心子网络 + Cytoscape 导出引导

---

## 四、新增功能

### 4.1 正向分析三种查询方式

Streamlit 正向模式支持实体类型选择器：

| 类型 | ID 选项 | 示例 |
|------|--------|------|
| 🌿 TCM | HVMID, cn_name, pinyin, en, latin | `柴胡`, `HVM0367` |
| 📜 Formula | HVPID 或 **方剂名（模糊搜索）** | `右归`, `小柴胡汤` |
| ⚗️ Monomer | 化合物名 | `deoxyshikonin`, `berberine` |

### 4.2 方剂名模糊搜索

- 使用 `str.contains(name, case=False)` 进行子串匹配
- 输入 `右归` → 匹配 `右归饮`
- 输入 `小柴胡` → 匹配 19 个含"小柴胡"的方剂
- 解析后显示确认信息，再执行分析

### 4.3 DeepSeek AI 辨证咨询

```
症状描述 → DeepSeek AI 辨证 → 结构化 JSON 输出
    ↓
自动交叉验证 JNUHERB 数据库
    ├── ✅ 标记在库草药/方剂（显示 HVMID/HVPID）
    └── ⚠️ 标记不在库建议
    ↓
一键按钮 → 将验证过的推荐送入 Forward 分析管道
```

**模块**：
- `src/jnuherb/chat.py`：DeepSeek API 调用 + JSON 解析 + 数据库验证
- `streamlit_app.py` 中的 `💬 AI Consultation` 标签页

**API Key 配置方式**：
1. 侧边栏 `⚙️ API Settings` 输入
2. 环境变量 `DEEPSEEK_API_KEY`
3. `.streamlit/secrets.toml`

---

## 五、部署方式

### 5.1 本地开发

```bash
cd jnuherb
pip install -e .
jnuherb --port 8501
# 或
./run.sh
```

### 5.2 Docker

```bash
docker build -t jnuherb .
docker run -p 8501:8501 jnuherb
# 持久化结果
docker run -p 8501:8501 -v $(pwd)/results:/home/jnuherb/results jnuherb
```

### 5.3 Docker Compose

```bash
docker compose up -d
docker compose logs -f
docker compose down
```

---

## 六、Streamlit 界面全览 (9 个标签页)

| 标签 | 功能 |
|------|------|
| 📋 **Overview** | 分析模式选择、快速指南、结果摘要、分析日志 |
| 💬 **AI Consultation** | DeepSeek 辨证分析、草药/方剂推荐、数据库交叉验证 |
| 🌿 **TCM Results** | 草药结果表、Importance Score 分布图、Top-15 条形图 |
| ⚗️ **Chemical Results** | 化合物结果表、评分分布图、Top-15 条形图 |
| 🧬 **Protein Results** | 靶点蛋白表、度统计摘要、方剂信息（正向模式） |
| 🔗 **Network Viz** | Plotly 核心子网络（Top-100 节点）、Cytoscape 导出引导 |
| 🔍 **Formula Filter** | 反向模式的方剂过滤子网络（边+节点表） |
| 📊 **Degree Stats** | 度分布箱线图、Top-20 节点条形图、完整度表 |
| ⚙️ **Optimization** | TCM/方剂组合优化结果、Boost 图表 |

---

## 七、测试覆盖

**23 个测试全部通过**，覆盖：

- **query**: 7 个查询函数 + 单字符串输入 + 去重 + score 过滤
- **score**: 概率联合评分 + 背包 DP + 组合优化
- **analysis**: 正向分析（TCM/Formula）+ 反向分析（Ensembl/gene_name）+ 错误处理
- **network**: 度计算 + 过滤 + 可视化 + 管道

---

## 八、依赖项

```toml
[project]
requires-python = ">=3.9"
dependencies = [
    "pandas>=1.5.0",
    "numpy>=1.24.0",
    "matplotlib>=3.7.0",
    "seaborn>=0.12.0",
    "plotly>=5.14.0",
    "streamlit>=1.28.0",
]
# DeepSeek AI: 使用 requests（已预装）调用 OpenAI 兼容 API
```

---

## 九、关键技术要点

1. **包名策略**: pip 包 `jnuherb`（小写），显示名 `JNUHERB`（大写）
2. **概率联合评分**: `1 - ∏(1 - sᵢ)` — 非标准 knapsack DP
3. **方剂自动检测**: ID 以 `HVP` 开头→ 方剂模式；否则→ TCM 模式
4. **DFS 过滤**: 4 层嵌套遍历确保只保留完整通路节点
5. **缓存策略**: `@st.cache_data` 加速分析；`importlib.resources` + dict 缓存 CSV
6. **AI JSON 解析**: 支持 clean JSON、```json  fences、嵌入式 JSON 三种格式
