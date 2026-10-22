# Streamlit 实施指南

## 1. 实施边界与现状

本目录只交付设计和说明。以下路径与函数均基于当前工作区阅读结果，行号可能随修改移动，应按名称定位。现有未提交修改不得覆盖。

当前入口链为 `app_entrypoint.py` → `src/jnuherb/data/streamlit_app.py`，`src/jnuherb/app.py` 负责启动。页面主文件超过三千行，包含主题 CSS、侧栏参数、运行调度、首页、11 个结果标签以及分子对接展示。核心逻辑已分布在 analysis.py、query.py、score.py、network.py、enrichment.py、chat.py、docking.py、output.py 等模块。

`.streamlit/config.toml` 已设置品牌绿；主页面 `_HIDE_STYLES` 已定义横幅和卡片。首页统计从 `src/jnuherb/data/dataset_manifest.json` 读取；应继续复用。网络图已有 D3.js，3D 展示已有 py3Dmol。因此这里的目标是保持 Python 主导的开发方式，而非禁止所有浏览器技术。

此次仅核对 UI 所需结构、依赖和部分相关源码，不宣称已逐个验证所有二进制文件、科学数据或历史设计报告中的结论。

## 2. 视觉规范

| 项目 | 建议值 | 用途 |
|---|---|---|
| 品牌绿 | #17684C | 主操作、选中状态、小面积强调 |
| 深绿 | #124D3B | 首页横幅 |
| 页面底色 | #F5F7F4 | 区分内容和导航 |
| 内容表面 | #FFFFFF | 卡片和表单 |
| 正文 | #213B32 | 主文本 |
| 次级文本 | #65796E | 说明、辅助标签 |
| 细分隔线 | #DFE7E0 | 轻量分组 |
| 字体 | Segoe UI / Microsoft YaHei / sans-serif | 不依赖远程字体 |
| 正文 / 标签 | 14–16px / 12–14px | 正式表单以可读性为先 |
| 主标题 | 约 30–36px | 每页只有一个主标题 |
| 卡片 / 输入圆角 | 10–12px / 7–8px | 统一控件语言 |
| 间距 | 8、12、16、24、32px | 避免每块独立随意设值 |

深色主题需单独定义背景与前景，不能只将白色替换成黑色。HTML 使用 `light-dark()` 自适应系统主题；Streamlit 应使用其主题机制。可将浅色作为第一阶段正式目标，再验收深色。

首页保留两列四模块。降低横幅高度，用短主标题、一段说明和一个主要按钮形成清晰主次；实体规模独立为一行，避免让插画比操作更突出。桌面侧栏约 205–240px；正式页面可使用 Streamlit 原生折叠侧栏，而不是强行固定宽度。内容区域保持适度最大宽度，超宽屏避免段落无限拉长。

原型在 760px 以下隐藏装饰侧栏，顶栏保留首页与工作台入口；正式站不能据此隐藏全部模块导航，需要保留原生侧栏展开按钮。窄屏卡片单列，统计三列换行；输入字体至少 16px，触屏按钮有效高度约 44px。

## 3. 推荐页面结构

首页：品牌横幅 → 六类数据规模 → 四个研究入口 → 示例入口。

开始分析：AI 辅助辨证 / 正向 / 反向。选定后才展示相关输入，避免首页展示完整参数。

探索结果：概览 / 实体（中药、成分、靶点）/ 网络 / 方剂筛选 / 度统计。

下游分析：富集 / 配伍优化 / 分子对接。

资源：数据来源与版本、使用指南、引用说明。若资源内容尚未准备好，先提供真实的简短页面，不保留无响应的假链接。

顶栏可在纯 Streamlit 阶段简化为品牌和当前位置；不必为了与 HTML 完全一致而自定义整套固定顶栏。原生导航与视觉稿允许存在合理差异，关键是布局比例、信息主次和流程一致。

## 4. 设计元素如何对应 Streamlit

| HTML 元素 | Streamlit 实现 | 说明 |
|---|---|---|
| 品牌横幅 | 自有 HTML/CSS 静态展示 + 原生按钮 | 文字可选择；操作不依赖 HTML onclick |
| 六项统计 | columns + metric 或自有静态 HTML | 读取 manifest，不重复硬编码 |
| 两列四卡片 | columns + container + button | 原生按钮绑定真实页面和状态 |
| 模块导航 | st.Page / st.navigation / st.switch_page | 先核对运行版本支持范围 |
| 表单 | st.form + 原生输入 + form_submit_button | 提交后才启动计算 |
| 高级参数 | st.expander | 保留原参数含义、数值与单位 |
| 运行反馈 | st.status / spinner / error | 展示实际执行阶段，不伪造百分比 |
| 表格 | st.dataframe + column_config | 清楚展示中文名、英文名、ID 与评分 |
| 图形 | 现有 Plotly / Matplotlib / D3 / py3Dmol | 统一字号、色板、工具栏位置 |
| 示例 | 现有 pending_prefill 机制 | 载入但不自动运行或调用模型 |

不要把完整 index.html 放到 components.html 中就宣称完成迁移。那样只得到一个隔离的原型，无法自动读取或更新 Python 状态。

## 5. 先核对版本，再选 API

在实际部署环境运行 `python -m pip show streamlit`，保存当前依赖快照。`pyproject.toml` 的 `streamlit>=1.28.0` 是下界，不代表服务器实际版本，也不能保证原生导航、新主题选项或组件 v2 可用。

实施前对照官方文档选择支持的版本，在开发环境验证后固定版本。不要只为装饰效果无约束升级整个科学计算环境。当前页面还导入 fcntl，并使用 Linux 风格分析锁目录；Windows 下直接运行可能失败。UI 开发可用现有容器环境，不应把平台兼容修改悄悄混入 UI 工作。

官方参考：
- https://docs.streamlit.io/develop/concepts/configuration/theming
- https://docs.streamlit.io/develop/api-reference/navigation/st.navigation
- https://docs.streamlit.io/develop/api-reference/execution-flow/st.form
- https://docs.streamlit.io/develop/concepts/custom-components/overview

## 6. 分阶段落地

### 阶段 A：只调整主题和首页

在独立开发分支或副本开始，先记录当前已验证案例的输入、参数和结果。保持启动入口不变，调整 config.toml 的基础主题；更丰富的圆角、字体与侧栏选项须按实际版本选择。

基础配置示意（替换前需与已有配置合并，不覆盖 server 等段）：

```toml
[theme]
primaryColor = "#17684C"
backgroundColor = "#F5F7F4"
secondaryBackgroundColor = "#FFFFFF"
textColor = "#213B32"
font = "sans serif"
```

定位 `_HIDE_STYLES`、`.jnu-hero`、`.jnu-card`。将首页样式移入单独资源文件，使用 `jnu-` 前缀，避免依赖 Streamlit 自动生成的哈希类名。局部静态 HTML 中不要插入未转义的用户输入；必要时用 html.escape。布局不要靠大量负 margin 和固定像素高度实现。

缩短横幅标题，统计移到横幅外；保留 manifest 读取与缺失时的明确提示。四卡每张只展示图标、标题、路径、一句说明和一个按钮。示例入口连接现有 `_example_dialog` / `_apply_pending_prefill`，不能保留原有 `onclick="return false"` 的装饰芯片作为正式操作。

### 阶段 B：导航与文件拆分

建议在包内新增以下结构，避免与启动器 `app.py` 混淆：

```text
src/jnuherb/
  ui/
    state.py
    theme.py
    home.py
    inputs.py
    results.py
    enrichment_view.py
    consultation_view.py
    docking_view.py
    assets/ui.css
  data/streamlit_app.py  # 保留为薄入口
```

这只是建议目录，当前没有创建这些源码。新资源若随 Python 包发布，需要在 pyproject.toml 中补充 package-data，并验证安装后的读取路径；不能只在源码目录中可用。

逐页提取渲染函数；先复用现有 `_run_forward`、`_run_reverse`、`_build_network_results` 与核心分析模块，不改评分逻辑。替换当前 `?tab=` 加 iframe 点击父页面 DOM 的跳转方式，用真实页面注册和导航。不要同时维护两套菜单路由。

如果暂不拆成多页，可先用一个原生导航选择值，只执行选中的渲染函数。不要仅因为将内容放进标签页就假设未选中的内容不会执行。

### 阶段 C：任务状态与表单

定义单一业务状态：当前模块、输入草稿、已提交参数、运行状态、有效结果、run_id 与错误信息。与控件 key 区分开，避免跨页时控件生命周期清理掉长期业务状态。

建议状态序列：idle → running → success / failed。页面显示应依据实际计算结果，而不是只依据按钮被点击。用新的输入指纹或参数快照关联结果；修改表单后明确标出结果来自哪次已提交任务。失败后不能把旧结果误标为当前成功。

把当前参数控件移到对应分析页，使用 form 一次提交。动态输入类型的选择可以放在 form 外，以便先切换表单结构；需要随输入即时变化的疾病解析不要错误放入无法即时回调的表单内部。高成本的疾病解析或模型调用必须有明确触发边界。

示例按原有 pending_prefill 思路，在控件创建前写入业务草稿；不在控件实例化后硬改其 session_state。三个示例分别映射方剂正向、疾病反向和症状输入，保留原有验证参数，不能只替换名称。

保留每会话输出隔离、锁、超时与 run_manifest.json 等现有机制。切换页面不可重新启动任务，载入示例不可自动启动计算。

### 阶段 D：结果界面

概览先展示本次输入、运行时间、参数摘要与真实结果数量。实体表统一排序、ID 列含义和下载命名。表格支持查看长名称，不用窄列省略掉关键 ID。

图形统一品牌绿及稳定的实体分类色；不同实体类别用颜色加文字或形状区分。浏览界面与导出科研图共享数据口径，导出时允许独立的出版字号、白底与分辨率配置。保留已有 PDF / CSV 导出。

富集页面保留 Disease / KEGG / GO 的参数和分组。配伍优化必须绑定已有候选结果；分子对接继续区分“从结果选择”和“上传结构”，其 box 参数和受体样式不能被通用表单取代。AI 页应保留多轮对话、候选及证据展示，不套用中药多行输入框。

## 7. 哪些地方值得加少量其他技术

第一阶段无需新增前端框架。静态横幅与装饰用有限 CSS 即可，输入、按钮、导航都保留原生组件。

若需要网络节点选择与 Python 筛选真正双向联动，再评估组件 v2；为它定义明确的输入、选中节点和事件协议，按版本选择 API，避免组件直接访问父页内部 DOM。当前 iframe 或 D3 图可继续复用，是否升级取决于真实交互需求。

React/Vue 全量界面只有在需要复杂拖拽画布、多窗口工作区、统一路由或独立产品前端时再考虑。它会引入 API、鉴权、状态同步和额外部署；不能把这些成本算成单纯美化。

## 8. 建议实施顺序

先确认首页视觉 → 完成主题与首页 → 替换导航 → 移动输入与状态 → 整理结果页 → 最后评估局部组件。每阶段保持可运行并可回退，不一次性重写主页面与核心算法。

验收目标不是与 HTML 每个像素完全相同，而是首屏能看清研究入口、输入层级明确、跨页状态可靠、科学结果和导出保持一致。
