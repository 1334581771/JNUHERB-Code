# JNUHERB Streamlit界面优化实施方案

日期：2026年10月21日。

目标：将已讨论的新版首页、导航、示例选择和分析衔接功能落实到现有Streamlit应用，同时整理运行状态、缓存和结果展示。当前工作目录中未发现网站应用入口、依赖配置或业务源码；本文件是可交付给开发者的实施方案，尚未改动或部署网站。界面视觉参考为同目录的JNUHERB_homepage_v2.png。

## 实施范围与对应工具

| 需求 | 推荐实现 | 验收方式 |
|---|---|---|
| 绿色横幅与插画 | 主题配置、列布局、本地图片；少量按组件限定的CSS | 标题和按钮为真实网页元素，图片中无操作文字 |
| 四个模块入口 | 带边框的容器、图片和按钮 | 每个按钮进入正确输入页面或模式 |
| 分组导航 | st.Page和st.navigation；需要精确菜单顺序时使用自定义原生导航 | 所有菜单选项均可访问，切换保留当前结果 |
| Try an example | st.dialog选择示例，写入待载入配置，再导航至目标页面 | 输入与参数正确预填，未自动调用LLM或运行对接 |
| 示例快捷按钮 | 与弹窗共用同一份示例配置及载入函数 | 两条入口的输入和参数完全一致 |
| 统一分析状态 | st.session_state中的业务对象 | 网络、富集、优化等页面使用同一批有效结果 |
| 参数提交 | st.form | 调整参数时不会反复执行昂贵分析 |
| 加载与错误反馈 | st.status、状态信息及异常处理 | 失败后不会显示成功，旧结果不被误标为新结果 |
| 性能 | 按页执行、数据缓存、局部重绘 | 未进入模块不运行对应分析，筛选不重复调用LLM |
| 数据口径 | 数据版本及统计清单作为单一数据来源 | 首页、导出和论文快照计数可对照 |

## 第一步 核对版本和现有结构

在网站实际使用的Python环境执行：

```powershell
python -m pip show streamlit
python -m pip freeze > requirements-before-ui-change.txt
```

记录网站启动命令、入口文件、当前导航方式和Python版本。建立独立开发分支或开发副本，使用开发端口预览。对照该版本的官方文档选用API，完成兼容性验证后再固定依赖；不直接假定设计稿中的全部交互已被旧版Streamlit支持。

现有核心算法通过适配层继续调用。界面优化不应顺带改变概率并集评分、疾病阈值、富集背景或对接准备路线。若需要算法调整，单独记录并重新核对案例结果。

## 第二步 拆分页面与公共组件

下面是建议结构，文件名可与现有项目保持一致：

```text
app.py
ui/
  navigation.py
  home.py
  cards.py
  example_dialog.py
  analysis_status.py
views/
  consultation.py
  forward.py
  reverse.py
  herbs.py
  ingredients.py
  targets.py
  network.py
  formula_filter.py
  degree.py
  enrichment.py
  optimization.py
  docking.py
services/
  backend_adapter.py
  loaders.py
state/
  analysis_state.py
config/
  dataset_manifest.json
  examples.json
assets/
  hero_graph.png
  consultation.png
  network.png
  optimization.png
  docking.png
.streamlit/
  config.toml
```

导航和示例载入应共用页面注册表。避免按钮、导航菜单和快捷案例分别维护不同的路径字符串。

## 第三步 实现首页和图片规范

保留两列四卡片。每张卡片仅包含标题、一句功能说明、一个短提示或示例、操作按钮和一张插画。长说明移入指南或模块帮助。

图片作为独立资产加入assets目录，使用st.image展示。不要把整张界面设计稿当作网页背景，也不要把按钮和数字烘焙到插画内。最终资产应单独制作或提取并逐项检查。

通用Disease统一使用疾病条目卡加医学十字图标，不用肺等特定器官。具体肺癌示例可以使用肺部图示。Formula、Herb、Ingredient、Protein、Disease保持稳定的视觉对应。

使用本地图片，按实际显示尺寸导出，控制文件体积。为说明性图片添加alt，装饰图使用空alt。桌面采用两列卡片，小屏幕检查是否需要单列。重点检查菜单、卡片按钮、长实体名和结果表，而不只是首页截图。

主题通过.streamlit/config.toml配置。优先使用原生列、容器、按钮和主题；必须补充CSS时，限定在自己命名的容器范围内，避免依赖整页元素顺序或框架内部自动生成类名。

## 第四步 导航与页面入口

固定导航结构：

- Overview：主页。
- Start analysis：AI-assisted consultation、Forward analysis、Reverse analysis。
- Explore：Herbs、Ingredients、Targets、Network visualization、Formula filter、Degree statistics。
- Enrichment：Disease enrichment、KEGG pathways、GO terms。
- Optimization：组合优化。
- Docking：分子对接。

优先用st.Page注册页面、st.navigation组织导航。若要求严格保持设计稿中“Overview—三个下拉菜单—Optimization—Docking”的排列，可把默认导航设为hidden，再用st.page_link和st.popover按顺序组合。所有目标页面仍必须事先注册。

首页按钮和菜单调用同一个导航方法。页面切换使用Streamlit导航组件或st.switch_page，不依赖指向内部页面的普通Markdown链接，以保留会话中的分析状态。

正向和反向页面可共用分析组件，但应有明确模式。进入优化或富集页面时，如缺少上游结果，应说明所需输入，并提供返回分析入口；不要只显示空白表格。

## 第五步 示例弹窗与预填

先配置三个示例：

| 示例 | 入口 | 预填内容 |
|---|---|---|
| Syndrome-guided consultation | Consultation | insomnia, dream-disturbed sleep, palpitations, forgetfulness, fatigue, pale tongue with thin coating |
| Gui Pi Decoction | Forward | HVP5811；成分–蛋白综合评分阈值900 |
| Lung cancer | Reverse | MONDO:0008903；关联评分下限0.1；前30靶点；成分–蛋白评分阈值900 |

如示例包含组合优化，再预填random_state=42和1,000次重启。输入草稿与预先保存的结果快照分别管理，不能把文稿数字作为新运行的实时输出。

交互顺序：点击Try an example → 打开st.dialog → 选择Load example → 保存待应用配置 → 退出弹窗并导航 → 目标页面应用预填 → 用户点击Run analysis。

弹窗只负责选择。可将pending_example写入Session State，再调用st.rerun，由主入口处理页面切换。目标页在创建输入控件之前应用预填值，避免在控件实例化后修改同名widget key。

为业务输入和widget使用不同的key。不同页面之间共享的是analysis_draft等业务对象，不单纯依赖临时控件值，避免页面离开后控件清理影响业务状态。显示Example标签，并提供恢复自定义输入的操作。

## 第六步 统一输入 运行与结果状态

建议至少保存以下逻辑字段：

```python
analysis_draft = {
    "mode": "forward",
    "input": "HVP5811",
    "parameters": {"interaction_score_threshold": 900},
    "example_id": "gui_pi",
}

analysis_run = {
    "run_id": None,
    "status": "idle",  # idle / running / success / failed
    "input_snapshot": None,
    "parameters_snapshot": None,
    "dataset_version": None,
    "result": None,
    "error": None,
}
```

这些是建议的数据结构，不是已接入的网站接口。实际字段应由适配层连接现有返回值。

输入变更只修改draft；点击运行后，校验输入，生成run_id和参数快照，再调用现有算法。成功后更新结果。参数变化但没有重新运行时，显示“输入已变更，当前显示上次结果”，并保留上次运行参数。

富集、优化和对接输出记录自己的上游run_id。上游改变后标记下游结果已过期，避免肺癌的网络结果与归脾汤的富集图被混合展示。

在各结果页顶部显示当前输入、模式、阈值、实体数量和运行状态。选择成分进入对接时传递稳定标识符、配体结构和靶点身份，而不只传递界面行号。

## 第七步 控制执行与等待反馈

LLM、网络构建、富集、优化和对接通过明确的运行按钮触发。多个参数放在st.form内统一提交。图表的局部过滤可用st.fragment控制刷新，切换颜色或筛选表格不应重新执行模型调用。

若保留st.tabs，检查是否在默认模式下一次执行了所有标签。当前官方文档提供选中状态跟踪和条件渲染，但需要开发者实际按选中状态限制计算；仅把代码放入标签并不自动解决耗时问题。改成多页面时，也应避免在公共入口提前运行各模块代码。

st.status可显示“解析输入—读取关联—构建网络—生成结果”等真实阶段。没有实际进度信息时使用阶段或加载状态，不伪造百分比。

耗时长的对接任务可逐步改为独立工作进程或任务队列。前端提交job_id并查询状态，输出文件按任务存放。st.status本身不是后台执行队列。

## 第八步 缓存与数据查询

- 公共数据表、GMT解析结果和确定性数据变换，可使用st.cache_data；缓存键包含数据版本及相关参数，并设置合理容量。
- 可安全共享且线程安全的资源可使用st.cache_resource。不要未经核查就让所有用户共用一个可变分析对象或不支持并发的数据库连接。
- 用户当前输入、实体选择和分析结果放在Session State；避免把用户结果写入未区分用户的公共可变变量。
- 共享缓存用于可复用的公共数据，按需查询关系数据。尤其注意大型DataFrame的缓存复制及内存占用，首页不需要读出全部蛋白–疾病关系来显示计数。
- 为靶点、成分及疾病标识符查询设置合适的数据库索引，首页统计读取预生成清单。

## 第九步 数据口径与结果导出

dataset_manifest.json集中保存数据版本、五类实体数、四类关系数及注释覆盖计数。示例数字应注明使用的快照，运行后的实时结果由算法返回。

按当前提供的表格，实体数为12,480、4,666、13,098、14,183、18,962；节点63,389，关系1,824,366。若疾病确为独立表，总记录数为1,887,755、共9张表。正式接入前与实际打包数据核对，不能只根据视觉稿定义后端结构。

每次下载至少包含实体ID、结果数据、输入、阈值、数据/软件版本及适用的随机种子。富集导出保留Overlap、原始P值和调整P值；组合导出保留最佳单项基线、组合分数及Boost定义；对接导出记录受体来源、准备流程、搜索盒、exhaustiveness、seed及真实构象文件。

## 实施顺序与验收

1. 盘点现有源码与Streamlit版本，保存三个案例的基线结果。
2. 完成导航、共享状态和业务适配层，验证跨页切换与旧结果处理。
3. 接入首页四卡片、独立插画、真实按钮和统计清单。
4. 接入示例弹窗、预填及首页快捷示例，验证二者参数一致。
5. 接通富集、优化和对接的上下游数据，补充错误与等待状态。
6. 优化缓存和查询，在桌面及窄屏预览，验证重复运行与多会话隔离。
7. 确认原有案例结果在相同数据和参数下保持一致，再更新部署依赖和发布。

跨页导航、示例预填、参数变更后结果过期、失败状态和会话隔离是需要重点验证的行为。对接只在真实准备路线及运行输出一致后认定完成，不以页面出现图标或评分框作为验收。

## 官方文档

- [页面导航](https://docs.streamlit.io/develop/concepts/multipage-apps/page-and-navigation)
- [st.navigation](https://docs.streamlit.io/develop/api-reference/navigation/st.navigation)
- [st.dialog](https://docs.streamlit.io/develop/api-reference/execution-flow/st.dialog)
- [st.switch_page](https://docs.streamlit.io/develop/api-reference/navigation/st.switch_page)
- [st.popover](https://docs.streamlit.io/develop/api-reference/layout/st.popover)
- [Session State](https://docs.streamlit.io/develop/api-reference/caching-and-state/st.session_state)
- [表单](https://docs.streamlit.io/develop/concepts/architecture/forms)
- [局部更新](https://docs.streamlit.io/develop/concepts/architecture/fragments)
- [标签页执行行为](https://docs.streamlit.io/develop/api-reference/layout/st.tabs)
- [数据缓存](https://docs.streamlit.io/develop/api-reference/caching-and-state/st.cache_data)
- [资源缓存](https://docs.streamlit.io/develop/api-reference/caching-and-state/st.cache_resource)
- [主题](https://docs.streamlit.io/develop/concepts/configuration/theming)
- [图片](https://docs.streamlit.io/develop/api-reference/media/st.image)
