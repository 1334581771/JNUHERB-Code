from pathlib import Path
import re

p = Path(__file__).with_name('ui-source.html')
s = p.read_text(encoding='utf-8')
s = s.replace('先试一个案例', '🧪 Verified example workflows')
s = s.replace('归脾汤 · 正向</button>', '🍵 Gui Pi Decoction</button>')
s = s.replace('肺癌 · 反向</button>', '🫁 Lung cancer</button>')
s = s.replace('失眠 · AI 辅助辨证</button>', '😴 Insomnia · AI</button>')
s = s.replace('</div></section>\n <section id="jh-work"', '</div><p class="small">💡 示例只预填输入与参数，不会自动启动分析。原项目案例作为入口参考，本原型未重新计算验证。</p></section>\n <section id="jh-work"')
icons = {'网络药理学':'🔗','AI 辅助辨证':'💬','配伍优化':'⚙️','分子对接':'🧪'}
for label, icon in icons.items():
    s = s.replace(f'<h3>{label}</h3>', f'<h3>{icon} {label}</h3>')
panel='''<section class="mode-panel" aria-label="Analysis Mode"><h3>Analysis Mode</h3><div class="direction-buttons"><button type="button" data-direction="0" aria-pressed="true">🚀 Forward</button><button type="button" data-direction="1" aria-pressed="false">🎯 Reverse</button></div><p id="jh-direction-help" aria-live="polite">🚀 Forward: TCMs / Formulas / Monomers → Proteins</p><p class="small">💡 先选择分析方向，再填写对应研究对象。参数将在工作台中设置。</p><button type="button" data-start="网络药理学">配置输入 →</button></section>'''
s=s.replace('<hr><div class="eyebrow">RESOURCES</div>', '<hr>'+panel+'<hr><div class="eyebrow">RESOURCES</div>')
s=s.replace('<label for="jh-mode">分析方向</label>', '<label for="jh-mode">Analysis Mode · 分析方向</label>')
s=s.replace('正向：中药 / 方剂 / 成分 → 靶点','🚀 Forward · 中药 / 方剂 / 单体 → 蛋白靶点')
s=s.replace('反向：疾病 / 靶点 → 中药与方剂','🎯 Reverse · 疾病 / 靶点 → 成分 → 中药 → 方剂')
s=s.replace('<label for="jh-type">输入类型</label>', '<p class="small" id="jh-work-direction">🚀 Forward: TCMs / Formulas / Monomers → Proteins</p><label for="jh-type">Query Entity Type · 输入类型</label>')
s=s.replace('<option>分子结构上传（仅示意）</option>', '<option>分子结构上传（仅示意）</option><option>Monomer · 单体名称</option><option>Target proteins · 靶点蛋白</option>')
s=s.replace('<label for="jh-input">研究对象</label>', '<p class="small">💡 中药、方剂和蛋白输入需选择匹配的 ID 类型；单体使用化合物名称。</p><label for="jh-input">研究对象 · one per line</label>')
s=s.replace('名称或 ID 的匹配结果应在正式界面中显示。', 'ℹ️ 每行输入一个名称或 ID。正式分析前应显示名称匹配结果；疾病匹配多个条目时需要选择具体疾病 ID。')
s=s.replace('<summary>高级参数 · 证据阈值与筛选条件</summary><p>在实际实现中保留现有参数、含义和默认值，按当前分析模式展示。</p>', '<summary>🎛️ Filters (optional) · 筛选与阈值说明</summary><p>中药 / 方剂正向分析默认综合评分阈值为 900，单体正向默认为 0。较低阈值保留更多关联，较高阈值侧重高置信度关联。</p><p>疾病反向可设置疾病–靶点证据评分与最大靶点数；这些参数与成分–蛋白的综合评分不是同一量纲。填写疾病时，现有逻辑优先使用疾病解析的靶点。</p><p>此处展示现有规则，不提供计算；正式实现保留全部模式专属参数及帮助。</p>')
s=s.replace('<button type="button" class="primary" id="jh-run">预览分析步骤 →</button>', '<button type="button" class="primary" id="jh-run">🚀 预览分析步骤 →</button>')
extra='''
#jh .mode-panel h3{margin:0 0 12px}#jh .direction-buttons{display:flex;gap:6px;flex-wrap:wrap}#jh .direction-buttons button{font-size:12px;padding:8px}#jh .direction-buttons button[aria-pressed=true]{background:var(--soft);border-color:var(--green);color:var(--green)}#jh .mode-panel p{font-size:12px;overflow-wrap:anywhere}#jh .examples{background:var(--soft);padding:16px;border:1px solid var(--line);border-radius:10px}#jh aside .mode-panel>button{margin-top:8px;font-size:12px}
'''
s=s.replace('</style>',extra+'</style>',1)
js='''
(()=>{const r=document.getElementById('jh');const m=r.querySelector('#jh-mode');const directions=['🚀 Forward: TCMs / Formulas / Monomers → Proteins','🎯 Reverse: Disease / Target Proteins → Chemicals → TCMs → Formulas'];function sync(){const i=m.selectedIndex;r.querySelector('#jh-direction-help').textContent=directions[i];r.querySelector('#jh-work-direction').textContent=directions[i];r.querySelectorAll('[data-direction]').forEach(b=>b.setAttribute('aria-pressed',String(Number(b.dataset.direction)===i)));}r.querySelectorAll('[data-direction]').forEach(b=>b.onclick=()=>{m.selectedIndex=Number(b.dataset.direction);m.dispatchEvent(new Event('change'));});m.addEventListener('change',sync);r.querySelectorAll('[data-example]').forEach(b=>b.addEventListener('click',()=>{if(b.dataset.example==='归脾汤'){r.querySelector('#jh-input').value='HVP5811';}sync();}));sync();})();
'''
s=s.replace('</script>',js+'</script>',1)
p.write_text(s,encoding='utf-8')
