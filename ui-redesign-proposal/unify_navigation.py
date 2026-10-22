from pathlib import Path
import re

p=Path(__file__).with_name('ui-source.html')
s=p.read_text(encoding='utf-8')
nav='''<aside><div class="eyebrow">RESEARCH PLATFORM</div>
<button type="button" class="link active" data-view="home"><span class="nav-icon" aria-hidden="true">🏠</span>研究首页</button>
<button type="button" class="link" data-start="网络药理学"><span class="nav-icon" aria-hidden="true">🔗</span>网络药理学</button>
<button type="button" class="link" data-start="AI 辅助辨证"><span class="nav-icon" aria-hidden="true">💬</span>AI 辅助辨证</button>
<button type="button" class="link" data-start="配伍优化"><span class="nav-icon" aria-hidden="true">⚙️</span>配伍优化</button>
<button type="button" class="link" data-start="分子对接"><span class="nav-icon" aria-hidden="true">🧪</span>分子对接</button>
<hr><div class="eyebrow">RESOURCES</div>
<div class="link resource-item"><span class="nav-icon" aria-hidden="true">📚</span>数据说明与来源</div>
<div class="link resource-item"><span class="nav-icon" aria-hidden="true">📖</span>使用指南 · 引用说明</div>
<hr><div class="small">证候与疾病线索<br>方药、网络与结构分析</div></aside>'''
s=re.sub(r'<aside>[\s\S]*?</aside>',nav,s,count=1)
css='''
#jh aside .link{display:flex;align-items:center;gap:11px;min-height:44px;padding:10px 12px;margin:4px 0;font-size:14px;line-height:1.5;border-radius:7px}
#jh aside .nav-icon{display:inline-flex;align-items:center;justify-content:center;flex:0 0 22px;width:22px;height:24px;font-size:18px;line-height:1}
#jh aside button.link:hover{background:var(--soft)}
#jh aside button.link.active{background:var(--soft);color:var(--green);font-weight:600}
#jh aside .resource-item{font-size:12px;color:var(--muted)}
@media(max-width:760px){#jh aside{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:4px;padding:12px 18px;border-right:0;border-bottom:1px solid var(--line)}#jh aside .eyebrow,#jh aside hr,#jh aside .small,#jh aside .resource-item{display:none}#jh aside .link{margin:0}}
'''
s=s.replace('</style>',css+'</style>',1)
js='''
(()=>{const r=document.getElementById('jh');function updateNav(){const home=!r.querySelector('#jh-home').hidden;const name=r.querySelector('#jh-worktitle').textContent;r.querySelectorAll('aside button.link').forEach(b=>{const selected=b.dataset.view==='home'?home:(!home&&b.dataset.start===name);b.classList.toggle('active',selected);if(selected)b.setAttribute('aria-current','page');else b.removeAttribute('aria-current');});}r.querySelectorAll('[data-view],[data-start],[data-example]').forEach(b=>b.addEventListener('click',updateNav));updateNav();})();
'''
s=s.replace('</script>',js+'</script>',1)
p.write_text(s,encoding='utf-8')
