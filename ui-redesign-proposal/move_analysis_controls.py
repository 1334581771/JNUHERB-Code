"""Move task-specific controls into the network pharmacology preview."""
from pathlib import Path
import re

p=Path(__file__).with_name('ui-source.html')
s=p.read_text(encoding='utf-8')
match=re.search(r'<section class="mode-panel"[\s\S]*?</section>',s)
if match:
    panel=match.group(0)
    panel=re.sub(r'<button type="button" data-start="网络药理学">配置输入 →</button>','',panel)
    panel=panel.replace('参数将在工作台中设置。','输入类型、名称和筛选条件在下方设置。')
    s=s[:match.start()]+s[match.end():]
    s=s.replace('<hr><hr>','<hr>')
    s=s.replace('<div class="workspace"><div class="steps">','<div class="workspace">'+panel+'<div class="steps">',1)
    s=s.replace('<label for="jh-mode">Analysis Mode · 分析方向</label>','<label for="jh-mode">分析方向详情</label>')
    s=s.replace('<select id="jh-mode">','<select id="jh-mode" aria-label="分析方向详情">')
    # The buttons already provide mode selection; keep a single visible control.
    s=s.replace('<label for="jh-mode">分析方向详情</label><select id="jh-mode" aria-label="分析方向详情">','<select id="jh-mode" aria-label="分析方向详情" hidden>')
    s=s.replace('<p class="small" id="jh-work-direction">','<p class="small" id="jh-work-direction" hidden>')
    s=s.replace('<h3>Analysis Mode</h3>','<h3>Analysis Mode · 分析方向</h3>')
    s=s.replace('<div class="link"><i data-lucide="network" aria-hidden="true"></i>网络药理学</div>','<button type="button" class="link" data-start="网络药理学">🔗 网络药理学</button>')
    s=s.replace('<div class="link active"><i data-lucide="layout-dashboard" aria-hidden="true"></i>研究导航</div>','<button type="button" class="link active" data-view="home">⌂ 研究首页</button>')
    css='''
#jh .workspace>.mode-panel{background:var(--soft);padding:18px;border-radius:9px;margin-bottom:24px}
#jh .workspace .direction-buttons button{font-size:14px;padding:10px 16px}
#jh .workspace .mode-panel p{margin-bottom:0}
#jh aside button.link{width:100%;text-align:left;border:0;background:transparent}
#jh aside button.link.active{background:var(--soft)}
'''
    s=s.replace('</style>',css+'</style>',1)
    js='''
(()=>{const r=document.getElementById('jh');function section(){const isNetwork=r.querySelector('#jh-worktitle').textContent==='网络药理学';r.querySelector('.mode-panel').hidden=!isNetwork;}r.querySelectorAll('[data-start],[data-example],[data-view]').forEach(b=>b.addEventListener('click',section));section();})();
'''
    s=s.replace('</script>',js+'</script>',1)
    p.write_text(s,encoding='utf-8')
print('Task controls moved from global sidebar to network pharmacology workspace.')
