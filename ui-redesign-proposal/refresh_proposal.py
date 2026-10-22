"""Refresh only the independent design proposal; never edits the Streamlit app."""
from pathlib import Path
import base64
import json
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
source = HERE / "ui-source.html"
html = source.read_text(encoding="utf-8")
html = html.replace('<i data-lucide="sprout" aria-hidden="true"></i>JNUHERB', '<span class="brand-herb" aria-hidden="true">🌿</span>JNUHERB')
html = html.replace('TRADITIONAL KNOWLEDGE · MOLECULAR INSIGHT', 'JNUHERB · RESEARCH PLATFORM')
html = html.replace('从本草线索，<br>走向分子机制。', 'An integrated platform for<br>traditional Chinese medicine research')
html = html.replace('连接证候、中药、成分、靶点与疾病，<br>为你的研究找到清晰的下一步。', 'AI-assisted syndrome differentiation, bidirectional network pharmacology,<br>formula optimization, and molecular docking in one platform.')
html = html.replace('将传统医药知识<br>与分子机制研究连接起来', '证候与疾病线索<br>方药、网络与结构分析')
html = html.replace('<i class="ornament" data-lucide="sprout" aria-hidden="true"></i>', '')
manifest = json.loads((ROOT / 'src/jnuherb/data/dataset_manifest.json').read_text(encoding='utf-8'))
entities = [
    ('syndrome_terms','tcm-syndrome-meridians.png','Syndromes','标准证候术语'),
    ('disease_terms','disease-medical-record.png','Diseases','疾病术语'),
    ('formulas','herbal-formula-transparent.png','Formulas','方剂'),
    ('herbs','medicinal-herb-plant.png','Herbs','中药'),
    ('ingredients','ingredient-molecule.png','Ingredients','独立成分'),
    ('proteins','protein-folded-chain.png','Proteins','蛋白 ID'),
]
tiles=[]
for key,filename,label,zh in entities:
    data=base64.b64encode((ROOT/'src/jnuherb/data/assets'/filename).read_bytes()).decode('ascii')
    count=manifest['entity_counts'][key]['unique']
    tiles.append(f'<div class="entity"><img src="data:image/png;base64,{data}" alt="{zh}"/><b>{count:,}</b><span class="entity-label">{label}</span><span>{zh}</span></div>')
stats='<div class="stats" aria-label="六类研究实体，统计来自当前 dataset_manifest.json">'+''.join(tiles)+'</div>'
html=re.sub(r'<div class="stats"[\s\S]*?(?=\s*<div class="section-title")',stats,html,count=1)
css='''
/* Revision: retain the original brand herb and all six entity illustrations. */
#jh .brand-herb{font-size:27px;line-height:1}
#jh h1{font-size:clamp(25px,2.5vw,38px);line-height:1.22;max-width:850px;letter-spacing:-.6px}
#jh .hero{padding:28px 30px}
#jh .hero p{max-width:850px;line-height:1.7;font-size:14px}
#jh .stats{background:var(--surface);margin-top:16px;border:1px solid var(--line);border-radius:12px;padding:18px 12px;gap:8px}
#jh .entity{display:flex;flex-direction:column;align-items:center;text-align:center;min-width:0;padding:0 5px}
#jh .entity img{width:100%;height:80px;object-fit:contain;margin-bottom:8px}
#jh .stats .entity-label{font-size:12px;color:var(--ink);margin:2px 0}
#jh .entity+.entity{border-left:1px solid var(--line)}
@media(max-width:760px){#jh .entity:nth-child(4){border-left:0}#jh .entity img{height:65px}#jh .stats{row-gap:18px}#jh .hero{padding:23px}#jh .hero p br{display:none}}
@media(max-width:440px){#jh .entity img{height:56px}#jh .stats{padding:14px 6px}#jh .stats b{font-size:17px}#jh .stats span{font-size:11px}#jh .hero h1 br{display:none}}
'''
if '/* Revision: retain' not in html:
    html=html.replace('</style>',css+'</style>',1)
html=html.replace('数据规模沿用截图 · 概念稿尚未接入数据库','数据规模来自当前统计清单 · 交互演示未接入计算')
source.write_text(html,encoding='utf-8')
standalone='''<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>JNUHERB UI 设计提案 · 原始图标保留版</title><style>:root{color-scheme:light dark}body{margin:0;padding:24px;background:light-dark(#e9eeea,#101813)}#jh{max-width:1440px;margin:auto}#jh .tile i,#jh aside i{display:none}@media(max-width:600px){body{padding:0}#jh{border-radius:0}}</style></head><body>'''+html+'</body></html>'
(HERE/'index.html').write_text(standalone,encoding='utf-8')
print('Updated independent proposal with six original PNG assets and current manifest counts.')
