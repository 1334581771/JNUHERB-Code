# JNUHERB 插画统一风格与生成 Prompt

用途：首页 hero 五实体图标 + 四张功能卡插画，论文 Figure 1 与网页共用一套资产。
生成方式：逐张单独生成（每个 prompt 自包含，不依赖上下文），建议模型 `google/gemini-3.1-flash-image-preview`，输出 PNG 后按实际显示尺寸压缩。

## 统一风格基线（每个 prompt 已内置，勿删改）

> flat vector line icon, uniform 2px stroke, rounded line caps and corners,
> single deep-green outline color #16634a, soft mint-green fill #e6f0ea,
> minimal geometric shapes, no gradients, no shadows, no text, no letters,
> centered composition with generous padding, isolated on solid white background,
> 1:1 square, clean scientific-platform aesthetic, consistent icon-family style

色彩令牌（与 `.streamlit/config.toml` 一致）：
- 主线条深绿 `#16634a`；软填充 `#e6f0ea`；强调金 `#aa7841`（仅允许在「组合优化」一张中作点缀）；辅助蓝 `#477b98`（仅蛋白一张可用）。

## 一、Hero 五实体图标（首页绿色横幅内并列展示）

### 1. 方剂 Formula — 药碗/汤剂
```
Flat vector line icon of a traditional Chinese herbal decoction bowl: a wide
ceramic bowl with three gentle steam wavy lines rising above it, a small
leaf resting on the rim, uniform 2px stroke, rounded line caps and corners,
single deep-green outline color #16634a, soft mint-green fill #e6f0ea on the
bowl body, minimal geometric shapes, no gradients, no shadows, no text,
centered composition with generous padding, isolated on solid white
background, 1:1 square, clean scientific-platform aesthetic
```

### 2. 药材 Herb — 草本植物
```
Flat vector line icon of a medicinal herb plant: one upright stem with three
pairs of oval serrated leaves and a small root tuft at the base, uniform 2px
stroke, rounded line caps and corners, single deep-green outline color
#16634a, soft mint-green fill #e6f0ea on the leaves, minimal geometric
shapes, no gradients, no shadows, no text, centered composition with
generous padding, isolated on solid white background, 1:1 square, clean
scientific-platform aesthetic
```

### 3. 成分 Ingredient — 分子结构
```
Flat vector line icon of a small organic molecule: a hexagonal benzene-like
ring with two short bond lines extending to small circle atoms, uniform 2px
stroke, rounded line caps and corners, single deep-green outline color
#16634a, soft mint-green fill #e6f0ea inside the ring, minimal geometric
shapes, no gradients, no shadows, no text, centered composition with
generous padding, isolated on solid white background, 1:1 square, clean
scientific-platform aesthetic
```

### 4. 蛋白 Protein — 折叠蛋白链
```
Flat vector line icon of a folded protein chain: a smooth ribbon-like
polypeptide strand with two gentle alpha-helix coils in the middle, drawn as
one continuous rounded line, uniform 2px stroke, rounded line caps and
corners, single deep-green outline color #16634a, soft mint-green fill
#e6f0ea between the coils, minimal geometric shapes, no gradients, no
shadows, no text, centered composition with generous padding, isolated on
solid white background, 1:1 square, clean scientific-platform aesthetic
```

### 5. 疾病 Disease — 病历卡 + 医疗十字
```
Flat vector line icon of a medical record card: a rounded rectangular card
with two horizontal lines suggesting text and a circular badge at the bottom
right corner containing a plus cross, uniform 2px stroke, rounded line caps
and corners, single deep-green outline color #16634a, soft mint-green fill
#e6f0ea on the card and badge, minimal geometric shapes, no gradients, no
shadows, no text, centered composition with generous padding, isolated on
solid white background, 1:1 square, clean scientific-platform aesthetic
```

## 二、功能卡插画（首页「选择工作流」四张卡片右侧）

### 6. 智能辨证咨询 AI consultation — 人物 + 对话气泡
```
Flat vector line icon of a person in consultation: a simple head-and-shoulder
figure on the left with a rounded speech bubble on the upper right containing
three small dots, uniform 2px stroke, rounded line caps and corners, single
deep-green outline color #16634a, soft mint-green fill #e6f0ea on the figure
and bubble, minimal geometric shapes, no gradients, no shadows, no text,
centered composition with generous padding, isolated on solid white
background, 1:1 square, clean scientific-platform aesthetic
```

### 7. 网络药理学 Network — 节点网络图
```
Flat vector line icon of a small network graph: six circular nodes connected
by thin straight edges, one slightly larger hub node near the center, uniform
2px stroke, rounded line caps and corners, single deep-green outline color
#16634a, soft mint-green fill #e6f0ea in every node, minimal geometric
shapes, no gradients, no shadows, no text, centered composition with
generous padding, isolated on solid white background, 1:1 square, clean
scientific-platform aesthetic
```

### 8. 组方优化 Combination optimization — 上升柱状图（唯一允许金色点缀）
```
Flat vector line icon of three ascending bar columns of increasing height
inside a thin baseline frame, the tallest bar carrying a small five-pointed
star on top, uniform 2px stroke, rounded line caps and corners, single
deep-green outline color #16634a, soft mint-green fill #e6f0ea on the bars,
the star filled with warm accent gold #aa7841, minimal geometric shapes, no
gradients, no shadows, no text, centered composition with generous padding,
isolated on solid white background, 1:1 square, clean
scientific-platform aesthetic
```

### 9. 分子对接 Molecular docking — 配体嵌入受体口袋
```
Flat vector line icon of molecular docking: a large rounded receptor pocket
shape formed by a partial ring of small circle atoms, with a small hexagonal
ligand molecule fitting into the opening, uniform 2px stroke, rounded line
caps and corners, single deep-green outline color #16634a, soft mint-green
fill #e6f0ea on both the pocket and the ligand, minimal geometric shapes, no
gradients, no shadows, no text, centered composition with generous padding,
isolated on solid white background, 1:1 square, clean scientific-platform
aesthetic
```

## 生成与验收清单

- [ ] 9 张全部生成后，平铺在同一张白底图上肉眼比对：线条粗细、圆角、留白必须一致，不像同一家的重出。
- [ ] 只允许第 8 张出现金色 `#aa7841`，其余 8 张只能深绿 + 薄荷绿。
- [ ] 生成结果位图如需上 hero 深绿横幅：用第 2 次编辑指令「把背景改为透明 / 白色图形置于透明背景」各出一份反白版本（白线条 + 深绿填充 `#16634a`），供深底使用。
- [ ] 保存路径建议：`src/jnuherb/data/assets/`（网页用）与 `figures/ui_icons/`（论文用）各存一份。
- [ ] 每张按实际显示尺寸导出（卡片插画 ≤ 400px，hero 图标 ≤ 200px），控制体积。
