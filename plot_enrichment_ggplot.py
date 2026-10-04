#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GO富集分析点图 - 按Ontology分类展示Top N
X轴: GeneRatio | 点大小: Count | 点颜色: -log10(P) 渐变色
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm
from matplotlib.lines import Line2D
import sys
import os


def plot_GO_dotplot(
    csv_path,
    output_path='GO_dotplot.pdf',
    top_n=3,
    figsize=(9, 7),
    dpi=300,
    title=None,
    x_label='GeneRatio',
    color_low='#8CC6ED',
    color_high='#F0A0A0',
    size_scale=18,
    category_order=None,
    category_labels=None,
):
    """
    绘制GO富集分析点图（ggplot2 风格）

    Parameters:
    -----------
    csv_path : str
        CSV文件路径，需包含 Term, Category, Count, GeneRatio, P-value 列
    output_path : str
        输出图片路径 (.pdf/.png/.svg)
    top_n : int
        每个Ontology展示Top N个条目（按P-value最显著选取）
    figsize : tuple
        图片尺寸 (宽, 高) inches
    dpi : int
        分辨率
    title : str or None
        图表标题
    x_label : str
        X轴标签
    color_low : str
        渐变色低端（-log10(P)较小）
    color_high : str
        渐变色高端（-log10(P)较大）
    size_scale : int
        点大小缩放系数
    category_order : list or None
        Ontology展示顺序（从上到下），默认 ['BP', 'CC', 'MF']
    category_labels : dict or None
        Ontology显示标签
    """

    # 默认值
    if category_order is None:
        category_order = ['BP', 'CC', 'MF']
    if category_labels is None:
        category_labels = {
            'BP': 'Biological Process (BP)',
            'CC': 'Cellular Component (CC)',
            'MF': 'Molecular Function (MF)',
        }

    # 1. 读取数据
    df = pd.read_csv(csv_path, sep=',', skipinitialspace=True)

    # 清理
    df['Term'] = df['Term'].astype(str).str.replace('"', '').str.strip()
    df['Category'] = df['Category'].astype(str).str.strip()
    df['Count'] = pd.to_numeric(df['Count'], errors='coerce')
    df['GeneRatio'] = pd.to_numeric(df['GeneRatio'], errors='coerce')
    df['P-value'] = pd.to_numeric(df['P-value'], errors='coerce')
    df = df.dropna(subset=['Term', 'Category', 'GeneRatio', 'P-value', 'Count'])

    # 2. 各Ontology取Top N（按P-value最显著），按 category_order 拼接
    selected = []
    for cat in category_order:
        subset = df[df['Category'] == cat].nsmallest(top_n, 'P-value').copy()
        if len(subset) > 0:
            subset['Display'] = subset['Term'].apply(
                lambda x: x.split('(GO:')[0].strip() if '(GO:' in x else x
            )
            # 每组内按GeneRatio降序（大的在上方）
            subset = subset.sort_values('GeneRatio', ascending=False).reset_index(drop=True)
            selected.append(subset)
        else:
            print(f"Warning: No data found for category '{cat}'")

    if not selected:
        print("Error: No data available for any category")
        sys.exit(1)

    df_plot = pd.concat(selected, ignore_index=True)

    labels = df_plot['Display'].values
    gene_ratios = df_plot['GeneRatio'].values
    counts = df_plot['Count'].values
    pvals = df_plot['P-value'].values
    cat_vals = df_plot['Category'].values
    neg_log_p = -np.log10(pvals)

    # 3. 渐变色映射
    norm = plt.Normalize(neg_log_p.min(), neg_log_p.max())
    cmap = mcolors.LinearSegmentedColormap.from_list('go_grad', [color_low, color_high])
    dot_colors = cmap(norm(neg_log_p))

    # 点大小
    dot_sizes = counts * size_scale

    # 4. 绘图
    fig, ax = plt.subplots(figsize=figsize)

    y_pos = np.arange(len(df_plot))

    # 散点
    ax.scatter(gene_ratios, y_pos, s=dot_sizes, c=dot_colors,
               edgecolors='black', linewidths=0.6, zorder=3, clip_on=False)

    # Y轴标签
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=10)

    # 反转y轴：category_order[0] 在最顶部
    ax.invert_yaxis()

    # 分类分隔线
    for i in range(len(df_plot) - 1):
        if cat_vals[i] != cat_vals[i + 1]:
            sep_y = (y_pos[i] + y_pos[i + 1]) / 2 + 0.5
            ax.axhline(y=sep_y, color='grey', linewidth=0.8, linestyle='-', alpha=0.5, zorder=2)

    # 连续渐变色图例 (colorbar: -log10(P))
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.5, aspect=18, pad=0.03, anchor=(0, 0.12))
    cbar.set_label(r'$-\log_{10}(P)$', fontsize=11, fontweight='bold')
    cbar.outline.set_edgecolor('black')
    cbar.outline.set_linewidth(0.8)
    cbar.ax.tick_params(width=0.6, length=3, color='black', labelsize=9)

    # 点大小图例 (Count)
    legend_counts = sorted(set(counts))
    if len(legend_counts) > 3:
        step = max(1, len(legend_counts) // 3)
        legend_counts = [legend_counts[0], legend_counts[step], legend_counts[-1]]
    legend_counts = sorted(set(legend_counts))

    legend_elements = []
    for c in legend_counts:
        s = c * size_scale
        legend_elements.append(Line2D([], [], marker='o', color='w',
                                       markerfacecolor='#999999',
                                       markeredgecolor='black', markeredgewidth=0.6,
                                       markersize=np.sqrt(s) / 2, label=f'{c}'))
    ax.legend(handles=legend_elements, loc='lower right',
              frameon=True, fancybox=False, edgecolor='black',
              fontsize=9, title='Count', title_fontsize=9,
              handletextpad=0.5, borderpad=0.8)

    # 坐标轴
    ax.set_xlabel(x_label, fontsize=12, fontweight='bold')
    if title is None:
        title = f'GO Enrichment Dot Plot (Top {top_n} per Ontology)'
    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)

    # theme_bw
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(0.6)
    ax.tick_params(axis='both', width=0.6, length=3, color='black', labelsize=10)
    ax.grid(False)

    plt.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches='tight')

    png_path = output_path.rsplit('.', 1)[0] + '.png'
    fig.savefig(png_path, dpi=dpi, bbox_inches='tight')

    print(f"Saved: {output_path}")
    print(f"Saved: {png_path}")
    plt.show()
    return fig, ax


# ============================================================================
# 主程序
# ============================================================================

if __name__ == '__main__':

    if len(sys.argv) > 1:
        csv_file = sys.argv[1]
    else:
        csv_file = 'GO_Enrichment.csv'

    if not os.path.exists(csv_file):
        print(f"Usage: python {sys.argv[0]} <csv_file>")
        print(f"Error: File '{csv_file}' not found")
        sys.exit(1)

    plot_GO_dotplot(
        csv_path=csv_file,
        output_path='GO_dotplot.pdf',
        top_n=3,
        figsize=(9, 7),
        dpi=300,
    )


# ============================================================================
# 使用说明:
# ---------
# 1. 命令行:
#    python plot_GO_dotplot.py your_GO.csv
#
# 2. 作为模块导入:
#    from plot_GO_dotplot import plot_GO_dotplot
#    plot_GO_dotplot('GO.csv', top_n=5)
#
# 3. CSV格式要求:
#    | Term | Category | Count | GeneRatio | P-value |
#    |------|----------|-------|-----------|---------|
#    | G protein-coupled receptor activity (GO:0004930) | MF | 126 | 0.1059 | 1.3e-59 |
#
# 关键特点:
#    - X轴: GeneRatio
#    - 排序: BP→CC→MF 从上到下，每组内按GeneRatio降序
#    - 点大小: Count
#    - 点颜色: -log10(P-value)，#8CC6ED → #F0A0A0 渐变
#    - 两个图例: colorbar(-log10P) + 大小(Count)
#    - 主题: theme_bw
# ============================================================================