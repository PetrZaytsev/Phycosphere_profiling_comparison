import os
import sys
import re

# Force non-interactive backend for headless server execution
import matplotlib
matplotlib.use('Agg')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

print("=== GENERATING CROSS-PLATFORM VERTICAL SPLIT HEATMAP ===")

# 1. Load unified RefSeq-aligned abundance matrices from data/
df_ngs_raw = pd.read_csv('NGS_TOTAL_CLEAN_bacterial_core_OI1-12_RefSeq.tsv', sep='\t', index_col=0)
df_wgs_raw = pd.read_csv('WMS_FINAL_abundance_matrix_RefSeq.txt', sep='\t', index_col=0)
df_ont_raw = pd.read_csv('ONT_16S_TOTAL_matrix_all_samples.tsv', sep='\t', index_col=0)

# Aggregate technical PCR replicates for NGS (map OI01_repl1/2 -> OI01)
ngs_mapped_cols = []
for col in df_ngs_raw.columns:
    num_match = re.search(r'OI\.?(\d+)', col)
    ngs_mapped_cols.append(f"OI{int(num_match.group(1)):02d}" if num_match else col)
df_ngs_raw.columns = ngs_mapped_cols
df_ngs = df_ngs_raw.T.groupby(df_ngs_raw.columns).sum().T

# Synchronize sample identifiers across WMS and ONT-16S
target_samples = [f"OI{i:02d}" for i in range(1, 13)]
if len(df_wgs_raw.columns) == 12:
    df_wgs_raw.columns = target_samples
if len(df_ont_raw.columns) == 12:
    df_ont_raw.columns = target_samples

# Experimental stage labels
centroids = [
    'Axenic Green', 
    'Axenic Red', 
    'Co-culture Green', 
    'Co-culture Red', 
    'Xenic Green', 
    'Xenic Red'
]

# Mapping of biological replicates to experimental stages
legacy_mapping = {
    'Axenic Green': ('OI01', 'OI02'),   
    'Axenic Red': ('OI03', 'OI04'),
    'Co-culture Green': ('OI05', 'OI06'), 
    'Co-culture Red': ('OI07', 'OI08'),
    'Xenic Green': ('OI09', 'OI10'),    
    'Xenic Red': ('OI11', 'OI12')
}

top_8_taxa = [
    'Arsenicitalea', 'Devosia', 'Methylobacillus', 'Microbacterium',
    'Hymenobacter', 'Acidovorax', 'Blastomonas', 'Sphingomonas'
]

def build_curated_matrix(df_raw, core_list, samples):
    df_clean = pd.DataFrame(0.0, index=core_list, columns=samples)
    for taxon_raw in df_raw.index:
        for biomarker in core_list:
            if re.search(biomarker, str(taxon_raw), re.IGNORECASE):
                df_clean.loc[biomarker, samples] += df_raw.loc[taxon_raw, samples].astype(float)
                break
    return df_clean

wgs_curated = build_curated_matrix(df_wgs_raw, top_8_taxa, target_samples)
ont_curated = build_curated_matrix(df_ont_raw, top_8_taxa, target_samples)
ngs_curated = build_curated_matrix(df_ngs, top_8_taxa, target_samples)

# Normalize within the core top-8 bacterial taxa pool (%)
wgs_rel = wgs_curated.div(wgs_curated.sum(axis=0) + 1e-9, axis=1) * 100
ont_rel = ont_curated.div(ont_curated.sum(axis=0) + 1e-9, axis=1) * 100
ngs_rel = ngs_curated.div(ngs_curated.sum(axis=0) + 1e-9, axis=1) * 100

# Assemble 18x16 master matrix
row_vectors = []
tech_labels = []

for center in centroids:
    s1, s2 = legacy_mapping[center]
    for tech, df_rel in [('W', wgs_rel), ('O', ont_rel), ('N', ngs_rel)]:
        row_line_expanded = []
        for taxon in top_8_taxa:
            row_line_expanded.append(df_rel.loc[taxon, s1])
            row_line_expanded.append(df_rel.loc[taxon, s2])
        row_vectors.append(row_line_expanded)
        tech_labels.append(tech)

x_rep_labels = ["1", "2"] * 8
df_vertical_master = pd.DataFrame(row_vectors, columns=x_rep_labels)

# Sequential white-to-red colormap with enhanced sensitivity for low abundances
color_positions = [
    (0.0, "#FFFFFF"),   # 0%   - Pure white for zero values
    (0.05, "#FFEBEE"),  # 5%   - Very soft blush for low abundance
    (0.25, "#FFCDD2"),  # 25%  - Light rose
    (0.55, "#EF5350"),  # 55%  - Coral red
    (1.0, "#B71C1C")    # 100% - Deep crimson for dominant taxa
]
custom_cmap = LinearSegmentedColormap.from_list("white_to_red", color_positions, N=256)

# Render canvas
fig, (ax_bio, ax_heat) = plt.subplots(1, 2, figsize=(22, 16), gridspec_kw={'width_ratios': [3.5, 10]})

ax_bio.set_ylim(0, 18)
ax_bio.invert_yaxis()
ax_bio.axis('off')

for idx, center in enumerate(centroids):
    ax_bio.text(0.90, idx * 3 + 1.5, center, fontsize=24, fontweight='bold', ha='right', va='center')

cbar_ax = fig.add_axes([0.94, 0.3, 0.02, 0.4])
sns.heatmap(df_vertical_master, cmap=custom_cmap, vmin=0, vmax=100,
            cbar_ax=cbar_ax, cbar_kws={'pad': 0.02}, 
            linewidths=0.5, linecolor="#B0BEC5", ax=ax_heat)

cbar_ax.tick_params(labelsize=20)
cbar_ax.set_ylabel('Relative Abundance (%)', fontsize=22, fontweight='bold', labelpad=15)

ax_heat.set_xticklabels([])
ax_heat.set_yticklabels(tech_labels, rotation=0, fontsize=22, fontweight='bold')

# Grid dividers for stages and taxa
for i in range(1, 6):
    ax_heat.axhline(y=i*3, color='black', linewidth=2.5)
for i in range(1, 16):
    if i % 2 != 0:
        ax_heat.axvline(x=i, color='#78909C', linewidth=1.0, linestyle='--')
for i in range(1, 8):
    ax_heat.axvline(x=i*2, color='black', linewidth=3.0)
ax_heat.axvline(x=0, color='black', linewidth=4.5)

# Top axis: bioreplicate labels (1, 2)
ax_top = ax_heat.twiny()
ax_top.set_xlim(ax_heat.get_xlim())
ax_top.xaxis.set_ticks_position('top')
ax_top.xaxis.set_label_position('top')
ax_top.set_xticks([i + 0.5 for i in range(16)])
ax_top.set_xticklabels(x_rep_labels, fontsize=22, fontweight='bold')
ax_top.spines['top'].set_visible(False)

# Bottom axis: bacterial genus names at 45-degree angle
ax_bottom = ax_heat.twiny()
ax_bottom.set_xlim(ax_heat.get_xlim())
ax_bottom.xaxis.set_ticks_position('bottom')
ax_bottom.xaxis.set_label_position('bottom')
ax_bottom.spines['bottom'].set_position(('outward', 15))
ax_bottom.set_xticks([idx * 2 + 1.0 for idx in range(8)])
ax_bottom.set_xticklabels(top_8_taxa, rotation=45, ha='right', fontsize=24, fontweight='bold', style='italic')

# Title
plt.suptitle("Cross-Platform Vertical Split Heatmap with Parallel Bioreplicates in Columns\n(W: Shotgun WMS | O: ONT Full-Length 16S | N: Illumina NGS V4)", 
             fontsize=14, fontweight='bold', y=0.96)

plt.subplots_adjust(left=0.05, right=0.92, top=0.90, bottom=0.22, wspace=0.02)

plt.savefig("heatmap_frequency.pdf", format='pdf', bbox_inches='tight')
plt.savefig("heatmap_frequency.png", dpi=300, bbox_inches='tight')
plt.close()

print("[+] Success! Saved heatmap_frequency.pdf and .png")
