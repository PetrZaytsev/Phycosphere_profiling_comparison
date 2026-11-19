import os
import sys
import re

# Force non-interactive backend for headless server execution
import matplotlib
matplotlib.use('Agg')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

try:
    from matplotlib_venn import venn3, venn3_circles
except ImportError:
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "matplotlib-venn", "--quiet"])
    from matplotlib_venn import venn3, venn3_circles

print("=== GENERATING FIGURE 2C: CROSS-PLATFORM VENN DIAGRAM ===")

# 1. Load unified RefSeq-aligned abundance matrices
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

target_samples = [f"OI{i:02d}" for i in range(1, 13)]
df_wgs_raw.columns = target_samples
if len(df_ont_raw.columns) == 12:
    df_ont_raw.columns = target_samples

def parse_genus_smart(name, file_type):
    s = str(name).strip().replace('"', '')
    if any(w in s.lower() for w in ['chloroplast', 'plastid', 'eukaryota', 'scenedesmus']):
        return None
    if re.search(r'Arsenic[i|e]talea', s, re.IGNORECASE):
        return 'Arsenicitalea'
        
    if file_type == 'NGS':
        # Support both legacy g__ prefixes and clean RefSeq genus names
        genus_match = re.search(r'g__([a-zA-Z0-9\-_]+)', s)
        if genus_match:
            g_name = genus_match.group(1).strip().replace('_', '')
        else:
            last_part = s.split(';')[-1].strip()
            g_name = re.sub(r'[^a-zA-Z0-9\-]', '', last_part)
        if len(g_name) > 2 and not any(w in g_name.lower() for w in ['unassigned', 'unknown', 'bacteria', 'incertae', 'taxon']):
            return g_name.capitalize()
                
    elif file_type == 'ONT':
        words = s.split()
        if len(words) > 0:
            g_name = re.sub(r'[^a-zA-Z]', '', words[0])
            if len(g_name) > 2 and g_name.lower() != 'taxon':
                return g_name.capitalize()
                
    elif file_type == 'WMS':
        parts = [p.strip() for p in re.split(r';|\|', s) if p.strip()]
        if len(parts) >= 5:
            last_part = parts[-1]
            words = last_part.split()
            if len(words) > 0:
                g_name = re.sub(r'[^a-zA-Z]', '', words[0])
                if g_name.lower().endswith('aceae') or g_name.lower().endswith('ales'):
                    return None
                if len(g_name) > 2 and g_name.lower() != 'bacteria':
                    return g_name.capitalize()
    return None

df_wgs_raw.index = [parse_genus_smart(i, 'WMS') for i in df_wgs_raw.index]
df_ont_raw.index = [parse_genus_smart(i, 'ONT') for i in df_ont_raw.index]
df_ngs.index = [parse_genus_smart(i, 'NGS') for i in df_ngs.index]

df_wgs_raw = df_wgs_raw.loc[df_wgs_raw.index.notnull()]
df_ont_raw = df_ont_raw.loc[df_ont_raw.index.notnull()]
df_ngs = df_ngs.loc[df_ngs.index.notnull()]

df_wgs_g = df_wgs_raw.groupby(df_wgs_raw.index).sum()
df_ont_g = df_ont_raw.groupby(df_ont_raw.index).sum()
df_ngs_g = df_ngs.groupby(df_ngs.index).sum()

all_global_taxa = set(df_wgs_g.index) | set(df_ont_g.index) | set(df_ngs_g.index)
all_global_taxa = sorted([t for t in all_global_taxa if t and t not in ['Bacteria', 'None']])

wgs_rel = df_wgs_g.reindex(all_global_taxa, fill_value=0.0).div(df_wgs_g.sum(axis=0) + 1e-9, axis=1) * 100
ont_rel = df_ont_g.reindex(all_global_taxa, fill_value=0.0).div(df_ont_g.sum(axis=0) + 1e-9, axis=1) * 100
ngs_rel = df_ngs_g.reindex(all_global_taxa, fill_value=0.0).div(df_ngs_g.sum(axis=0) + 1e-9, axis=1) * 100

wgs_set = set([t for t in all_global_taxa if wgs_rel.loc[t].mean() > 0.005])
ont_set = set([t for t in all_global_taxa if ont_rel.loc[t].mean() > 0.005])
ngs_set = set([t for t in all_global_taxa if ngs_rel.loc[t].mean() > 0.005])

plt.figure(figsize=(10, 10))
v = venn3(subsets=[wgs_set, ont_set, ngs_set], 
          set_labels=('WMS-ONT', '16S-ONT', '16S-NGS'),
          set_colors=('#B3E5FC', '#D32F2F', '#F8F9FA'), alpha=0.75)

for text in v.set_labels:
    if text: 
        text.set_fontsize(16)
        text.set_fontweight('bold')

for text in v.subset_labels:
    if text: 
        text.set_fontsize(14)
        text.set_fontweight('bold')

c = venn3_circles(subsets=[wgs_set, ont_set, ngs_set], linestyle='solid', linewidth=1.2, color='black')

plt.tight_layout()
plt.savefig("venn_diagram.png", dpi=300, bbox_inches='tight')
plt.savefig("venn_diagram.pdf", format='pdf', bbox_inches='tight')
plt.close()

print("[+] Success! Saved venn_diagram.png and .pdf")
