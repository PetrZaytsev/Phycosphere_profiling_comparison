import os
import sys
import re
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform

try:
    import skbio
except ImportError:
    import subprocess
    print("Package scikit-bio not found. Installing in background...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "scikit-bio", "--quiet"])
    import skbio

from skbio.stats.distance import mantel
from scipy.spatial import procrustes

print("\n=== CALCULATING MULTI-OMICS CONGRUENCE METRICS FOR TABLE 3 ===")

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

# Synchronize sample identifiers across WMS and ONT-16S
target_samples = [f"OI{i:02d}" for i in range(1, 13)]
if len(df_wgs_raw.columns) == 12:
    df_wgs_raw.columns = target_samples
if len(df_ont_raw.columns) == 12:
    df_ont_raw.columns = target_samples

# Core consortium bacterial biomarkers for cross-platform benchmarking
core_biomarkers = ['Arsenicitalea', 'Devosia', 'Methylobacillus', 'Microbacterium', 'Hymenobacter']

def build_curated_matrix(df_raw, core_list, samples):
    df_clean = pd.DataFrame(0.0, index=core_list, columns=samples)
    for taxon_raw in df_raw.index:
        for biomarker in core_list:
            if re.search(biomarker, str(taxon_raw), re.IGNORECASE):
                df_clean.loc[biomarker, samples] += df_raw.loc[taxon_raw, samples].astype(float)
                break
    for col in samples:
        if df_clean[col].sum() == 0:
            df_clean[col] = 1e-5
    return df_clean

# Build aligned biomarker matrices
wgs_curated = build_curated_matrix(df_wgs_raw, core_biomarkers, target_samples)
ont_curated = build_curated_matrix(df_ont_raw, core_biomarkers, target_samples)
ngs_curated = build_curated_matrix(df_ngs, core_biomarkers, target_samples)

print(f"-> Matrices aligned across {len(core_biomarkers)} core phycosphere bacterial genera.")

# Convert read counts to relative abundances (%)
wgs_rel = wgs_curated.div(wgs_curated.sum(axis=0), axis=1) * 100
ont_rel = ont_curated.div(ont_curated.sum(axis=0), axis=1) * 100
ngs_rel = ngs_curated.div(ngs_curated.sum(axis=0), axis=1) * 100

# 1. Sample-wise Bray-Curtis dissimilarity
bc_wgs_vs_ont = [np.sum(np.abs(wgs_rel[col] - ont_rel[col])) / np.sum(wgs_rel[col] + ont_rel[col]) for col in target_samples]
bc_wgs_vs_ngs = [np.sum(np.abs(wgs_rel[col] - ngs_rel[col])) / np.sum(wgs_rel[col] + ngs_rel[col]) for col in target_samples]
bc_ont_vs_ngs = [np.sum(np.abs(ont_rel[col] - ngs_rel[col])) / np.sum(ont_rel[col] + ngs_rel[col]) for col in target_samples]

# 2. Distance matrices for Mantel test (Spearman correlation, 999 permutations)
dm_wgs = squareform(pdist(wgs_rel.T, metric='braycurtis'))
dm_ont = squareform(pdist(ont_rel.T, metric='braycurtis'))
dm_ngs = squareform(pdist(ngs_rel.T, metric='braycurtis'))

r_mantel_w_o, p_mantel_w_o, _ = mantel(dm_wgs, dm_ont, method='spearman', permutations=999)
r_mantel_w_n, p_mantel_w_n, _ = mantel(dm_wgs, dm_ngs, method='spearman', permutations=999)
r_mantel_o_n, p_mantel_o_n, _ = mantel(dm_ont, dm_ngs, method='spearman', permutations=999)

# 3. Procrustes analysis (M2 disparity)
_, _, m2_w_o = procrustes(wgs_rel.values.T, ont_rel.values.T)
_, _, m2_w_n = procrustes(wgs_rel.values.T, ngs_rel.values.T)
_, _, m2_o_n = procrustes(ont_rel.values.T, ngs_rel.values.T)

# Assemble Table 3
table_3_data = {
    'Comparison Pair': [
        'WMS-ONT vs. 16S-ONT', 
        'WMS-ONT vs. 16S-NGS',
        '16S-ONT vs. 16S-NGS'
    ],
    'Mean Bray-Curtis Dissimilarity': [
        round(np.mean(bc_wgs_vs_ont), 4), 
        round(np.mean(bc_wgs_vs_ngs), 4),
        round(np.mean(bc_ont_vs_ngs), 4)
    ],
    'Mantel Test (r)': [
        round(r_mantel_w_o, 4), 
        round(r_mantel_w_n, 4),
        round(r_mantel_o_n, 4)
    ],
    'Mantel Test (p-value)': [p_mantel_w_o, p_mantel_w_n, p_mantel_o_n],
    'Procrustes Analysis (M2 error)': [
        round(float(m2_w_o), 4), 
        round(float(m2_w_n), 4),
        round(float(m2_o_n), 4)
    ]
}

df_table_3 = pd.DataFrame(table_3_data)
df_table_3.to_csv("Table3_congruence_metrics.tsv", sep='\t', index=False)

print("\n--- UPDATED CONGRUENCE METRICS FOR TABLE 3 ---")
print(df_table_3.to_string(index=False))
print("==================================================")
