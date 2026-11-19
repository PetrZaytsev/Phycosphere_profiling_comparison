import os
import sys
import numpy as np
import pandas as pd

# Force non-interactive backend for headless server execution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

print("=== GENERATING FIGURE 3: PAIRWISE CONGRUENCE DIAGONAL HEATMAPS ===")

# 1. Dynamically load calculated metrics from Table 3 TSV
metrics_file = "Table3_congruence_metrics.tsv"
if not os.path.exists(metrics_file):
    raise FileNotFoundError(f"Metrics file '{metrics_file}' not found. Run calculate_metrics.py first.")

df_metrics = pd.read_csv(metrics_file, sep='\t').set_index('Comparison Pair')

# Helper function to safely retrieve pairwise values regardless of pair order
def get_pair_val(t1, t2, col_name, diag_val):
    if t1 == t2:
        return diag_val
    p1 = f"{t1} vs. {t2}"
    p2 = f"{t2} vs. {t1}"
    if p1 in df_metrics.index:
        return float(df_metrics.loc[p1, col_name])
    elif p2 in df_metrics.index:
        return float(df_metrics.loc[p2, col_name])
    else:
        raise KeyError(f"Pair '{p1}' or '{p2}' not found in {metrics_file}")

tech_labels = ["16S-ONT", "16S-NGS", "WMS-ONT"]

def build_sym_matrix(col_name, diag_val=0.0):
    mat = np.zeros((3, 3))
    for i, t1 in enumerate(tech_labels):
        for j, t2 in enumerate(tech_labels):
            mat[i, j] = get_pair_val(t1, t2, col_name, diag_val)
    return mat

bc_data = build_sym_matrix('Mean Bray-Curtis Dissimilarity', diag_val=0.0)
pm_data = build_sym_matrix('Procrustes Analysis (M2 error)', diag_val=0.0)
mr_data = build_sym_matrix('Mantel Test (r)', diag_val=1.0)
mp_data = build_sym_matrix('Mantel Test (p-value)', diag_val=0.0)

matrix_layout = [
    {"title": "Bray-Curtis\nDissimilarity", "data": bc_data, "cmap": "YlOrRd", "is_corr": False},
    {"title": "Procrustes Analysis\nResidual Error (M2)", "data": pm_data, "cmap": "OrRd", "is_corr": False},
    {"title": "Mantel Test\nCorrelation (r)", "data": mr_data, "cmap": "coolwarm", "is_corr": True},
    {"title": "Mantel Test\nSignificance (p)", "data": mp_data, "cmap": "BuPu", "is_corr": False}
]

fig, axes = plt.subplots(2, 2, figsize=(13, 12))
axes = axes.flatten()

# Mask upper triangle above main diagonal
mask = np.triu(np.ones((3, 3), dtype=bool), k=1)

for idx, item in enumerate(matrix_layout):
    ax = axes[idx]
    df = pd.DataFrame(item["data"], index=tech_labels, columns=tech_labels)
    
    vmin, vmax = (-1.0, 1.0) if item["is_corr"] else (0.0, 1.0)
        
    sns.heatmap(df, mask=mask, cmap=item["cmap"], vmin=vmin, vmax=vmax, annot=False,
                cbar_kws={"shrink": 0.7, "pad": 0.04},
                linewidths=1.5, linecolor="white", square=True, ax=ax)
    
    for i in range(3):
        for j in range(3):
            if j <= i:
                val = item["data"][i, j]
                ax.text(j + 0.5, i + 0.5, f"{val:.3f}", 
                        fontsize=13, fontweight="bold", color="black",
                        ha="center", va="center")
    
    ax.set_title(item["title"], fontsize=13, fontweight="bold", pad=12)
    ax.set_xticklabels(tech_labels, rotation=45, ha="right", fontsize=11, fontweight="bold")
    ax.set_yticklabels(tech_labels, rotation=0, fontsize=11, fontweight="bold")

plt.suptitle("Pairwise Cross-Platform Congruence Topology and Statistical Metrics", 
             fontsize=15, fontweight="bold", y=0.98)

plt.tight_layout()
plt.savefig("heatmaps_congruence.pdf", format="pdf", bbox_inches="tight")
plt.savefig("heatmaps_congruence.png", dpi=300, bbox_inches='tight')
plt.close()

print("[+] Success! Saved updated heatmaps_congruence.pdf and .png")
