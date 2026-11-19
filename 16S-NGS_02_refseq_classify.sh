#!/bin/bash
# ==============================================================================
# QIIME2 Downstream Taxonomy Classification and Matrix Extraction Pipeline
# Reference Database: Custom Trained NCBI RefSeq + rrnDB v5.6
# ==============================================================================
set -e

# Define input arguments with defaults
CLASSIFIER=${1:-"custom_refseq_v5.6_classifier.qza"}
OUT_DIR=${2:-"qiime2_refseq_results"}

mkdir -p "$OUT_DIR"

echo "=================================================================="
echo "STARTING SKLEARN TAXONOMY CLASSIFICATION VIA CUSTOM REFseq v5.6"
echo "=================================================================="

# 1. Run Sklearn taxonomic classification for both sequencing runs
echo "[*] Classifying Representative Sequences for Run 1..."
qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads rep-seqs-run1.qza \
  --o-classification "${OUT_DIR}/taxonomy-run1-refseq.qza"

echo "[*] Classifying Representative Sequences for Run 2..."
qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads rep-seqs-run2.qza \
  --o-classification "${OUT_DIR}/taxonomy-run2-refseq.qza"

# 2. Collapse feature tables to Genus level (Level 6)
echo "[*] Collapsing Feature Tables to Genus level (Level 6)..."
qiime taxa collapse \
  --i-table table-run1.qza \
  --i-taxonomy "${OUT_DIR}/taxonomy-run1-refseq.qza" \
  --p-level 6 \
  --o-collapsed-table "${OUT_DIR}/table-run1-genus.qza"

qiime taxa collapse \
  --i-table table-run2.qza \
  --i-taxonomy "${OUT_DIR}/taxonomy-run2-refseq.qza" \
  --p-level 6 \
  --o-collapsed-table "${OUT_DIR}/table-run2-genus.qza"

# 3. Export BIOM tables and convert to standard TSV matrices
echo "[*] Exporting Genus tables to TSV format..."
qiime tools export --input-path "${OUT_DIR}/table-run1-genus.qza" --output-path "${OUT_DIR}/exported_table1"
qiime tools export --input-path "${OUT_DIR}/table-run2-genus.qza" --output-path "${OUT_DIR}/exported_table2"

biom convert -i "${OUT_DIR}/exported_table1/feature-table.biom" -o "${OUT_DIR}/matrix_run1_raw.tsv" --to-tsv
biom convert -i "${OUT_DIR}/exported_table2/feature-table.biom" -o "${OUT_DIR}/matrix_run2_raw.tsv" --to-tsv

# 4. Perform matrix merging, header parsing, and strict host decontamination via Python
echo "[*] Running final matrix merge, header resolution, and host decontamination..."
python3 -c "
import pandas as pd
import re

# Read raw quantified abundance tables
r1 = pd.read_csv('${OUT_DIR}/matrix_run1_raw.tsv', sep='\t', skiprows=1, index_col=0)
r2 = pd.read_csv('${OUT_DIR}/matrix_run2_raw.tsv', sep='\t', skiprows=1, index_col=0)

r1.columns = [f'{c}_run1' for c in r1.columns]
r2.columns = [f'{c}_run2' for c in r2.columns]
combined = pd.concat([r1, r2], axis=1).fillna(0.0)

# Normalize multiplexed sample headers into target nomenclature (OI01_repl1 ... OI12_repl2)
unified = []
for col in combined.columns:
    num_match = re.search(r'OI\.?(\d+)', str(col), re.IGNORECASE)
    repl_match = re.search(r'repl\.?(\d+)', str(col), re.IGNORECASE)
    if num_match and repl_match:
        unified.append(f'OI{int(num_match.group(1)):02d}_repl{int(repl_match.group(1))}')
    else:
        unified.append(str(col))
combined.columns = unified

# Slice strict 24 core experimental columns
target_24 = [f'OI{i:02d}_repl{r}' for i in range(1, 13) for r in (1, 2)]
combined = combined.T.groupby(combined.columns).sum().T
combined = combined[[c for c in target_24 if c in combined.columns]]

# Strict filtration and isolation of eukaryotic microalgal host rRNA and organelle noise
rows_clean, rows_host = [], []
for taxon in combined.index:
    t_str = str(taxon).strip()
    genus = t_str.split(';')[-1].strip()
    
    # Filter out empty paths, unassigned slots, chloroplasts (Stanieria), and host artifacts
    if genus == '' or genus.startswith('__') or any(w in t_str.lower() for w in ['chloroplast', 'plastid', 'eukaryota', 'scenedesmus', 'tetradesmus', 'acetonema', 'stanieria', 'unassigned']):
        rows_clean.append('DROP_MEMBER')
        rows_host.append('Host_Scenedesmus_rRNA_and_Organelles')
    else:
        rows_clean.append(genus)
        rows_host.append(genus)

# Export Matrix 1: Purified Bacterial Core Profile (For Cross-Platform Benchmarking)
df_clean = combined.copy()
df_clean.index = rows_clean
mat_clean = df_clean.groupby(df_clean.index).sum().drop('DROP_MEMBER', errors='ignore')
mat_clean.to_csv('NGS_TOTAL_RefSeq_aligned_matrix.tsv', sep='\t')

# Export Matrix 2: Isolated Non-Target Host Fraction Profile (For Quality Control Validation)
df_host = combined.copy()
df_host.index = rows_host
mat_host = df_host.groupby(df_host.index).sum().drop('DROP_MEMBER', errors='ignore')
mat_host.to_csv('NGS_TOTAL_RefSeq_WITH_HOST_matrix.tsv', sep='\t')

print('[+] Downstream matrix assembly and decontamination finished successfully.')
"

echo "=================================================================="
echo "SUCCESS: ALL METADATA TABLES GENERATED VIA REFSEQ v5.6 RE-INDEX"
echo "=================================================================="
