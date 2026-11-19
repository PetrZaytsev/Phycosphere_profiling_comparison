#!/bin/bash
# ==============================================================================
# QIIME2 Pipeline for Raw Amplicon Processing (Illumina NGS 16S V4 Single-End)
# Optimized for variable run lengths: Run 1 (130 bp) and Run 2 (220 bp)
# ==============================================================================
# Usage: ./run_qiime2_pipeline.sh <manifest_r1> <manifest_r2> <out_dir>
# ==============================================================================

set -e

MANIFEST_R1=$1
MANIFEST_R2=$2
OUT_DIR=$3
THREADS=16

if [ -z "$MANIFEST_R1" ] || [ -z "$MANIFEST_R2" ] || [ -z "$OUT_DIR" ]; then
    echo "Usage: $0 <manifest_run1.txt> <manifest_run2.txt> <output_directory>"
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "[*] Importing raw sequencing reads into QIIME2 artifacts..."
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "$MANIFEST_R1" \
  --output-path "${OUT_DIR}/single-end-run1.qza" \
  --input-format SingleEndFastqManifestPhred33V2

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "$MANIFEST_R2" \
  --output-path "${OUT_DIR}/single-end-run2.qza" \
  --input-format SingleEndFastqManifestPhred33V2

echo "[*] Running DADA2 denoising, quality filtering, and chimera removal..."

echo "    [+] Processing Run 1 (Truncation length: 130 bp)..."
qiime dada2 denoise-single \
  --i-demultiplexed-seqs "${OUT_DIR}/single-end-run1.qza" \
  --p-trim-left 0 \
  --p-trunc-len 130 \
  --p-n-threads $THREADS \
  --o-table "${OUT_DIR}/table-run1.qza" \
  --o-representative-sequences "${OUT_DIR}/rep-seqs-run1.qza" \
  --o-denoising-stats "${OUT_DIR}/stats-run1.qza"

echo "    [+] Processing Run 2 (Truncation length: 220 bp)..."
qiime dada2 denoise-single \
  --i-demultiplexed-seqs "${OUT_DIR}/single-end-run2.qza" \
  --p-trim-left 0 \
  --p-trunc-len 220 \
  --p-n-threads $THREADS \
  --o-table "${OUT_DIR}/table-run2.qza" \
  --o-representative-sequences "${OUT_DIR}/rep-seqs-run2.qza" \
  --o-denoising-stats "${OUT_DIR}/stats-run2.qza"

echo "[*] Exporting unique representative sequences (ASVs) to FASTA format..."
qiime tools export --input-path "${OUT_DIR}/rep-seqs-run1.qza" --output-path "${OUT_DIR}/exported_run1"
qiime tools export --input-path "${OUT_DIR}/rep-seqs-run2.qza" --output-path "${OUT_DIR}/exported_run2"

cat "${OUT_DIR}/exported_run1/dna-sequences.fasta" "${OUT_DIR}/exported_run2/dna-sequences.fasta" > "${OUT_DIR}/all_ngs_asv.fasta"

echo "[*] Exporting feature abundance matrices to TSV format..."
qiime tools export --input-path "${OUT_DIR}/table-run1.qza" --output-path "${OUT_DIR}/exported_table1"
qiime tools export --input-path "${OUT_DIR}/table-run2.qza" --output-path "${OUT_DIR}/exported_table2"

# Convert BIOM to standard TSV and place them in the root working directory for downstream script access
biom convert -i "${OUT_DIR}/exported_table1/feature-table.biom" -o "NGS_RUN1_raw_abundance_matrix.tsv" --to-tsv
biom convert -i "${OUT_DIR}/exported_table2/feature-table.biom" -o "NGS_RUN2_raw_abundance_matrix.tsv" --to-tsv

echo "=============================================================================="
echo "QIIME2 STAGE COMPLETED SUCCESSFULLY!"
echo " Raw matrices generated in working directory: NGS_RUN1_raw_abundance_matrix.tsv"
echo "                                             NGS_RUN2_raw_abundance_matrix.tsv"
echo " Combined ASV FASTA file located at:         ${OUT_DIR}/all_ngs_asv.fasta"
echo "=============================================================================="
