#!/usr/bin/env bash
# ==============================================================================
# Step 3: Two-step Metaxa2 WMS profiling pipeline.
#         1) Native metaxa2 extraction (--classifications F)
#         2) Custom database classification via metaxa2_c
# Usage:  ./03_run_metaxa2_classification.sh -i <fastq_dir> -d <custom_db_dir> -o <results_dir> [-c <threads>]
# ==============================================================================
set -euo pipefail

THREADS=16
FASTQ_DIR=""
DB_DIR=""
OUT_DIR=""

while getopts "i:d:o:c:h" opt; do
    case $opt in
        i) FASTQ_DIR="$OPTARG" ;;
        d) DB_DIR="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        c) THREADS="$OPTARG" ;;
        h)
            echo "Usage: $0 -i <chunked_fastq_dir> -d <compiled_metaxa_db_dir> -o <output_dir> [-c <cpu_threads>]"
            exit 0
            ;;
        *) exit 1 ;;
    esac
done

if [[ -z "$FASTQ_DIR" || -z "$DB_DIR" || -z "$OUT_DIR" ]]; then
    echo "Error: Missing required arguments (-i, -d, -o)." >&2
    exit 1
fi

# Append /blast prefix required by metaxa2_c index lookup
BLAST_DB_PREFIX="${DB_DIR%/}/blast"
EXTRACT_DIR="${OUT_DIR}/01_extracted_rrna"
CLASS_DIR="${OUT_DIR}/02_classified_taxonomy"

mkdir -p "$EXTRACT_DIR" "$CLASS_DIR"

shopt -s nullglob
FQ_FILES=("${FASTQ_DIR}"/*.fastq)

if [[ ${#FQ_FILES[@]} -eq 0 ]]; then
    echo "Error: No .fastq files found in ${FASTQ_DIR}" >&2
    exit 1
fi

for fq_path in "${FQ_FILES[@]}"; do
    sample_id=$(basename "$fq_path" .fastq)
    sample_id="${sample_id%_chunked}"

    echo "[*] Extracting 16S rRNA sequences for ${sample_id}..."
    metaxa2 \
        -i "$fq_path" \
        -o "${EXTRACT_DIR}/${sample_id}" \
        --classifications F \
        --cpu "$THREADS"

    bac_fasta="${EXTRACT_DIR}/${sample_id}.bacteria.fasta"

    if [[ -s "$bac_fasta" ]]; then
        echo "[*] Classifying bacterial 16S reads for ${sample_id} against custom DB..."
        # Explicit -T vector bypasses empty cutoff file array subscript bug in metaxa2_c
        metaxa2_c \
            -i "$bac_fasta" \
            -d "$BLAST_DB_PREFIX" \
            -g E \
            -T 0,60,70,75,85,90,97 \
            -o "${CLASS_DIR}/${sample_id}" \
            --cpu "$THREADS"
    else
        echo "[_] No bacterial 16S reads detected for ${sample_id}. Creating empty taxonomy file."
        touch "${CLASS_DIR}/${sample_id}.taxonomy.txt"
    fi
done

echo "[+] Classification complete. Outputs saved in: ${CLASS_DIR}"
