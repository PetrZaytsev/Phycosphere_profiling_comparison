#!/usr/bin/env bash
# ==============================================================================
# Step 2: Convert SUP BAM files to FASTQ and chunk ultra-long reads (>90kb)
#         into overlapping 50kb windows to prevent HMMER3 100K pipeline overflow.
# Usage:  ./02_convert_and_chunk_sup.sh -i <bam_dir> -o <chunked_fastq_dir> [-c <threads>]
# ==============================================================================
set -euo pipefail

THREADS=16
BAM_DIR=""
OUT_DIR=""

while getopts "i:o:c:h" opt; do
    case $opt in
        i) BAM_DIR="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        c) THREADS="$OPTARG" ;;
        h)
            echo "Usage: $0 -i <input_bam_dir> -o <output_fastq_dir> [-c <cpu_threads>]"
            exit 0
            ;;
        *) exit 1 ;;
    esac
done

if [[ -z "$BAM_DIR" || -z "$OUT_DIR" ]]; then
    echo "Error: Both -i <bam_dir> and -o <out_dir> are required." >&2
    exit 1
fi

TMP_FQ_DIR="${OUT_DIR}/tmp_raw_fastq"
mkdir -p "$TMP_FQ_DIR" "$OUT_DIR"

shopt -s nullglob
BAM_FILES=("${BAM_DIR}"/*.bam)

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "Error: No .bam files found in ${BAM_DIR}" >&2
    exit 1
fi

for bam_path in "${BAM_FILES[@]}"; do
    base_name=$(basename "$bam_path" .bam)
    raw_fq="${TMP_FQ_DIR}/${base_name}.fastq"
    safe_fq="${OUT_DIR}/${base_name}_chunked.fastq"

    echo "[*] Converting ${base_name}.bam to FASTQ..."
    samtools fastq -@ "$THREADS" "$bam_path" > "$raw_fq"

    echo "[*] Chunking ultra-long reads (>90kb) for ${base_name}..."
    python3 - "$raw_fq" "$safe_fq" << 'EOF'
import sys

fin_path, fout_path = sys.argv[1:3]
MAX_LEN = 90000
WIN_SIZE = 50000
STEP_SIZE = 45000

with open(fin_path, 'r') as fin, open(fout_path, 'w') as fout:
    while True:
        header = fin.readline()
        if not header:
            break
        seq = fin.readline().strip()
        plus = fin.readline()
        qual = fin.readline().strip()

        if len(seq) <= MAX_LEN:
            fout.write(f"{header}{seq}\n{plus}{qual}\n")
        else:
            clean_id = header.strip().split()[0]
            for idx, start in enumerate(range(0, len(seq), STEP_SIZE)):
                sub_seq = seq[start:start + WIN_SIZE]
                sub_qual = qual[start:start + WIN_SIZE]
                if len(sub_seq) > 500:
                    fout.write(f"{clean_id}_chunk{idx}\n{sub_seq}\n+\n{sub_qual}\n")
EOF
    rm -f "$raw_fq"
done

rmdir "$TMP_FQ_DIR" 2>/dev/null || true
echo "[+] All samples converted and chunked into: ${OUT_DIR}"
