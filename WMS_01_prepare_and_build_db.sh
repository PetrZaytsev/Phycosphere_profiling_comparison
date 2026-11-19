#!/usr/bin/env bash
# ==============================================================================
# Step 1: Normalize NCBI RefSeq + rrnDB taxonomy to strict 7-level hierarchy
#         and compile a custom Metaxa2 database.
# Usage:  ./01_prepare_and_build_db.sh -f <input.fasta> -t <taxonomy.tsv> -o <out_db_dir> [-c <threads>]
# ==============================================================================
set -euo pipefail

THREADS=16
FASTA_IN=""
TAX_IN=""
DB_OUT=""

while getopts "f:t:o:c:h" opt; do
    case $opt in
        f) FASTA_IN="$OPTARG" ;;
        t) TAX_IN="$OPTARG" ;;
        o) DB_OUT="$OPTARG" ;;
        c) THREADS="$OPTARG" ;;
        h)
            echo "Usage: $0 -f <species_taxid.fasta> -t <taxonomy.tsv> -o <output_db_name> [-c <cpu_threads>]"
            exit 0
            ;;
        *) exit 1 ;;
    esac
done

if [[ -z "$FASTA_IN" || -z "$TAX_IN" || -z "$DB_OUT" ]]; then
    echo "Error: Missing required arguments (-f, -t, -o)." >&2
    exit 1
fi

WORK_DIR=$(dirname "$DB_OUT")
mkdir -p "$WORK_DIR"

CLEAN_FASTA="${WORK_DIR}/cleaned_sequences.fasta"
CLEAN_TAX="${WORK_DIR}/normalized_7level_taxonomy.tsv"

echo "[*] Normalizing FASTA headers and enforcing 7-level taxonomy..."
python3 - "$FASTA_IN" "$TAX_IN" "$CLEAN_FASTA" "$CLEAN_TAX" << 'EOF'
import sys
import pandas as pd

fasta_in, tax_in, fasta_out, tax_out = sys.argv[1:5]

df_tax = pd.read_csv(tax_in, sep='\t')
ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

tax2path = {}
for _, row in df_tax.iterrows():
    levels = []
    for r in ranks:
        val = str(row[r]).strip() if pd.notnull(row[r]) else ''
        if not val or val.lower() in {'nan', 'none', 'null'}:
            prev = levels[-1].split('_')[0] if levels else 'Bacteria'
            val = f"{prev}_unclassified"
        levels.append(val)
    tax2path[str(row['tax_id']).strip()] = ';'.join(levels)

with open(fasta_in, 'r') as fin, open(fasta_out, 'w') as fout, open(tax_out, 'w') as tout:
    counter = 0
    for line in fin:
        if line.startswith('>'):
            counter += 1
            raw_id = line[1:].strip().split()[0]
            tax_id = raw_id.split(':')[0].split('_')[0]
            new_id = f"{tax_id}_refseq_{counter}"
            lineage = tax2path.get(
                tax_id,
                "Bacteria;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified"
            )
            fout.write(f">{new_id}\n")
            tout.write(f"{new_id}\t{lineage}\n")
        else:
            fout.write(line)
EOF

echo "[*] Compiling custom Metaxa2 database using ${THREADS} threads..."
metaxa2_dbb \
    -i "$CLEAN_FASTA" \
    -t "$CLEAN_TAX" \
    -o "$DB_OUT" \
    --description "Custom NCBI RefSeq plus rrnDB 7-level database" \
    --mode conserved \
    --cpu "$THREADS"

echo "[+] Database compilation finished: ${DB_OUT}"
