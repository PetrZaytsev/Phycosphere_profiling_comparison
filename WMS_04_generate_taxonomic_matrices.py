#!/usr/bin/env python3
"""
Step 4: Aggregate Metaxa2 classification outputs (*.taxonomy.txt) into clean,
decontaminated abundance matrices across all taxonomic ranks.

Usage:
    python3 04_generate_taxonomic_matrices.py -i <classified_taxonomy_dir> -o <matrices_out_dir>
"""

import argparse
import glob
import os
import pandas as pd


DEFAULT_NOISE_KEYWORDS = [
    "chloroplast",
    "plastid",
    "eukaryota",
    "scenedesmus",
    "tetradesmus",
    "acetonema",
    "stanieria",
    "unassigned",
]

RANKS_MAP = {
    "0_Full_Lineage": -1,
    "1_Phylum": 1,
    "2_Class": 2,
    "3_Order": 3,
    "4_Family": 4,
    "5_Genus": 5,
    "6_Species": 6,
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate multi-rank taxonomic abundance matrices from Metaxa2 outputs."
    )
    parser.add_argument(
        "-i", "--input-dir", required=True, help="Directory containing *.taxonomy.txt files"
    )
    parser.add_argument(
        "-o", "--output-dir", required=True, help="Output directory for TSV abundance matrices"
    )
    parser.add_argument(
        "--exclude",
        nargs="*",
        default=DEFAULT_NOISE_KEYWORDS,
        help="Case-insensitive substrings to filter out host/organelle noise",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    tax_files = sorted(glob.glob(os.path.join(args.input_dir, "*.taxonomy.txt")))
    if not tax_files:
        raise FileNotFoundError(f"No *.taxonomy.txt files found in {args.input_dir}")

    data_by_level = {lvl: {} for lvl in RANKS_MAP}
    noise_lower = [w.lower() for w in args.exclude]

    for fpath in tax_files:
        sample_name = os.path.basename(fpath).replace(".taxonomy.txt", "")
        for lvl in RANKS_MAP:
            data_by_level[lvl][sample_name] = {}

        if os.path.getsize(fpath) == 0:
            continue

        df = pd.read_csv(
            fpath, sep="\t", header=None, names=["id", "tax", "ident", "len", "score"]
        )

        for tax_path in df["tax"].dropna():
            tax_str = str(tax_path).strip()
            if not tax_str or any(w in tax_str.lower() for w in noise_lower):
                continue

            parts = [p.strip() for p in tax_str.split(";") if p.strip()]

            for lvl_name, idx in RANKS_MAP.items():
                if idx == -1:
                    taxon = tax_str
                elif len(parts) > idx:
                    taxon = parts[idx]
                else:
                    continue

                if taxon and taxon.lower() not in {"nan", "none", "null", "unknown", ""}:
                    data_by_level[lvl_name][sample_name][taxon] = (
                        data_by_level[lvl_name][sample_name].get(taxon, 0.0) + 1.0
                    )

    for lvl_name, samples_dict in data_by_level.items():
        df_out = pd.DataFrame(samples_dict).fillna(0.0)
        df_out = df_out.reindex(sorted(df_out.columns), axis=1)
        df_out["Total_Reads"] = df_out.sum(axis=1)
        df_out = df_out.sort_values(by="Total_Reads", ascending=False)

        out_path = os.path.join(args.output_dir, f"WMS_abundance_{lvl_name}.tsv")
        df_out.to_csv(out_path, sep="\t")
        print(f"[+] Saved {out_path} ({len(df_out)} unique taxa)")


if __name__ == "__main__":
    main()
