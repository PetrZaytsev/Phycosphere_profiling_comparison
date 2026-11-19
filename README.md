# Microalgal Phycosphere Cross-Platform Sequencing Benchmark

Reproducible bioinformatic workflows and statistical benchmarking scripts 
comparing short-read amplicon sequencing (**16S-NGS**, Illumina V4), full-length 
amplicon sequencing (**16S-ONT**, Oxford Nanopore via standard Emu pipeline), and 
assembly-free shotgun metagenomics (**WMS-ONT**, Super-Accurate basecalling) 
across microalgal phycosphere consortia.

To eliminate database-induced taxonomic discrepancies across platforms, all 
sequencing datasets were harmonized against a unified, custom-compiled 
**RefSeq+rrnDBv5.6** reference database.

---

## Prerequisites & Dependencies

* **QIIME 2** (`v2023.9` or higher)  for 16S-NGS quality filtering, denoising, and 
taxonomic classification
* **Metaxa2** (`v2.2.3`)  includes `metaxa2`, `metaxa2_c`, and `metaxa2_dbb`
* **HMMER3** (`v3.1b2` or higher) & **NCBI BLAST+** (`v2.9.0` or higher)  backend 
engines for rRNA extraction and classification
* **SAMtools** (`v1.10` or higher)  for converting SUP-basecalled BAM files to 
FASTQ
* **Python** (`v3.8+`) with packages: `pandas`, `numpy`, `scipy`, `matplotlib`, 
`seaborn`, `scikit-bio`, `matplotlib-venn`

---

## Repository Contents

### 1. Short-Read Amplicon Workflow (`16S-NGS`)
* `16S-NGS_01_qiime2_pipeline.sh`  Imports single-end Illumina R1 reads from two 
independent sequencing runs (handling variable read lengths: 130 bp and 220 bp), 
performs quality filtering (Phred score >= 30), and denoises reads into ASV 
feature tables and representative sequences.
* `16S-NGS_02_refseq_classify.sh`  Classifies denoised representative sequences 
against a custom-trained **RefSeq+rrnDBv5.6** Naive Bayes classifier (`*.qza`), 
collapses feature tables to the genus level, pools technical PCR replicates, and 
filters out host/plastid sequences.

### 2. Shotgun Metagenomic Workflow (`WMS-ONT`)
* `WMS_01_prepare_and_build_db.sh`  Normalizes NCBI **RefSeq+rrnDBv5.6** taxonomy 
into a strict 7-level hierarchy and compiles the custom Metaxa2 reference 
database.
* `WMS_02_convert_and_chunk_sup.sh`  Converts demultiplexed Super-Accurate (SUP) 
BAM files to FASTQ and splits ultra-long ONT reads (`> 90 kb`) into overlapping 50 
kb windows.
* `WMS_03_run_metaxa2_classification.sh`  Executes decoupled two-step prokaryotic 
16S rRNA extraction (`metaxa2 --classifications F`) and classification 
(`metaxa2_c`) against **RefSeq+rrnDBv5.6**.
* `WMS_04_generate_taxonomic_matrices.py`  Parses per-read taxonomic assignments, 
removes host/organelle noise, and outputs decontaminated abundance matrices across 
all taxonomic ranks.

### 3. Cross-Platform Statistical Benchmarking & Reference Matrices
* `NGS_TOTAL_CLEAN_bacterial_core_OI1-12_RefSeq.tsv`  Pooled Illumina V4 
genus-level read counts across technical replicates (`OI01_repl1``OI12_repl2`).
* `ONT_16S_TOTAL_matrix_all_samples.tsv`  Full-length ONT 16S relative abundances 
(`AS_1``AS_12`).
* `WMS_FINAL_abundance_matrix_RefSeq.txt`  Decontaminated ONT shotgun 16S rRNA 
read counts (`barcode01``barcode12`).
* `calculate_metrics.py`  Computes pairwise sample-wise Bray-Curtis dissimilarity, 
Mantel test correlations (Spearman, 999 permutations), and Procrustes disparity 
($M^2$) across core consortium genera (`Table3_congruence_metrics.tsv`).
* `plot_heatmap_frequency.py`  Generates a cross-platform vertical split heatmap 
comparing relative abundances of the top-8 core bacterial genera across biological 
replicates and experimental stages (`heatmap_frequency.pdf/.png`).
* `plot_heatmaps_congruence.py`  Generates a 2x2 grid of lower-triangular 
congruence heatmaps summarizing Bray-Curtis, Procrustes, and Mantel metrics 
(`heatmaps_congruence.pdf/.png`).
* `plot_venn_diagram.py`  Generates a 3-way Venn diagram illustrating genus-level 
taxonomic overlap across **WMS-ONT**, **16S-ONT**, and **16S-NGS** 
(`venn_diagram.pdf/.png`).

---

## Execution Guide

Grant execution permissions to all shell scripts before running:
```bash
chmod +x *.sh
```

### Stage A: Short-Read 16S-NGS Processing (QIIME 2)
```bash
# 1. Import, quality-filter, and denoise single-end R1 reads from Run 1 and Run 2
./16S-NGS_01_qiime2_pipeline.sh manifest_run1.txt manifest_run2.txt qiime2_output/

# 2. Classify representative sequences against custom RefSeq+rrnDBv5.6 classifier 
and aggregate replicates
cd qiime2_output/
../16S-NGS_02_refseq_classify.sh /path/to/custom_refseq_v5.6_classifier.qza 
refseq_results/
```

### Stage B: Shotgun WMS-ONT Processing (Metaxa2)
```bash
# 1. Compile custom 7-level database (pass directory prefix without /blast)
./WMS_01_prepare_and_build_db.sh -f raw_db/species_taxid.fasta -t 
raw_db/taxonomy.tsv -o compiled_db/refseq_rrndb_v5.6 -c 16

# 2. Convert SUP BAMs and chunk ultra-long reads (>90 kb)
./WMS_02_convert_and_chunk_sup.sh -i raw_sup_bams/ -o chunked_fastq/ -c 16

# 3. Extract and classify bacterial 16S rRNA reads
./WMS_03_run_metaxa2_classification.sh -i chunked_fastq/ -d 
compiled_db/refseq_rrndb_v5.6 -o metaxa_results/ -c 16

# 4. Assemble multi-rank decontaminated matrices
python3 WMS_04_generate_taxonomic_matrices.py -i 
metaxa_results/02_classified_taxonomy/ -o final_matrices/
```

### Stage C: Cross-Platform Comparison & Visualization
Run sequentially from the repository root where the three harmonized `.tsv`/`.txt` 
matrices are located:

```bash
python3 calculate_metrics.py
python3 plot_heatmap_frequency.py
python3 plot_heatmaps_congruence.py
python3 plot_venn_diagram.py
```

---

## Methodological Rationale & Technical Workarounds

1. **Variable-Length Single-End Denoising (`16S-NGS_01_qiime2_pipeline.sh`):**
   Illumina V4 forward reads (R1) from independent sequencing trials with 
differing read lengths (130 bp and 220 bp) were processed separately to prevent 
truncation artifacts and subsequently pooled per biological sample to eliminate 
run-specific batch effects.
2. **Strict 7-Level Taxonomy Normalization (`WMS_01_prepare_and_build_db.sh`):**
   Metaxa2 requires every reference sequence to adhere to a rigid 7-rank hierarchy 
(`Superkingdom;Phylum;Class;Order;Family;Genus;Species`). Missing intermediate 
ranks in NCBI taxonomy are padded with `<ParentRank>_unclassified` placeholders to 
prevent array out-of-bounds errors (`subscript -1`) in the Perl parser without 
discarding valid reference sequences.
3. **Sliding-Window Chunking of Ultra-Long Reads 
(`WMS_02_convert_and_chunk_sup.sh`):**
   SUP basecalling recovers genomic reads exceeding 100 kb, which triggers fatal 
sequence-length overflows in the underlying `HMMER3` engine (`p7_pipeline.c` 100K 
limit). Reads longer than 90,000 bp are partitioned into 50,000 bp windows with a 
45,000 bp step size (5,000 bp overlap), ensuring complete ~1,500 bp 16S rRNA genes 
are preserved intact across window boundaries.
4. **Decoupled Extraction & Classification 
(`WMS_03_run_metaxa2_classification.sh`):**
   Prokaryotic 16S rRNA sequences are isolated first via the native `metaxa2` 
wrapper (`--classifications F`) and subsequently classified via `metaxa2_c` 
directly against the compiled `/blast` index (`-g E`) with an explicit identity 
cutoff vector (`-T 0,60,70,75,85,90,97`).
5. **Non-Phylogenetic Beta-Diversity (`calculate_metrics.py`):**
   Bray-Curtis dissimilarity was selected over UniFrac because cross-platform 
benchmarking integrates short V4 amplicons (~250 bp), full-length 16S amplicons 
(~1,500 bp), and variable-length shotgun 16S fragments, avoiding phylogenetic 
tree-building artifacts caused by unequal sequence lengths.

