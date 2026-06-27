# STR-TSS enrichment analysis

This repository contains the scripts and output files for a genome-wide analysis of short tandem repeats (STRs) relative to transcript-level transcription start sites (TSS) in the human hg38/GRCh38 genome.

The workflow uses STR coordinates from the HipSTR hg38 reference catalog and transcript annotations from the GENCODE v49 basic annotation. Transcript-level TSS positions are extracted in a strand-aware way and expanded to ±2 kb windows. STR overlap, nearest-TSS distances, and chromosome-aware enrichment are calculated with BEDTools. Final summary tables and figures are generated with Python.

## Repository structure

```text
scripts/       Bash and Python scripts for the analysis
data/          Notes on required external input files
results/       Final summary tables
figures/       Final figures
intermediate/  Temporary files generated locally
```

## Requirements

The workflow was run in Ubuntu via Windows Subsystem for Linux (WSL).

Required tools:

* BEDTools
* awk/gawk
* gzip
* coreutils
* wget
* Python 3 with pandas, numpy and matplotlib

Example installation:

```bash
sudo apt update
sudo apt install -y bedtools gawk gzip unzip coreutils wget python3 python3-pip
pip3 install pandas numpy matplotlib
```

## Input files

The large reference input files are not included in this repository. Place the following files in the project root directory before running the workflow:

```text
gencode.v49.basic.annotation.gtf.gz
hg38.hipstr_reference.bed.gz
```

Already decompressed versions are also accepted:

```text
gencode.v49.basic.annotation.gtf
hg38.hipstr_reference.bed
```

The UCSC hg38 chromosome-size file is downloaded automatically by the workflow.

## Running the workflow

Run the complete workflow with:

```bash
bash scripts/run_all.sh
```

Or run the steps separately:

```bash
bash scripts/00_check_environment.sh
bash scripts/01_prepare_inputs_and_tss_windows.sh
bash scripts/02_overlap_and_distance_analysis.sh
bash scripts/03_permutation_enrichment.sh 100
bash scripts/04_coding_noncoding_analysis.sh
python3 scripts/05_make_figures.py
```

## Main outputs

Important result tables:

* `results/input_summary.tsv`
* `results/str_overlap_summary.tsv`
* `results/tss_window_str_count_summary.tsv`
* `results/distance_threshold_counts.tsv`
* `results/distance_bins_250bp.tsv`
* `results/permutation_overlap_counts_100.tsv`
* `results/enrichment_summary_100.tsv`
* `results/coding_noncoding_summary.tsv`
* `results/coding_vs_noncoding_tss_distance_100bp.tsv`

Important figures:

* `figures/figure_inside_vs_outside_tss_pm2kb.png`
* `figures/figure_distance_250bp_bins.png`
* `figures/figure_coding_vs_noncoding_tss_2kb_100bp_normalized.png`
* `figures/figure_observed_vs_expected_str_overlap_near_tss.png`

## Notes

TSS-proximal windows are generated with `bedtools slop` using the full UCSC hg38 chromosome-size file. For the permutation analysis, a chromosome-filtered version of the same file is used with `bedtools shuffle -chrom` so that randomized STR intervals remain on chromosomes present in the STR dataset.
