#!/usr/bin/env bash
set -euo pipefail

bash scripts/00_check_environment.sh
bash scripts/01_prepare_inputs_and_tss_windows.sh
bash scripts/02_overlap_and_distance_analysis.sh
bash scripts/03_permutation_enrichment.sh 100
bash scripts/04_coding_noncoding_analysis.sh
python3 scripts/05_make_figures.py
