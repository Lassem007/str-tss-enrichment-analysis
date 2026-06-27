#!/usr/bin/env bash
set -euo pipefail

mkdir -p results intermediate/permutations_100

STR_SORTED="intermediate/hg38.hipstr_reference.sorted.bed"
TSS_WIN="intermediate/tss_pm2kb.sorted.bed"
CHROMS="intermediate/hg38.str.chrom.sizes"
N=${1:-100}

for f in "$STR_SORTED" "$TSS_WIN" "$CHROMS"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: Missing $f. Run scripts/01_prepare_inputs_and_tss_windows.sh first." >&2
    exit 1
  fi
done

OBS=$(bedtools intersect \
  -a "$STR_SORTED" \
  -b "$TSS_WIN" \
  -u -sorted | wc -l)

echo "Observed overlap: $OBS"
rm -f results/permutation_overlap_counts_${N}.tsv

for i in $(seq 1 "$N")
do
  echo "Permutation $i / $N"
  bedtools shuffle \
    -i "$STR_SORTED" \
    -g "$CHROMS" \
    -chrom \
    -seed "$i" \
  | sort -k1,1 -k2,2n \
  > intermediate/permutations_100/tmp_shuffle.bed

  count=$(bedtools intersect \
    -a intermediate/permutations_100/tmp_shuffle.bed \
    -b "$TSS_WIN" \
    -u -sorted | wc -l)

  echo -e "${i}\t${count}" >> results/permutation_overlap_counts_${N}.tsv
done

python3 - "$OBS" "$N" <<'PY'
import sys
import pandas as pd

observed = int(sys.argv[1])
n = int(sys.argv[2])
path = f"results/permutation_overlap_counts_{n}.tsv"
perm = pd.read_csv(path, sep="\t", header=None, names=["permutation", "overlap_count"])
expected_mean = perm["overlap_count"].mean()
expected_sd = perm["overlap_count"].std(ddof=1)
fold_enrichment = observed / expected_mean
z_score = (observed - expected_mean) / expected_sd
n_extreme = (perm["overlap_count"] >= observed).sum()
empirical_p = (n_extreme + 1) / (len(perm) + 1)

summary = pd.DataFrame({
    "metric": [
        "observed_overlap",
        "expected_mean",
        "expected_sd",
        "fold_enrichment",
        "z_score",
        "empirical_p_value",
        "number_of_permutations",
    ],
    "value": [
        observed,
        expected_mean,
        expected_sd,
        fold_enrichment,
        z_score,
        empirical_p,
        len(perm),
    ],
})
out = f"results/enrichment_summary_{n}.tsv"
summary.to_csv(out, sep="\t", index=False)
print(summary)
PY
