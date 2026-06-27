#!/usr/bin/env bash
set -euo pipefail

mkdir -p results intermediate

STR_SORTED="intermediate/hg38.hipstr_reference.sorted.bed"
TSS="intermediate/tss.sorted.bed"
TSS_WIN="intermediate/tss_pm2kb.sorted.bed"

for f in "$STR_SORTED" "$TSS" "$TSS_WIN"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: Missing $f. Run scripts/01_prepare_inputs_and_tss_windows.sh first." >&2
    exit 1
  fi
done

printf "Finding unique STR loci overlapping +/- 2 kb TSS windows...\n"
bedtools intersect \
  -a "$STR_SORTED" \
  -b "$TSS_WIN" \
  -u -sorted \
> results/str_in_tss_pm2kb_unique.bed

printf "Finding STR loci outside +/- 2 kb TSS windows...\n"
bedtools intersect \
  -a "$STR_SORTED" \
  -b "$TSS_WIN" \
  -v -sorted \
> results/str_outside_tss_pm2kb.bed

printf "Saving full STR-TSS-window overlap table...\n"
bedtools intersect \
  -a "$STR_SORTED" \
  -b "$TSS_WIN" \
  -wa -wb -sorted \
> results/str_in_tss_pm2kb.tsv

printf "Counting STRs per TSS window...\n"
bedtools intersect \
  -a "$TSS_WIN" \
  -b "$STR_SORTED" \
  -c -sorted \
> results/tss_window_str_counts.tsv

printf "Calculating nearest TSS distance for each STR record...\n"
bedtools closest \
  -a "$STR_SORTED" \
  -b "$TSS" \
  -d -sorted \
> results/str_to_nearest_tss.tsv

TOTAL_STR=$(wc -l < "$STR_SORTED")
INSIDE=$(wc -l < results/str_in_tss_pm2kb_unique.bed)
OUTSIDE=$(wc -l < results/str_outside_tss_pm2kb.bed)
FRACTION=$(awk -v i="$INSIDE" -v t="$TOTAL_STR" 'BEGIN{printf "%.6f", i/t}')

{
  echo -e "metric\tvalue"
  echo -e "total_str_loci\t${TOTAL_STR}"
  echo -e "str_loci_inside_pm2kb_tss_windows\t${INSIDE}"
  echo -e "str_loci_outside_pm2kb_tss_windows\t${OUTSIDE}"
  echo -e "fraction_inside_pm2kb_tss_windows\t${FRACTION}"
} > results/str_overlap_summary.tsv

awk 'BEGIN{FS=OFS="\t"}
NR==1 { }
{
  sum += $NF;
  if ($NF > 0) with_str++;
  if ($NF > max) max = $NF;
}
END {
  print "metric", "value";
  print "tss_windows", NR;
  print "tss_windows_with_at_least_one_str", with_str;
  print "fraction_tss_windows_with_at_least_one_str", with_str/NR;
  print "total_str_tss_window_overlaps", sum;
  print "mean_strs_per_tss_window", sum/NR;
  print "maximum_strs_in_one_tss_window", max;
}' results/tss_window_str_counts.tsv > results/tss_window_str_count_summary.tsv

awk 'BEGIN{FS=OFS="\t"; c500=0; c1000=0; c2000=0}
{
  d=$NF;
  if(d<0)d=-d;
  if(d<=500)c500++;
  if(d<=1000)c1000++;
  if(d<=2000)c2000++;
}
END{
  print "distance_threshold", "str_associated_output_lines";
  print "<=500_bp", c500;
  print "<=1000_bp", c1000;
  print "<=2000_bp", c2000;
}' results/str_to_nearest_tss.tsv > results/distance_threshold_counts.tsv

cat results/str_overlap_summary.tsv
cat results/tss_window_str_count_summary.tsv
cat results/distance_threshold_counts.tsv
