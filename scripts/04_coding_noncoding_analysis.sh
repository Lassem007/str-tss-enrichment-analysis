#!/usr/bin/env bash
set -euo pipefail

mkdir -p results intermediate

GTF="gencode.v49.basic.annotation.gtf"
STR_SORTED="intermediate/hg38.hipstr_reference.sorted.bed"
TSS="intermediate/tss.sorted.bed"

for f in "$GTF" "$STR_SORTED" "$TSS"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: Missing $f. Run previous scripts first and make sure input files are in project root." >&2
    exit 1
  fi
done

printf "Extracting CDS intervals from GENCODE GTF...\n"
awk 'BEGIN{FS=OFS="\t"}
$3=="CDS" {
  start=$4-1;
  end=$5;
  if(start<0) start=0;
  print $1, start, end;
}' "$GTF" \
| sort -k1,1 -k2,2n \
> intermediate/cds.sorted.bed

printf "Classifying STRs as coding or non-coding based on CDS overlap...\n"
bedtools intersect \
  -a "$STR_SORTED" \
  -b intermediate/cds.sorted.bed \
  -u -sorted \
> results/str_coding.bed

bedtools intersect \
  -a "$STR_SORTED" \
  -b intermediate/cds.sorted.bed \
  -v -sorted \
> results/str_noncoding.bed

printf "Calculating signed distance from coding/non-coding STRs to nearest TSS...\n"
bedtools closest \
  -a results/str_coding.bed \
  -b "$TSS" \
  -D b -sorted \
> results/coding_to_nearest_tss_signed.tsv

bedtools closest \
  -a results/str_noncoding.bed \
  -b "$TSS" \
  -D b -sorted \
> results/noncoding_to_nearest_tss_signed.tsv

{
  echo -e "metric\tvalue"
  echo -e "coding_str_loci\t$(wc -l < results/str_coding.bed)"
  echo -e "noncoding_str_loci\t$(wc -l < results/str_noncoding.bed)"
} > results/coding_noncoding_summary.tsv

cat results/coding_noncoding_summary.tsv
