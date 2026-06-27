#!/usr/bin/env bash
set -euo pipefail

mkdir -p intermediate results figures

STR_GZ="hg38.hipstr_reference.bed.gz"
STR_BED="hg38.hipstr_reference.bed"
GTF_GZ="gencode.v49.basic.annotation.gtf.gz"
GTF="gencode.v49.basic.annotation.gtf"

if [[ ! -f "$STR_BED" ]]; then
  if [[ -f "$STR_GZ" ]]; then
    gunzip -c "$STR_GZ" > "$STR_BED"
  else
    echo "ERROR: Missing $STR_BED or $STR_GZ in project root" >&2
    exit 1
  fi
fi

if [[ ! -f "$GTF" ]]; then
  if [[ -f "$GTF_GZ" ]]; then
    gunzip -c "$GTF_GZ" > "$GTF"
  else
    echo "ERROR: Missing $GTF or $GTF_GZ in project root" >&2
    exit 1
  fi
fi

printf "Inspecting input files...\n"
head "$STR_BED" > results/input_preview_str.txt
head "$GTF" > results/input_preview_gtf.txt

printf "Sorting STR BED file...\n"
sort -k1,1 -k2,2n "$STR_BED" > intermediate/hg38.hipstr_reference.sorted.bed

printf "Extracting transcript-level TSS positions in a strand-aware manner...\n"
awk 'BEGIN{FS=OFS="\t"}
$3=="transcript" {
  tid="."; gid="."; gname=".";
  match($9,/transcript_id "([^"]+)"/,a); if(a[1]!="") tid=a[1];
  match($9,/gene_id "([^"]+)"/,b); if(b[1]!="") gid=b[1];
  match($9,/gene_name "([^"]+)"/,c); if(c[1]!="") gname=c[1];

  if($7=="+") {
    start=$4-1; end=$4;
  } else if($7=="-") {
    start=$5-1; end=$5;
  } else {
    next;
  }

  if(start<0) start=0;
  print $1, start, end, tid "|" gname "|" gid, ".", $7;
}' "$GTF" \
| sort -k1,1 -k2,2n \
> intermediate/tss.sorted.bed

printf "Downloading hg38 chromosome sizes...\n"
wget -O intermediate/hg38.chrom.sizes \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

printf "Restricting chromosome sizes to chromosomes present in STR file...\n"
cut -f1 intermediate/hg38.hipstr_reference.sorted.bed | sort -u > intermediate/str.chroms.txt
grep -Fwf intermediate/str.chroms.txt intermediate/hg38.chrom.sizes > intermediate/hg38.str.chrom.sizes

printf "Generating +/- 2 kb TSS windows with chromosome-boundary clipping...\n"
bedtools slop \
  -i intermediate/tss.sorted.bed \
  -g intermediate/hg38.chrom.sizes \
  -b 2000 \
| sort -k1,1 -k2,2n \
> intermediate/tss_pm2kb.sorted.bed

{
  echo -e "metric\tvalue"
  echo -e "transcript_level_tss_records\t$(wc -l < intermediate/tss.sorted.bed)"
  echo -e "unique_genomic_tss_positions\t$(cut -f1-3 intermediate/tss.sorted.bed | sort -u | wc -l)"
  echo -e "tss_pm2kb_windows\t$(wc -l < intermediate/tss_pm2kb.sorted.bed)"
} > results/input_tss_summary.tsv

cat results/input_tss_summary.tsv
