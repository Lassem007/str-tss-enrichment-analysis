# Input data

This folder documents the external input files used for the analysis. The large reference files are not included in this repository because they are public datasets and can be downloaded from their original sources.

Before running the workflow, place the following files in the project root directory:

```text
gencode.v49.basic.annotation.gtf.gz
hg38.hipstr_reference.bed.gz

```

Already decompressed files can also be used:

```text
gencode.v49.basic.annotation.gtf
hg38.hipstr_reference.bed
```

The GENCODE annotation is used to extract transcript-level TSS positions. The HipSTR reference file provides the STR coordinates. The UCSC hg38 chromosome-size file is downloaded automatically by the workflow and is used for chromosome boundaries.
