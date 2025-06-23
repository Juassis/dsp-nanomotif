#!/bin/bash

# Set paths
RUN1="/mnt/nanomotif/Sharing/NP260_run1_Peter/pass"
RUN3="/mnt/nanomotif/Sharing/NP260_run3_Peter/no_sample_id/20250307_1539_MD-101172_FAY38103_cb646fcb"
OUT="/mnt/nanomotif/Sharing/Merged"

# Make output dir
mkdir -p "$OUT"

# Get all unique barcodes from both FASTq/fastq_pass and BAM/bam_pass
barcodes=$(find "$RUN1/FASTq" "$RUN1/BAM" "$RUN3/fastq_pass" "$RUN3/bam_pass" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort -u)

# Loop through each barcode
for bc in $barcodes; do
  echo "Processing $bc..."

  # Create target dirs
  mkdir -p "$OUT/$bc/FASTq"
  mkdir -p "$OUT/$bc/BAM"

  # Copy FASTq files
  find "$RUN1/FASTq/$bc" "$RUN3/fastq_pass/$bc" -type f -name "*.fastq.gz" -exec cp {} "$OUT/$bc/FASTq/" \; 2>/dev/null

  # Copy BAM files
  find "$RUN1/BAM/$bc" "$RUN3/bam_pass/$bc" -type f -name "*.bam" -exec cp {} "$OUT/$bc/BAM/" \; 2>/dev/null
done

