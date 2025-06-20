#!/bin/bash
set -e

# Initialize conda and activate the environment
eval "$(conda shell.bash hook)"
conda activate bam2fastq

# Define base output directory
base_out="/home/azureuser/nanomotif/02_Working/Merging"

# Loop through barcodes, skipping 62–64
for i in 61 {65..80}; do
  input_dir="/mnt/nanomotif/Sharing/Merged/barcode${i}/BAM"
  out_dir="${base_out}/barcode${i}"
  mkdir -p "${out_dir}"

  merged_bam="${out_dir}/ALL_bam${i}.bam"
  sorted_bam="${out_dir}/ALL_bam${i}.sorted.bam"

  echo "Processing barcode${i}..."

  # Merge all BAMs in the input directory
  samtools merge -f "${merged_bam}" "${input_dir}"/*.bam

  # Sort the merged BAM
  samtools sort -o "${sorted_bam}" "${merged_bam}"

  # Index the sorted BAM
  samtools index "${sorted_bam}"

  # Optional: remove unsorted BAM
  # rm -f "${merged_bam}"

  echo "Done with barcode${i}"
done

