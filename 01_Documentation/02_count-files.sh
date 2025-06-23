#!/bin/bash

# Set the base directory containing all barcodes
BASE_DIR="/mnt/nanomotif/Sharing/Merged/"

# Output header
echo -e "Barcode\tFASTQ_files_count"

# Loop through each barcode directory
for barcode in "$BASE_DIR"/barcode*; do
    if [ -d "$barcode/FASTq" ]; then
        count=$(find "$barcode/FASTq" -type f -name "*.fastq.gz" | wc -l)
        echo -e "$(basename "$barcode")\t$count"
    fi
done
