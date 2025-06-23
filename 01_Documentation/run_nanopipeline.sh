#!/bin/bash
set -e

# Initialize conda and activate the environment
eval "$(conda shell.bash hook)"
conda activate bam2fastq

# Base paths
MERGING_DIR="/home/azureuser/nanomotif/02_Working/Merging"
ASSEMBLY_BASE="/home/azureuser/nanomotif/02_Working"
OUTPUT_DIR="/home/azureuser/nanomotif/02_Working/Nanomotif"

# Make output dir if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through each barcode directory in Merging
for BARCODE_DIR in "$MERGING_DIR"/barcode*; do
    BARCODE=$(basename "$BARCODE_DIR")  # e.g., barcode65
    echo "Processing $BARCODE..."

#    BAM_SORTED="$BARCODE_DIR/ALL_${BARCODE}.sorted.bam"
 BAM_SORTED="$BARCODE_DIR/ALL_bam${BARCODE#barcode}.sorted.bam"

 ASSEMBLY="$ASSEMBLY_BASE/$BARCODE/assembly.fasta"

    # Check if necessary files exist
    if [[ ! -f "$BAM_SORTED" ]]; then
        echo "Missing BAM file for $BARCODE: $BAM_SORTED"
        continue
    fi
    if [[ ! -f "$ASSEMBLY" ]]; then
        echo "Missing assembly.fasta for $BARCODE: $ASSEMBLY"
        continue
    fi

    # Define output paths
    MAPPED_BAM="$OUTPUT_DIR/${BARCODE}_mapping.bam"
    PILEUP_FILE="$OUTPUT_DIR/${BARCODE}_pileup.bed"

    # Run modcall2.sh (assumes it runs in place or is global)
#    modcall2.sh

    # Convert, align, sort
    samtools fastq -T MM,ML "$BAM_SORTED" | \
        minimap2 -ax map-ont -y "$ASSEMBLY" - | \
        samtools view -bS - | \
        samtools sort -o "$MAPPED_BAM"

    # Index
    samtools index "$MAPPED_BAM"

    # Pileup
    modkit pileup --only-tabs "$MAPPED_BAM" "$PILEUP_FILE"

    echo "Done with $BARCODE"
done

echo "All done!"

