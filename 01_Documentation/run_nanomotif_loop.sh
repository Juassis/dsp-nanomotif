#!/bin/bash

for barcode in {61..80}; do
    # Skip barcodes 62, 63, and 64
    if [[ "$barcode" == "62" || "$barcode" == "63" || "$barcode" == "64" ]]; then
        continue
    fi

    echo "Running nanomotif for barcode${barcode}"

    nanomotif motif_discovery \
        /home/azureuser/nanomotif/02_Working/barcode${barcode}/assembly.fasta \
        ../Pileup/barcode${barcode}_pileup.bed \
        ../Bins/contig_bin_barcode${barcode}.tsv \
        --out barcode${barcode}
done

