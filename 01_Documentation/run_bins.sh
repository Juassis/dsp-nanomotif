#!/bin/bash

for dir in barcode*; do
    BINFILENAME="${dir}/assembly.fasta"
    BINNAME="${dir}"
    OUT="contig_bin_${BINNAME}.tsv"

    if [ -f "$BINFILENAME" ]; then
        grep ">" "$BINFILENAME" | sed 's/>//' | awk -v bin="$BINNAME" '{print $1 "\t" bin}' > "$OUT"
    else
        echo "Warning: $BINFILENAME not found, skipping." >&2
    fi
done

