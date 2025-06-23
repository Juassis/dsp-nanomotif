#!/bin/bash

for f in barcode*/FASTq/*.fastq.gz; do
  if [[ -f "$f" ]]; then
    echo -n "$f: "
    zcat "$f" | paste - - - - | awk '{print length($2)}' | sort -n | tail -1
  fi
done
