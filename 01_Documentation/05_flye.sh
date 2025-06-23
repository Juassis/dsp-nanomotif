conda activate flye
for i in {69..80}; do
  flye \
    --nano-raw /mnt/nanomotif/Sharing/Merged/barcode${i}/FASTq/ALL_fastq_barcode${i}.fastq.gz \
    --out-dir barcode${i} \
    --threads 3 \
    --asm-coverage 50 \
    --genome-size 400m
done
