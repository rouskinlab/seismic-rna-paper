#!/bin/bash
#SBATCH -c 8

# Download the single-end FASTQ files from the NCBI Short Read Archive. 

set -eux -o pipefail

for RUN in $@; do
    if [[ ! ( -f "${RUN}.fastq" || -f "${RUN}.fastq.gz" ) ]]; then
        prefetch $RUN --max-size=30g
        fasterq-dump -v --threads=8 --mem=1024MB $RUN
    fi
    if [[ -d $RUN ]]; then
        rm -r $RUN
    fi
    FILE=${RUN}.fastq
    if [[ -f $FILE ]]; then
        bgzip -@ 8 $FILE
    fi
done

