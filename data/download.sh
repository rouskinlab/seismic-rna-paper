#!/bin/bash
#SBATCH -c 8

# Download the FASTQ files from the NCBI Short Read Archive. 

set -eux -o pipefail

for RUN in $@; do
    if [[ ! ( ( -f "${RUN}_1.fastq" || -f "${RUN}_1.fastq.gz" ) && ( -f "${RUN}_2.fastq" || -f "${RUN}_2.fastq.gz" ) ) ]]; then
        prefetch $RUN --max-size=30g
        fasterq-dump -v --threads=8 --mem=1024MB $RUN
    fi
    if [[ -d $RUN ]]; then
        rm -r $RUN
    fi
    for MATE in 1 2; do
        FILE=${RUN}_${MATE}.fastq
        if [[ -f $FILE ]]; then
            bgzip -@ 8 $FILE
        fi
    done
done

