#!/bin/bash

set -eux -o pipefail

# The FASTQ files are in NCBI: GEO GSE262014, BioProject PRJNA1090038
for SRR in {28398659..28398691}; do
    ../download.sh SRR${SRR}
done

# Make an alias for each file to describe the plate(s) and treatment(s).
fq_alias () {
    for MATE in 1 2; do
        SRR=$1
        PLATE_NAME=$2
        DMS_NAME=dms-$3
        ALIAS_DIR=$PLATE_NAME
        ALIAS_DIR=${PLATE_NAME}_small  #FIXME
        mkdir -p $ALIAS_DIR
        FQ_ALIAS=$ALIAS_DIR/${PLATE_NAME}_${DMS_NAME}_${MATE}.fastq.gz
        FQ_ALIAS=$ALIAS_DIR/${PLATE_NAME}_${DMS_NAME}_${MATE}.fastq  #FIXME
        if [[ ! -f $FQ_ALIAS ]]; then
            FQ_SRR=${SRR}_${MATE}.fastq.gz
            #ln $FQ_SRR $FQ_ALIAS  #FIXME
            zcat $FQ_SRR | head -n 1000000 > $FQ_ALIAS || [[ $? -eq 141 ]]
            echo "Wrote $FQ_ALIAS"
        fi
    done
}

fq_alias SRR28398659 mRNA_plate-8 ut
fq_alias SRR28398660 mRNA_plate-7 ut
fq_alias SRR28398661 mRNA_plate-6 ut
fq_alias SRR28398662 mRNA_plate-5 ut
fq_alias SRR28398662 mRNA_plate-12 b
fq_alias SRR28398663 mRNA_plate-4 ut
fq_alias SRR28398664 mRNA_plate-3 ut
fq_alias SRR28398665 mRNA_plate-2 ut
fq_alias SRR28398666 mRNA_plate-1 ut
fq_alias SRR28398667 mRNA_plate-12 ut
fq_alias SRR28398668 mRNA_plate-11 ut
fq_alias SRR28398669 mRNA_plate-10 ut
fq_alias SRR28398669 mRNA_plate-3 b
fq_alias SRR28398670 mRNA_plate-9 ut
fq_alias SRR28398671 mRNA_plate-9 a
fq_alias SRR28398672 mRNA_plate-8 a
fq_alias SRR28398673 mRNA_plate-7 a
fq_alias SRR28398674 mRNA_plate-5 b
fq_alias SRR28398674 mRNA_plate-6 b
fq_alias SRR28398675 mRNA_plate-5 a
fq_alias SRR28398676 pri-miRNA_plate-3 b
fq_alias SRR28398677 pri-miRNA_plate-3 a
fq_alias SRR28398678 mRNA_plate-3 a
fq_alias SRR28398679 mRNA_plate-2 b
fq_alias SRR28398680 pri-miRNA_plate-2 b
fq_alias SRR28398681 mRNA_plate-2 a
fq_alias SRR28398682 pri-miRNA_plate-2 a
fq_alias SRR28398683 mRNA_plate-1 b
fq_alias SRR28398684 pri-miRNA_plate-1 b
fq_alias SRR28398685 mRNA_plate-1 a
fq_alias SRR28398686 pri-miRNA_plate-1 a
fq_alias SRR28398687 mRNA_plate-12 a
fq_alias SRR28398688 mRNA_plate-11 b
fq_alias SRR28398689 mRNA_plate-11 a
fq_alias SRR28398690 mRNA_plate-6 a
fq_alias SRR28398690 mRNA_plate-10 b
fq_alias SRR28398691 mRNA_plate-4 b
fq_alias SRR28398691 mRNA_plate-10 a

