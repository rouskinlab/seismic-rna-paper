#!/bin/bash
#SBATCH -p short
#SBATCH -t 6:00:00
#SBATCH -c 6
#SBATCH --mem 72G

# Process the Manfredonia et al. (2020) dataset

set -eux -o pipefail

CPUS=6
REF=sars-cov-2
FQ_FULL_DIR=../data/manfredonia-2020
FA_FILE=$FQ_FULL_DIR/$REF.fa
OUT_DIR=out-manfredonia-2020
POOL=pool
VITRO_NAI_1O=SRR11859187
VITRO_NAI_1E=SRR11859188
VITRO_NAI_1=vitro-nai-rep1
VITRO_NAI_2O=SRR11859189
VITRO_NAI_2E=SRR11859190
VITRO_NAI_2=vitro-nai-rep2
VITRO_NAI_POOL=vitro-nai-$POOL
VIVO_NAI_1=SRR13020886
VIVO_NAI_2=SRR13020887
VIVO_NAI_POOL=vivo-nai-$POOL
MUTDIST=mutdist

# Non-default options (all steps):
# --num-cpus $CPUS: On a cluster, the default, os.cpu_count(), may return the total number of cores in the node, not the number requested using -c.

# Single-end
for SAMPLE in $VITRO_NAI_1O $VITRO_NAI_1E $VITRO_NAI_2O $VITRO_NAI_2E; do
    seismic align --num-cpus $CPUS -z $FQ_FULL_DIR/${SAMPLE}.fastq.gz -o $OUT_DIR $FA_FILE
done
# Paired-end
for SAMPLE in $VIVO_NAI_1 $VIVO_NAI_2; do
    seismic align --num-cpus $CPUS -x $FQ_FULL_DIR/${SAMPLE}_1.fastq.gz -x $FQ_FULL_DIR/${SAMPLE}_2.fastq.gz -o $OUT_DIR $FA_FILE
done
seismic relate --num-cpus $CPUS --batch-size 4000000 -o $OUT_DIR $FA_FILE $OUT_DIR
seismic pool --num-cpus $CPUS --pooled $VITRO_NAI_1 $OUT_DIR/$VITRO_NAI_1O $OUT_DIR/$VITRO_NAI_1E
seismic pool --num-cpus $CPUS --pooled $VITRO_NAI_2 $OUT_DIR/$VITRO_NAI_2O $OUT_DIR/$VITRO_NAI_2E
seismic graph scatter --num-cpus $CPUS -o $OUT_DIR $OUT_DIR/$VITRO_NAI_1/relate/$REF $OUT_DIR/$VITRO_NAI_2/relate/$REF
seismic graph scatter --num-cpus $CPUS -o $OUT_DIR $OUT_DIR/$VIVO_NAI_1/relate/$REF $OUT_DIR/$VIVO_NAI_2/relate/$REF
seismic pool --num-cpus $CPUS --pooled $VITRO_NAI_POOL $OUT_DIR/$VITRO_NAI_1/ $OUT_DIR/$VITRO_NAI_2
seismic pool --num-cpus $CPUS --pooled $VIVO_NAI_POOL $OUT_DIR/$VIVO_NAI_1/ $OUT_DIR/$VIVO_NAI_2
# Non-default options:
# -b mutdist: Indicate that these results should be used only for seismic graph mutdist, not for an accurate mutational profile.
# -i manfredonia-primers.csv: The in vitro (but not in vivo) samples were RT-PCRed as amplicons with these primers.
# --min-ncov-read 110: Because the amplicons overlap by 45-101 nt, reads that come from adjacent amplicons must be excluded; thus require reads to have 110 nt in the amplicon -- more than the maximum that could come from an adjacent amplicon, but less than the read length.
# --max-mask-iter 2: Because the only positions that are masked out are those with 0 coverage, no reads should be masked out on iteration 2, hence masking should be complete within 2 iterations.
# --keep-gu: These data are from the SHAPE reagent NAI, which modifies all four nucleotides.
# --keep-del --keep-ins: Indels are not typically counted because doing so increases the mutation rates of bases where indels can be unambiguous but not where indels can be ambiguous (i.e. repetitive sequences), causing a bias in the mutational profile. Calculating distances between mutations does not require an accurate mutational profile but should include any NAI-induced mutations, including indels.
# --mask-polya 0: Poly(A) sequences are not typically counted because their mutation rates are lower due to an artifact during reverse transcription, causing a bias in the mutational profile. Calculating distances between mutations does not require an accurate mutational profile but should include any DMS-induced mutations, including within poly(A) sequences.
# --min-mut-gap 0: To calculate the observer bias, reads must be kept regardless of the distance between mutations.
# --min-ninfo-pos 1: Low-coverage positions are typically masked out because there is too much uncertainty in their mutation rates. Calculating distances between mutations does not require an accurate mutational profile but should include any DMS-induced mutations, including at low-coverage positions.
# --no-mask-read-table: Calculating the table of relationships per read is not necessary for seismic graph mutdist and causes a crash from needing >32G memory.
#seismic -vvvv mask --num-cpus $CPUS -b mutdist -i manfredonia-primers.csv --min-ncov-read 110 --max-mask-iter 2 --keep-gu --keep-del --keep-ins --mask-polya 0 --min-mut-gap 0 --min-ninfo-pos 1 --no-mask-pos-table --no-mask-read-table $OUT_DIR/$VITRO_NAI_POOL
#seismic join --num-cpus $CPUS --joined $JOINED --no-mask-read-table $OUT_DIR/$VITRO_NAI_POOL/mask_mutdist/$REF/amplicon-*
seismic -vvvv mask --num-cpus $CPUS -b $MUTDIST --max-mask-iter 2 --keep-gu --keep-del --keep-ins --mask-polya 0 --min-mut-gap 0 --min-ninfo-pos 1 --no-mask-read-table $OUT_DIR/*$POOL
seismic -vvvv mask --num-cpus $CPUS -b ${MUTDIST}-sub --max-mask-iter 2 --keep-gu --mask-polya 0 --min-mut-gap 0 --min-ninfo-pos 1 --no-mask-read-table $OUT_DIR/*$POOL
seismic -vvvv graph profile --num-cpus $CPUS -r n --use-count $OUT_DIR/*$POOL/mask_${MUTDIST}*/$REF/full
seismic -vvvv graph profile --num-cpus $CPUS $OUT_DIR/*$POOL/mask_${MUTDIST}*/$REF/full
seismic -vvvv graph mutdist --num-cpus $CPUS $OUT_DIR/*$POOL/mask_${MUTDIST}*/$REF/full

