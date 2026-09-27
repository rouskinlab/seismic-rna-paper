#!/bin/bash
#SBATCH -p short
#SBATCH -t 1:00:00
#SBATCH -c 6
#SBATCH --mem 12G

# Process the Manfredonia et al. (2020) dataset

set -eux -o pipefail

CPUS=6
REF=CP009262.1
FQ_FULL_DIR=../data/olson-2022
FA_FILE=$FQ_FULL_DIR/$REF.fa
OUT_DIR=out-olson-2022
POOL=pool
REP1=SRR15560843
REP2=SRR15560851
MUTDIST=mutdist

# Non-default options (all steps):
# --num-cpus $CPUS: On a cluster, the default, os.cpu_count(), may return the total number of cores in the node, not the number requested using -c.

seismic align --num-cpus $CPUS --no-fastp-detect-adapter-for-pe --bt2-mixed -x $FQ_FULL_DIR -o $OUT_DIR $FA_FILE
seismic relate --num-cpus $CPUS -o $OUT_DIR $FA_FILE $OUT_DIR
seismic graph scatter --num-cpus $CPUS -o $OUT_DIR $OUT_DIR/$REP1/relate/$REF $OUT_DIR/$REP2/relate/$REF
seismic pool --num-cpus $CPUS --pooled $POOL $OUT_DIR/$REP1/ $OUT_DIR/$REP2
# Non-default options:
# -b mutdist: Indicate that these results should be used only for seismic graph mutdist, not for an accurate mutational profile.
# --max-mask-iter 2: Because the only positions that are masked out are those with 0 coverage, no reads should be masked out on iteration 2, hence masking should be complete within 2 iterations.
# --keep-gu: These data are from the SHAPE reagent NAI, which modifies all four nucleotides.
# --keep-del --keep-ins: Indels are not typically counted because doing so increases the mutation rates of bases where indels can be unambiguous but not where indels can be ambiguous (i.e. repetitive sequences), causing a bias in the mutational profile. Calculating distances between mutations does not require an accurate mutational profile but should include any NAI-induced mutations, including indels.
# --mask-polya 0: Poly(A) sequences are not typically counted because their mutation rates are lower due to an artifact during reverse transcription, causing a bias in the mutational profile. Calculating distances between mutations does not require an accurate mutational profile but should include any DMS-induced mutations, including within poly(A) sequences.
# --min-mut-gap 0: To calculate the observer bias, reads must be kept regardless of the distance between mutations.
# --min-ninfo-pos 1: Low-coverage positions are typically masked out because there is too much uncertainty in their mutation rates. Calculating distances between mutations does not require an accurate mutational profile but should include any DMS-induced mutations, including at low-coverage positions.
# --no-mask-read-table: Calculating the table of relationships per read is not necessary for seismic graph mutdist and causes a crash from needing >32G memory.
seismic -vvvv mask --num-cpus $CPUS -b $MUTDIST -i olson-primers.csv --max-mask-iter 2 --keep-gu --keep-del --keep-ins --mask-polya 0 --min-mut-gap 0 --min-ninfo-pos 1 --no-mask-read-table $OUT_DIR/*$POOL
seismic -vvvv mask --num-cpus $CPUS -b ${MUTDIST}-sub -i olson-primers.csv --max-mask-iter 2 --keep-gu --mask-polya 0 --min-mut-gap 0 --min-ninfo-pos 1 --no-mask-read-table $OUT_DIR/*$POOL
seismic -vvvv graph profile --num-cpus $CPUS -r n --use-count $OUT_DIR/*$POOL/mask_${MUTDIST}*/$REF
seismic -vvvv graph profile --num-cpus $CPUS $OUT_DIR/*$POOL/mask_${MUTDIST}*/$REF
seismic -vvvv graph mutdist --num-cpus $CPUS $OUT_DIR/*$POOL/mask_${MUTDIST}*/$REF

