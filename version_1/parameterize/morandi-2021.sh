#!/bin/bash
#SBATCH -p short
#SBATCH -t 8:00:00
#SBATCH -c 8
#SBATCH --mem 32G

# Process the Morandi et al. (2021) dataset

set -eux -o pipefail

CPUS=8
REF=sars-cov-2
FQ_FULL_DIR=../data/morandi-2021
FA_FILE=$FQ_FULL_DIR/$REF.fa
OUT_DIR=out-morandi-2021
REP1=SRR12653367
REP2=SRR12653368
POOL=pooled

# Non-default options (all steps):
# --num-cpus 8: On a cluster, the default, os.cpu_count(), may return the total number of cores in the node, not the number requested using -c (which is 8).

seismic align --num-cpus $CPUS -x $FQ_FULL_DIR -o $OUT_DIR $FA_FILE
seismic relate --num-cpus $CPUS -o $OUT_DIR $FA_FILE $OUT_DIR
seismic graph scatter --num-cpus $CPUS -o $OUT_DIR $OUT_DIR/$REP1/relate/$REF $OUT_DIR/$REP2/relate/$REF
seismic pool --num-cpus $CPUS --pooled $POOL --no-relate-pos-table $OUT_DIR
# Non-default options:
# -b mutdist: Indicate that these results should be used only for seismic graph mutdist, not for an accurate mutational profile.
# --max-mask-iter 2: Because the only positions that are masked out are those with 0 coverage, no reads should be masked out on iteration 2, hence masking should be complete within 2 iterations.
# --keep-del --keep-ins: Indels are not typically counted because doing so increases the mutation rates of bases where indels can be unambiguous but not where indels can be ambiguous (i.e. repetitive sequences), causing a bias in the mutational profile. Calculating distances between mutations does not require an accurate mutational profile but should include any DMS-induced mutations, including indels.
# --mask-polya 0: Poly(A) sequences are not typically counted because their mutation rates are lower due to an artifact during reverse transcription, causing a bias in the mutational profile. Calculating distances between mutations does not require an accurate mutational profile but should include any DMS-induced mutations, including within poly(A) sequences.
# --mask-pos $REF 12 --mask-pos $REF 15101: These positions mutate (but not >98% to the same base) >50% of the time, which are probably (mostly) not real DMS-induced mutations.
# --min-mut-gap 0: To calculate the observer bias, reads must be kept regardless of the distance between mutations.
# --min-ninfo-pos 1: Low-coverage positions are typically masked out because there is too much uncertainty in their mutation rates. Calculating distances between mutations does not require an accurate mutational profile but should include any DMS-induced mutations, including at low-coverage positions.
# --no-mask-read-table: Calculating the table of relationships per read is not necessary for seismic graph mutdist and causes a crash from needing >32G memory.
seismic -vvvv mask --num-cpus $CPUS -b mutdist --max-mask-iter 2 --keep-del --keep-ins --mask-polya 0 --mask-pos $REF 12 --mask-pos $REF 15101 --min-mut-gap 0 --min-ninfo-pos 1 --no-mask-read-table $OUT_DIR/$POOL
seismic -vvvv graph profile --num-cpus $CPUS $OUT_DIR/$POOL/mask_mutdist
seismic -vvvv graph profile --num-cpus $CPUS -r n --use-count $OUT_DIR/$POOL/mask_mutdist
seismic -vvvv graph mutdist --num-cpus $CPUS $OUT_DIR/$POOL/mask_mutdist

