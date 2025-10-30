#!/bin/bash
#SBATCH -p medium
#SBATCH -t 5-00:00:00
#SBATCH -c 20
#SBATCH --mem 64G

# Process the Morandi et al. (2021) dataset

set -eux -o pipefail

CPUS=20
REF=sars-cov-2
FQ_FULL_DIR=../../data/morandi-2021
FA_FILE=$FQ_FULL_DIR/$REF.fa
OUT_DIR=out-seismic
POOL=pooled

seismic -vv align --num-cpus 40 --sep-strands -x $FQ_FULL_DIR -o $OUT_DIR $FA_FILE
seismic -vv relate --num-cpus $CPUS --sep-strands --batch-size 4194304 -o $OUT_DIR $FA_FILE $OUT_DIR/*/align/$REF.bam
seismic pool --num-cpus $CPUS --pooled $POOL --no-relate-pos-table $OUT_DIR/*/relate/$REF/
# Non-default options:
# --mask-pos $REF 12 --mask-pos $REF 15101: These positions mutate (but not >98% to the same base) >50% of the time, which are probably (mostly) not real DMS-induced mutations.
# --no-mask-read-table: Calculating the table of relationships per read is not necessary and causes a crash from needing >32G memory.
#seismic -vvvv ensembles --num-cpus $CPUS --mask-pos $REF 12 --mask-pos $REF 15101 --no-mask-read-table --no-jackpot $OUT_DIR/$POOL/relate/
seismic -vvvv ensembles --num-cpus $CPUS -b fdr01 --mask-pos $REF 12 --mask-pos $REF 15101 --no-mask-read-table --pair-fdr 0.01 --no-jackpot $OUT_DIR/$POOL/relate/
seismic -vvvv ensembles --num-cpus $CPUS -b fdr001 --mask-pos $REF 12 --mask-pos $REF 15101 --no-mask-read-table --pair-fdr 0.001 --no-jackpot $OUT_DIR/$POOL/relate/
seismic -vv fold --num-cpus $CPUS --fold-table --fold-mfe $OUT_DIR/$POOL/cluster*/
seismic -vv graph profile --num-cpus $CPUS $OUT_DIR/$POOL/cluster*/
seismic -vv graph aucroll --num-cpus $CPUS --fold-table $OUT_DIR/$POOL/cluster*/

rm -r log/

