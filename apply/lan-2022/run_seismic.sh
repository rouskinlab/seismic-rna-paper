#!/bin/bash
#SBATCH -p medium
#SBATCH -t 5-00:00:00
#SBATCH -c 20
#SBATCH --mem 64G

# Process the Lan et al. (2022) dataset

set -eux -o pipefail

CPUS=20
REF=sars-cov-2
FQ_FULL_DIR=../../data/lan-2022
FA_FILE=$FQ_FULL_DIR/$REF.fa
OUT_DIR=out-seismic
CTRL=SRR12153162
REPS=SRR139589
POOL=pooled

seismic -vvvv align --num-cpus 40 --min-mapq 20 --sep-strands --bt2-no-un -x $FQ_FULL_DIR -o $OUT_DIR $FA_FILE
seismic -vvv relate --num-cpus $CPUS --min-mapq 20 --batch-size 1048576 -o $OUT_DIR $FA_FILE $OUT_DIR/*/align/$REF.bam
seismic -v list --min-ninfo-pos 1 --max-fmut-pos 0.07 $OUT_DIR/$CTRL
seismic -vv pool --num-cpus $CPUS --pooled $POOL $OUT_DIR/$REPS*
seismic -v graph profile --num-cpus $CPUS $OUT_DIR/$POOL/relate/
seismic -v graph profile --num-cpus $CPUS -racgt $OUT_DIR/$POOL/relate/
seismic -vvv ensembles --num-cpus $CPUS --mask-pos-file $OUT_DIR/$CTRL/relate/$REF/relate-position-list.csv --no-mask-read-table --no-jackpot $OUT_DIR/$POOL
seismic -vvv fold --num-cpus $CPUS --fold-table --fold-mfe $OUT_DIR/$POOL/cluster/
seismic -v graph profile --num-cpus $CPUS $OUT_DIR/$POOL/cluster/
seismic -v graph aucroll --num-cpus $CPUS --fold-table $OUT_DIR/$POOL/cluster/

rm -r log/

