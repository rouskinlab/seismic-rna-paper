#!/bin/bash
#SBATCH -p short
#SBATCH -t 1:00:00
#SBATCH -c 10
#SBATCH --mem 20G

set -eux -o pipefail

CPUS=10

OUT=out-seismic
FASTA=all-rnas.fa

# --no-fastp: Do not waste time running FASTP because all simulated reads are already high-quality.
# --min-mapq 1: For the 140 nt reference, using the default minimum mapping quality discards about 40% of reads, especially those with more mutations.
seismic -vv align --num-cpus $CPUS --no-fastp --min-mapq 1 -o $OUT -X sim/samples $FASTA
# --min-mapq 1: For the 140 nt reference, using the default minimum mapping quality discards about 40% of reads, especially those with more mutations.
seismic -vv relate --num-cpus $CPUS --min-mapq 1 --brotli-level 1 --no-relate-pos-table -o $OUT $FASTA $OUT/*/align*/
seismic -vv table --num-cpus $CPUS $OUT/sample-biased-200000/relate/
# --mask-regions-file: Exclude primers from the 140 nt and 280 nt references.
# --keep-del: Deletions are included in the simulations, so test how well SEISMIC-RNA can handle them.
# --no-mask-pos-table: The ensemble average mutation rates are not used while analyzing the results.
seismic -vv mask --num-cpus $CPUS --mask-regions-file regions.csv --keep-del --brotli-level 1 --no-mask-pos-table $OUT/*/relate*/
seismic -vv table --num-cpus $CPUS $OUT/sample-biased-200000/mask/
# --min-mut-gap 0: Turn off the dropout bias correction to analyze the effects of that bias and its correction.
seismic -vv mask --num-cpus $CPUS -b gap-0 --min-mut-gap 0 --mask-regions-file regions.csv --keep-del --brotli-level 1 --no-mask-pos-table $OUT/sample-*bias*/relate*/
# --no-jackpot: Do not waste time checking for jackpotting because the samples contain no PCR-like bias.
seismic -vv cluster --num-cpus $CPUS --no-jackpot --brotli-level 1 $OUT/*bias*/mask*/
# --jackpot: The jackpot samples do contain simulated jackpotting bias.
# --max-jackpot-quotient 1000000: Permit virtually any amount of jackpotting.
seismic -vv cluster --num-cpus $CPUS --jackpot --max-jackpot-quotient 1000000 --brotli-level 1 --no-cluster-pos-table --no-cluster-abundance-table $OUT/*jackpot*/mask*/
seismic -v graph histread --num-cpus $CPUS --use-count --no-html $OUT/*bias*/mask*/
seismic -v graph abundance --num-cpus $CPUS --no-html $OUT/*bias*/cluster*/
seismic -v graph profile --num-cpus $CPUS --no-html $OUT/*bias*/mask*/
seismic -v graph profile --num-cpus $CPUS --rels n --use-count --no-html $OUT/*bias*/relate*/
seismic -v graph profile --num-cpus $CPUS --no-html $OUT/*bias*/cluster*/
seismic -v graph mutdist --num-cpus $CPUS --no-html $OUT/*bias*/mask_gap-0/

rm -r log/

