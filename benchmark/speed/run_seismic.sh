#!/bin/bash
#SBATCH -p short
#SBATCH -t 60:00
#SBATCH -c 1
#SBATCH --mem 2G

set -eux -o pipefail

# Define general parameters.
CPUS=1
OUT=$PWD/out-seismic-$CPUS-cpu
mkdir -p $OUT

# Count the samples and references.
SIM=sim
REGIONS=regions.csv
SAMPLES_DIR=$SIM/samples
SAMPLE_DIRS=( $(ls -d $SAMPLES_DIR/*-biased-200000) )
REFS_DIR=$SIM/refs
REF_FILES=( $(ls $REFS_DIR/rna-280-4-*.fa) )
NSAMP=${#SAMPLE_DIRS[@]}
NREF=${#REF_FILES[@]}
NTASK=$(( NSAMP * NREF ))

# Loop over all possible tasks.
for T in $(seq 0 $(( NTASK - 1 ))); do
    # Determine the parameters for this task.
    SAMPLE_IDX=$(( T / NREF ))
    SAMPLE_DIR=${SAMPLE_DIRS[$SAMPLE_IDX]}
    SAMPLE=$(basename $SAMPLE_DIR)
    REF_IDX=$(( T % NREF ))
    REF_FILE=$(basename ${REF_FILES[$REF_IDX]})
    REF="${REF_FILE%.*}"
    SAMPLE_REF=${SAMPLE}_${REF}
    TIME_FILE_A=$OUT/time_seismic-align_${SAMPLE_REF}.txt
    TIME_FILE_R=$OUT/time_seismic-relate_${SAMPLE_REF}.txt
    TIME_FILE_M=$OUT/time_seismic-mask_${SAMPLE_REF}.txt
    TIME_FILE_C=$OUT/time_seismic-cluster_${SAMPLE_REF}.txt

    # Check if this script was assigned one task.
    if [[ -v TASK ]]; then
        # If so, then continue until the current task (T) equals the assigned task (TASK).
        if [ $T -ne $TASK ]; then
            continue
        fi

        echo "TASK $TASK, SAMPLE $SAMPLE, REF $REF"

        # Run SEISMIC.
        if [[ ! -f $TIME_FILE_A ]]; then
            /usr/bin/time -f "%e" -o $TIME_FILE_A seismic -vv align --num-cpus $CPUS -X $SAMPLES_DIR/$SAMPLE/${REF}_R1.fq -X $SAMPLES_DIR/$SAMPLE/${REF}_R2.fq --no-fastp-detect-adapter-for-pe --fastp-adapter-1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --fastp-adapter-2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --min-mapq 1 -o $OUT $REFS_DIR/$REF_FILE --force
        fi
        if [[ ! -f $TIME_FILE_R ]]; then
            /usr/bin/time -f "%e" -o $TIME_FILE_R seismic -vv relate --num-cpus $CPUS --min-mapq 1 -o $OUT $REFS_DIR/$REF_FILE $OUT/$SAMPLE/align/${REF}.bam --force
        fi
        if [[ ! -f $TIME_FILE_M ]]; then
            /usr/bin/time -f "%e" -o $TIME_FILE_M seismic -vv mask --num-cpus $CPUS --mask-regions-file regions.csv --keep-del $OUT/$SAMPLE/relate/$REF --force
        fi
        if [[ ! -f $TIME_FILE_C ]]; then
            /usr/bin/time -f "%e" -o $TIME_FILE_C seismic -vv cluster --num-cpus $CPUS $OUT/$SAMPLE/mask/$REF --force
        fi
    else
        # Otherwise, submit a job for this task if any output files do not exist.
        if [[ ! -f $TIME_FILE_C ]]; then
            sbatch --export=ALL,TASK="$T" run_seismic.sh
        fi
    fi
done

