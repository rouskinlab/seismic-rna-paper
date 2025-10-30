#!/bin/bash
#SBATCH -p short
#SBATCH -t 12:00:00
#SBATCH -c 1
#SBATCH --mem 64G

set -eux -o pipefail

# Define general parameters.
CPUS=1
OUT=$PWD/out-seismic
mkdir -p $OUT

# Count the samples and references.
SIM=sim
REGIONS=regions.csv
SAMPLES_DIR=$SIM/samples
SAMPLE_DIRS=( $(ls -d $SAMPLES_DIR/sample-biased-*) )
REFS_DIR=$SIM/refs
REF_FILES=( $(ls $REFS_DIR/long-rna-??.fa) )
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
    TIME_FILE_E=$OUT/time_seismic-ensembles_${SAMPLE_REF}.txt
    TIME_FILE_P=$OUT/time_seismic-profile_${SAMPLE_REF}.txt

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
        if [[ ! -f $TIME_FILE_E ]]; then
            /usr/bin/time -f "%e" -o $TIME_FILE_E seismic -vv ensembles --num-cpus $CPUS --gap-mode insert $OUT/$SAMPLE/relate/$REF --force
        fi
        python sliding_cluster.py -vvvv --num-cpus $CPUS -b sliding -L 100 -O 0.5 --joined module $OUT/$SAMPLE/relate/$REF --no-verify-times
        seismic -v graph profile --num-cpus $CPUS --no-html $OUT/$SAMPLE/cluster_sliding/$REF
        if [[ ! -f $TIME_FILE_P ]]; then
            /usr/bin/time -f "%e" -o $TIME_FILE_P seismic -v graph profile --num-cpus $CPUS --no-html $OUT/$SAMPLE/cluster/$REF --force
        fi
    else
        # Otherwise, submit a job for this task if any output files do not exist.
        if [[ ! -f $TIME_FILE_P ]]; then
            sbatch --export=ALL,TASK="$T" run_seismic.sh
        fi
    fi
done

