#!/bin/bash
#SBATCH -p medium
#SBATCH -t 5-0:00:00
#SBATCH -c 1
#SBATCH --mem 32G

set -eux -o pipefail

# Define general parameters.
DREEM_CODE=$HOME/dreem/code
DREEM_PY=Run_DREEM.py
DREEM_OUT=$PWD/out-dreem
mkdir -p $DREEM_OUT 

# Count the samples and references.
SIM=sim
REGIONS=regions.csv
SAMPLES_DIR=$SIM/samples
SAMPLE_DIRS=( $(ls -d $SAMPLES_DIR/*-biased-200000) )
REFS_DIR=$SIM/refs
REF_FILES=( $(ls $REFS_DIR/rna-140-*-*.fa) )
NSAMP=${#SAMPLE_DIRS[@]}
NREF=${#REF_FILES[@]}
NTASK=$(( NSAMP * NREF ))

# Get the array task ID and use it to determine which sample and reference to process.

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
    OUT=$DREEM_OUT/$SAMPLE_REF
    TMP=$PWD/tmp-dreem-$SAMPLE_REF

    # Determine the 5' and 3' ends.
    read END5 END3 <<< "$(awk -F',' -v REF=$REF 'BEGIN { found = 0 } $2 == REF { print $3, $4; found = 1; exit } END { if (!found) { split(REF, a, "-"); reflen = a[2]; print 1, reflen } }' $REGIONS)"
    
    # Determine the expected output files.
    LOG_FILE=$OUT/EM_Clustering/${SAMPLE_REF}_${END5}_${END3}/log.txt

    # Check if this script was assigned one task.
    if [[ -v TASK ]]; then
        # If so, then continue until the current task (T) equals the assigned task (TASK).
        if [ $T -ne $TASK ]; then
            continue
        fi

        echo "TASK $TASK, SAMPLE $SAMPLE, REF $REF"

        # Run DREEM.
        if [[ ! -f $LOG_FILE ]]; then
            mkdir -p $OUT
            mkdir $TMP
            ln $SAMPLES_DIR/$SAMPLE/${REF}_R1.fq $TMP/${SAMPLE}_mate1.fastq
            ln $SAMPLES_DIR/$SAMPLE/${REF}_R2.fq $TMP/${SAMPLE}_mate2.fastq
            ln $REFS_DIR/$REF_FILE $TMP/${REF}.fasta
            bowtie2-build --quiet $TMP/${REF}.fasta $TMP/$REF
            cd $DREEM_CODE
            python $DREEM_PY $TMP $OUT $SAMPLE $REF $END5 $END3 --fastq
            cd -
            rm -r $TMP
        fi
    else
        # Otherwise, submit a job for this task if any output files do not exist.
        if [[ ! -f $LOG_FILE ]]; then
            sbatch --export=ALL,TASK="$T" run_dreem.sh
        fi
    fi
done

