#!/bin/bash
#SBATCH -p short
#SBATCH -t 3:00:00
#SBATCH -c 1
#SBATCH --mem 8G

set -eux -o pipefail

# Define general parameters.
CPUS=1
DM_PY=$HOME/DanceMapper/DanceMapper.py
OUT=$PWD/out-dance-$CPUS-cpu
ALL_PRIMERS=sm-primers.txt
OUT_SM=$OUT/shapemapper
OUT_DM=$OUT/dancemapper
mkdir -p $OUT $OUT_SM $OUT_DM

# Count the samples and references.
SIM=sim
SAMPLES_DIR=$SIM/samples
SAMPLE_DIRS=( $(ls -d $SAMPLES_DIR/*-biased-200000) )
REFS_DIR=$SIM/refs
REF_FILES=( $(ls $REFS_DIR/rna-280-2-*) )
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
    TMP=tmp-dance-$SAMPLE_REF
    TIME_FILE_SM=$OUT/time_dance-sm_${SAMPLE_REF}.txt
    TIME_FILE_DM=$OUT/time_dance-dm_${SAMPLE_REF}.txt
    
    # Determine the expected output files.
    PROFILE=$OUT_SM/${SAMPLE_REF}_profile.txt
    PARSED_MUT=$OUT_SM/${SAMPLE}_Modified_${REF}_parsed.mut
    DM_PREFIX=$OUT_DM/$SAMPLE_REF
    REACTS=$DM_PREFIX-reactivities.txt

    # Check if this script was assigned one task.
    if [[ -v TASK ]]; then
        echo "TASK $TASK, SAMPLE $SAMPLE, REF $REF"
        # If so, then continue until the current task (T) equals the assigned task (TASK).
        if [ $T -ne $TASK ]; then
            continue
        fi

        # Run shapemapper.
        if [[ ! -f $PROFILE || ! -f $PARSED_MUT ]]; then
            rm -f $PARSED_MUT
            if grep -q "^>$REF$" $ALL_PRIMERS; then
                # This reference has primers.
                REF_PRIMERS=sm-primers-${SAMPLE_REF}.txt
                grep -A 1 "^>$REF$" $ALL_PRIMERS > $REF_PRIMERS
                /usr/bin/time -f "%e" -o $TIME_FILE_SM shapemapper --name $SAMPLE --nproc $CPUS --min-depth 1000 --min-mapq 1 --min-mutation-separation 4 --modified --R1 $SAMPLES_DIR/$SAMPLE/${REF}_R1.fq --R2 $SAMPLES_DIR/$SAMPLE/${REF}_R2.fq --target $REFS_DIR/$REF_FILE --output-parsed-mutations --out $OUT_SM --temp $TMP --primers $REF_PRIMERS
                rm $REF_PRIMERS
            else
                # This reference does not have primers.
                /usr/bin/time -f "%e" -o $TIME_FILE_SM shapemapper --name $SAMPLE --nproc $CPUS --min-depth 1000 --min-mapq 1 --min-mutation-separation 4 --modified --R1 $SAMPLES_DIR/$SAMPLE/${REF}_R1.fq --R2 $SAMPLES_DIR/$SAMPLE/${REF}_R2.fq --target $REFS_DIR/$REF_FILE --output-parsed-mutations --out $OUT_SM --temp $TMP
            fi
        fi
        
        # Run DanceMapper.
        if [[ ! -f $REACTS ]]; then
            /usr/bin/time -f "%e" -o $TIME_FILE_DM python $DM_PY --mod $PARSED_MUT --prof $PROFILE --outputprefix $DM_PREFIX --fit --maxcomponents 5 --maskG --maskU --mincoverage 1
        fi
    else
        # Otherwise, submit a job for this task if any output files do not exist.
        if [[ ! -f $REACTS ]]; then
            sbatch --export=ALL,TASK="$T" run_dance.sh
        fi
    fi
done

