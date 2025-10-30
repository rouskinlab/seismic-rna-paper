#!/bin/bash
#SBATCH -p short
#SBATCH -t 12:00:00
#SBATCH -c 1
#SBATCH --mem 6G

set -eux -o pipefail

# Define general parameters.
CPUS=1
ABSWINLEN=100
OUT=$PWD/out-draco
mkdir -p $OUT

# Count the samples and references.
SIM=sim
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
    TMP=tmp-draco-$SAMPLE_REF
    TIME_FILE_MAP=$OUT/time_draco-map_${SAMPLE_REF}.txt
    TIME_FILE_COU=$OUT/time_draco-count_${SAMPLE_REF}.txt
    TIME_FILE_CLU=$OUT/time_draco-cluster-${ABSWINLEN}_${SAMPLE_REF}.txt
    TIME_FILE_JAS=$OUT/time_draco-json2rc-${ABSWINLEN}_${SAMPLE_REF}.txt
    TIME_FILE_NOR=$OUT/time_draco-norm-${ABSWINLEN}_${SAMPLE_REF}.txt
    
    # Determine the expected output files.
    OUT_SAMPLE=$OUT/$SAMPLE
    OUT_REF=$OUT_SAMPLE/$REF

    OUT_MAP=$OUT_REF/rf_map
    BAM=$OUT_MAP/${REF}.bam

    OUT_COUNT=$OUT_REF/rf_count
    RC_COUNT=$OUT_COUNT/${REF}.rc
    MM=$OUT_COUNT/${REF}.mm
    
    OUT_DRACO=$OUT_REF/draco-$ABSWINLEN
    JSON=$OUT_DRACO/${REF}.json

    OUT_JSON2RC=$OUT_REF/rf_json2rc-$ABSWINLEN
    RC_DRACO=$OUT_JSON2RC/${REF}.rc

    OUT_NORM=$OUT_REF/rf_norm-$ABSWINLEN

    # Check if this script was assigned one task.
    if [[ -v TASK ]]; then
        # If so, then continue until the current task (T) equals the assigned task (TASK).
        if [ $T -ne $TASK ]; then
            continue
        fi

        mkdir $TMP
        mkdir -p $OUT_SAMPLE $OUT_REF

        # rf-map
        if [[ ! -f $BAM ]]; then
            FQ1=$SAMPLES_DIR/$SAMPLE/${REF}_R1.fq
            FQ2=$SAMPLES_DIR/$SAMPLE/${REF}_R2.fq
            FQ_PFX=$SAMPLES_DIR/$SAMPLE/${REF}
            FQ12=$FQ_PFX.fq
            $HOME/pear-0.9.11-linux-x86_64/bin/pear -f $FQ1 -r $FQ2 -o $FQ_PFX
            mv $FQ_PFX.assembled.fastq $FQ12
            rm $FQ_PFX.discarded.fastq $FQ_PFX.unassembled.forward.fastq $FQ_PFX.unassembled.reverse.fastq
            INDEX=$TMP/$REF
            bowtie2-build --quiet $REFS_DIR/$REF_FILE $INDEX
            rm -rf $OUT_MAP
            /usr/bin/time -f "%e" -o $TIME_FILE_MAP rf-map -p 1 --working-threads $CPUS --bowtie2 --bowtie-softclip -t $TMP/rf_map -o $OUT_MAP --bowtie-index $INDEX $FQ12
            rm $FQ12 ${INDEX}.*.bt2
        fi

        # rf-count
        if [[ ! -f $MM ]]; then
            rm -rf $OUT_COUNT
            /usr/bin/time -f "%e" -o $TIME_FILE_COU rf-count -p 1 --working-threads $CPUS --fast -m --mutation-map -o $OUT_COUNT -f $REFS_DIR/$REF_FILE $BAM
        fi

        rmdir $TMP

        # draco
        if [[ ! -f $JSON ]]; then
            rm -rf $OUT_DRACO
            mkdir $OUT_DRACO
            # Write to a temporary file, then rename so that if DRACO aborts in the middle of processing, the output file is not confused for a finished output.
            JSON_TMP=${JSON}.tmp
            rm -f $JSON_TMP
            /usr/bin/time -f "%e" -o $TIME_FILE_CLU draco --processors $CPUS --absWinLen $ABSWINLEN --reportNonInformative --allNonInformativeToOne --mm $MM --output $JSON_TMP
            mv $JSON_TMP $JSON
        fi

        # rf-json2rc
        if [[ ! -f $RC_DRACO ]]; then
            rm -rf $OUT_JSON2RC
            /usr/bin/time -f "%e" -o $TIME_FILE_JAS rf-json2rc --median-pre-cov 1 --median-cov 1 --min-confs 1 -j $JSON -r $RC_COUNT -o $OUT_JSON2RC
        fi

        # rf-norm
        if compgen -G "$OUT_NORM/*.xml" > /dev/null; then
            echo "Already finished $OUT_NORM"
        else
            rm -rf $OUT_NORM
            /usr/bin/time -f "%e" -o $TIME_FILE_NOR rf-norm -p $CPUS --scoring-method 4 --raw --reactive-bases AC -D 6 -t $RC_DRACO -o $OUT_NORM
        fi
    else
        # Otherwise, submit a job for this task if any output files do not exist.
        if compgen -G "$OUT_NORM/*.xml" > /dev/null; then
            echo "Already finished $OUT_NORM"
        else
            sbatch --export=ALL,TASK="$T" run_draco.sh
        fi
    fi
done
