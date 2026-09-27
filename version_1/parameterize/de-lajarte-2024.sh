#!/bin/bash
#SBATCH -p short
#SBATCH -t 1:00:00
#SBATCH -c 10
#SBATCH --mem 10G

# Process the de Lajarte et al. (2024) dataset

set -eux -o pipefail

CPUS=10
DATA_DIR=../data/de-lajarte-2024
RNANDRIA_FILE=$DATA_DIR/rnandria_structs.csv
USE_RNANDRIA=0
OUT_DIR=out-de-lajarte-2024

# Non-default options (all steps):
# --num-cpus $CPUS: On a cluster, the default, os.cpu_count(), may return the total number of cores in the node, not the number requested using -c (which is 20).

for PLATE in {1..12}; do
    PLATE_NAME=mRNA_plate-$PLATE
    FA_FILE=$DATA_DIR/${PLATE_NAME}.fa
    SAMPLE_PREFIX=${PLATE_NAME}_dms-
    seismic align --num-cpus $CPUS -x $DATA_DIR/$PLATE_NAME -o $OUT_DIR $FA_FILE
    seismic relate --num-cpus $CPUS -o $OUT_DIR $FA_FILE $OUT_DIR/${SAMPLE_PREFIX}*
    # List positions with excessive DMS reactivities (or insufficient coverage to determine) in the untreated sample.
    seismic -v list --num-cpus $CPUS --min-ninfo-pos 500 --max-fmut-pos 0.01 $OUT_DIR/${SAMPLE_PREFIX}ut
    # Graph profiles of the DMS-treated replicates and the untreated control.
    seismic -v graph profile --num-cpus $CPUS --no-html $OUT_DIR/${SAMPLE_PREFIX}*/relate
    seismic -v graph profile --num-cpus $CPUS -r n --use-count --no-html $OUT_DIR/${SAMPLE_PREFIX}*/relate
    # Compare the DMS-treated replicates A and B (if both exist).
    seismic -v graph scatter --num-cpus $CPUS --no-html -o $OUT_DIR $OUT_DIR/${SAMPLE_PREFIX}[ab]/relate
    # Pool the DMS-treated replicates A and B for references with untreated controls.
    POOLED_SAMPLE=${SAMPLE_PREFIX}pool
    for REF_DIR in $OUT_DIR/${SAMPLE_PREFIX}ut/relate/*; do
        REF=$(basename $REF_DIR)
        if [[ ! -f $OUT_DIR/$POOLED_SAMPLE/relate/$REF/relate-report.json ]]; then
            # Confirm that replicates A and B exist and have sufficient correlation.
            SCATTER=$OUT_DIR/${SAMPLE_PREFIX}a_VS_${SAMPLE_PREFIX}b/graph/$REF/full/scatter_all_m-ratio-q0.csv
            if [[ -f $SCATTER ]]; then
                SUFFICIENT=$(python check_pearson.py $SCATTER 0.9)
                if [[ $SUFFICIENT -eq 1 ]]; then
                    echo $SUFFICIENT
                    seismic pool --num-cpus $CPUS --pooled $POOLED_SAMPLE --no-relate-pos-table $OUT_DIR/${SAMPLE_PREFIX}[ab]/relate/$REF
                fi
            fi
        fi
    done
    if [[ -d $OUT_DIR/$POOLED_SAMPLE ]]; then
        # Non-default options:
        # -b mutdist: Indicate that these results should be used only for seismic graph mutdist, not for an accurate mutational profile.
        # --max-mask-iter 2: Because the only positions that are masked out are those with 0 coverage, no reads should be masked out on iteration 2, hence masking should be complete within 2 iterations.
        # --keep-del --keep-ins: Indels are not typically counted because doing so increases the mutation rates of bases where indels can be unambiguous but not where indels can be ambiguous (i.e. repetitive sequences), causing a bias in the mutational profile. Calculating distances between mutations does not require an accurate mutational profile but should include any DMS-induced mutations, including indels.
        # --mask-polya 0: Poly(A) sequences are not typically counted because their mutation rates are lower due to an artifact during reverse transcription, causing a bias in the mutational profile. Calculating distances between mutations does not require an accurate mutational profile but should include any DMS-induced mutations, including within poly(A) sequences.
        # --mask-pos-file $OUT_DIR/${PlATE_NAME}_dms-ut: Mask out these positions that were too highly mutated in the untreated control.
        # --min-mut-gap 0: To calculate the observer bias, reads must be kept regardless of the distance between mutations.
        # --min-ninfo-pos 1: Low-coverage positions are typically masked out because there is too much uncertainty in their mutation rates. Calculating distances between mutations does not require an accurate mutational profile but should include any DMS-induced mutations, including at low-coverage positions.
        # --no-mask-read-table: Calculating the table of relationships per read is not necessary for seismic graph mutdist and causes a crash from needing >32G memory.
        seismic -vv mask --num-cpus $CPUS -b mutdist --max-mask-iter 2 --keep-del --keep-ins --mask-polya 0 --mask-pos-file $OUT_DIR/${SAMPLE_PREFIX}ut --min-mut-gap 0 --min-ninfo-pos 1 --no-mask-read-table $OUT_DIR/$POOLED_SAMPLE
        # Non-default options:
        # -b profile: Indicate that these results should be used only for seismic graph profile, not for distances between mutations.
        # --keep-del: Include deletions so that the frequency of deletions can be estimated by seismic sim abstract.
        # --mask-pos-file $OUT_DIR/${PlATE_NAME}_dms-ut: Mask out these positions that were too highly mutated in the untreated control.
        # --no-mask-read-table: Calculating the table of relationships per read is not necessary for seismic graph mutdist and causes a crash from needing >32G memory.
        seismic -vv mask --num-cpus $CPUS -b profile --keep-del --mask-pos-file $OUT_DIR/${SAMPLE_PREFIX}ut --no-mask-read-table $OUT_DIR/$POOLED_SAMPLE
        if [[ $USE_RNANDRIA -eq 1 ]]; then
            # Create DB and CT files using the structure in the RNAndria database.
            FOLD_DIR=$OUT_DIR/$POOLED_SAMPLE/fold_rnandria
            mkdir -p $FOLD_DIR
            for REF_DIR in $OUT_DIR/$POOLED_SAMPLE/mask_profile/*; do
                REF=$(basename $REF_DIR)
                # Check if the reference occurs in the RNAndria database.
                COUNT=$(awk "/$REF/{n++} END{print n+0}" $RNANDRIA_FILE)
                if [[ $COUNT -eq 1 ]]; then
                    FOLD_REF_DIR=$FOLD_DIR/$REF
                    mkdir -p $FOLD_REF_DIR
                    FOLD_REF_REG_DIR=$FOLD_REF_DIR/full
                    mkdir -p $FOLD_REF_REG_DIR
                    DB_FILE=$FOLD_REF_REG_DIR/rnandria.db
                    if [[ ! -f $DB_FILE ]]; then
                        # Write the dot-bracket file line by line.
                        echo ">$REF" > $DB_FILE
                        # Use the sequence from the FASTA file (the RNAndria sequence may differ by several bases).
                        # Convert it from a DNA to an RNA sequence.
                        grep -A 1 "^>$REF$" $FA_FILE | tail -n 1 | sed "s/T/U/g" >> $DB_FILE
                        # Use the structure from the RNAndria database.
                        grep "^$REF" $RNANDRIA_FILE | cut -d "," -f 3 >> $DB_FILE
                    fi
                fi
            done
            seismic db2ct $FOLD_DIR
        else
            # Model the structures from the DMS data.
            seismic -vv fold --num-cpus $CPUS --fold-mfe $OUT_DIR/$POOLED_SAMPLE/mask_profile/*/full/mask-position-table.csv
        fi
    fi
done

seismic -v graph profile --num-cpus $CPUS --no-html $OUT_DIR/*pool/mask_profile
seismic -v graph mutdist --num-cpus $CPUS --no-html $OUT_DIR/*pool/mask_mutdist
bash run-calc-enrichment.sh

# Calculate the parameters for simulation from these mutation rates and structures.
seismic -v sim abstract --num-cpus $CPUS --struct-file $OUT_DIR $OUT_DIR/*/mask_profile/*/full

