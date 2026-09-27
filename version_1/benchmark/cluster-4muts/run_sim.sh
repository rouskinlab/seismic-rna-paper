#!/bin/bash
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH -c 20
#SBATCH --mem 20G

set -eux -o pipefail

CPUS=20
python run_sim.py $CPUS
python generate_mask_files.py rf-mask.txt sm-primers.txt all-rnas.fa regions.csv

