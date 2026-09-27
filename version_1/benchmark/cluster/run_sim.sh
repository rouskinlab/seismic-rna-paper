#!/bin/bash
#SBATCH -p short
#SBATCH -t 1:00:00
#SBATCH -c 16
#SBATCH --mem 16G

set -eux -o pipefail

CPUS=16
python run_sim.py $CPUS
python generate_mask_files.py rf-mask.txt sm-primers.txt all-rnas.fa regions.csv

