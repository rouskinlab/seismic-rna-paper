#!/bin/bash
#SBATCH -p short
#SBATCH -t 6:00:00
#SBATCH -c 10
#SBATCH --mem 40G

set -eux -o pipefail

CPUS=10
python run_sim.py $CPUS

