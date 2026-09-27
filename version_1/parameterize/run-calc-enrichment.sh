#!/bin/bash
#SBATCH -p short
#SBATCH -c 10
#SBATCH -t 2:00:00
#SBATCH --mem 128G

set -eux -o pipefail

python calc_enrichment.py out-de-lajarte-2024/ out-de-lajarte-2024/mu-enrich.csv out-de-lajarte-2024/mu-enrich.pdf 10

