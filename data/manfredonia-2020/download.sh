#!/bin/bash

# The FASTQ files are in NCBI GEO: Series GSE151327
# NAI in vitro, odd  set, replicate 1: Sample GSM4572511, Experiment SRX8409268, Run SRR11859187
# NAI in vitro, even set, replicate 1: Sample GSM4572512, Experiment SRX8409269, Run SRR11859188
# NAI in vitro, odd  set, replicate 2: Sample GSM4572513, Experiment SRX8409270, Run SRR11859189
# NAI in vitro, even set, replicate 2: Sample GSM4572514, Experiment SRX8409271, Run SRR11859190
# NAI in vivo, replicate 1: Sample GSM4892124, Experiment SRX9471129, Run SRR13020886
# NAI in vivo, replicate 2: Sample GSM4892125, Experiment SRX9471130, Run SRR13020887

../download-single.sh SRR11859187 SRR11859188 SRR11859189 SRR11859190
../download.sh SRR13020886 SRR13020887

