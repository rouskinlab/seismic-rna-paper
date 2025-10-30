#!/bin/bash

# The FASTQ files are in NCBI GEO: Series GSE153851
# Control (DMS-untreated):                       Sample GSM4656061, Experiment SRX8672926, Run SRR12153162
# Biological Replicate 1, Technical Replicate 1: Sample GSM4656059, Experiment SRX8672924, Run SRR13958908
# Biological Replicate 2, Technical Replicate 1: Sample GSM4656060, Experiment SRX8672925, Run SRR13958909
# Biological Replicate 2, Technical Replicate 2: Sample GSM4656060, Experiment SRX8672925, Run SRR13958910
../download.sh SRR12153162 SRR13958908 SRR13958909 SRR13958910

