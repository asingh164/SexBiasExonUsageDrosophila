#!/usr/bin/env bash

###########################
####  Program details:  ###
###########################

## Author: Amardeep Singh -- amardeep.singh[at]utoronto.ca
## RNAseq data obtained from reposity associated with Osada et al. Genetics DOI: 10.1534/genetics.117.201459/-/DC1.1

### Program Details ###
#  This is a script to run FastQC on the fastq files for the DGRP lines that were sequenced  ####
#  Run this script in the folder containing the fastq files that you want to run through FastQC
#  After running FastQC it combines data from all runs using MultiQC
#  NOTE: I don't know if this is founded, but I have a feeling that FastQC works much quicker on uncompressed files
#        so I unzipped all files using bzip2 -df prior to running this bash script


####################################################################################################

# Decompress all .bz2 files
parallel -j 30 bzip2 ::: *.bz2

# First, we will create an output directory for fastqc outputs one directory above the current directory assuming it doesn't already exist
mkdir ../fastqc.output.files

# Run fastqc on files. I set the number of threads to 30
fastqc *.fastq -t 36 -o ../fastqc.output.files

# Run MultiQC to aggregate FastQC outputs for all samples
# Navigate over to fastqc output directory
cd ../fastqc.output.files

# Run MultiQC to aggregate fastqc outputs
multiqc .

## END
