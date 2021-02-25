#!/usr/bin/env bash

###########################
####  Program details:  ###
###########################

## Author: Amardeep Singh -- amardeep.singh[at]utoronto.ca
## RNAseq data obtained from reposity associated with Osada et al. Genetics DOI: 10.1534/genetics.117.201459/-/DC1.1

### Program Details ###
# This is a program that I wrote to deal with adapter content in the RNAseq fastq files from the DGRP identified as an issue after running fastqc
# NOTE: the path to trimmomatic has been added as an alias in my bashrc. Users will need to do the same or modify code.
# NOTE: This file uses 40 cores, may need to update based on resources available to you

mkdir ../fastq.trimmed.files

ls RAL-786*.fastq | nice -n 15 parallel -j 40 java -jar /plas1/amardeep.singh/apps/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 {} /plas1/amardeep.singh/RNA.Seq.Data/fastq.trimmed.files/trimmed.{} ILLUMINACLIP:/plas1/amardeep.singh/apps/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10

echo DONE

# END
