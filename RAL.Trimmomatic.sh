#!/usr/bin/env bash

###########################
####  Program details:  ###
###########################

## Author: Amardeep Singh -- amardeep.singh[at]utoronto.ca
## RNAseq data obtained from reposity associated with Osada et al. Genetics DOI: 10.1534/genetics.117.201459/-/DC1.1

# This is a program that I wrote to deal with adapter content in the RNAseq fastq files from the DGRP identified as an issue after running fastqc
# NOTE: the path to trimmomatic has been added as an alias in my bashrc. Users will need to do the same or modify code.

trimmomatic SE -phred33



# END
