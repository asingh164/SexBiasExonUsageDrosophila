#!/bin/bash

# This is a script to run FastQC on the fastq files for the DGRP lines that were sequenced
# After running FastQC it combines data from all runs using MultiQC

# Move to directory containing fastq files
cd /plas1/amardeep.singh/RNA.Seq.Data/fastq.files

# Run fastqc on files. I set the number of threads to 30
fastqc RAL*.fastq.bz2 -t 36 -o /plas1/amardeep.singh/RNA.Seq.Data/fastqc.output.files




## END
