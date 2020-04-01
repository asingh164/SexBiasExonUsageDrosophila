#!/bin/bash

# This is a script to run FastQC on the fastq files for the DGRP lines sequenced

# Move to directory containing fastq files
cd ~/plas1/amardeep.singh/RNA.Seq.Data/Fastq.files

# Run fastqc on files. I set the number of threads to 30
fastqc -t 30 -o ~/plas1/amardeep.singh/RNA.Seq.Data/Fastq.files/fastqc.output.files




##
