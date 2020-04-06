#!/usr/bin/env bash

#########################################
###### STAR Genome Index Generator  #####
#########################################

## Author: Amardeep Singh -- amardeep.singh[at]utoronto.ca
## This program makes use of public releases of the Drosophila genome from Flybase (dmel_r6.32_FB2020_01)
## To retrieve the genome assembly file yourself:
# wget ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.32_FB2020_01/

# This is a program that I wrote to generate a genome index file for use in an RNAseq alignment
##

# Make a directory in the directory containing the genome files downloaded from Ensembl that will contain
#   the genome index that STAR generates

mkdir STAR.genome.index

# STAR command to generate a genome index file
# NOTE: Based on the genome size reported by STAR, it is suggested that I scale the --genomeSAindexNbases to 13 instead of 14 (Not sure why to be honest)
STAR --runThreadN 30 --runMode genomeGenerate --genomeDir STAR.genome.index --genomeFastaFiles /plas1/amardeep.singh/Flybase.Dmel.Genome.Release/uncompressed.files/dmel-all-chromosome-r6.32.fasta --sjdbGTFfile /plas1/amardeep.singh/Flybase.Dmel.Genome.Release/uncompressed.files/dmel-all-r6.32.gtf --sjdbOverhang 99 --genomeSAindexNbases 13

echo Done!












#
