#!/usr/bin/env bash

# STAR alignment pipeline
# This is a program that I wrote to align trimmed reads from

# STAR uses a lot of resources. The current settings (running in 30 jobs in parallel, each job using 10 threads = 300 threads + 32gb per job = 900gb of memory)

ls *.fastq | parallel -j 30 STAR --genomeDir /plas1/amardeep.singh/Flybase.Dmel.Genome.Release/STAR.genome.index/ --runThreadN 10 --readFilesIn {} --outFileNamePrefix /plas1/amardeep.singh/RNA.Seq.Data/fastq.trimmed.files/mapped.files/{}_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard


# Alignment quality control
# Quality control of RNAseq alignment done with


# Sorting BAM files by read name instead of coordinate for downstream analysis with DESeq
#This is a directory made to hold temprary files generated by samtools to sort the bam files
mkdir samtools.sort.out
ls *_Aligned.sortedByCoord.out.bam | parallel -j 30 samtools sort -n -@ 10 -T samtools.sort.out -o SAM.sorted.{} {}




#
