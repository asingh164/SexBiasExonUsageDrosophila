#
#####################################################################
###### JunctionSeq Script for Differential Exon Usage Analysis  #####
#####################################################################

## Author: Amardeep Singh -- amardeep.singh[at]utoronto.ca
## This program makes use of public releases of the Drosophila genome from Ensembl (BDGP6.28)
## To retrieve the genome assembly file yourself:
# Fasta file
# wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-46/fasta/drosophila_melanogaster/dna//Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa.gz
# GTF file
# wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-46/fasta/drosophila_melanogaster/

### Script Details ###
# This script is written primarily in R, Bash commands are denoted by #---- BASH ---- and end with # ---- ----


# ---- Linux ----
# Make a flat gff file from the genome annotation provided by Ensembl
QoRTs makeFlatGff --stranded /plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/gtf/Drosophila_melanogaster.BDGP6.28.99.gtf.gz /plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/Dmel.JunctionSeq.flat.gff.gz

## Processing bam files
# First, run QoRTs to perform QC checks on aligned bam files. This needs to be done so that the size factor for normalization can be done
# NOTE: QoRTs seems to be pretty memory hungry. I've allocated Java 10gb of memory here, you may need to change this depending on available resources
ls *.bam | parallel -j 50 "mkdir {.} && java -Xmx10G -jar /plas1/amardeep.singh/apps/QoRTs-STABLE.jar --stranded --maxReadLength 101 --singleEnded {} /plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/gtf/Drosophila_melanogaster.BDGP6.28.99.gtf.gz {.}/"

# Make a 'decoder' file that will use the sample names as the directories to find outputs
# I had thes in the following folder, but this might change depending on where and what is in the directory.
ls /plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/RAL.mapped.files.to.ensembl/bam.files/samtools.sort.out/qorts.output/  > QoRTs.decoder.file.txt

# ----  ----

# Read QoRTs outputs into R to generate a sizeFactor by which to normalize read counts
# Tell QoRTs which directory to look into for directory that contain outputs. Each directory should be named to correspond to a sample (which should be decoded by the decoder file) and give it the location of your decoder file
qorts.results = read.qc.results.data("/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/RAL.mapped.files.to.ensembl/bam.files/samtools.sort.out/qorts.output/", decoder.files = "/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/RAL.mapped.files.to.ensembl/bam.files/samtools.sort.out/qorts.output/QoRTs.decoder.file.txt", calc.DESeq2 = TRUE)

# Save size factors for each sample in a text file
get.size.factors(qorts.results, outfile = "/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/RAL.size.Factors.GEO.txt");

#












#
