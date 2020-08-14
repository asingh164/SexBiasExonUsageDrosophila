################################################################
###### DEXSeq Script for Differential Exon Usage Analysis  #####
################################################################


#### D E P R E C A T E D -- S W I T C H E D   T O   U S I N G    J U C T I O N S E Q ####


## Author: Amardeep Singh -- amardeep.singh[at]utoronto.ca
## This program makes use of public releases of the Drosophila genome from Flybase (dmel_r6.32_FB2020_01)
## To retrieve the genome assembly file yourself:
# wget ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.32_FB2020_01/

### Script Details ###
# This script is written primarily in R, Bash commands are denoted by #---- BASH ---- and end with # ---- ----

# Python scripts:

# You will need to have installed DEXSeq in R which will automatically download two Python scripts that converts the GTF file into a
# format necessary for DEXSeq
# Location of Python scripts
pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )
list.files(pythonScriptsDir)
system.file( "python_scripts", package="DEXSeq", mustWork=TRUE ) # This is the directory where the python scripts have been saved
#

# There are two steps that need to be done to prepare bam files for DEXSeq and both require you to run provided Python scripts
# 1. Prepare genome annotation file (GTF) by converting it to a flat GFF file by running Python script provided by DEXSeq to convert a gtf file to "flat" format required for DEXSeq
#---- BASH ----

# NOTE: There seems to be some error in the .gtf file that I'm not sure how to rectify.
Running the following script with the original .gtf file supplied by FlyBase I get the following error:

#    "Traceback (most recent call last): File "/plas1/amardeep.singh/R/x86_64-redhat-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_prepare_annotation.py", line 127, in <module>
#    assert l[i].iv.end <= l[i+1].iv.start, str(l[i+1]) + " starts too early"  AssertionError: <GenomicFeature: exonic_part 'FBgn0267649+FBgn0266176+FBgn0261839+FBgn0261845+FBgn0266178+FBgn0267652+FBgn0266174+FBgn0261837+FBgn0267648+FBgn0261843+FBgn0266170+FBgn0266172+FBgn0261841+FBgn0261838+FBgn0261844+FBgn0266177+FBgn0261840+FBgn0266173+FBgn0002781+FBgn0261842+FBgn0267651+FBgn0267650+FBgn0266171+FBgn0266175' at 3R: 21351242 -> 21350844 (strand '-')> starts too early"

# These FlyBase IDs appear to be related to a single gene: mod(mdg4) (FlyBaseId: FBgn0002781). When I grep this out:
grep -v "mod(mdg4)" dmel-all-r6.32.gtf > dmel-all-r6.32.filtered.gtf

# And then run the python script on this filtered file, I get no such error.
python /plas1/amardeep.singh/R/x86_64-redhat-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_prepare_annotation.py \
/plas1/amardeep.singh/Flybase.Dmel.Genome.Release/gtf/dmel-all-r6.32.filtered.gtf \
/plas1/amardeep.singh/RNA.Seq.Data/DEXSeq.Files/dmel-all-r6.32.filtered.flat.gff

# 2. Counting reads in "exon bins" for each BAM file that will form the basis for analysis by DEXSeq
ls SAM*.bam | parallel -j 40 python /plas1/amardeep.singh/R/x86_64-redhat-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -r name -f bam -a 20 /plas1/amardeep.singh/RNA.Seq.Data/DEXSeq.Files/dmel-all-r6.32.filtered.flat.gff {} /plas1/amardeep.singh/RNA.Seq.Data/DEXSeq.Files/DEXSeq.bam.count.files/{}.read.counts.txt

#---- ----

# DEXSeq Script

## In this first pass, I'm only going to be concerned with sex differences uing whole body tissue
require(DEXSeq) # Main package used for differential exon usage
require(BiocParallel) # This is a package used for paralellization of jobs
require(doBy)

# Set working directory to the directory with samples
setwd("/plas1/amardeep.singh/RNA.Seq.Data/DEXSeq.Files/DEXSeq.bam.count.files/")

# Generate a table that contains information on the sample runs
# Vectors containing file names for the samples
file.names.for.samples = dir("/plas1/amardeep.singh/RNA.Seq.Data/DEXSeq.Files/DEXSeq.bam.count.files")
# grep file names for whole body tissue only
#file.names.for.samples = Filter(function(x) grepl("body", x), file.names.for.samples)

# Vectors containing sample names (cleaned up a bit -- Probably a more elegant way of doing this)
sample.names.tmp = gsub(".fastq_Aligned.*","",file.names.for.samples)
#sample.names.tmp = gsub(".replicate.1.*","",sample.names.tmp)
sample.names = gsub("SAM.sorted.trimmed.","",sample.names.tmp)


# Vector containing sex of sample
sex.tmp = gsub("RAL.....", "", sample.names)
sex = gsub("\\..*", "", sex.tmp)

# Vector containing tissue of sample
#tissue.tmp = gsub("RAL.*male.", "", sample.names)
#tissue = gsub("replicate..*", "", tissue.tmp)

# Make a dataframe for samples
sample.table = data.frame(row.names = file.names.for.samples, sample.names = sample.names, sex = sex)

# Constructing the DEXSeq data set
data.table = DEXSeqDataSetFromHTSeq(file.names.for.samples,
                                    sampleData = sample.table,
                                    design = ~sample + exon + exon:sex,
                                    flattenedfile = "/plas1/amardeep.singh/RNA.Seq.Data/DEXSeq.Files/dmel-all-r6.32.filtered.flat.gff")

full.model = ~sample + exon + sex:exon
#reduced.model = ~sample + exon

# Normalization of read counts
data.table = estimateSizeFactors(data.table)

# Estimate dispersion of read counts per exon bin ## NOTE: I've set this to use 50 cores, you may need to adjust depending on available resources
BPPARAM = MulticoreParam(50)
data.table = estimateDispersions(data.table, BPPARAM = MulticoreParam(50), formula = full.model)

# Filtering exons based on a average read count ## For now
# Calculate average read count per exon
exon.read.counts = as.data.frame(counts(data.table))
exon.read.counts$mean.reads.mapped = apply(exon.read.counts, 1, mean)

# Remove exons with fewer than an average of 50 reads mapping to it
exon.read.counts.subset = exon.read.counts[exon.read.counts$mean.reads.mapped > 25,]

# Count number of exons per gene
gene.IDs = as.data.frame(substring(row.names(exon.read.counts.subset), 1, 11))
colnames(gene.IDs) = "gene.ID"
gene.counts.unique = summaryBy(gene.ID ~ gene.ID, data = gene.IDs, FUN=c(length))

# Remove genes that only have a single exon
multi.exon.genes = gene.counts.unique[gene.counts.unique$gene.ID.length > 1,]

gene.counts.unique = as.data.frame(unique(gene.IDs$gene.ID))
colnames(gene.counts.unique) = "gene.ID"

# Perform differential exon usage analysis
data.table = testForDEU(data.table, BPPARAM = MulticoreParam(50), fullModel = full.model)

# Calculating Exon Fold Change
data.table = estimateExonFoldChanges(data.table, fitExpToVar = "sex", BPPARAM = MulticoreParam(50))


# Piping differential exon usage results into a new dataframe
DEU.results = as.data.frame(DEXSeqResults(data.table))
DEU.results.filtered = DEU.results[as.character(DEU.results$groupID) %in% as.character(multi.exon.genes$gene.ID),]

# Table of results -- number of genes that show signnificant exon * sex interaction
table(tapply(DEU.results.filtered$padj < 0.05, DEU.results.filtered$groupID, any))


test.gene.data = DEU.results.filtered[DEU.results.filtered$groupID == "FBgn0000008",]

female.data = test.gene.data[,c(grep("female",names(test.gene.data))), drop=F]
female.data$mean = apply(female.data, 1, function(x) mean(x))
female.data$se= apply(female.data, 1, function(x) sd(x)/sqrt(18))
female.exon.usage = female.data[,19:20]

male.data = test.gene.data[,c(grep("\\.male",names(test.gene.data))), drop=F]
male.data$mean = apply(male.data, 1, function(x) mean(x))
male.data$se= apply(male.data, 1, function(x) sd(x)/sqrt(18))
male.exon.usage = male.data[,19:20]

#
