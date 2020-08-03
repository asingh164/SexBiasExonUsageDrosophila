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

# This script requires R package which can be downloaded through Bioconductor

### Script Details ###
# This script is written primarily in R, Bash commands are denoted by #---- BASH ---- and end with # ---- ----


# ---- Linux ----
# Make a flat gff file from the genome annotation provided by Ensembl
QoRTs makeFlatGff --stranded \
  /plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/gtf/Drosophila_melanogaster.BDGP6.28.99.gtf.gz \
  /plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/Dmel.JunctionSeq.flat.gff.gz

## Processing bam files
# First, run QoRTs to perform QC checks on aligned bam files. This needs to be done so that the size factor for normalization can be done
# NOTE: QoRTs seems to be pretty memory hungry. I've allocated Java 10gb of memory here, you may need to change this depending on available resources
ls *.bam | parallel -j 50 \
"mkdir {.} && java -Xmx10G -jar /plas1/amardeep.singh/apps/QoRTs-STABLE.jar \
                    --stranded --maxReadLength 101 \
                    --singleEnded {} /plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/gtf/Drosophila_melanogaster.BDGP6.28.99.gtf.gz {.}/"

# Make a 'decoder' file that will use the sample names as the directories to find outputs
# I had thes in the following folder, but this might change depending on where and what is in the directory.
ls /plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/RAL.mapped.files.to.ensembl/bam.files/samtools.sort.out/qorts.output/  > QoRTs.decoder.file.txt
mv QoRTs.decoder.file.txt /plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.file

# ----  ----

# ---- R Code ---
require(QoRTs)

# Read QoRTs outputs into R to generate a sizeFactor by which to normalize read counts
# Tell QoRTs which directory to look into for directory that contain outputs. Each directory should be named to correspond to a sample (which should be decoded by the decoder file) and give it the location of your decoder file
qorts.results <- read.qc.results.data("/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/RAL.mapped.files.to.ensembl/bam.files/samtools.sort.out/qorts.output/",
                                      decoder.files = "/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.file/QoRTs.decoder.file.txt",
                                      calc.DESeq2 = TRUE)

# Save size factors for each sample in a text file
get.size.factors(qorts.results, outfile = "/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/RAL.size.Factors.GEO.txt");

# ---- Linux ----
# Create novel junction splice sites
java -Xmx10G -jar /plas1/amardeep.singh/apps/QoRTs-STABLE.jar mergeNovelSplices --minCount 20 --stranded \
/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/RAL.mapped.files.to.ensembl/bam.files/samtools.sort.out/qorts.output \
/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/RAL.size.Factors.GEO.txt \
/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/gtf/Drosophila_melanogaster.BDGP6.28.99.gtf.gz \
/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/count.files

# ---- R Code ----
### NOTE: Junctonseq (as of July 2020) will require that you install an older version of DESeq2 (versio 1.10.1)
##                https://www.bioconductor.org/packages/3.2/bioc/src/contrib/DESeq2_1.10.1.tar.gz
## NOTE: This same script was used to independently examine differential exon usage in both whole body and brain tissue but the code only refers to a single tissue.

### JunctionSeq Script for Body Tissue

rm(list=ls())
require(DESeq2)
require(JunctionSeq)
require(BiocParallel) # This is a package used for paralellization of jobs

# Load in decoder files and add fields for conditions

# Run this only once to edit the decoder file
#decoder.file =  read.delim("/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/QoRTs.decoder.file.txt", header = FALSE, stringsAsFactors=FALSE)
#colnames(decoder.file) = "unique.ID"
#sex.tmp = gsub("SAM.sorted.trimmed.RAL.....","", decoder.file$unique.ID)
#decoder.file$sex =  gsub("\\..*","", sex.tmp)
#write.table(decoder.file, "/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/QoRTs.decoder.file.for.JunctionSeq.txt", quote = F, row.names = FALSE, col.names = TRUE, sep = "\t")

# Loading in decoder file
decoder <- read.table("/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/QoRTs.decoder.file.for.JunctionSeq.txt", header = TRUE, stringsAsFactors = FALSE)
# Here I am subsetting only the samples that I want to use (for this analysis I started just with body tissue and only a single female replicate)
decoder.for.junctionseq <- decoder[!(grepl("head|replicate.2", decoder$unique.ID)),]

#Providing the directory for the count files:
countFiles <- paste0("/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/count.files/body.only.replicate.1/", decoder.for.junctionseq$unique.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")

###  Run the differential exon usage (DEU) analysis:

# Creast a design dataframe which has a column for the condition you are testng (This will need to change when adding factors not sure how yet)
design.df <- data.frame(condition = factor(decoder.for.junctionseq$sex))

# Building the count set object that JunctionSeq will analyze and add to it all of the parameters of the analysis
count.set.object <- readJunctionSeqCounts(countfiles = countFiles,
                                          samplenames = decoder.for.junctionseq$unique.ID,
                                          design = design.df,
                                          flat.gff.file = "/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/count.files/withNovel.forJunctionSeq.gff.gz",
                                          nCores = 70,
                                          verbose = TRUE)

# Generate size factors for normalization and load them into the count.set.object
count.set.object <- estimateJunctionSeqSizeFactors(count.set.object)

# Generate test specific dispersion estimates and load into count.set.object
count.set.object <- estimateJunctionSeqDispersions(count.set.object, nCores = 70)

# Fit the observed dispersions to a regression to create a fitted dispersion
count.set.object <- fitJunctionSeqDispersionFunction(count.set.object)

# Perform the hypothesis tests to test for differential splice junction/exon usage (DEU)
count.set.object <- testForDiffUsage(count.set.object, nCores = 60, )

# Calculate effect sizes and parameter estimates
count.set.object <- estimateEffectSizes(count.set.object)

# Save output to file
writeCompleteResults(count.set.object, "/plas1/amardeep.singh/RNA.Seq.Data/JunctionSeq.Files/BodyOutput/Aug1.Body.Only",
                    gzip.output = TRUE,
                    FDR.threshold = 0.01,
                    save.allGenes = TRUE, save.sigGenes = TRUE,
                    save.fit = FALSE, save.VST = FALSE,
                    save.bedTracks = TRUE,
                    save.jscs = TRUE,
                    bedtrack.format = c("BED", "GTF", "GFF3"),
                    verbose = TRUE)


##### JunctionSeq Script for Head tissue

rm(list=ls())
require(DESeq2)
require(JunctionSeq)
require(BiocParallel) # This is a package used for paralellization of jobs

# Load in decoder files and add fields for conditions

# Run this only once to edit the decoder file
#decoder.file =  read.delim("/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/QoRTs.decoder.file.txt", header = FALSE, stringsAsFactors=FALSE)
#colnames(decoder.file) = "unique.ID"
#sex.tmp = gsub("SAM.sorted.trimmed.RAL.....","", decoder.file$unique.ID)
#decoder.file$sex =  gsub("\\..*","", sex.tmp)
#write.table(decoder.file, "/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/QoRTs.decoder.file.for.JunctionSeq.txt", quote = F, row.names = FALSE, col.names = TRUE, sep = "\t")

# Loading in decoder file
decoder <- read.table("/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/QoRTs.decoder.file.for.JunctionSeq.txt", header = TRUE, stringsAsFactors = FALSE)
# Here I am subsetting only the samples that I want to use (for this analysis I started just with body tissue and only a single female replicate)
decoder.for.junctionseq <- decoder[!(grepl("body|replicate.2", decoder$unique.ID)),]

#Providing the directory for the count files:
countFiles <- paste0("/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/count.files/head.only.replicate.1/", decoder.for.junctionseq$unique.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")

###  Run the differential exon usage (DEU) analysis:

# Creast a design dataframe which has a column for the condition you are testng (This will need to change when adding factors not sure how yet)
design.df <- data.frame(condition = factor(decoder.for.junctionseq$sex))

# Building the count set object that JunctionSeq will analyze and add to it all of the parameters of the analysis
count.set.object <- readJunctionSeqCounts(countfiles = countFiles,
                                          samplenames = decoder.for.junctionseq$unique.ID,
                                          design = design.df,
                                          flat.gff.file = "/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/count.files/withNovel.forJunctionSeq.gff.gz",
                                          nCores = 70,
                                          verbose = TRUE)

# Generate size factors for normalization and load them into the count.set.object
count.set.object <- estimateJunctionSeqSizeFactors(count.set.object)

# Generate test specific dispersion estimates and load into count.set.object
count.set.object <- estimateJunctionSeqDispersions(count.set.object, nCores = 70)

# Fit the observed dispersions to a regression to create a fitted dispersion
count.set.object <- fitJunctionSeqDispersionFunction(count.set.object)

# Perform the hypothesis tests to test for differential splice junction/exon usage (DEU)
count.set.object <- testForDiffUsage(count.set.object, nCores = 70, )

# Calculate effect sizes and parameter estimates
count.set.object <- estimateEffectSizes(count.set.object)

# Save output to file
writeCompleteResults(count.set.object, "/plas1/amardeep.singh/RNA.Seq.Data/JunctionSeq.Files/HeadOutput/Aug1.Head.Only",
                    gzip.output = TRUE,
                    FDR.threshold = 0.01,
                    save.allGenes = TRUE, save.sigGenes = TRUE,
                    save.fit = FALSE, save.VST = FALSE,
                    save.bedTracks = TRUE,
                    save.jscs = TRUE,
                    bedtrack.format = c("BED", "GTF", "GFF3"),
                    verbose = TRUE)


## Analyzing junctionseq output

jseq.allgenes.body = gzfile('/plas1/amardeep.singh/RNA.Seq.Data/JunctionSeq.Files/BodyOutput/Aug1.Body.OnlyallGenes.results.txt.gz', 'rt')
jseq.allgenes.head = gzfile('/plas1/amardeep.singh/RNA.Seq.Data/JunctionSeq.Files/HeadOutput/Aug1.Head.OnlyallGenes.results.txt.gz', 'rt')

jseq.allgenes.body = read.delim(jseq.allgenes.body, header = TRUE)
jseq.allgenes.head  = read.delim(jseq.allgenes.head, header = TRUE)

length(unique((jseq.allgenes.body[jseq.allgenes.body$geneWisePadj < 0.01,])$geneID))
length(unique((jseq.allgenes.head[jseq.allgenes.head$geneWisePadj < 0.01,])$geneID))


# Remove all sites where there are fewer than 5 reads mapping to a site
jseq.allgenes.body.remove.low.expression = jseq.allgenes.body[jseq.allgenes.body$expr_female > 20 | jseq.allgenes.body$expr_male > 20, ]
length(unique((jseq.allgenes.body.remove.low.expression[jseq.allgenes.body.remove.low.expression$geneWisePadj < 0.01,])$geneID))

jseq.allgenes.head.remove.low.expression = jseq.allgenes.head[jseq.allgenes.head$expr_female > 20 | jseq.allgenes.head$expr_male > 20, ]
length(unique((jseq.allgenes.head.remove.low.expression[jseq.allgenes.head.remove.low.expression$geneWisePadj < 0.01,])$geneID))



#
