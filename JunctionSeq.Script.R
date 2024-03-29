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

# This script requires the 'JunctionSeq' R packages which can be downloaded through Bioconductor

### Script Details ###
# This script is written primarily in R, Bash commands are denoted by #---- Linux ---- and end with # ---- ----


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
count.set.object <- testForDiffUsage(count.set.object, nCores = 60)

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
require(VennDiagram)
require(grDevices)

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


#########################################
##    Analyzing junctionseq output    ###
#########################################
require(VennDiagram)
require(grDevices)

#jseq.allgenes.body = gzfile('/plas1/amardeep.singh/RNA.Seq.Data/JunctionSeq.Files/BodyOutput/Aug1.Body.OnlyallGenes.results.txt.gz', 'rt')
#jseq.allgenes.head = gzfile('/plas1/amardeep.singh/RNA.Seq.Data/JunctionSeq.Files/HeadOutput/Aug1.Head.OnlyallGenes.results.txt.gz', 'rt')
jseq.allgenes.body = read.table("/plas1/amardeep.singh/RNA.Seq.Data/JunctionSeq.Files/BodyOutput/Aug1.Body.OnlyallGenes.results.txt", sep = "\t", header = TRUE)
jseq.allgenes.head = read.table("/plas1/amardeep.singh/RNA.Seq.Data/JunctionSeq.Files/HeadOutput/Aug1.Head.OnlyallGenes.results.txt", sep = "\t", header = TRUE)

# Clean up files and remove any exons that were not tested
jseq.allgenes.body = jseq.allgenes.body[!(is.na(jseq.allgenes.body$expr_female)),]
jseq.allgenes.body = jseq.allgenes.body[!(is.na(jseq.allgenes.body$expr_male)),]
jseq.allgenes.head = jseq.allgenes.head[!(is.na(jseq.allgenes.head$expr_female)),]
jseq.allgenes.head = jseq.allgenes.head[!(is.na(jseq.allgenes.head$expr_male)),]

# Obtain a list of significant and non-significant genes and output into seperate txt files
# first, remove all genes with fewer than 50 reads mapping to exons in males and females
jseq.allgenes.body = jseq.allgenes.body[jseq.allgenes.body$expr_female > 50 & jseq.allgenes.body$expr_male > 50,]
jseq.allgenes.head = jseq.allgenes.head[jseq.allgenes.head$expr_female > 50 & jseq.allgenes.head$expr_male > 50,]

# Unique genes with significant DEU
sig.genes.body = unique(jseq.allgenes.body$geneID[jseq.allgenes.body$geneWisePadj < 0.01])
sig.genes.head = unique(jseq.allgenes.head$geneID[jseq.allgenes.head$geneWisePadj < 0.01])
nonsig.genes.body = unique(jseq.allgenes.body$geneID[jseq.allgenes.body$geneWisePadj > 0.01])
nonsig.genes.head = unique(jseq.allgenes.head$geneID[jseq.allgenes.head$geneWisePadj > 0.01])


write.table(sig.genes.body, file = "/plas1/amardeep.singh/RNA.Seq.Data/JunctionSeq.Files/body.significant.genes.Aug18.txt", sep = "\t", row.name = FALSE, col.names = FALSE, quote = FALSE)
write.table(sig.genes.head, file = "/plas1/amardeep.singh/RNA.Seq.Data/JunctionSeq.Files/head.significant.genes.Aug18.txt", sep = "\t", row.name = FALSE, col.names = FALSE, quote = FALSE)
write.table(nonsig.genes.body, file = "/plas1/amardeep.singh/RNA.Seq.Data/JunctionSeq.Files/body.non.significant.genes.Aug18.txt", sep = "\t", row.name = FALSE, col.names = FALSE, quote = FALSE)
write.table(nonsig.genes.head, file = "/plas1/amardeep.singh/RNA.Seq.Data/JunctionSeq.Files/head.non.significant.genes.Aug18.txt", sep = "\t", row.name = FALSE, col.names = FALSE, quote = FALSE)


# Remove all sites where there are fewer than 50 reads mapping to a site
jseq.allgenes.body.remove.low.expression = jseq.allgenes.body[jseq.allgenes.body$expr_female > 50 | jseq.allgenes.body$expr_male > 50, ]
length(unique((jseq.allgenes.body.remove.low.expression[jseq.allgenes.body.remove.low.expression$geneWisePadj < 0.01,])$geneID))

jseq.allgenes.head.remove.low.expression = jseq.allgenes.head[jseq.allgenes.head$expr_female > 50 | jseq.allgenes.head$expr_male > 50, ]
length(unique((jseq.allgenes.head.remove.low.expression[jseq.allgenes.head.remove.low.expression$geneWisePadj < 0.01,])$geneID))


# Venn Diagrams to present tissue specific data

venndiagram = venn.diagram(x = list(sig.genes.body, sig.genes.head, nonsig.genes.body, nonsig.genes.head),
                           category.names = c("Body Sig", "Head Sig", "Body NonSig", "Head NonSig"),
                           filename = NULL,
                            # Circles
                            lwd = 2, lty = 'blank', fill = c("#7294d4", "#72D481", "#D4B272", "#D472C5"),
                            # Numbers
                            cex = 1.5, fontface = "bold", fontfamily = "Helvetica",
                            )
pdf(file="/plas1/amardeep.singh/tmp/VennDiagramHeadBody.pdf", height = 10, width = 10)
    grid.draw(venndiagram)
dev.off()


# Plotting fraction of genes by chromosome
# Load in JunctionSeq outputs

# First, I generated a list of genes and which chromosome they are in using the Flybase D. mel gtf file
#---- Bash ----
# First pull out a list of every gene in the gtf
cat dmel-all-r6.32.filtered.gtf | cut -f9 | cut -f2 | cut -c10-20 > list.of.genes.txt
# Next, pull out a list of each chromosome location
cat dmel-all-r6.32.filtered.gtf | cut -f1 > list.of.chromosomes.txt
# Merge the two files together
paste list.of.chromosomes.txt list.of.genes.txt > gene.locations.txt

# ----/----

# --- R Code ---
rm(list = ls())
require(doBy)
require(ggplot2)

## Loading in data files
# Read in JunctionSeq results and gene location data file
junctionseq.results.body = read.table("/plas1/amardeep.singh/RNA.Seq.Data/JunctionSeq.Files/BodyOutput/Aug1.Body.OnlyallGenes.results.txt", header = TRUE, sep = "\t")
junctionseq.results.head = read.table("/plas1/amardeep.singh/RNA.Seq.Data/JunctionSeq.Files/HeadOutput/Aug1.Head.OnlyallGenes.results.txt", header = TRUE, sep = "\t")

## Cleaning up junctionseq files
# Remove any sites that were not tested
junctionseq.results.body = junctionseq.results.body[!(is.na(junctionseq.results.body$pvalue)),]
junctionseq.results.head = junctionseq.results.head[!(is.na(junctionseq.results.head$pvalue)),]
# Remove any rows where the expr in males or females is less than 50
junctionseq.results.body.filtered = junctionseq.results.body[junctionseq.results.body$expr_male > 50 & junctionseq.results.body$expr_female > 50 ,]
junctionseq.results.head.filtered = junctionseq.results.head[junctionseq.results.head$expr_male > 50 & junctionseq.results.head$expr_female > 50 ,]
# Keep only columns that I need
junctionseq.results.body.filtered = junctionseq.results.body.filtered[,c(2,14,25)]
junctionseq.results.head.filtered = junctionseq.results.head.filtered[,c(2,14,25)]

# Remove duplicates
junctionseq.results.body.filtered = junctionseq.results.body.filtered[!duplicated(junctionseq.results.body.filtered[1:3]),]
junctionseq.results.head.filtered = junctionseq.results.head.filtered[!duplicated(junctionseq.results.head.filtered[1:3]),]

# Assign SDIU and non-SDIU genes
junctionseq.results.body.filtered$sig.hit[junctionseq.results.body.filtered$geneWisePadj <= 0.01] = 1
junctionseq.results.body.filtered$sig.hit[junctionseq.results.body.filtered$geneWisePadj > 0.01] = 0
junctionseq.results.head.filtered$sig.hit[junctionseq.results.head.filtered$geneWisePadj <= 0.01] = 1
junctionseq.results.head.filtered$sig.hit[junctionseq.results.head.filtered$geneWisePadj > 0.01] = 0

# Calcualte fractions of genes that are SDIU and non-SDIU
sdiu.body.fraction = nrow(junctionseq.results.body.filtered[junctionseq.results.body.filtered$sig.hit == 1 & junctionseq.results.body.filtered$chr=="X",]) / nrow(junctionseq.results.body.filtered[junctionseq.results.body.filtered$sig.hit == 1,])
non.sdiu.body.fraction = nrow(junctionseq.results.body.filtered[junctionseq.results.body.filtered$sig.hit == 0 & junctionseq.results.body.filtered$chr=="X",]) / nrow(junctionseq.results.body.filtered[junctionseq.results.body.filtered$sig.hit == 0,])
sdiu.head.fraction = nrow(junctionseq.results.head.filtered[junctionseq.results.head.filtered$sig.hit == 1 & junctionseq.results.head.filtered$chr=="X",]) / nrow(junctionseq.results.head.filtered[junctionseq.results.head.filtered$sig.hit == 1,])
non.sdiu.head.fraction = nrow(junctionseq.results.head.filtered[junctionseq.results.head.filtered$sig.hit == 0 & junctionseq.results.head.filtered$chr=="X",]) / nrow(junctionseq.results.head.filtered[junctionseq.results.head.filtered$sig.hit == 0,])

# Resampling
# Subset genes in SDIU and non-SDIU genes
sdiu.body.subset = junctionseq.results.body.filtered[junctionseq.results.body.filtered$sig.hit == 1,]
non.sdiu.body.subset = junctionseq.results.body.filtered[junctionseq.results.body.filtered$sig.hit == 0,]
sdiu.head.subset = junctionseq.results.head.filtered[junctionseq.results.head.filtered$sig.hit == 1,]
non.sdiu.head.subset = junctionseq.results.head.filtered[junctionseq.results.head.filtered$sig.hit == 0,]

sdiu.body.fraction.resample = vector(mode = "numeric", length = 10000)
non.sdiu.body.fraction.resample = vector(mode = "numeric", length = 10000)
sdiu.head.fraction.resample = vector(mode = "numeric", length = 10000)
non.sdiu.head.fraction.resample = vector(mode = "numeric", length = 10000)


# Resampling
for (i in 1:10000){
  sample.body.sdiu.fraction = sample(sdiu.body.subset$chr, nrow(sdiu.body.subset), replace = TRUE)
  sample.body.non.sdiu.fraction = sample(non.sdiu.body.subset$chr, nrow(non.sdiu.body.subset), replace = TRUE)
  sample.head.sdiu.fraction = sample(sdiu.head.subset$chr, nrow(sdiu.head.subset), replace = TRUE)
  sample.head.non.sdiu.fraction = sample(non.sdiu.head.subset$chr, nrow(non.sdiu.head.subset), replace = TRUE)

  sdiu.body.fraction.resample[i] = length(sample.body.sdiu.fraction[sample.body.sdiu.fraction == "X"]) / length(sample.body.sdiu.fraction)
  non.sdiu.body.fraction.resample[i] = length(sample.body.non.sdiu.fraction[sample.body.non.sdiu.fraction == "X"]) / length(sample.body.non.sdiu.fraction)
  sdiu.head.fraction.resample[i] = length(sample.head.sdiu.fraction[sample.head.sdiu.fraction == "X"]) / length(sample.head.sdiu.fraction)
  non.sdiu.head.fraction.resample[i] = length(sample.head.non.sdiu.fraction[sample.head.non.sdiu.fraction == "X"]) / length(sample.head.non.sdiu.fraction)
print(i)
}

body.lower.CI = as.vector( c( quantile(sdiu.body.fraction.resample, 0.025)[1],
                              quantile(non.sdiu.body.fraction.resample, 0.025)[1]))

body.upper.CI = as.vector( c( quantile(sdiu.body.fraction.resample, 0.975)[1],
                              quantile(non.sdiu.body.fraction.resample, 0.975)[1]))

head.lower.CI = as.vector( c( quantile(sdiu.head.fraction.resample, 0.025)[1],
                              quantile(non.sdiu.head.fraction.resample, 0.025)[1]))

head.upper.CI = as.vector( c( quantile(sdiu.head.fraction.resample, 0.975)[1],
                              quantile(non.sdiu.head.fraction.resample, 0.975)[1]))

# Make dataframe for plotting
body.df.tmp = as.data.frame(c(sdiu.body.fraction,non.sdiu.body.fraction))
colnames(body.df.tmp) = "fraction"
body.df.tmp$tissue = "body"
body.df.tmp$sdiu = c(1,0)
body.df.tmp$upper.CI = body.upper.CI
body.df.tmp$lower.CI = body.lower.CI

head.df.tmp = as.data.frame(c(sdiu.head.fraction,non.sdiu.head.fraction))
colnames(head.df.tmp) = "fraction"
head.df.tmp$tissue = "head"
head.df.tmp$sdiu = c(1,0)
head.df.tmp$upper.CI = head.upper.CI
head.df.tmp$lower.CI = head.lower.CI

fraction.df = rbind(body.df.tmp,head.df.tmp)

# Plotting
fraction.of.SDIU.X.chromosome = ggplot(fraction.df, aes(y=fraction, x = as.factor(tissue), colour = as.factor(sdiu))) +
                                  geom_point(aes(y=fraction, x = as.factor(tissue), colour = as.factor(sdiu)), size = 8,
                                  position = position_dodge(width = 0.5)) + ylim(0,0.5) +
                                  geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI, x = as.factor(tissue)),
                                  width = 0, position = position_dodge(width = 0.5)) +
                                  theme_bw()  + scale_colour_manual(values = c("#800080", "#65c86e")) +
                                  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                                      axis.text = element_text(face="bold", color="black",size=15, family = "Helvetica"),
                                      axis.title=element_blank(), legend.position = "none")
                                  #facet_grid(tissue~., scales = "free")
pdf("/plas1/amardeep.singh/tmp/sdiu.fraction.on.x.chromosome.sept16.pdf",height = 10, width = 5)
fraction.of.SDIU.X.chromosome
dev.off()




#
