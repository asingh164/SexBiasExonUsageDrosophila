########################################################################################
###### JunctionSeq Script for Differential Exon Usage Analysis -- Resample script  #####
########################################################################################

## Author: Amardeep Singh -- amardeep.singh[at]utoronto.ca
## This program makes use of public releases of the Drosophila genome from Ensembl (BDGP6.28)
## To retrieve the genome assembly file yourself:
# Fasta file
# wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-46/fasta/drosophila_melanogaster/dna//Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa.gz
# GTF file
# wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-46/fasta/drosophila_melanogaster/

# This script requires R package which can be downloaded through Bioconductor

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

rm(list=ls())
# Loading in decoder file
decoder <- read.table("/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/QoRTs.decoder.file.for.JunctionSeq.txt", header = TRUE, stringsAsFactors = FALSE)
# Here I am subsetting only the samples that I want to use (for this analysis I started just with body tissue and only a single female replicate)
decoder.for.junctionseq <- decoder[!(grepl("head|replicate.2", decoder$unique.ID)),]
rownames(decoder.for.junctionseq) = c(1:nrow(decoder.for.junctionseq))

# Resample the condition factor (i.e., sex)
# Add an index to the decoder file
decoder.for.junctionseq$index <- seq(from = 1, to = nrow(decoder.for.junctionseq), by = 1)
male.rows <- (decoder.for.junctionseq[decoder.for.junctionseq$sex == "male", ])$index
female.rows <- (decoder.for.junctionseq[decoder.for.junctionseq$sex == "female", ])$index

#Providing the directory for the count files:
countFiles <- paste0("/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/count.files/body.only.replicate.1/", decoder.for.junctionseq$unique.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")


for (i in 1:10) {
  male.rows.sample <- sample(x = (decoder.for.junctionseq[decoder.for.junctionseq$sex == "male", ])$index,,
                            size = (nrow(decoder.for.junctionseq) / 4),
                            replace = FALSE)

  female.rows.sample <- male.rows.sample - 1

  decoder.for.junctionseq.resample <- decoder.for.junctionseq
  decoder.for.junctionseq.resample$sex.resample <- as.character(decoder.for.junctionseq.resample$sex)
  decoder.for.junctionseq.resample$sex.resample[decoder.for.junctionseq.resample$index %in% male.rows.sample] <- "female"
  decoder.for.junctionseq.resample$sex.resample[decoder.for.junctionseq.resample$index %in% female.rows.sample] <- "male"

  ###  Run the differential exon usage (DEU) analysis:

  # Creast a design dataframe which has a column for the condition you are testng (This will need to change when adding factors not sure how yet)
  design.df.resample <- data.frame(condition = factor(decoder.for.junctionseq.resample$sex.resample))

  # Run Junctionseq
  count.set.object.resample <- readJunctionSeqCounts(countfiles = countFiles,
                                                    samplenames = decoder.for.junctionseq.resample$unique.ID,
                                                    design = design.df.resample,
                                                    flat.gff.file = "/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/count.files/body.only.replicate.1/withNovel.forJunctionSeq.gff.gz",
                                                    nCores = 50,
                                                    verbose = TRUE)

  # Generate size factors for normalization and load them into the count.set.object
  count.set.object.resample <- estimateJunctionSeqSizeFactors(count.set.object.resample)


  # Generate test specific dispersion estimates and load into count.set.object
  count.set.object.resample <- estimateJunctionSeqDispersions(count.set.object.resample, nCores = 50)

  # Fit the observed dispersions to a regression to create a fitted dispersion
  count.set.object.resample <- fitJunctionSeqDispersionFunction(count.set.object.resample)

  # Perform the hypothesis tests to test or differential splice junction/exon usage (DEU)
  count.set.object.resample  <- testForDiffUsage(count.set.object.resample , nCores = 50)

  # Calculate effect sizes and parameter estimates
  count.set.object.resample <- estimateEffectSizes(count.set.object.resample)

  # Save output to file
  setwd("/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/ResampledData")
  outputprefix = paste0("resample", i)
  writeCompleteResults(count.set.object.resample,
                      outfile.prefix = outputprefix,
                      gzip.output = TRUE,
                      FDR.threshold = 0.01,
                      save.allGenes = TRUE, save.sigGenes = TRUE,
                      save.fit = FALSE, save.VST = FALSE,
                      save.bedTracks = TRUE,
                      save.jscs = TRUE,
                      bedtrack.format = c("BED", "GTF", "GFF3"),
                      verbose = TRUE)
  print(i)
}
