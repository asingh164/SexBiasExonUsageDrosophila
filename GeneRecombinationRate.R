# This script was used to assign a recombination rate to genes in the Drosophila
rm(list = ls())

# Source the R script to find which recombination 'window' each gene falls into. Named 'Gene.Locations.from.GTF.R'
source("/plas1/amardeep.singh/RNA.Seq.Data/RecombinationRateData/Gene.Locations.from.GTF.R")

# Load in gtf and recombination data files
gtf = read.delim("/plas1/amardeep.singh/Flybase.Dmel.Genome.Release/release5.26_gtf/dmel-all-r5.36.gtf", header = FALSE)
recombinationRateList = read.delim("/plas1/amardeep.singh/RNA.Seq.Data/RecombinationRateData/Comeron_tables/Comeron_100kb_allchr.txt", header = FALSE)

# Clean gtf of the Y chromosome and mitochondria
gtf = gtf[!(gtf$V1 == "dmel_mitochondrion_genome" | gtf$V1 == "Y" | gtf$V1 == "YHet" | gtf$V1 == "U" | gtf$V1 == "4"),]

# Get recombination rate for each gene
GeneRecombinationRate = GetLocalRecombinationRatePerGene(gtf, recombinationRateList)

# Write the output file to disk
write.table(GeneRecombinationRate, file = "/plas1/amardeep.singh/RNA.Seq.Data/RecombinationRateData/Gene.Recombination.Rates.txt", sep = "\t", col.name = T, row.names = F, quote = F)

# Remove all objects
rm(list=ls())

## Recomombination rate and SDIU/SBGE relationship
# --- R Code ---

rm(list = ls())
require(dplyr)
require(doBy)
require(ggplot2)
require(permuco)
require(car)

# Load in Recombination data and intersect with JunctionSeq output

## Loading in data files
# Read in JunctionSeq results
junctionseq.results.body = read.table("/plas1/amardeep.singh/RNA.Seq.Data/JunctionSeq.Files/BodyOutput/Aug1.Body.OnlyallGenes.results.txt", header = TRUE, sep = "\t")
junctionseq.results.head = read.table("/plas1/amardeep.singh/RNA.Seq.Data/JunctionSeq.Files/HeadOutput/Aug1.Head.OnlyallGenes.results.txt", header = TRUE, sep = "\t")
# Read in Tajimas D per gene results
recombination.rate.per.gene = read.table("/plas1/amardeep.singh/RNA.Seq.Data/RecombinationRateData/Gene.Recombination.Rates.txt", header = TRUE, sep = "\t")
## Load in the differential gene expression data
DGE.data.body = read.table("/plas1/amardeep.singh/RNA.Seq.Data/GeneExpression/RAL.DifferentialGeneExpression.body.txt", header = TRUE, sep = "\t")
DGE.data.head = read.table("/plas1/amardeep.singh/RNA.Seq.Data/GeneExpression/RAL.DifferentialGeneExpression.head.txt", header = TRUE, sep = "\t")

# Cleaning up files

# Remove any sites that were not tested
junctionseq.results.body = junctionseq.results.body[!(is.na(junctionseq.results.body$pvalue)),]
junctionseq.results.head = junctionseq.results.head[!(is.na(junctionseq.results.head$pvalue)),]
# Lets remove any rows where the expr in males or females is less than 10 in the junctionseq data
junctionseq.results.filtered.body = junctionseq.results.body[junctionseq.results.body$expr_male > 50 & junctionseq.results.body$expr_female > 50,]
junctionseq.results.filtered.head = junctionseq.results.head[junctionseq.results.head$expr_male > 50 & junctionseq.results.head$expr_female > 50,]
junctionseq.results.filtered.head = junctionseq.results.filtered.head[complete.cases(junctionseq.results.filtered.head[ ,1]),]
# Merge the junctionseq and Recombination rate data files and clean them up
junctionseq.results.body.merged = merge(junctionseq.results.filtered.body, recombination.rate.per.gene, by.x = "geneID", by.y = "geneID", sort = FALSE)
junctionseq.results.body.merged = junctionseq.results.body.merged[, c(1,14,25,29)]
junctionseq.results.body.merged = junctionseq.results.body.merged[!duplicated(junctionseq.results.body.merged[1:3]),]

junctionseq.results.head.merged = merge(junctionseq.results.filtered.head, recombination.rate.per.gene, by.x = "geneID", by.y = "geneID", sort = FALSE)
junctionseq.results.head.merged = junctionseq.results.head.merged[, c(1,14,25,29)]
junctionseq.results.head.merged = junctionseq.results.head.merged[!duplicated(junctionseq.results.head.merged[1:3]),]

# Filter out genes with fewer than 30 4fold sites
#junctionseq.results.body.merged = junctionseq.results.body.merged[junctionseq.results.body.merged$numberOf4FoldSites >= 30, ]
#junctionseq.results.head.merged = junctionseq.results.head.merged[junctionseq.results.head.merged$numberOf4FoldSites >= 30, ]

# Assign significant hits
junctionseq.results.body.merged$sig.hit = NA
junctionseq.results.body.merged$sig.hit[junctionseq.results.body.merged$geneWisePadj <= 0.01] = 1
junctionseq.results.body.merged$sig.hit[!(junctionseq.results.body.merged$geneWisePadj <= 0.01)] = 0

junctionseq.results.head.merged$sig.hit = NA
junctionseq.results.head.merged$sig.hit[junctionseq.results.head.merged$geneWisePadj <= 0.01] = 1
junctionseq.results.head.merged$sig.hit[!(junctionseq.results.head.merged$geneWisePadj <= 0.01)] = 0

# Collapse duplicates
junctionseq.results.body.merged.unique = unique(junctionseq.results.body.merged)
junctionseq.results.head.merged.unique = unique(junctionseq.results.head.merged)

# Merge Junctionseq/TajimasD data with DGE data
junctionseq.results.body.merged.unique = merge(junctionseq.results.body.merged.unique, DGE.data.body, by.x = "geneID", by.y = "FlyBaseID")
junctionseq.results.head.merged.unique = merge(junctionseq.results.head.merged.unique, DGE.data.head, by.x = "geneID", by.y = "FlyBaseID")

# Subset out the columns of interest (i.e., gene ID, TajimasD, significant differences in exon usage, and log2FC)
expression.data.body = junctionseq.results.body.merged.unique[, c(1:8)]
expression.data.head = junctionseq.results.head.merged.unique[, c(1:8)]

# Assign quantiles for sex-averaged mean expression
expression.data.body = expression.data.body %>% mutate(sex.averaged.expression.quantile = ntile(baseMean, 3))
expression.data.head = expression.data.head %>% mutate(sex.averaged.expression.quantile = ntile(baseMean, 3))

# Seperate out genes into male and female biased gene expression
male.biased.body = expression.data.body[expression.data.body$log2FoldChange > 0,]
female.biased.body = expression.data.body[expression.data.body$log2FoldChange < 0,]
male.biased.head = expression.data.head[expression.data.head$log2FoldChange > 0,]
female.biased.head = expression.data.head[expression.data.head$log2FoldChange < 0,]

# Remove any rows that have an NA added to them
male.biased.body = male.biased.body[!(is.na(male.biased.body$log2FoldChange)),]
female.biased.body = female.biased.body[!(is.na(female.biased.body$log2FoldChange)),]
male.biased.head = male.biased.head[!(is.na(male.biased.head$log2FoldChange)),]
female.biased.head = female.biased.head[!(is.na(female.biased.head$log2FoldChange)),]

# Assign quantiles for sex-averaged mean expression
#male.biased.body = male.biased.body %>% mutate(sex.averaged.expression.quantile = ntile(baseMean, 3))
#male.biased.head = male.biased.body %>% mutate(sex.averaged.expression.quantile = ntile(baseMean, 3))
#female.biased.body = female.biased.body %>% mutate(sex.averaged.expression.quantile = ntile(baseMean, 3))
#female.biased.head = female.biased.head %>% mutate(sex.averaged.expression.quantile = ntile(baseMean, 3))

# Assign quantile for both MBG and FBG within each sex-averaged expression quantile ## Careful in using this variable, it is only relevant when parsed by sex-averaged gene expression!
male.biased.body.lowexp = male.biased.body[male.biased.body$sex.averaged.expression.quantile == 1,] %>% mutate(sex.averaged.log2FC.quantile = ntile(log2FoldChange, 3))
male.biased.body.intermediateexp = male.biased.body[male.biased.body$sex.averaged.expression.quantile == 2,] %>% mutate(sex.averaged.log2FC.quantile = ntile(log2FoldChange, 3))
male.biased.body.highexp = male.biased.body[male.biased.body$sex.averaged.expression.quantile == 3,] %>% mutate(sex.averaged.log2FC.quantile = ntile(log2FoldChange, 3))
male.biased.head.lowexp = male.biased.head[male.biased.head$sex.averaged.expression.quantile == 1,] %>% mutate(sex.averaged.log2FC.quantile = ntile(log2FoldChange, 3))
male.biased.head.intermediateexp = male.biased.head[male.biased.head$sex.averaged.expression.quantile == 2,] %>% mutate(sex.averaged.log2FC.quantile = ntile(log2FoldChange, 3))
male.biased.head.highexp = male.biased.head[male.biased.head$sex.averaged.expression.quantile == 3,] %>% mutate(sex.averaged.log2FC.quantile = ntile(log2FoldChange, 3))

female.biased.body.lowexp = female.biased.body[female.biased.body$sex.averaged.expression.quantile == 1,] %>% mutate(sex.averaged.log2FC.quantile = ntile(log2FoldChange, 3))
female.biased.body.intermediateexp = female.biased.body[female.biased.body$sex.averaged.expression.quantile == 2,] %>% mutate(sex.averaged.log2FC.quantile = ntile(log2FoldChange, 3))
female.biased.body.highexp = female.biased.body[female.biased.body$sex.averaged.expression.quantile == 3,] %>% mutate(sex.averaged.log2FC.quantile = ntile(log2FoldChange, 3))
female.biased.head.lowexp = female.biased.head[female.biased.head$sex.averaged.expression.quantile == 1,] %>% mutate(sex.averaged.log2FC.quantile = ntile(log2FoldChange, 3))
female.biased.head.intermediateexp = female.biased.head[female.biased.head$sex.averaged.expression.quantile == 2,] %>% mutate(sex.averaged.log2FC.quantile = ntile(log2FoldChange, 3))
female.biased.head.highexp = female.biased.head[female.biased.head$sex.averaged.expression.quantile == 3,] %>% mutate(sex.averaged.log2FC.quantile = ntile(log2FoldChange, 3))

# For the male genes, lets add 3 to each quantile
male.biased.body.lowexp$sex.averaged.log2FC.quantile = male.biased.body.lowexp$sex.averaged.log2FC.quantile + 3
male.biased.body.intermediateexp$sex.averaged.log2FC.quantile = male.biased.body.intermediateexp$sex.averaged.log2FC.quantile + 3
male.biased.body.highexp$sex.averaged.log2FC.quantile = male.biased.body.highexp$sex.averaged.log2FC.quantile + 3
male.biased.head.lowexp$sex.averaged.log2FC.quantile = male.biased.head.lowexp$sex.averaged.log2FC.quantile + 3
male.biased.head.intermediateexp$sex.averaged.log2FC.quantile = male.biased.head.intermediateexp$sex.averaged.log2FC.quantile + 3
male.biased.head.highexp$sex.averaged.log2FC.quantile = male.biased.head.highexp$sex.averaged.log2FC.quantile + 3

# Combine back into male-biased and female-biased dataframes
male.biased.body = rbind(male.biased.body.lowexp, male.biased.body.intermediateexp,male.biased.body.highexp)
male.biased.head = rbind(male.biased.head.lowexp, male.biased.head.intermediateexp,male.biased.head.highexp)
female.biased.body = rbind(female.biased.body.lowexp, female.biased.body.intermediateexp,female.biased.body.highexp)
female.biased.head = rbind(female.biased.head.lowexp, female.biased.head.intermediateexp,female.biased.head.highexp)

# Assign quantile for both MBG and FBG overall
male.biased.body = male.biased.body %>% mutate(log2FC.quantile = ntile(log2FoldChange, 3))
male.biased.head = male.biased.head %>% mutate(log2FC.quantile = ntile(log2FoldChange, 3))
female.biased.body = female.biased.body %>% mutate(log2FC.quantile = ntile(log2FoldChange, 3))
female.biased.head = female.biased.head %>% mutate(log2FC.quantile = ntile(log2FoldChange, 3))

# For the male genes, lets add 3 to each quantile
male.biased.body$log2FC.quantile = male.biased.body$log2FC.quantile + 3
male.biased.head$log2FC.quantile = male.biased.head$log2FC.quantile + 3

# Merge data back
expression.data.body = rbind(male.biased.body,female.biased.body)
expression.data.head = rbind(male.biased.head,female.biased.head)
expression.data.body$tissue= "body"
expression.data.head$tissue= "head"

# Assign genes to predefined 'bins' of log2FC in gene expression
# Add upper and lower CIs for log2FC means
expression.data.body$log2FC.upperCI = expression.data.body$log2FoldChange + expression.data.body$lfcSE * 2
expression.data.body$log2FC.lowerCI = expression.data.body$log2FoldChange - expression.data.body$lfcSE * 2
expression.data.head$log2FC.upperCI = expression.data.head$log2FoldChange + expression.data.head$lfcSE * 2
expression.data.head$log2FC.lowerCI = expression.data.head$log2FoldChange - expression.data.head$lfcSE * 2

expression.data.body$log2FC.bins[expression.data.body$log2FC.upperCI < -2 ] = 1
expression.data.body$log2FC.bins[expression.data.body$log2FC.lowerCI >= -2 & expression.data.body$log2FC.upperCI <= -0.5 ] = 2
expression.data.body$log2FC.bins[expression.data.body$log2FC.lowerCI > -0.5 & expression.data.body$log2FC.upperCI < 0.5 ] = 3
expression.data.body$log2FC.bins[expression.data.body$log2FC.lowerCI >= 0.5 & expression.data.body$log2FC.upperCI <= 2 ] = 4
expression.data.body$log2FC.bins[expression.data.body$log2FC.lowerCI > 2 ] = 5

expression.data.head$log2FC.bins[expression.data.head$log2FC.upperCI < -2 ] = 1
expression.data.head$log2FC.bins[expression.data.head$log2FC.lowerCI >= -2 & expression.data.head$log2FC.upperCI <= -0.5 ] = 2
expression.data.head$log2FC.bins[expression.data.head$log2FC.lowerCI > -0.5 & expression.data.head$log2FC.upperCI < 0.5 ] = 3
expression.data.head$log2FC.bins[expression.data.head$log2FC.lowerCI >= 0.5 & expression.data.head$log2FC.upperCI <= 2 ] = 4
expression.data.head$log2FC.bins[expression.data.head$log2FC.lowerCI > 2 ] = 5

expression.data.recombination = rbind(expression.data.head, expression.data.body)

# Remove X-linked genes from recomnbination file for plotting
#recombination.rate.per.gene = recombination.rate.per.gene[complete.cases(recombination.rate.per.gene[,1]),]
recombination.rate.per.gene = expression.data.recombination[!(expression.data.recombination$chr.x=="X"),]

# For model fitting retain x and assign each gene to the X or autosome
expression.data.recombination$XorAutosome =  NA
# Assign a 1 to autosomal genes and a 0 to X-linked genes
expression.data.recombination$XorAutosome[!(expression.data.recombination$chr.x == "X")] = 1
expression.data.recombination$XorAutosome[expression.data.recombination$chr.x == "X"] = 0

###    Model Fitting    ###
# Fitting with lm()
recombination.model.body = lm(recombination.rate ~ baseMean +abs(log2FoldChange) + sig.hit, data = expression.data.recombination[expression.data.recombination$tissue == "body", ])
recombination.model.head = lm(recombination.rate ~ baseMean + abs(log2FoldChange) + sig.hit, data = expression.data.recombination[expression.data.recombination$tissue == "head", ])
# Model fitting with permuco
# This model should read 1,000,000 iterations for each effect w/interaction effects
#recombination.body.interaction.model = lmperm(recombination.rate ~ baseMean + abs(log2FoldChange) + sig.hit + abs(log2FoldChange)*sig.hit, data = expression.data.recombination[expression.data.recombination$tissue == "body", ], np=10000)
#recombination.head.interaction.model = lmperm(recombination.rate ~ baseMean + abs(log2FoldChange) + sig.hit + abs(log2FoldChange)*sig.hit, data = expression.data.recombination[expression.data.recombination$tissue == "head", ], np=10000)

recombination.body.model = lmperm(recombination.rate ~ baseMean + abs(log2FoldChange) + sig.hit + XorAutosome , data = expression.data.recombination[expression.data.recombination$tissue == "body", ], np=10000)
recombination.head.model = lmperm(recombination.rate ~ baseMean + abs(log2FoldChange) + sig.hit + XorAutosome, data = expression.data.recombination[expression.data.recombination$tissue == "head", ], np=10000)
recombination.body.model
recombination.head.model

#################
### Plotting #### -- Three different plots
#################

# 1. First plot is just averaging Tajima's D among significant and non-significant genes
recombinationRate.summary = summaryBy(recombination.rate ~ sig.hit + tissue, data = expression.data.recombination)

resample.data = as.data.frame(cbind(vector(mode="numeric", length = 10000 * nrow(recombinationRate.summary)),
                                    vector(mode="numeric", length = 10000 * nrow(recombinationRate.summary)),
                                    vector(mode="numeric", length = 10000 * nrow(recombinationRate.summary)),
                                    vector(mode="numeric", length = 10000 * nrow(recombinationRate.summary))))

output.row.start = 1
output.row.end = 4
for (i in 1:10000){

    resampled.df = rbind( sample_n(expression.data.recombination[expression.data.recombination$sig.hit == 0 & expression.data.recombination$tissue == "body", ], nrow(expression.data.recombination[expression.data.recombination$sig.hit == 0 & expression.data.recombination$tissue == "body", ]), replace = TRUE),
                         sample_n(expression.data.recombination[expression.data.recombination$sig.hit == 1 & expression.data.recombination$tissue == "body", ], nrow(expression.data.recombination[expression.data.recombination$sig.hit == 1 & expression.data.recombination$tissue == "body", ]), replace = TRUE),
                         sample_n(expression.data.recombination[expression.data.recombination$sig.hit == 0 & expression.data.recombination$tissue == "head", ], nrow(expression.data.recombination[expression.data.recombination$sig.hit == 0 & expression.data.recombination$tissue == "head", ]), replace = TRUE),
                         sample_n(expression.data.recombination[expression.data.recombination$sig.hit == 1 & expression.data.recombination$tissue == "head", ], nrow(expression.data.recombination[expression.data.recombination$sig.hit == 1 & expression.data.recombination$tissue == "head", ]), replace = TRUE))

                  # Output summary of resampled data to resample dataframe
                  resample.data[output.row.start:output.row.end, ] = summaryBy(recombination.rate ~ sig.hit + tissue, FUN=c(mean, length), data = resampled.df)
                  output.row.start = output.row.end  + 1
                  output.row.end = output.row.end + 4
                  print(i)

}
# Save outputs

# Rename columns of resampled data
colnames(resample.data) = c("sig.hit","tissue","RecombinationRate.mean","Number.of.Genes")

RecombinationRate.data.summary = summaryBy(RecombinationRate.mean ~ sig.hit + tissue, FUN = function(x) c(Mean = mean(x), lower.ci = quantile(x, probs = 0.05), upper.ci=quantile(x, probs = 0.95)), data = resample.data)
colnames(RecombinationRate.data.summary) = c("sig.hit","tissue","RecombinationRate.Mean","RecombinationRate.lowerCI","RecombinationRate.upperCI")

# Replace resample.data.summary "mean" with true observed mean
RecombinationRate.data.summary$RecombinationRate.Mean = recombinationRate.summary$recombination.rate.mean

# Plotting
recombination.summary.plot = ggplot(RecombinationRate.data.summary, aes(y = RecombinationRate.Mean, x = as.factor(tissue), colour = as.factor(sig.hit))) +
                                      geom_point(aes(y = RecombinationRate.Mean, x = as.factor(tissue), colour = as.factor(sig.hit)),size = 10,
                                      position = position_dodge(width = 0.5)) + scale_y_continuous(limits = c(1.75,2.5)) +
                                      geom_errorbar(aes(ymin = RecombinationRate.lowerCI, ymax = RecombinationRate.upperCI, x = as.factor(tissue), colour=as.factor(sig.hit)),
                                      width = 0, position = position_dodge(width = 0.5)) +
                                      theme_bw()  + scale_colour_manual(values = c("#800080", "#65c86e")) +
                                      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                                            panel.border = element_rect(colour = "black", fill=NA, size=1),
                                            axis.text = element_text(face="bold", color="black",size=20, family = "Helvetica"),
                                            axis.title=element_blank(), legend.position = "none")
pdf("/plas1/amardeep.singh/tmp/RecombinationRate.summary.plot.April6.pdf", height = 10, width = 5)
recombination.summary.plot
dev.off()


##
## 2B. Plotting Recombination Rate for all genes against log2FC without parsing out by sex-averaged gene expression   ######
##

# Remove unbinned genes
expression.data.recombination.plotting = expression.data.recombination[!(is.na(expression.data.recombination$log2FC.bins)),]

# Summary for means
expression.data.recombination.summary = summaryBy(recombination.rate ~ log2FC.bins + tissue, FUN=c(mean, length), data = expression.data.recombination.plotting)

# Bootstrap means for each subset
# Data frame to hold resampled results
resample.data = as.data.frame(cbind(vector(mode="numeric", length = 10000 * nrow(expression.data.recombination.summary)),
                                    vector(mode="numeric", length = 10000 * nrow(expression.data.recombination.summary)),
                                    vector(mode="numeric", length = 10000 * nrow(expression.data.recombination.summary)),
                                    vector(mode="numeric", length = 10000 * nrow(expression.data.recombination.summary)),
                                    vector(mode="numeric", length = 10000 * nrow(expression.data.recombination.summary))))

## Loop to resample TajimasD data
output.row.start = 1
output.row.end = 10
for (i in 1:10000){
  resampled.df=rbind( sample_n(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "body" & expression.data.recombination.plotting$log2FC.bins == 1, ], nrow(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "body" & expression.data.recombination.plotting$log2FC.bins == 1, ]), replace = TRUE),
                      sample_n(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "body" & expression.data.recombination.plotting$log2FC.bins == 2, ], nrow(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "body" & expression.data.recombination.plotting$log2FC.bins == 2, ]), replace = TRUE),
                      sample_n(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "body" & expression.data.recombination.plotting$log2FC.bins == 3, ], nrow(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "body" & expression.data.recombination.plotting$log2FC.bins == 3, ]), replace = TRUE),
                      sample_n(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "body" & expression.data.recombination.plotting$log2FC.bins == 4, ], nrow(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "body" & expression.data.recombination.plotting$log2FC.bins == 4, ]), replace = TRUE),
                      sample_n(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "body" & expression.data.recombination.plotting$log2FC.bins == 5, ], nrow(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "body" & expression.data.recombination.plotting$log2FC.bins == 5, ]), replace = TRUE),

                      sample_n(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "head" & expression.data.recombination.plotting$log2FC.bins == 1, ], nrow(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "head" & expression.data.recombination.plotting$log2FC.bins == 1, ]), replace = TRUE),
                      sample_n(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "head" & expression.data.recombination.plotting$log2FC.bins == 2, ], nrow(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "head" & expression.data.recombination.plotting$log2FC.bins == 2, ]), replace = TRUE),
                      sample_n(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "head" & expression.data.recombination.plotting$log2FC.bins == 3, ], nrow(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "head" & expression.data.recombination.plotting$log2FC.bins == 3, ]), replace = TRUE),
                      sample_n(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "head" & expression.data.recombination.plotting$log2FC.bins == 4, ], nrow(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "head" & expression.data.recombination.plotting$log2FC.bins == 4, ]), replace = TRUE),
                      sample_n(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "head" & expression.data.recombination.plotting$log2FC.bins == 5, ], nrow(expression.data.recombination.plotting[expression.data.recombination.plotting$tissue == "head" & expression.data.recombination.plotting$log2FC.bins == 5, ]), replace = TRUE))

                    # Output summary of resampled data to resample dataframe
                    resample.data[output.row.start:output.row.end, ] = summaryBy(recombination.rate ~ log2FC.bins + tissue, FUN=c(mean, length), data = resampled.df)
                    output.row.start = output.row.end  + 1
                    output.row.end = output.row.end + 10
                    print(i)
}

resample.data=resample.data[,1:4]
colnames(resample.data) = c("log2FC.bins","tissue","recombinationRate","number.of.genes")
write.table(resample.data, file = "/plas1/amardeep.singh/RNA.Seq.Data/RecombinationRateData/ResampleData/Recombination.resampled.By.log2FC.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

resample.data = read.delim("/plas1/amardeep.singh/RNA.Seq.Data/RecombinationRateData/ResampleData/Recombination.resampled.By.log2FC.txt", header = TRUE)

# Rename columns of resampled data
recombination.data.summary = summaryBy(recombinationRate ~ log2FC.bins + tissue, FUN = function(x) c(Mean = mean(x), lower.ci = quantile(x, probs = 0.05), upper.ci=quantile(x, probs = 0.95)), data = resample.data)

# Replace resample.data.summary "mean" with true observed mean
colnames(recombination.data.summary) = c("log2FC.bins","tissue","recombinationRate.Mean","recombinationRate.lowerCI","recombinationRate.upperCI")
recombination.data.summary$recombinationRate.Mean = expression.data.recombination.summary$recombination.rate.mean
# Add new column to denote MB, FB and unbiased genes
recombination.data.summary$sexbias[recombination.data.summary$log2FC.bins == 1 & recombination.data.summary$tissue == "body"] = "FB"
recombination.data.summary$sexbias[recombination.data.summary$log2FC.bins == 1 & recombination.data.summary$tissue == "head"] = "FB.1"
recombination.data.summary$sexbias[recombination.data.summary$log2FC.bins == 2] = "FB"
recombination.data.summary$sexbias[recombination.data.summary$log2FC.bins == 3] = "UB"
recombination.data.summary$sexbias[recombination.data.summary$log2FC.bins == 4] = "MB"
recombination.data.summary$sexbias[recombination.data.summary$log2FC.bins == 5 & recombination.data.summary$tissue == "body"] = "MB"
recombination.data.summary$sexbias[recombination.data.summary$log2FC.bins == 5 & recombination.data.summary$tissue == "head"] = "MB.1"

# Colours: darkred: #b30000; light red: #e09999; Dark Blue: #0000b3; light blue: #9999e0

recombination.expression.plot = ggplot(recombination.data.summary, aes(y = recombinationRate.Mean, x = as.factor(log2FC.bins), colour = as.factor(sexbias))) +
                                      geom_point(aes(y = recombinationRate.Mean, x = as.factor(log2FC.bins), colour = as.factor(sexbias)),size = 10,
                                      position = position_dodge(width = 0.5)) + #ylim(1,3.5) +
                                      geom_errorbar(aes(ymin = recombinationRate.lowerCI, ymax = recombinationRate.upperCI, x = as.factor(log2FC.bins), colour=as.factor(sexbias)),
                                      width = 0, position = position_dodge(width = 0.5)) +
                                      theme_bw()  + scale_colour_manual(values = c("#b30000","#e09999", "#0000b3","#9999e0", "#000000")) +
                                      #theme_bw()  + scale_colour_manual(values = c("#b30000","#e09999", "#0000b3", "#000000")) +
                                      #theme_bw()  + scale_colour_manual(values = c("#b30000", "#0000b3", "#000000")) +
                                      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                                            panel.border = element_rect(colour = "black", fill=NA, size=1),
                                            axis.text = element_text(face="bold", color="black",size=20, family = "Helvetica"),
                                            axis.title=element_blank(), legend.position = "none") +
                                      facet_grid(tissue ~., scales = "free")

pdf("/plas1/amardeep.singh/tmp/recombination.by.log2FC.expression.April6.pdf", height = 10, width = 10)
recombination.expression.plot
dev.off()



###
#3A.  Plotting Recombination Rate for SDIU and non-SDIU genes against sex-averaged gene expression without parsing by log2FC  ######
###

# Summary for means
expression.data.recombination.summary = summaryBy(recombination.rate ~ sex.averaged.expression.quantile + sig.hit + tissue, FUN=c(mean, length), data = expression.data.recombination)

# Bootstrap means for each subset
# Data frame to hold resampled results
resample.data = as.data.frame(cbind(vector(mode="numeric", length = 10000 * nrow(expression.data.recombination.summary)),
                                    vector(mode="numeric", length = 10000 * nrow(expression.data.recombination.summary)),
                                    vector(mode="numeric", length = 10000 * nrow(expression.data.recombination.summary)),
                                    vector(mode="numeric", length = 10000 * nrow(expression.data.recombination.summary)),
                                    vector(mode="numeric", length = 10000 * nrow(expression.data.recombination.summary))))

## Loop to resample TajimasD data
output.row.start = 1
output.row.end = 12
for (i in 1:10000){
  resampled.df=rbind( sample_n(expression.data.recombination[expression.data.recombination$sig.hit == 1 & expression.data.recombination$tissue == "body" & expression.data.recombination$sex.averaged.expression.quantile == 1, ], nrow(expression.data.recombination[expression.data.recombination$sig.hit == 1 & expression.data.recombination$tissue == "body" & expression.data.recombination$sex.averaged.expression.quantile == 1, ]), replace = TRUE),
                    sample_n(expression.data.recombination[expression.data.recombination$sig.hit == 0 & expression.data.recombination$tissue == "body" & expression.data.recombination$sex.averaged.expression.quantile == 1, ], nrow(expression.data.recombination[expression.data.recombination$sig.hit == 0 & expression.data.recombination$tissue == "body" & expression.data.recombination$sex.averaged.expression.quantile == 1, ]), replace = TRUE),
                    sample_n(expression.data.recombination[expression.data.recombination$sig.hit == 1 & expression.data.recombination$tissue == "body" & expression.data.recombination$sex.averaged.expression.quantile == 2, ], nrow(expression.data.recombination[expression.data.recombination$sig.hit == 1 & expression.data.recombination$tissue == "body" & expression.data.recombination$sex.averaged.expression.quantile == 2, ]), replace = TRUE),
                    sample_n(expression.data.recombination[expression.data.recombination$sig.hit == 0 & expression.data.recombination$tissue == "body" & expression.data.recombination$sex.averaged.expression.quantile == 2, ], nrow(expression.data.recombination[expression.data.recombination$sig.hit == 0 & expression.data.recombination$tissue == "body" & expression.data.recombination$sex.averaged.expression.quantile == 2, ]), replace = TRUE),
                    sample_n(expression.data.recombination[expression.data.recombination$sig.hit == 1 & expression.data.recombination$tissue == "body" & expression.data.recombination$sex.averaged.expression.quantile == 3, ], nrow(expression.data.recombination[expression.data.recombination$sig.hit == 1 & expression.data.recombination$tissue == "body" & expression.data.recombination$sex.averaged.expression.quantile == 3, ]), replace = TRUE),
                    sample_n(expression.data.recombination[expression.data.recombination$sig.hit == 0 & expression.data.recombination$tissue == "body" & expression.data.recombination$sex.averaged.expression.quantile == 3, ], nrow(expression.data.recombination[expression.data.recombination$sig.hit == 0 & expression.data.recombination$tissue == "body" & expression.data.recombination$sex.averaged.expression.quantile == 3, ]), replace = TRUE),

                    sample_n(expression.data.recombination[expression.data.recombination$sig.hit == 1 & expression.data.recombination$tissue == "head" & expression.data.recombination$sex.averaged.expression.quantile == 1, ], nrow(expression.data.recombination[expression.data.recombination$sig.hit == 1 & expression.data.recombination$tissue == "head" & expression.data.recombination$sex.averaged.expression.quantile == 1, ]), replace = TRUE),
                    sample_n(expression.data.recombination[expression.data.recombination$sig.hit == 0 & expression.data.recombination$tissue == "head" & expression.data.recombination$sex.averaged.expression.quantile == 1, ], nrow(expression.data.recombination[expression.data.recombination$sig.hit == 0 & expression.data.recombination$tissue == "head" & expression.data.recombination$sex.averaged.expression.quantile == 1, ]), replace = TRUE),
                    sample_n(expression.data.recombination[expression.data.recombination$sig.hit == 1 & expression.data.recombination$tissue == "head" & expression.data.recombination$sex.averaged.expression.quantile == 2, ], nrow(expression.data.recombination[expression.data.recombination$sig.hit == 1 & expression.data.recombination$tissue == "head" & expression.data.recombination$sex.averaged.expression.quantile == 2, ]), replace = TRUE),
                    sample_n(expression.data.recombination[expression.data.recombination$sig.hit == 0 & expression.data.recombination$tissue == "head" & expression.data.recombination$sex.averaged.expression.quantile == 2, ], nrow(expression.data.recombination[expression.data.recombination$sig.hit == 0 & expression.data.recombination$tissue == "head" & expression.data.recombination$sex.averaged.expression.quantile == 2, ]), replace = TRUE),
                    sample_n(expression.data.recombination[expression.data.recombination$sig.hit == 1 & expression.data.recombination$tissue == "head" & expression.data.recombination$sex.averaged.expression.quantile == 3, ], nrow(expression.data.recombination[expression.data.recombination$sig.hit == 1 & expression.data.recombination$tissue == "head" & expression.data.recombination$sex.averaged.expression.quantile == 3, ]), replace = TRUE),
                    sample_n(expression.data.recombination[expression.data.recombination$sig.hit == 0 & expression.data.recombination$tissue == "head" & expression.data.recombination$sex.averaged.expression.quantile == 3, ], nrow(expression.data.recombination[expression.data.recombination$sig.hit == 0 & expression.data.recombination$tissue == "head" & expression.data.recombination$sex.averaged.expression.quantile == 3, ]), replace = TRUE))

                    # Output summary of resampled data to resample dataframe
                    resample.data[output.row.start:output.row.end, ] = summaryBy(recombination.rate ~ sex.averaged.expression.quantile + sig.hit + tissue, FUN=c(mean, length), data = resampled.df)
                    output.row.start = output.row.end  + 1
                    output.row.end = output.row.end + 24
                    print(i)
}

colnames(resample.data) = c("sex.averaged.expression.quantile","sig.hit","tissue","RecombinationRate","number.of.genes")
write.table(resample.data, file = "/plas1/amardeep.singh/RNA.Seq.Data/RecombinationRateData/ResampleData/Recombination.resampled.By.sex.averaged.expression.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

# Load in Resample Data
resample.data = read.delim("/plas1/amardeep.singh/RNA.Seq.Data/RecombinationRateData/ResampleData/Recombination.resampled.By.sex.averaged.expression.txt", header = TRUE)

# Rename columns of resampled data
recombination.data.summary = summaryBy(RecombinationRate ~ sex.averaged.expression.quantile + sig.hit + tissue, FUN = function(x) c(Mean = mean(x), lower.ci = quantile(x, probs = 0.05), upper.ci=quantile(x, probs = 0.95)), data = resample.data)

# Replace resample.data.summary "mean" with true observed mean
colnames(recombination.data.summary) = c("sex.averaged.expression.quantile","sig.hit","tissue","recombinationRate.Mean","recombinationRate.lowerCI","recombinationRate.upperCI")
recombination.data.summary$recombinationRate.Mean = expression.data.recombination.summary$recombination.rate.mean

recombination.sex.averaged.expression.plot = ggplot(recombination.data.summary, aes(y = recombinationRate.Mean, x = as.factor(sex.averaged.expression.quantile), colour = as.factor(sig.hit))) +
                                      geom_point(aes(y = recombinationRate.Mean, x = as.factor(sex.averaged.expression.quantile), colour = as.factor(sig.hit)),size = 10,
                                      position = position_dodge(width = 0.5)) + #ylim(-0.05, 0.05) +
                                      geom_errorbar(aes(ymin = recombinationRate.lowerCI, ymax = recombinationRate.upperCI, x = as.factor(sex.averaged.expression.quantile), colour=as.factor(sig.hit)),
                                      width = 0, position = position_dodge(width = 0.5)) +
                                      theme_bw()  + scale_colour_manual(values = c("#800080", "#65c86e")) +
                                      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                                            panel.border = element_rect(colour = "black", fill=NA, size=1),
                                            axis.text = element_text(face="bold", color="black",size=20, family = "Helvetica"),
                                            axis.title=element_blank(), legend.position = "none") +
                                      facet_grid(tissue ~., scales = "free" )

pdf("/plas1/amardeep.singh/tmp/recombinationRate.by.sex.averaged.expression.quantile.filtered.April6.pdf", height = 10, width = 10)
recombination.sex.averaged.expression.plot
dev.off()
##

##


##
