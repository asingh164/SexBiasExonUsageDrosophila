# Pull coordinates for genes in the D. mel reference genome

# Load packages

# Take old gtf and update the coordinates to the old genome



# GTF File location
gtf = read.delim("/plas1/amardeep.singh/Flybase.Dmel.Genome.Release/release5.26_gtf/dmel-all-r5.36.gtf", header = FALSE)
recombinationRateList = read.delim("/plas1/amardeep.singh/apps/Comeron_tables/Comeron_100kb_allchr.txt", header = FALSE)

GetGeneLocations<-function(gtf,vectorOfGenes){
  # Load necessary packages
  require(stringr)

  #Split gtf information into seperate columns
    gtf.gene.names = substring((str_split_fixed(gtf$V9, ";", 3)[,2]), 10,20)
    gtf.file = cbind(gtf[,c(1,3:5)], gtf.gene.names)
    gtf.file$V1 = paste("chr", gtf.file$V1, sep = "")
    colnames(gtf.file) = c("chr","feature","start.pos","end.pos","geneID")
    vectorOfGenes = unique(gtf.file$geneID)

    #set up output object
    output = as.data.frame(matrix(NA, ncol = 4, nrow = length(gtf.gene.names)))
    colnames(output) = c("geneID", "start.pos", "end.pos", "coordinates.formatted")

    # Add end coordinates to recombination rate list
    recombinationRateList$V4 = recombinationRateList$V2 + 99999
counter = 0
counter2 = 0
    for (i in 1:length(vectorOfGenes)){
      gene = vectorOfGenes[i]
      gtf.sub = gtf.file[gtf.file$geneID == gene,]
      gene.range = c(gtf.sub$chr[1], c(range(gtf.sub$start.pos)[1],range(gtf.sub$end.pos)[2]) )

      recombinationRateList.sub = recombinationRateList[(recombinationRateList$V1 == gene.range[1] & recombinationRateList$V2 <= as.numeric(gene.range[2]) & recombinationRateList$V2 >= (as.numeric(gene.range[2]) - 100000) | recombinationRateList$V1 == gene.range[1] & recombinationRateList$V2 <= as.numeric(gene.range[3]) & recombinationRateList$V4 >= as.numeric((gene.range[3]))), ]

      if (nrow(recombinationRateList.sub) > 1) {
        counter = counter + 1} else
        { counter2 = counter2 + 1
        }
        print(i)
}
      # Find 100kb recombination range


      # output gene range data to output object
      output[i,] = c(gene,gene.range[1],gene.range[2],paste())

    }

    output.df =

  gtf = gtf[,]

}



# Load in list of genes to query
