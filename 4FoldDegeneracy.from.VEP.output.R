##########################################################################
##########################################################################
#####               Four Fold Degenercy from VEP output              #####
##########################################################################
##########################################################################

# This program will take the output of vep and tell you which sites are four fold degenerate versus not
# Requires VEP output in standard VEP format (i.e., not a vcf)

## Identifying four-fold degenerate sites
#---- R CODE ----
rm(list=ls())

# Load in VEP outputs
vep.output.unmodified = read.delim("/plas1/amardeep.singh/RNA.Seq.Data/piNpiS.Analysis/DGRP.vep.output", header = TRUE)

## Clean up data file
# Remove duplicate sites in vep output
vep.output = vep.output.unmodified[!duplicated(vep.output.unmodified[c(2:4,7,12)]), ]
# Give the
# I only care about synonymous sites because the non-synonymous sites are obviously not 4-fold degenerate, so subset out only the synonymous sites
vep.output.syn = vep.output[ with(vep.output, grepl("synonymous_variant", Consequence)), ]
vep.output.nonsyn = vep.output[ with(vep.output, !(grepl("synonymous_variant", Consequence))), ]

# Subset out list of codons, convert to upper case and keep only the first three elements
vep.output.syn$ref.codon = substring( vep.output.syn$Codons, 1, 3)
vep.output.syn$fourfold.degenerate = NA
# Dictionary of fourfold degenerate codons
fourfold.sites.list = c("TCT","TCC","TCG","TCA","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG",
                   "CGT","CGC","CGA","CGG","ACT","ACC","ACA","ACG","GTT","GTC","GTA","GTG",
                   "GCT","GCC","GCA","GCG","GGT","GGC","GGA","GGG")

# First thing we want to do is figure out which position in a codon we are interested in based on which letter is uppercase
for (i in 1:nrow(vep.output.syn)){
  if ( (substring(vep.output.syn[i,"ref.codon"],3,3)) != "A" | (substring(vep.output.syn[i,"ref.codon"],3,3)) != "C" | (substring(vep.output.syn[i,"ref.codon"],3,3)) != "T" | (substring(vep.output.syn[i,"ref.codon"],3,3)) != "G"){
    if (as.character(toupper(vep.output.syn$ref.codon[i])) %in% fourfold.sites.list) {
      vep.output.syn$fourfold.degenerate[i] = 1
    } else {
      vep.output.syn$fourfold.degenerate[i] = 0
    }
  } else {
    vep.output.syn$fourfold.degenerate[i] = 0
  }
  print(i)
}

# Redirect all rows with four fold degenerate sites
fourfold.sites = vep.output.syn[vep.output.syn$fourfold.degenerate == 1,]

# Clean up file, join back to original file with new column signifying whether that site is four fold degenerate or not
vep.output.unmodified$Fourfold_degenerate = 0
vep.output.unmodified$Fourfold_degenerate[vep.output.unmodified$Location %in% fourfold.sites$Location] = 1

# Write file out
write.table(vep.output.unmodified, file = "/plas1/amardeep.singh/RNA.Seq.Data/piNpiS.Analysis/DGRP.vep.output.with.4fold.degeneracy",quote = F, sep = "\t", col.names = T, row.names = F)


#
