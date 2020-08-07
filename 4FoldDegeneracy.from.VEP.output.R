##########################################################################
##########################################################################
#####               Four Fold Degenercy from VEP output              #####
##########################################################################
##########################################################################

# This program will take the output of vep and tell you which sites are four fold degenerate versus not
# Requires VEP output in standard VEP format (i.e., not a vcf)

## Identifying 4fold and 0fold degenerate sites
#---- R CODE ----
rm(list=ls())

# Load in VEP outputs
vep.output.unmodified = read.delim("/plas1/amardeep.singh/RNA.Seq.Data/piNpiS.Analysis/DGRP.vep.output", header = TRUE)

## Clean up data file
# Remove duplicate sites in vep output
vep.output = vep.output.unmodified[!duplicated(vep.output.unmodified[c(2:4,7,12)]), ]

# Seperate out synonymous and non-synonymous variants -- I know that any site that has a synonymous variants is NOT 0fold
#   degenerate and any site with a non-synonymous variant is not 4fold degenerate
vep.output.syn = vep.output[ with(vep.output, grepl("synonymous_variant", Consequence)), ]
vep.output.nonsyn = vep.output[ with(vep.output, !(grepl("synonymous_variant", Consequence))), ]

# Subset out list of codons, convert to upper case and keep only the first three elements
vep.output.syn$ref.codon = substring( vep.output.syn$Codons, 1, 3)
vep.output.syn$fourfold.degenerate = NA

vep.output.nonsyn$ref.codon = substring( vep.output.nonsyn$Codons, 1, 3)
vep.output.nonsyn$zerofold.firstposition.degenerate = NA
vep.output.nonsyn$zerofold.secondposition.degenerate = NA

# Dictionary of 4fold degenerate codons and 0fold degenerate sites at position 1 and 2
third.site.4fold = c("TCT","TCC","TCG","TCA","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG",
                   "CGT","CGC","CGA","CGG","ACT","ACC","ACA","ACG","GTT","GTC","GTA","GTG",
                   "GCT","GCC","GCA","GCG","GGT","GGC","GGA","GGG")

# If first site is variable the following codons are 0-fold degenerate
first.site.0fold = c("TTT","TTC","CTT","CTC","ATT","ATC","ATG","GTT","GTC",
                    "TCT","TCC","TCA","TCG","CCT","CCC","CCA","CCG","ACT",
                    "ACC","ACA","ACG","GCT","GCC","GCA","GCG","TAT","CAT",
                    "AAT","GAT","TAC","CAC","AAC","GAC","TAA","CAA","AAA",
                    "GAA","TAG","CAG","AAG","GAG","TGT","CGT","AGT","GGT",
                    "TGC","CGC","AGC","GGC")

# If second site is variable the following codons are 0-fold degenerate
second.site.0fold = c("TTT","TCT","TAT","TGT","TTC","TCC","TAC","TGC","TTG",
                      "TCG","TAG","TGG","CTT","CCT","CAT","CGT","CTC","CCC",
                      "CAC","CGC","CTA","CCA","CAA","CGA","CTG","CCG","CAG",
                      "CGG","ATT","ACT","AAT","AGT","ATC","ACC","AAC","AGC",
                      "ATA","ACA","AAA","AGA","ATG","ACG","AAG","AGG","GTT",
                      "GCT","GAT","GGT","GTC","GCC","GAC","GGC","GTA","GCA",
                      "GAA","GGA","GTG","GCG","GAG","GGG")


# Lets deal with identifying 4fold degenerate sites first
# First thing we want to do is figure out which position in a codon we are interested in based on which letter is uppercase
for (i in 1:nrow(vep.output.syn)){
  if ( (substring(vep.output.syn[i,"ref.codon"],3,3)) == "A" | (substring(vep.output.syn[i,"ref.codon"],3,3)) == "C" | (substring(vep.output.syn[i,"ref.codon"],3,3)) == "T" | (substring(vep.output.syn[i,"ref.codon"],3,3)) == "G"){
    if (as.character(toupper(vep.output.syn$ref.codon[i])) %in% third.site.4fold) {
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
rm(i)

# Next, lets deal with identifying 0fold degenerate sites
# First thing we want to do is figure out which position in a codon we are interested in based on which letter is uppercase
for (i in 1:nrow(vep.output.nonsyn)){
  # Start by looking at position 1
  if ( (substring(vep.output.nonsyn[i,"ref.codon"],1,1)) == "A" | (substring(vep.output.nonsyn[i,"ref.codon"],1,1)) == "C" | (substring(vep.output.nonsyn[i,"ref.codon"],1,1)) == "T" | (substring(vep.output.nonsyn[i,"ref.codon"],1,1)) == "G" ){
    if (as.character(toupper(vep.output.nonsyn$ref.codon[i])) %in% first.site.0fold) {
      vep.output.nonsyn$zerofold.firstposition.degenerate[i] = 1
    } else {
      vep.output.nonsyn$zerofold.firstposition.degenerate[i] = 0
    }
  } else if ( (substring(vep.output.nonsyn[i,"ref.codon"],2,2)) == "A" | (substring(vep.output.nonsyn[i,"ref.codon"],2,2)) == "C" | (substring(vep.output.nonsyn[i,"ref.codon"],2,2)) == "T" | (substring(vep.output.nonsyn[i,"ref.codon"],2,2)) == "G" ){
    if (as.character(toupper(vep.output.nonsyn$ref.codon[i])) %in% second.site.0fold) {
      vep.output.nonsyn$zerofold.secondposition.degenerate[i] = 1
    } else {
      vep.output.nonsyn$zerofold.secondposition.degenerate[i] = 0
    }
  } else {}
    print(i)
}

# Redirect all rows with four fold degenerate sites
fourfold.sites = vep.output.syn[vep.output.syn$fourfold.degenerate == 1,]
zerofold.sites = vep.output.nonsyn[vep.output.nonsyn$zerofold.firstposition.degenerate == 1 | vep.output.nonsyn$zerofold.secondposition.degenerate == 1,]
# Remove rows with NA
zerofold.sites = zerofold.sites[!(is.na(zerofold.sites$Location)),]

# Clean up file, join back to original file with new column signifying whether that site is four fold degenerate or not
vep.output.unmodified$Degeneracy =  NA
vep.output.unmodified$Fourfold_degenerate

vep.output.unmodified$Degeneracy[vep.output.unmodified$Location %in% fourfold.sites$Location] = 4
vep.output.unmodified$Degeneracy[vep.output.unmodified$Location %in% zerofold.sites$Location] = 0

# Write file out
write.table(vep.output.unmodified, file = "/plas1/amardeep.singh/RNA.Seq.Data/piNpiS.Analysis/DGRP.vep.output.with.degeneracy",quote = F, sep = "\t", col.names = T, row.names = F)


#

## 0 fold Sites







##
