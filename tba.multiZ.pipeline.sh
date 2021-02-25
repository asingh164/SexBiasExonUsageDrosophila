# tba/multiZ pipeline

## 1. Prepare the phylogenetic tree by converting from standard newick format to format necessary for tba by running the
#     bash program I wrote for this
bash PreparePhyloTree.sh -i ~/plas1/amardeep.singh/RNA.Seq.Data/HomologyData/drosophila-25spec-time-tree.tre

## 2. Use lastZ to make pairwise alignment files between the three target species
# Melanogaster - Simulans alignment
bash /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/Fasta.files.for.tba/LastZ.alignment.sh \
  -a /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/Fasta.files.for.tba/D.melanogaster \
  -b /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/Fasta.files.for.tba/D.simulans \
  -o /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/Fasta.files.for.tba/mel.sim.alignment &

# Melanogaster - Yakuba alignment
bash /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/Fasta.files.for.tba/LastZ.alignment.sh \
  -a /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/Fasta.files.for.tba/D.melanogaster \
  -b /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/Fasta.files.for.tba/D.yakuba \
  -o /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/Fasta.files.for.tba/mel.yak.alignment &

# Simulans - Yakuba alignment
bash /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/Fasta.files.for.tba/LastZ.alignment.sh \
  -a /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/Fasta.files.for.tba/D.simulans \
  -b /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/Fasta.files.for.tba/D.yakuba \
  -o /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/Fasta.files.for.tba/sim.yak.alignment &


  
tba "((Dsimulans Dyakuba))((Dmelanogaster))" *.*.maf tba.maf >&tba.log
