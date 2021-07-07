# This is a pipeline to use Mercator to generate 1-to-1 orthologs
# This pipeline assumes you have Mercator and BLAT installed and that you have a fasta and gff for your references

# Export file paths
export PATH_TO_Dmel_GFF=/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/gtf/Drosophila_melanogaster.BDGP6.32.51.chr.gff3.gz
export PATH_TO_Dsim_GFF=/plas1/amardeep.singh/Ensembl.Dsim.Genome.Release/gff/Drosophila_simulans.ASM75419v3.51.chr.gff3.gz
export PATH_TO_Dmel_FASTA=/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/fasta/Drosophila_melanogaster.BDGP6.32.dna_rm.toplevel.fa.gz
export PATH_TO_Dsim_FASTA=/plas1/amardeep.singh/Ensembl.Dsim.Genome.Release/fasta/Drosophila_simulans.ASM75419v3_.dna_rm.toplevel.fa.gz
# Specify Mercator utility paths
#export makeMercatorInput=/plas1/apps/cndsrc-2017.11.30/bin/makeMercatorInput.py
export makeMercatorInput_AmardeepsEdit=/plas1/amardeep.singh/apps/cndsrc-2017.11.30/bin/makeMercatorInput_AmardeepEdit.py
#export blat2hits="python2 /plas1/apps/cndsrc-2017.11.30/bin/blat2hits.py"
export blat2hits_AmardeepsEdit="python2 /plas1/amardeep.singh/apps/cndsrc-2017.11.30/bin/blat2hits_AmardeepEdit.py"
export phits2constraints=/plas1/amardeep.singh/apps/cndsrc-2017.11.30/bin/phits2constraints.py
export makeAlignmentInput=/plas1/amardeep.singh/apps/cndsrc-2017.11.30/bin/makeAlignmentInput
## Create input files for Mercator

# 1. Generate a .sbd index file for each fasta file to pass to mercator
mkdir -p ./MercatorAnalysis/FastaIndexFiles
mkdir -p ./MercatorAnalysis/GFFs
mkdir -p ./MercatorAnalysis/MercatorInputFiles
mkdir -p ./MercatorAnalysis/MercatorOutputFiles

# Export paths
zcat ${PATH_TO_Dmel_FASTA} | faReformat --strip-title | fa2sdb -c ./MercatorAnalysis/FastaIndexFiles/Dmel.sdb
zcat ${PATH_TO_Dsim_FASTA} | faReformat --strip-title | fa2sdb -c ./MercatorAnalysis/FastaIndexFiles/Dsim.sdb

# 2. Generate genome annotation file  (gff/gtfs)
#   Already done and paths listed above
# Move files into files with names to match fasta index files
cp ${PATH_TO_Dmel_GFF} ./MercatorAnalysis/GFFs/Dmel.gff.gz && gunzip ./MercatorAnalysis/GFFs/Dmel.gff.gz
cp ${PATH_TO_Dsim_GFF} ./MercatorAnalysis/GFFs/Dsim.gff.gz && gunzip ./MercatorAnalysis/GFFs/Dsim.gff.gz

# 3. Generate Mercator input files
cd MercatorAnalysis/
python2 ${makeMercatorInput_AmardeepsEdit} --genome-dir=./FastaIndexFiles --gff-dir=./GFFs --out-dir=./MercatorInputFiles Dmel Dsim
# The script seems to have an issue making the '.hits' file. I think I can generate it myself using the .blat file as follows
cat Dmel-Dsim.blat | awk 'BEGIN {OFS = FS = "\t"} {print $1 "\t" $2 "\t" $3 "\t" $11}' > Dmel-Dsim.hits

# 4. Generate constraints file with Mercator
python2 ${phits2constraints} --input-dir=/plas1/amardeep.singh/RNA.Seq.Data/HomologyData/MercatorAnalysis/MercatorInputFiles < /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/MercatorAnalysis/MercatorInputFiles/Dmel-Dsim.hits > constraints

#5. Generate alignment inputs with Mercaotor
makeAlignmentInput . ALIGNMENTDIR

#6. Local alignment with MAVID
python2 /plas1/amardeep.singh/apps/cndsrc-2017.11.30/bin/mavidAlignDirs.py  --skip-completed --init-dir=ALIGNMENTDIR


# 6. Local alignments of syntenic regions with fsa
fsa --fast --log 1 --mercator /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/MercatorAnalysis/MercatorOutputFiles/constraints /plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/fasta/Dmel.fa /plas1/amardeep.singh/Ensembl.Dsim.Genome.Release/fasta/Dsim.fa > DmelDsim_Alignment.mfa

fsa --fast --log 1 --logcopy fsa.test.log --mercator /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/MercatorAnalysis/MercatorOutputFiles/map /plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/fasta/Dmel.fa /plas1/amardeep.singh/Ensembl.Dsim.Genome.Release/fasta/Dsim.fa > DmelDsim_test_Alignment.mfa


vcf_subset=vcf[!(vcf$V5=="*" && vcf$V10== 1,)]



#
