# tba/multiZ pipeline

# This pipeline will take fasta files from the species we are interested in aligning (focusing on Dmel and Dsim) and
# producing MAF files and then converting MAF to VCf files

## 1. Prepare the FASTA references for each species
sed -i 's/^>/>Dmel./' D.melanogaster
sed -i 's/^>/>Dsim./' D.simulans

# Prep D. mel fasta file
# Lets pull the line numbers for each chromosome arm from the fast
export Chr2L_LINES_TMP=`cat /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/FastaReferenceSequences/D.melanogaster | sed -n -e '/^>Dmel.2L type/,/^>Dmel./ p' | wc -l | cut -f1 -d' ' `
export Chr2R_LINES_TMP=`cat /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/FastaReferenceSequences/D.melanogaster | sed -n -e '/^>Dmel.2R type/,/^>Dmel./ p' | wc -l | cut -f1 -d' ' `
export Chr3L_LINES_TMP=`cat /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/FastaReferenceSequences/D.melanogaster | sed -n -e '/^>Dmel.3L type/,/^>Dmel./ p' | wc -l | cut -f1 -d' ' `
export Chr3R_LINES_TMP=`cat /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/FastaReferenceSequences/D.melanogaster | sed -n -e '/^>Dmel.3R type/,/^>Dmel./ p' | wc -l | cut -f1 -d' ' `
export ChrX_LINES_TMP=`cat /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/FastaReferenceSequences/D.melanogaster | sed -n -e '/^>Dmel.X type/,/^>Dmel./ p' | wc -l | cut -f1 -d' ' `

# Assign chromosome line numbers to global variable
export Chr2L_LINES=`expr ${Chr2L_LINES_TMP} - 1`
export Chr2R_LINES=`expr ${Chr2R_LINES_TMP} - 1`
export Chr3L_LINES=`expr ${Chr3L_LINES_TMP} - 1`
export Chr3R_LINES=`expr ${Chr3R_LINES_TMP} - 1`
export ChrX_LINES=`expr ${ChrX_LINES_TMP} - 1`

# Pull each chromosome arm and assign to seperate file
cat /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/FastaReferenceSequences/D.melanogaster | sed -n -e '/^>Dmel.2L type/,/^>Dmel./ p' | head -n $Chr2L_LINES > Dmel_Chr2L.fasta
cat /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/FastaReferenceSequences/D.melanogaster | sed -n -e '/^>Dmel.2R type/,/^>Dmel./ p' | head -n $Chr2R_LINES > Dmel_Chr2R.fasta
cat /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/FastaReferenceSequences/D.melanogaster | sed -n -e '/^>Dmel.3L type/,/^>Dmel./ p' | head -n $Chr3L_LINES > Dmel_Chr3L.fasta
cat /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/FastaReferenceSequences/D.melanogaster | sed -n -e '/^>Dmel.3R type/,/^>Dmel./ p' | head -n $Chr3R_LINES > Dmel_Chr3R.fasta
cat /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/FastaReferenceSequences/D.melanogaster | sed -n -e '/^>Dmel.X type/,/^>Dmel./ p' | head -n $ChrX_LINES > Dmel_ChrX.fasta

# Concatenate all of the chromosome arms into a final fasta file
cat Dmel_Chr2L.fasta Dmel_Chr2R.fasta Dmel_Chr3L.fasta Dmel_Chr3R.fasta Dmel_ChrX.fasta > Dmel.modified.fasta

# 2. Use Santiago's version of lastZ to make pairwise alignment files between the three target species
export DMEL_FASTA_FILE=/plas1/amardeep.singh/RNA.Seq.Data/HomologyData/FastaReferenceSequences/D.melanogaster
export DSIM_FASTA_FILE=/plas1/amardeep.singh/RNA.Seq.Data/HomologyData/FastaReferenceSequences/D.simulans
export DMEL_DSIM_OUTPUT=/plas1/amardeep.singh/RNA.Seq.Data/HomologyData/Fasta.files.for.tba/output/Dmel.Dsim.May2.maf
# Run LastZ
/plas1/amardeep.singh/RNA.Seq.Data/HomologyData/Fasta.files.for.tba/lastzWrapper -sp1 ${DMEL_FASTA_FILE} -sp2 ${DSIM_FASTA_FILE} -threads 60 > ${DMEL_DSIM_OUTPUT}

# 3. Generate a vcf file from the MAF files using MAF FILTER
# Add this script to an 'option_file' and run the command maffilter param=option_file
DATA=/plas1/amardeep.singh/RNA.Seq.Data/HomologyData/Fasta.files.for.tba/output/Dmel.Dsim.May2_alt //A user-defined variable, representing the input maf file, without extension
input.file=$(DATA).maf.gz  //Input maf file, gzipped
input.file.compression=gzip
output.log=$(DATA).maffilter.log //Output log file
maf.filter=VcfOutput(file=./DmelDsim.alt.vcf.gz,compression=gzip,reference=Dmel,genotypes=(Dmel,Dsim),all=yes) //A coma separated list of filters, here empty (the program only read the input file and exits)

#4. Filter MAF file for any alignment scores below 11793 (75th quantile)
maf.filter=QualFilter(file = ./test.gz, compression=gzip,species=(Dmel,Dsim), window.size=10, window.step=1,                      \
        min.qual=0.8,                       \
        file=data.trash_qual.maf.gz,        \
        compression=gzip),                  \
    [...]




############## D E P R E C A T E D     C O D E    ##############
# maffilter param=option_file

# Pull MAF file alignment scores and alignment lengths
# Scores
zcat Dmel.Dsim.May2_alt.maf.gz | grep '^a' | cut -f2 | cut -c 9- > Dmel.Dsim.May2_alt.maf.scores
# Alignment lengths in Dmel
zcat Dmel.Dsim.May2_alt.maf.gz | grep 'Dmel' | sed 's/ /:/g' | tr -s ':' | sed 's/:/   /g'> Dmel.Alignment.Lengths.maf
# Alignment lengths in Dsim
zcat Dmel.Dsim.May2_alt.maf.gz | grep 'Dsim' | sed 's/ /:/g' | tr -s ':' | sed 's/:/   /g'> Dsim.Alignment.Lengths.maf
# combine files into a single file to analyze
paste Dmel.Dsim.May2_alt.maf.scores Dmel.Alignment.Lengths.maf Dsim.Alignment.Lengths.maf > DmelDsim_Alignment_Scores_and_Length.txt




## 2. Use lastZ to make pairwise alignment files between the three target species
# Melanogaster - Simulans alignment
bash /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/Fasta.files.for.tba/LastZ.alignment.sh \
  -a /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/FastaReferenceSequences/D.melanogaster \
  -b /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/FastaReferenceSequences/D.simulans \
  -o /plas1/amardeep.singh/RNA.Seq.Data/HomologyData/Fasta.files.for.tba/output/mel.sim.alignment &

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

## Use MAFFILTER to convert the .maf file to .vcf

tba "((Dsimulans Dyakuba))((Dmelanogaster))" *.*.maf tba.maf >&tba.log
