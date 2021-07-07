#!/usr/bin/env bash

# This is a program to use lastZ for pair-wise sequence alignments. This program assumes that the first inputted fasta
#   file is the 'reference' file that you want both sequences to be aligned to and that the fasta file contains multiple
#   sequences.

while getopts "a:b:o:" FLAG
do
	case $FLAG in
		a) fastaA=${OPTARG} ;;		# This should be the target
		b) fastaB=${OPTARG} ;;
    o) outputFile=${OPTARG} ;;
		\?)
		 echo "Invalid option: -$OPTARG"

		 exit
			;;
	esac
done

#lastz $fastaA[multiple] $fastaA --notransition --step=20 --nogapped --format=axt > $outputFile.maf

lastz ${fastaA}[multiple] ${fastaB} --gapped --format=axt --step=19 --hspthresh=2200 --inner=2000 --ydrop=3400 --gappedthresh=10000 --chain > ${outputFile}.axt

#lastz_cmd = ('lastz_latest ' +
#             target_2bit + '[multiple][nameparse=darkspace] ' +
#             query_2bit + '[nameparse=darkspace] '
#             '--format=axt '
#             '--step=19 '
#             '--hspthresh=2200 '
#             '--inner=2000 '
#             '--ydrop=3400 '
#             '--gappedthresh=10000 '
#             '--scores=/fastdata/bop15hjb/bird_alignments/UCSC_pipeline/Scores/HoxD55 '
#             '--chain '
