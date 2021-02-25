#!/usr/bin/env bash

# This is a program to use lastZ for pair-wise sequence alignments. This program assumes that the first inputted fasta
#   file is the 'reference' file that you want both sequences to be aligned to and that the fasta file contains multiple
#   sequences.

while getopts "a:b:o:" FLAG
do
	case $FLAG in
		a) fastaA=${OPTARG} ;;
		b) fastaB=${OPTARG} ;;
    o) outputFile=${OPTARG} ;;
		\?)
		 echo "Invalid option: -$OPTARG"

		 exit
			;;
	esac
done

lastz $fastaA[multiple] $fastaA --notransition --step=20 --nogapped --format=maf > $outputFile.maf
