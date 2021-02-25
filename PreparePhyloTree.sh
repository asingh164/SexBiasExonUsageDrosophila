#!/usr/bin/env bash

# This is program that will take a phylogentic tree in a traditiaonl Newick format and convert it to a format that can
# be passed to tba/multiz
# Written by Amardeep Singh -- Jan. 2021
# amardeep.singh [at] utexas.edu

# Pass the phylogenetic tree in Newick form to the -i flag. Currently only handles one file at a time


while getopts i: flag
do
    case "${flag}" in
        i) inputfile=${OPTARG};;
    esac
done

cat $inputfile | sed 's/[0-9]/ /g' | sed 's/://g' | sed 's/D\'.'/Drosophila/g'  | sed 's/\.//g' |
  sed 's/ //g' | sed 's/,/ /g' | sed 's/Drosophila/D./g' > $inputfile.tba.format.txt
