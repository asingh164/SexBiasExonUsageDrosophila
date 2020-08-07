#!/bin/bash
## Originally written by Tyler Kent 2017

#take gff, choose annotation type, output bed-like file
#only argument is annotation type as a string

sed -e 's/;/\t/g' -e 's/:/\t/g' $1 | awk -v var=$2 'BEGIN {OFS="\t"}{if ($3==var) print $1,($4-1),$5,$7,$8,$10}' |\
sort -k1,1 -k2,2n
