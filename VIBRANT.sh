#!/bin/bash


dir='file_list'
for i in $dir;
do

python ./VIBRANT_run.py -i ${i}.fasta -t 20 -folder ${i}_vibrant

done

## then remove viral contigs from the ${i}.fasta and generate the ${i}_no_phage.fasta file.
