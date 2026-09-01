#!/bin/bash

dir='file_list'

mkdir metawrap


for i in $dir
do

#assembly
metawrap assembly -1 ${i}_clean_1.fastq -2 ${i}_clena_1.fastq -t 80 -o metawrap/${i}_assembly -m 1000

#binning and refinement
metawrap binning -o metawrap/${i}_binning -m 1000 -t 64 -a metawrap/${i}_assembly/final_assembly.fasta --metabat2 --maxbin2 --concoct ${i}_clean_1.fastq ${i}_clena_1.fastq
metawrap bin_refinement -m 1000 -o metawrap/${i}_refine -t 60 -m 1000 -A metawrap/${i}_binning/metabat2_bins/ -B metawrap/${i}_binning/maxbin2_bins/ -C metawrap/${i}_binning/concoct_bins/ -c 50 -x 10

done
