#!/bin/bash
dir='file_list'
for i in $dir;
do
EukRep -i ${i}.fasta -o ${i}_Eukaryote_output --prokarya ${i}_no_euk.fasta
done
