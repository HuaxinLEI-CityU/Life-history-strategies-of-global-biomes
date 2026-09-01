#!/bin/bash
dir='file_list'
for i in $dir;
do
EukRep -i ${i}_marine.fasta -o ${i}_marine_Eukaryote_output --prokarya ${i}_marine_no_euk.fasta
done
