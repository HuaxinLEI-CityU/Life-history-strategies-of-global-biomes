#!/bin/bash

for i in $dir
do

python /gpfs1/home/huaxinlei2/scratch/VIBRANT/VIBRANT_run.py -i ${i}_final_assembly.fasta -t 20 -folder ${i}_vibrant


done
#!/bin/bash