#!/bin/bash
#PBS -l nodes=5:ppn=8
#PBS -o /disk/rdisk10/huaxinlei2/air_sample2/18-19-kraken2/data_txz/pier/vibrant_o
#PBS -e /disk/rdisk10/huaxinlei2/air_sample2/18-19-kraken2/data_txz/pier/vibrant_e
#PBS -N vibrant_pier
#PBS -l walltime=96:00:00
#PBS -m bea
#PBS -M huaxin_l@163.com

source deactivate
conda deactivate
conda activate vibrant

cd /disk/rdisk10/huaxinlei2/air_sample2/18-19-kraken2/data_txz/pier/
mkdir no_phage

python /disk/rdisk10/huaxinlei2/tools/VIBRANT/VIBRANT_run.py -i /disk/rdisk10/huaxinlei2/air_sample2/18-19-kraken2/data_txz/pier/contig/all_pier.fasta -t 96 -folder no_phage/pier_vibrant
