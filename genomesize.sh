#!/bin/bash
##PBS -l select=1:ncpus=16:host=compute-2-4+1:ncpus=12:host=compute-3-2+1:ncpus=13:host=compute-4-0+1:ncpus=16:host=compute-2-4
#PBS -l nodes=4:ppn=8
#PBS -o /disk/rdisk10/huaxinlei2/air_sample2/18-19-kraken2/data_txz/pier/genomesize_res
#PBS -e /disk/rdisk10/huaxinlei2/air_sample2/18-19-kraken2/data_txz/pier/genomesize_err
#PBS -N ags_pier
#PBS -l walltime=96:00:00
#PBS -m bea
#PBS -M huaxin_l@163.com
#pigz /disk/rdisk10/huaxinlei2/air_sample2/18-19-kraken2/decontam_18/2018_clean_reads/*fastq -p 80
#source deactivate
#conda deactivate
#conda activate base
#
##assembly based on metaspades, the files pathways must use the PWD
##time 


dir='
SL310995
SL342522
SL342528
'
cd /disk/rdisk10/huaxinlei2/air_sample2/18-19-kraken2/data_txz/all_clean_reads/

for i in $dir
do
python /disk/rdisk10/huaxinlei2/tools/MicrobeCensus/scripts/run_microbe_census.py -n 100000000 -t 90 ${i}_clean_paired_rm_euk_phage_1.fastq.gz,${i}_clean_paired_rm_euk_phage_1.fastq.gz ${i}_ags.txt
done


