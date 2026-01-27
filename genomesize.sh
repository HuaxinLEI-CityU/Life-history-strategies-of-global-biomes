#!/bin/bash

dir='file_list'
cd $pathway

for i in $dir
do
python ./run_microbe_census.py -n 100000000 -t 90 ${i}_no_euk_phage.gasta ${i}_ags.txt
done
