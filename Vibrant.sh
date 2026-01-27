#!/bin/bash

source deactivate
conda deactivate
conda activate vibrant

mkdir no_phage

python ./VIBRANT_run.py -i ${i}_assembly.fasta -t 96 -folder no_phage/${i}_vibrant
