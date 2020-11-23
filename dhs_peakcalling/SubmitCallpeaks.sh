#!/bin/bash


mkdir -p logs


for bam in data/bams/*.bam;
do

FILE=$(basename -- "$bam")
BASE="${FILE%.*}"


sbatch --out="logs/${BASE}.out" --err="logs/${BASE}.err" Callpeaks.sh ${bam}

done
