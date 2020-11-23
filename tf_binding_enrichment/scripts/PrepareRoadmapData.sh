#!/bin/bash

ponce=../data/summary_statistics/ponce_summary_stats_b37.bed
roadmap=../data/roadmap/coreMarks
output=../output

for file in $(ls $roadmap/*.filtered.bed);
do

cell=$(basename $file | sed 's/\_15\_coreMarks\_mnemonics.bed.filtered.bed//g')

echo "[INFO] Prepping $cell"

mkdir -p $output/$cell/clumping
mkdir -p $output/$cell/reli_input

sbatch PrepareReliDataDHS.sh \
$ponce \
$output/$cell/clumping \
$file

done
