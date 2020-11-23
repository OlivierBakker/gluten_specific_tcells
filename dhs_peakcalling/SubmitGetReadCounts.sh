#!/bin/bash


for bam in bams/*.bam;
do
	echo "$bam"
	sbatch GetReadCountsFromBam.sh $bam
done
