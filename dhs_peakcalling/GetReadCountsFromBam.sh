#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1

ml SAMtools
ml BEDTools

base=$(basename $1)
base="${base%%.*}"

samtools view -q 30 $1 -o bams/${base}.filtered.bam

bedtools intersect \
-b bams/${base}.filtered.bam \
-a merged.collapsed.sorted.hotspots.bed \
-bed -wa | \
uniq -c | \
awk '{print $2":"$3"-"$4"\t"$2"\t"$3"\t"$4"\t"$1}' > \
counts/${base}.counts

