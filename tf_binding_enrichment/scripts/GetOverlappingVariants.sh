#!/bin/bash


ml BEDTools

reli_output=$1
locifile=$2


for statsfile in ${reli_output}/*.stats;
do

	filename=$(basename $statsfile)
	celltype="${filename%.*}"
	celltype="${celltype%.*}"
	echo "[INFO] prepping $celltype"

	awk '{print $1"\t"$2"\t"$3}'  /groups/umcg-wijmenga/tmp04/umcg-obbakker/tools/RELI/data/ChIP-seq/${celltype} | uniq > ${reli_output}/${celltype}.chipseq.tmp
	bedtools intersect -a ${reli_output}/${celltype}.chipseq.tmp -b ${locifile} -wo | awk -v var="$celltype" '{print var"\t"$0}' > ${reli_output}/${celltype}.RELI.snps

	rm ${reli_output}/${celltype}.chipseq.tmp
done
