#!/bin/bash



celltypeInfo=../parameter_files/immune_cells_sorted.txt
reliDir=$1


cell=$(basename $reliDir)
cd $reliDir

echo "Id" > ids.tmp
ls -1 reli_output/*.stats | sed 's/reli_output\///g' | sed 's/.RELI.stats//g' >> ids.tmp

#head -n 1 reli_output/hg19_1268.RELI.stats > merged.tmp

echo -e "FormalPhenotype	Ancestry	Source	Cell	FormalCell	Label	Intersect	Total	Ratio	Mean	Std	Zscore	RelativeRisk	Pval	CorrectedPval	NullModel	Species" > merged.tmp
tail -q -n +2  reli_output/*.stats >> merged.tmp

paste ids.tmp merged.tmp > output_merged.tmp


echo "BroadType" > broad_type.tmp
awk '{print $7}' $celltypeInfo >> broad_type.tmp

paste output_merged.tmp broad_type.tmp > reli_output_merged.tmp


# Reformat
awk 'BEGIN{OFS="\t"} { print $1,$19,$5,$6,$4,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$18,$2,$3}'  reli_output_merged.tmp >  ${cell}_reli_output_merged.tsv




rm *.tmp





