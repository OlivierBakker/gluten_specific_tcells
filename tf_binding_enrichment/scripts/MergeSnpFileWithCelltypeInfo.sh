#!/bin/bash



rm $2.tmp
for id in $(awk '{print $1}' $1);
do
	grep $id ../parameter_files/immune_cells.txt >> $2.tmp

done

paste $1 $2.tmp > $2.withCellInfo

rm $2.withSnp
for snp in $(awk '{print $8}' $2.withCellInfo);
do
	tmp="$(grep "\\s"$snp"\\s" ../data/summary_statistics/ponce_summary_stats_b37.bed)"
	ld_block="$(grep -P "\\s"$snp"\\s{0,1}" $3 | awk '{print $1}' | sed 's/://g')"

	if [ -z "$tmp" ]
	then
		tmp="NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA"
	fi
	echo -e "$tmp	$ld_block" >> $2.withSnp

done

paste $2.withCellInfo $2.withSnp > $2

rm $2.tmp $2.withCellInfo $2.withSnp
