#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=21G
#SBATCH --time=00:30:00


set -e 

# INPUT
MEM=20000
DATA_DIR=/groups/umcg-wijmenga/tmp04/umcg-obbakker/data/
REFERENCE=$DATA_DIR/reference/1000G/filtered/phase3/EUR/no_hla/1000G_phase3_snps_only_no_hla
SUMMARY=$1

# OUTPUT
CLUMP_DIR=$2

module load plink

FILE=$(basename "$SUMMARY")
FILE="${FILE%.*}"
echo "-------------------------------------------------------------"
echo "[INFO]	Proccessing file $FILE"
echo "-------------------------------------------------------------"

for SNP_THRESH in 5e-8
do
	FILE_2=$(basename "$SNP_THRESH")
	FILE_2="${FILE}_clumped_${FILE_2%.*}"

	echo "-------------------------------------------------------------"
	echo "[INFO]	Clump SNPS based on LD"
	echo "-------------------------------------------------------------"

	echo "[INFO] Reference population: " $(basename $REFERENCE)
	echo "[INFO] Calculating clumps for $FILE and threshold $FILE_2"

	mkdir -p $CLUMP_DIR
	cd $CLUMP_DIR
	

	plink -bfile $REFERENCE --memory $MEM --out $FILE_2 --clump $SUMMARY --clump-p1 $SNP_THRESH --clump-p2 $SNP_THRESH --clump-r2 0 --clump-kb 1000 --clump-field pvalue --clump-snp-field rsId37
	awk '{print $3}' < $FILE_2.clumped | tail -n +2 > $FILE_2.snps


	awk '{print $12}' < $FILE_2.clumped  | tail -n +2 | sed "s/,/\\n/g" | sed "s/(.*)//g" | grep -v NONE > $FILE_2.allClumps
	cat $FILE_2.snps >> $FILE_2.allClumps


	#plink -bfile $REFERENCE --memory $MEM --extract $FILE_2.allClumps --out $FILE_2 --show-tags $FILE_2.snps --tag-r2 0.8 --tag-kb 1000 --list-all
        plink -bfile $REFERENCE --memory $MEM --out $FILE_2 --show-tags $FILE_2.snps --tag-r2 0.8 --tag-kb 1000 --list-all


	awk '{print $1":\t"$1"|"$8}' $FILE_2.tags.list | tail -n +2 | sed 's/|/\t/g' | sed 's/NONE//g' > $FILE_2.ld


	plink -bfile $REFERENCE --memory $MEM --extract $FILE_2.tags --make-bed --out ${FILE_2}_tmp

	awk '{print "chr"$1"\t"$4"\t"$4"\t"$2}' ${FILE_2}_tmp.bim > ${FILE_2}.snps

	mkdir -p ../reli_input
	mv ${FILE_2}.ld ${FILE_2}.snps ../reli_input

done

echo "-------------------------------------------------------------"
echo "[INFO]	Done"
echo "-------------------------------------------------------------"

