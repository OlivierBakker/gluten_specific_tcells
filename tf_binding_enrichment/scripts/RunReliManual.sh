#!/bin/bash
#SBATCH --mem=10G
#SBATCH --time=05:49:00
#SBATCH --nodes=1
#SBATCH --output=reli_logs/%j


ml GSL
ml GCC

RELI_DATA_DIR=/groups/umcg-wijmenga/tmp04/umcg-obbakker/tools/RELI/

for cell in $(seq 1 389);
do

celltype="$(awk '{print $1}' ../parameter_files/immune_cells.txt | tail -n $cell | head -n 1)"
echo "[INFO] Processing" $celltype

if [ -f "$3/$celltype.RELI.stats" ];
then
        echo "[INFO] $celltype already present, skipping"
        continue;
fi

/groups/umcg-wijmenga/tmp04/umcg-obbakker/tools/RELI/RELI \
-snp $1 \
-ld $2 \
-index $RELI_DATA_DIR/data/ChIPseq.index \
-data $RELI_DATA_DIR/data/ChIP-seq \
-target $celltype \
-build $RELI_DATA_DIR/data/GenomeBuild/hg19.txt \
-null $RELI_DATA_DIR/data/Null/CommonSNP_MAFmatch \
-dbsnp $RELI_DATA_DIR/data/SNPtable/SNPtable \
-out $3 \
-match \
-rep 2000 \
-corr 389 \
-phenotype celiac_disease \
-ancestry EU 

done
