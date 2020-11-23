#!/bin/bash
#SBATCH --time=05:59:00
#SBATCH --mem=16G


set -e

ml RPlus
ml SAMtools
ml Python/3.6.3-foss-2015b
source macs2_venv/bin/activate


INBAM=$1
OUTDIR=./output
FILE=$(basename -- "$INBAM")
BASE="${FILE%.*}"


# Filter non-primary and duplicated alignments as well as lower quality's then 30
samtools view -h -F 1280 -q 30 -b ${INBAM} > ${OUTDIR}/${BASE}.filtered.bam

CMD="macs2 callpeak \
--treatment ${OUTDIR}/${BASE}.filtered.bam \
--name ${BASE}.filtered \
--outdir ${OUTDIR} \
--format BAM \
--gsize hs \
--broad \
--nomodel \
--shift -125 \
--extsize 250 \
--bdg"

echo $CMD
eval $CMD

gzip output/${BASE}.filtered_control_lambda.bdg
gzip output/${BASE}.filtered_treat_pileup.bdg

Rscript output/${BASE}.sorted_model.r
