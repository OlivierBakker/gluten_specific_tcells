#!/bin/bash

# INPUT
MEM=20000
DATA_DIR=/groups/umcg-wijmenga/tmp04/umcg-obbakker/data
REFERENCE=$DATA_DIR/reference/1000G/1KG_EUR_1000G_all

# OUTPUT
CLUMP_DIR=./

module load plink
module load R
set -e


plink -bfile $REFERENCE --out ced_tagging --show-tags $1 --tag-r2 0.8 --tag-kb 500 --list-all
awk '{print $1":\t"$8}' ced_tagging.tags.list | tail -n +2 | sed 's/|/\t/g' | sed 's/NONE//g' > ced_tagging.ld

