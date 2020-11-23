#!/bin/bash



ponce=/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/wd_gs_tcell/tf_binding_enrichment/data/ponce_summary_stats_b37.bed
roadmap=/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/wd_gs_tcell/tf_binding_enrichment/data/roadmap/coreMarks
output=/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/wd_gs_tcell/tf_binding_enrichment/output

while read cell;
do

echo "[INFO] Submitting $cell"

mkdir -p $output/$cell/reli_output

sbatch RunReli.sh \
$output/$cell/reli_input/ponce_summary_stats_b37_clumped_5e-8.snps \
$output/$cell/reli_input/ponce_summary_stats_b37_clumped_5e-8.ld \
$output/$cell/reli_output/

done < ../current_batch.txt
