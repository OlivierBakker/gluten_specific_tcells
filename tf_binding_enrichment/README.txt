# Readme
Created by:		umcg-obakker (Olivier Bakker)
Created on:		wo apr 29 10:59:04 CEST 2020
Contact at:		work.olivierbakker@gmail.com
Last updated  on:	
Principal investigator:		
Metadata:		
Source data:		

# Folders
/groups/umcg-wijmenga/tmp04/umcg-obbakker/projects/pr_gs_tcell/tf_binding_enrichment
├── data		Input data
├── output		Output from RELI analysis
├── r			R scripts and code for plotting
├── parameter_files	Files containing info on the celltypes to use for RELI analysis
└── scripts		Bash scripts to run RELI analysis


# General info
This folder contains the analysis of TF binding (based on chip seq) in CeD disease loci. 
Enrichment was performed in several senario's using Regulatory Trait Locus Intersection (RELI).
https://github.com/WeirauchLab/RELI, paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6022759/

RELI intersects a chip seq signal with genetic markers (and their LD proxies) and compares
that to a MAF matched genetic background signal based on permutations. If a observed #
of intersections is higher then random (zscore) it can be considered significantly enriched.

Here RELI is run with the Ponce meta-analysis for CeD. First top snps are determined,
by clumping in 1mb window using 0 ld threshold. Then LD proxies (r2 >0.8) are determined
and given as an input to RELI. This analysis is called "ponce_top"

A second analysis strategy involves first intersecting the genetic signal with
DHS data for gluten specific tcells, or core marks from the epigenome roadmap
project, to see if these likely enhancer regions have different enrichments.
these are either "dhs_..." or "E..." in the output folder.

Final RELI output can be found in the files *_reli_output_merged.tsv
