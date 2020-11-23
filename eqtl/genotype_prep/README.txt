
# General info
MDS plots of gstcell genotypes.

source_files contains the MDS results. Its a bit messy, sorry.
Probably it is best to re-check this PCA/MDS whenever doing something
more serious with this data. 

MDS_plot_updated_genotypes > based on newly genotypes samples
MDS_with_evind_samples     > based on the data from Evind's genotypes
genotype_pca               > based on some other version of these genotypes?


# MDS calculation
MDS plots were made using MDS components calculated using plink on the pruned GStcell genotypes,
excluding the HLA region (chr6 25mb-26mb) 1000G EUR data was used and matched to provide reference.
Pruning probably done on 0.2 r2 as the LD threshold.
