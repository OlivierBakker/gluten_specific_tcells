# Analysis on gluten specific T-cell clones

This repository contains all code and scripts relating to the manuscript on gluten specific T cell clones (Bakker & Ramirez-Sanchez 2021)
The paper: https://www.nature.com/articles/s41598-021-86612-5

## Abstract
Celiac disease is an auto-immune disease in which an immune response to dietary gluten leads to inflammation and subsequent atrophy of small intestinal villi, causing severe bowel discomfort and malabsorption of nutrients. The major instigating factor for the immune response in celiac disease is the activation of gluten-specific CD4+ T cells expressing T cell receptors that recognize gluten peptides presented in the context of HLA-DQ2 and DQ8. Here we provide an in-depth characterization of 28 gluten-specific T cell clones. We assess their transcriptional and epigenetic response to T cell receptor stimulation and link this to genetic factors associated with celiac disease. Gluten-specific T cells have a distinct transcriptional profile that mostly resembles that of Th1 cells but also express cytokines characteristic of other types of T-helper cells. This transcriptional response appears not to be regulated by changes in chromatin state, but rather by early upregulation of transcription factors and non-coding RNAs that likely orchestrate the subsequent activation of genes that play a role in immune pathways. Finally, integration of chromatin and transcription factor binding profiles suggest that genes activated by T cell receptor stimulation of gluten‑specific T cells may be impacted by genetic variation at several genetic loci associated with celiac disease.


# Repo structure

The code in this repo is divided per analysis into a folder, appart from plotting_function.r which is sourced by several downstream R scripts. Each analysis is mostly done in R and bash. This repo is not intended as a pipeline, but represents bespoke code for the analysis of the data. Most scripts should be supplied with inline documentation, if not I appologise. Non R analyses should have a readme detailling the scripts.
