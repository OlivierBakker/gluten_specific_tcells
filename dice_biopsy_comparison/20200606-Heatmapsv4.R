library(plyr)
library(dplyr)
library(tidyr)
library(DOSE)
library(AnnotationDbi)
library(GO.db)
library(org.Hs.eg.db)
library(GSEABase)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(biomaRt)
library(stringi)
library(ComplexHeatmap)
library(gtools)
library(circlize)
library(scico)

#setting working directory
setwd("/Users/R_projects/gsTcells")

#reading data from gluten-specific T cells
df.de.gstcells <- read.csv(file = "/Users/R_projects/gsTcells/data/de_analysis/output/gs_tcells_deseq2_batch_and_gender_corrected_leiden_fdr_005_replication_olso.tsv", sep = "\t", header = T, stringsAsFactors = F)

#making sure is a data frame
df.de.gstcells <- as.data.frame(df.de.gstcells)

#formatting columng to ENSEMBL for further work
names(df.de.gstcells)[4] <- "ENSEMBL"

ensembl <- read.table(file = "/Users/R_projects/gsTcells/data/OliverRNA_DEA/ensembl_gene_position_b37_v75.txt", sep = "\t", header = T, stringsAsFactors = F)

ensembl <- ensembl[,c(1,3,5)]
names(ensembl) <- c("ENSEMBL", "ENTREZID", "SYMBOL")
ensembl <- ensembl %>% distinct(ENSEMBL, .keep_all = TRUE)

#including annotation of genes
df.de.gstcells <- dplyr::left_join(df.de.gstcells, ensembl, by = "ENSEMBL")

#making a short df with ENSEMBL L2FC and condition
df.de.gstcells.short <- df.de.gstcells[c("ENSEMBL", "log2FoldChange", "condition")]

#melting df toget a matrix shape
df.de.gstcells.short <- tidyr::spread(df.de.gstcells.short, condition, log2FoldChange)
df.de.gstcells.short[duplicated(df.de.gstcells.short$ENSEMBL),]
rownames(df.de.gstcells.short) <- df.de.gstcells.short$ENSEMBL

#df.de.gstcells.short$ENSEMBL <- NULL

t10 <- read.csv(file = "/Users/R_projects/gsTcells/data/OlivierRNA2020/leiden_t10_raw.csv", header = T, stringsAsFactors = F)

t30 <- read.csv(file = "/Users/R_projects/gsTcells/data/OlivierRNA2020/leiden_t30_raw.csv", header = T, stringsAsFactors = F)

t180 <- read.csv(file = "/Users/R_projects/gsTcells/data/OlivierRNA2020/leiden_t180_raw.csv", header = T, stringsAsFactors = F)

t10 <- t10[c("gene", "log2FoldChange")]
t30 <- t30[c("gene", "log2FoldChange")]
t180 <- t180[c("gene", "log2FoldChange")]

names(t10) <- c("ENSEMBL", "t10")
names(t30) <- c("ENSEMBL", "t30")
names(t180) <- c("ENSEMBL", "t180")

df.de.gstcells.short.2 <- dplyr::left_join(df.de.gstcells.short, t10, by = "ENSEMBL")
df.de.gstcells.short.2 <- dplyr::left_join(df.de.gstcells.short.2, t30, by = "ENSEMBL")
df.de.gstcells.short.2 <- dplyr::left_join(df.de.gstcells.short.2, t180, by = "ENSEMBL")

rownames(df.de.gstcells.short.2) <- df.de.gstcells.short.2$ENSEMBL
df.de.gstcells.short.2$t0.t10 <- NULL
df.de.gstcells.short.2$t10.t30 <- NULL
df.de.gstcells.short.2$t30.t180 <- NULL
df.de.gstcells.short.2$ENSEMBL <- NULL
######Making Heatmaps######

#Creating matrix of gene expression using averages per cell type and category
#Creating matrix
all.deg.L2FC.m <- as.matrix(df.de.gstcells.short.2)
####all DEG, clustering
mat <- all.deg.L2FC.m
#mat[is.na(mat)] <- 0
mat.gene.order <- rownames(mat)

#getting annotation data from 6 clusters
clusters10 <- read.csv(file = "/Users/R_projects/gsTcells/data/cluster_memberships_kmeans_nclust6.tsv", header = F, stringsAsFactors = F, sep = "\t")

#renaming colnames to match with environment
names(clusters10) <- c("ENSEMBL", "Cluster", "SYMBOL", "Function")

#making sure Cluster column is a string vector to proper use mtach() and paste()
clusters10$Cluster <- as.character(clusters10$Cluster)

#reordering dataframe of annotation to match with current matrix
clusters10 <- clusters10[match(mat.gene.order, clusters10$ENSEMBL), ]

#reassigning cluster number based on discussion with Iris and Olivier
clusters10$Cluster <- paste(clusters10$Cluster, "A", sep = "")
clusters10$Cluster <- gsub("4A","1", clusters10$Cluster)
clusters10$Cluster <- gsub("6A","2", clusters10$Cluster)
clusters10$Cluster <- gsub("5A","3", clusters10$Cluster)
clusters10$Cluster <- gsub("2A","4", clusters10$Cluster)
clusters10$Cluster <- gsub("1A","5", clusters10$Cluster)
clusters10$Cluster <- gsub("3A","6", clusters10$Cluster)

# non-coding = includes 3prime_overlapping_ncrna, antisense, lincRNA, miRNA, processed_transcript, sense_overlapping, sense_intronic, misc_RNA
# coding = includes protein_coding, pseudogene, IG_C_pseudogene, IG_V_gene
# r/sn/snoRNAs = blank includes snRNA, rRNA, snoRNA (the last three are clearly defined ncRNAs with specific function, thus not of interest to this class)

clusters10$Function <- gsub("protein_coding", "coding RNA", clusters10$Function)
clusters10$Function <- gsub("pseudogene", "coding RNA", clusters10$Function)
clusters10$Function <- gsub("IG_C_pseudogene", "coding RNA", clusters10$Function)
clusters10$Function <- gsub("IG_V_gene", "coding RNA", clusters10$Function)
clusters10$Function <- gsub("IG_C_coding RNA", "coding RNA", clusters10$Function)
clusters10$Function <- gsub("polymorphic_coding RNA", "coding RNA", clusters10$Function)

clusters10$Function <- gsub("3prime_overlapping_ncrna", "non-coding RNA", clusters10$Function)
clusters10$Function <- gsub("antisense", "non-coding RNA", clusters10$Function)
clusters10$Function <- gsub("lincRNA", "non-coding RNA", clusters10$Function)
clusters10$Function <- gsub("miRNA", "non-coding RNA", clusters10$Function)
clusters10$Function <- gsub("processed_transcript", "non-coding RNA", clusters10$Function)
clusters10$Function <- gsub("sense_overlapping", "non-coding RNA", clusters10$Function)
clusters10$Function <- gsub("sense_intronic", "non-coding RNA", clusters10$Function)
clusters10$Function <- gsub("misc_RNA", "non-coding RNA", clusters10$Function)

clusters10$Function <- gsub("snoRNA", "non-coding RNA", clusters10$Function)
clusters10$Function <- gsub("snRNA", "non-coding RNA", clusters10$Function)
clusters10$Function <- gsub("rRNA", "non-coding RNA", clusters10$Function)


#Set annotation
#ann <- data.frame(metadata$Type, metadata$Type2)

type = gsub("s\\d+_", "", colnames(mat))
ha = HeatmapAnnotation(type = type, annotation_name_side = "left")


#colnames(ann) <- c("Type", "Type2")
colours <- list("Type"=c("M"="red2","N"="royalblue"), "Type2"=c("AT"="limegreen","PT"="gold"))
#colAnn <- HeatmapAnnotation(df=ann, which="col", col=colours, annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))

Dark2colors <- brewer.pal(n = 6, name = "Dark2")
#plot heatmap

#3 colors palette for Gene function
#tricolors <- brewer.pal(n = 3, name = "Dark2")
tricolors <- c("#EC4F4B", "#2F9486")
#Custom palette for clusters
custom.colors <- c("#3BB273", "#6B6174", "#EF6461", "#FABC2A", "#8576B6", "#2274A5")

#Heatmap for all deg splitted by cluster
m <- Heatmap(mat, 
        name = "L2FC", 
        col = colorRamp2(c(6, 2, 0, -2, -6), brewer.pal(n = 5, name = "RdBu")),
        na_col = "white",
        left_annotation =  rowAnnotation(Function = clusters10$Function, 
                                         col = list(Function = c("coding RNA" = tricolors[1], 
                                                                 "non-coding RNA" = tricolors[2])),
                                         annotation_name_gp = gpar(fontsize = 11, fontfamily = "ArialMT")),
        right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 2:7, 
                                                                     col = NA, 
                                                                    fill = custom.colors,
                                                                    fontsize = 11, 
                                                                    fontfamily = "ArialMT"))),
        column_title_rot = 0,
        row_title_rot = 0,
        row_title_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
        row_names_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
        column_names_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
        show_row_names = FALSE,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        row_split = clusters10$Cluster,
        cluster_row_slices = FALSE,
        row_gap = unit(2, "mm"),
        use_raster = F
)

pdf(file="Results/heatmap_perclusterv3.pdf", width=4, height=6,  family = "ArialMT", paper="a4")
m
dev.off()

#Heatmap for all deg hierarchical clustered 
m2 <- Heatmap(mat, 
             name = "L2FC", 
             col = colorRamp2(c(6, 2, 0, -2, -6), brewer.pal(n = 5, name = "RdBu")),
             na_col = "white",
             left_annotation =  rowAnnotation(Function = clusters10$Function, 
                                              col = list(Function = c("coding RNA" = tricolors[1], 
                                                                      "non-coding RNA" = tricolors[2])),
                                              annotation_name_gp = gpar(fontsize = 11, fontfamily = "ArialMT")),
             right_annotation = rowAnnotation(Cluster = clusters10$Cluster,
                                              col = list(Cluster = c("1" = custom.colors[1], 
                                                                     "2" = custom.colors[2],
                                                                     "3" = custom.colors[3],
                                                                     "4" = custom.colors[4],
                                                                     "5" = custom.colors[5],
                                                                     "6" = custom.colors[6]) ), 
                                              annotation_name_gp = gpar(fontsize = 11, fontfamily = "ArialMT")),
             column_title_rot = 0,
             row_title_rot = 0,
             row_title_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
             row_names_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
             column_names_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
             show_row_names = FALSE,
             cluster_columns = FALSE,
             cluster_rows = TRUE,
             #row_split = clusters10$Cluster,
             cluster_row_slices = FALSE,
             row_gap = unit(2, "mm"),
             use_raster = F
)

pdf(file="Results/heatmap_hierarchicalclusterv2.pdf", width=4, height=6,  family = "ArialMT", paper="a4")
m2
dev.off()

#Heatmaps for eQTLs

###eQTL genes
#getting reference for eQTL genes with evidence
eQTL.genes.ref.evidence <- read.table(file = "/Users/R_projects/gsTcells/data/2020Apr-eQTLgenes-128prior-Adriaan.csv", sep = ",", header = TRUE)
#getting all associated eQTL
eQTL.genes.ref.all <- read.table(file = "/Users/R_projects/gsTcells/data/2020Apr-eQTLgenes-all-Adriaan.csv", sep = ",", header = TRUE)

#custom L2FC breaks
l2fc.breaks <- c(seq(from=-5 , to=-0.5, length.out = 5),-1, -0.50 ,0, 0.5 ,1, seq(from=1.1 , to=5, length.out = 5))
#Custom L2FC colors
l2fc.colors <- c(rev(brewer.pal(name = "Blues", n=9)[5:9]), "#525252","#252525", "#252525", "#525252", brewer.pal(name = "Reds", n=9)[5:9])

#intersection of all deg with eQTLs
intersect.alldeg_eQTLref <- intersect(eQTL.genes.ref.all$ENSEMBL, rownames(all.deg.L2FC.m))

#intersection of all prioritized eQTL
deg_eQTLref_evidence <- intersect(eQTL.genes.ref.evidence$ENSEMBL, rownames(all.deg.L2FC.m))

#making matrix
matrix.deg.eQTLref <- all.deg.L2FC.m[intersect.alldeg_eQTLref, ]

#making metadata
deg.eQTLref <- eQTL.genes.ref.all[eQTL.genes.ref.all$ENSEMBL %in% intersect.alldeg_eQTLref, ]
deg.eQTLref <- left_join(deg.eQTLref, clusters10, by = "ENSEMBL")

#Heatmap with all eqTLgenes
Heatmap(matrix.deg.eQTLref, 
        name = "expression", 
        col = colorRamp2(c(6, 2, 0, -2, -6), brewer.pal(n = 5, name = "RdBu")),
        na_col = "white", 
        rect_gp = gpar(col = "gray", lwd = 0.5),
        column_title_rot = 0,
        row_title_rot = 0,
        row_labels = deg.eQTLref$SYMBOL.x,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        right_annotation =  rowAnnotation(Cluster = deg.eQTLref$Cluster, 
                                          col = list(Cluster = c("1" = custom.colors[1], 
                                                                 "2" = custom.colors[2], 
                                                                 "3" = custom.colors[3], 
                                                                 "4" = custom.colors[4], 
                                                                 "5" = custom.colors[5], 
                                                                 "6" = custom.colors[6]))),
        row_title_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
        column_title_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
        column_names_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
        row_names_gp = gpar(fontsize = 8, fontfamily = "ArialMT",
                            col = ifelse(rownames(matrix.deg.eQTLref) %in% deg_eQTLref_evidence, 
                                         "black", NA)),
        row_split = factor(deg.eQTLref$region, 
                           levels = mixedsort(as.character(unique(deg.eQTLref$region)))),
        column_split = gsub(
          gsub(colnames(matrix.deg.eQTLref), 
               pattern = '_.*', 
               replacement = ''), 
          pattern = ' ', 
          replacement = '_'),
        row_gap = unit(1, "mm"),
        column_gap = unit(2, "mm")
)


####Prioritized genes
#making matrix, noticed that is overwritting terms
matrix.deg.eQTLref <- all.deg.L2FC.m[deg_eQTLref_evidence, ]

#making metadata, overwritting previous one
deg.eQTLref <- eQTL.genes.ref.evidence[eQTL.genes.ref.evidence$ENSEMBL %in% deg_eQTLref_evidence, ]
deg.eQTLref <- left_join(deg.eQTLref, clusters10, by = "ENSEMBL")

#heatmap with prioritized genes
m <- Heatmap(matrix.deg.eQTLref, 
        name = "L2FC", 
        col = colorRamp2(c(6, 2, 0, -2, -6), brewer.pal(n = 5, name = "RdBu")),
        na_col = "white", 
        rect_gp = gpar(col = "gray", lwd = 0.5),
        column_title_rot = 0,
        row_title_rot = 0,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        right_annotation =  rowAnnotation(Cluster = deg.eQTLref$Cluster, 
                                          col = list(Cluster = c("1" = custom.colors[1], 
                                                                 "2" = custom.colors[2], 
                                                                 "3" = custom.colors[3], 
                                                                 "4" = custom.colors[4], 
                                                                 "5" = custom.colors[5], 
                                                                 "6" = custom.colors[6]))),
        row_labels = deg.eQTLref$SYMBOL.x,
        row_title_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
        column_title_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
        row_names_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
        column_names_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
        row_split = factor(deg.eQTLref$region, 
                           levels = mixedsort(as.character(unique(deg.eQTLref$region)))),
        column_split = gsub(
          gsub(colnames(matrix.deg.eQTLref), 
               pattern = '_.*', 
               replacement = ''), 
          pattern = ' ', 
          replacement = '_'),
        row_gap = unit(2, "mm"),
        column_gap = unit(2, "mm")
)

pdf(file="Results/heatmap_prioritizedgenes.pdf", width=6, height=6,  family = "ArialMT", paper="a4")
m
dev.off()
