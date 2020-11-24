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
library(ggpubr)
library(grid)
library(ggplotify)
library(UpSetR)

## ------------------------------------------------------------------------
# Simple plotting theme for ggplot using arial family font
theme.plain <- function(p, base_size = 11, base_family = "ArialMT") {
  p <- p + theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black", size=0.75),
          axis.ticks = element_line(size=0.75),
          axis.text = element_text(size=base_size, family="ArialMT", face="plain"),
          strip.background = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size=base_size, family = "ArialMT", face="plain"),
          complete = TRUE,
          plot.title = element_text(hjust=0.5))
  return(p)
}

theme.plain.2 <- function(p, base_size, base_family = "ArialMT") {
  p <- p + theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black", size=0.75),
          axis.ticks = element_line(size=0.75),
          axis.text = element_text(size=base_size, family="ArialMT", face="plain"),
          strip.background = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size=base_size, family = "ArialMT", face="plain"),
          complete = TRUE,
          plot.title = element_text(hjust=0.5))
  return(p)
}
#setting working directory
setwd("//Users/R_projects/gsTcells")

#reading data from gluten-specific T cells
de.leiden <- read.csv(file = "/Users/R_projects/gsTcells/data/OliverRNA_DEA/leiden_t0t180_raw.csv", header = T, stringsAsFactors = F)

#making sure is a data frame
de.leiden <- as.data.frame(de.leiden)

#formatting columng to ENSEMBL for further work
names(de.leiden)[1] <- "ENSEMBL"

#Creating an extra column with the overlpa category of the genes
de.leiden$Direction <- paste(de.leiden$log2FoldChange > 0)
de.leiden$Direction <- gsub("TRUE", "Upregulated", de.leiden$Direction)
de.leiden$Direction <- gsub("FALSE", "Downregulated", de.leiden$Direction)

de.leiden <- de.leiden[order(de.leiden[,1], de.leiden[,7], -de.leiden[,3]), ]
de.leiden <- de.leiden %>% distinct(ENSEMBL, .keep_all = TRUE)
de.leiden <- de.leiden[-1,]
de.leiden <- de.leiden[-1,]
de.leiden <- de.leiden[-1,]

de.leiden.unfiltered <- de.leiden
de.leiden <- de.leiden.unfiltered
de.leiden <- de.leiden[!is.na(de.leiden$padj),]
de.leiden <- de.leiden[de.leiden$padj < 0.05,]
de.leiden <- de.leiden[!is.na(de.leiden$log2FoldChange),]
de.leiden <- de.leiden[abs(de.leiden$log2FoldChange) > 1,]
de.leiden <- de.leiden[order(de.leiden[,1], -de.leiden[,3]), ]
de.leiden <- de.leiden %>% distinct(ENSEMBL, .keep_all = TRUE)

####Upset plots

reference.biopsy <- read.table(file = "/Users/R_projects/gsTcells/Biopsy_data/DEseq_DonatellaBiopsies_21072017.txt", sep = "\t", header = T, stringsAsFactors = F)
reference.Tcell.stim <- read.csv(file = "/Users/R_projects/gsTcells/Pattern_genes/CD4_STIM.Vs.CD4_NAIVE.csv", header = T, stringsAsFactors = F)
reference.solid.healthy <- read.csv(file = "/Users/R_projects/gsTcells/data/DEG_solid_healthy_comparison.csv", header = T, stringsAsFactors = F)
reference.solid.tet <- read.csv(file = "/Users/R_projects/gsTcells/data/DEG_solid_tet_comparison.csv", header = T, stringsAsFactors = F)


reference.solid.healthy$comparison <- "healthy"
reference.solid.tet$comparison <- "tet"
reference.solid <- rbind(reference.solid.tet, reference.solid.healthy)
reference.solid.2 <- full_join(reference.solid.healthy, reference.solid.tet, by = "ENSEMBL")

reference.solid <- reference.solid %>% distinct(ENSEMBL, .keep_all = TRUE)

#plotting concordant solid data
annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("Q3", "Q4", "Q2","Q1"),
  hjustvar = c(-1, -1, 2, 2) ,
  vjustvar = c(-1, 2, -1, 2)) #<- adjust

p.solid.only <- ggplot(data = reference.solid.2, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
  geom_point(alpha = 0.5, size = 2) + 
  #geom_point(data = subset(de.gstcell.dice, Overlap == "DICE"), aes(colour = Overlap), alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #geom_text(data = subset(de.gstcell.dice, Gene == "IL21"), aes(label = Gene),hjust=0,vjust=0) +
  ylab("L2FC in tet+ vs tet-") +
  xlab("L2FC in tet+ vs healthy") +
  geom_text(data=annotations,aes(x = xpos, y = ypos,hjust=hjustvar,vjust=vjustvar,label= annotateText), size = 4) +
  xlim(-8, 18) +
  ylim(-8, 18) +
  coord_fixed() +
  #scale_color_manual(values=c("#8CED3C", "#4B4BCD", "#30BCBC")) +
  ggtitle("")

p.solid.only <- theme.plain(p.solid.only)

#formatting columng to ENSEMBL for further work
names(reference.biopsy)[1] <- "ENSEMBL"
names(reference.Tcell.stim)[7] <- "ENSEMBL"

names(reference.Tcell.stim)[6] <- "padj"
names(reference.Tcell.stim)[5] <- "log2FoldChange"

reference.biopsy <- subset(reference.biopsy, 
                           abs(reference.biopsy$log2FoldChange) > 1) 
reference.biopsy <- subset(reference.biopsy, 
                           reference.biopsy$padj < 0.05)
reference.biopsy <- reference.biopsy %>% distinct(ENSEMBL, .keep_all = TRUE)

reference.Tcell.stim <- subset(reference.Tcell.stim, 
                               abs(reference.Tcell.stim$log2FoldChange) > 1) 
reference.Tcell.stim <- subset(reference.Tcell.stim, 
                               reference.Tcell.stim$padj < 0.05)
reference.Tcell.stim <- reference.Tcell.stim %>% distinct(ENSEMBL, .keep_all = TRUE)

#UPSET
#Creating function to create matrix with the DEG 
list_of_genes <- list_to_matrix(list(gsTcells = as.character(de.leiden$ENSEMBL),
                                     Biopsy = as.character(reference.biopsy$ENSEMBL),
                                     DICE = as.character(reference.Tcell.stim$ENSEMBL),
                                     biopsy_derived_gsTcells = as.character(reference.solid$ENSEMBL)))

upsetcolors <- c("#8CED3C", "#4B4BCD", "#30BCBC")
#Creating matrix of combinations for T cells
#upsetcolors <- c("#8CED3C", "#4B4BCD", "#30BCBC")
m.upset <- make_comb_mat(list_of_genes)

#if using L2FC >2
#upset.p <- UpSet(m.upset, 
#                 row_title_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
#                 column_title_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
#                 column_names_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
#                 row_names_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
#                 comb_col = c(upsetcolors[2],
#                              upsetcolors[1],
#                              upsetcolors[1],
#                              upsetcolors[3],
#                              upsetcolors[1],
#                              upsetcolors[3],
#                              upsetcolors[3],
#                              upsetcolors[1],
#                              upsetcolors[3],
#                              upsetcolors[1],
#                              upsetcolors[3],
#                              upsetcolors[1])
#)

#if using L2FC >1

ss <- set_size(m.upset)
cs <- comb_size(m.upset)

upset.p <- UpSet(m.upset, 
                 set_order = order(ss),
                 comb_order = order(comb_degree(m.upset), -cs),
                 row_title_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
                 column_title_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
                 column_names_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
                 row_names_gp = gpar(fontsize = 11, fontfamily = "ArialMT"),
                 comb_col = c(upsetcolors[2],
                              upsetcolors[1],
                              upsetcolors[1],
                              upsetcolors[3],
                              upsetcolors[3],
                              upsetcolors[1],
                              upsetcolors[1],
                              upsetcolors[3],
                              upsetcolors[1],
                              upsetcolors[3],
                              upsetcolors[1],
                              upsetcolors[3],
                              upsetcolors[3],
                              upsetcolors[1],
                              upsetcolors[3])
)

colnames(list_of_genes)
test.combination <- make_comb_mat(list_of_genes)
comb_size(test.combination)

upset.p
ht <- draw(upset.p)
od <- column_order(ht)
cs <- comb_size(test.combination)
decorate_annotation("Intersection\nsize", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))})


shared.all <- extract_comb(test.combination, "1111")
shared.tet <- extract_comb(test.combination, "1011")
shared.ctrl <- extract_comb(test.combination, "1101")
shared.gst <- extract_comb(test.combination, "1001")
specific <- extract_comb(test.combination, "1000")
specific.2 <- extract_comb(test.combination, "1100")


de.gstcells <- read.csv(file = "/Users/R_projects/gsTcells/Results/de.gstcell.dice.v3.csv", header = T, stringsAsFactors = F)
shared.de.genes <- de.gstcells[de.gstcells$ENSEMBL %in% c(shared.all, 
                                                          shared.tet, 
                                                          shared.ctrl, 
                                                          shared.gst),]



specific.de.genes <- de.gstcells[de.gstcells$ENSEMBL %in% c(specific,
                                                            specific.2
                                                            ),]
specific.specific.de.genes <- specific.de.genes[specific.de.genes$Overlap == "gsTcells",]

ensembl <- read.table(file = "/Users/R_projects/gsTcells/data/OliverRNA_DEA/ensembl_gene_position_b37_v75.txt", sep = "\t", header = T, stringsAsFactors = F)

ensembl <- ensembl[,c(1,3,5)]
names(ensembl) <- c("ENSEMBL", "ENTREZID", "SYMBOL")
shared.de.genes <- dplyr::left_join(shared.de.genes, ensembl, by = "ENSEMBL")
specific.specific.de.genes <- dplyr::left_join(specific.specific.de.genes, ensembl, by = "ENSEMBL")

shared.de.genes <- shared.de.genes %>% distinct(ENSEMBL, .keep_all = TRUE)
specific.specific.de.genes <- specific.specific.de.genes %>% distinct(ENSEMBL, .keep_all = TRUE)

shared.de.genes$Overlap.external.gstcell <- "Shared"
specific.specific.de.genes$Overlap.external.gstcell <- "Unique"

de.overlap.gstcells <- rbind(shared.de.genes, specific.specific.de.genes)
de.overlap.gstcells <- de.overlap.gstcells %>% distinct(ENSEMBL, .keep_all = TRUE)
#write.csv(de.overlap.gstcells, file = "./Results/de.overlap.gstcells.solid.csv")

#New format for oing enrichment, sorting genes based on L2FC
formatting.gene.list.second <- function(gene.list.2columns){
  df <- gene.list.2columns
  geneList = df$log2FoldChange
  names(geneList) = as.character(df$ENTREZID)
  geneList = sort(geneList, decreasing = TRUE)
  return(geneList)
}
#New format to only include gene names already sorted
gene.list.only.names <- function(gene.list.vector){
  df <- gene.list.vector
  return(names(df))
}


### Specific genes reactome enrichment gstcells compared with DICE
### gstcells
de.overlap.gstcells.2 <- de.overlap.gstcells[,c("ENSEMBL", "ENTREZID.x", "SYMBOL.x", "log2FoldChange.x", "padj.x", "Overlap", "Direction", "Overlap.external.gstcell")]
names(de.overlap.gstcells.2) <- c("ENSEMBL", "ENTREZID", "SYMBOL", "log2FoldChange", "padj", "Overlap", "Direction", "Overlap.external.gstcell")
de.overlap.gstcells.2 <- de.overlap.gstcells.2[abs(de.overlap.gstcells.2$log2FoldChange) > 2 & de.overlap.gstcells.2$padj < 0.05,]

de.overlap.gstcells.2 <- de.overlap.gstcells.2 %>% distinct(ENSEMBL, .keep_all = TRUE)
de.overlap.gstcells.2 <- de.overlap.gstcells.2[de.overlap.gstcells.2$Overlap.external.gstcell == "Shared",]

short.de.overlap.gstcells <- de.overlap.gstcells.2[,c("ENTREZID", "log2FoldChange", "Direction")]

list.de.overlap.gstcells <- list(Upregulated = subset(short.de.overlap.gstcells, Direction == "Upregulated"), 
                                 Downregulated = subset(short.de.overlap.gstcells, Direction == "Downregulated"))
list.de.overlap.gstcells <- lapply(list.de.overlap.gstcells, function(Q){
  Q$Direction <- NULL
  return(Q)
})

#list with vector names as genes and L2FC
list.de.overlap.gstcells.reactomeformatted <- lapply(list.de.overlap.gstcells , formatting.gene.list.second)
#list with only gene names
list.de.overlap.gstcells.reactomeformatted.2 <- lapply(list.de.overlap.gstcells.reactomeformatted, gene.list.only.names)

res.df.gst <- compareCluster(list.de.overlap.gstcells.reactomeformatted.2, fun="enrichPathway", readable=T)

res.df.gst@compareClusterResult$p.adjust <- -log10(res.df.gst@compareClusterResult$p.adjust)


p.react.gst.shared <- theme.plain(dotplot(res.df.gst, font.size = 8, showCategory = 10), base_size = 8) +
  scale_color_gradient(low="lightgrey", high = "red", name="-log10(p-adjust)") 
#write.csv(res.df.gst, file = "./Results/reactome.de.gstcell.solid.shared.csv")



### Specific genes reactome enrichment specific gstcells 
### gstcells
de.overlap.gstcells.2 <- de.overlap.gstcells[,c("ENSEMBL", "ENTREZID.x", "SYMBOL.x", "log2FoldChange.x", "padj.x", "Overlap", "Direction", "Overlap.external.gstcell")]
names(de.overlap.gstcells.2) <- c("ENSEMBL", "ENTREZID", "SYMBOL", "log2FoldChange", "padj", "Overlap", "Direction", "Overlap.external.gstcell")
de.overlap.gstcells.2 <- de.overlap.gstcells.2[abs(de.overlap.gstcells.2$log2FoldChange) > 2 & de.overlap.gstcells.2$padj < 0.05,]

de.overlap.gstcells.2 <- de.overlap.gstcells.2 %>% distinct(ENSEMBL, .keep_all = TRUE)
de.overlap.gstcells.2 <- de.overlap.gstcells.2[de.overlap.gstcells.2$Overlap.external.gstcell == "Unique",]

short.de.overlap.gstcells <- de.overlap.gstcells.2[,c("ENTREZID", "log2FoldChange", "Direction")]

list.de.overlap.gstcells <- list(Upregulated = subset(short.de.overlap.gstcells, Direction == "Upregulated"), 
                                 Downregulated = subset(short.de.overlap.gstcells, Direction == "Downregulated"))
list.de.overlap.gstcells <- lapply(list.de.overlap.gstcells, function(Q){
  Q$Direction <- NULL
  return(Q)
})

#list with vector names as genes and L2FC
list.de.overlap.gstcells.reactomeformatted <- lapply(list.de.overlap.gstcells , formatting.gene.list.second)
#list with only gene names
list.de.overlap.gstcells.reactomeformatted.2 <- lapply(list.de.overlap.gstcells.reactomeformatted, gene.list.only.names)

res.df.gst <- compareCluster(list.de.overlap.gstcells.reactomeformatted.2, fun="enrichPathway", readable=T)

res.df.gst@compareClusterResult$p.adjust <- -log10(res.df.gst@compareClusterResult$p.adjust)


p.react.gst.specific <- theme.plain(dotplot(res.df.gst, font.size = 8, showCategory = 10), base_size = 8) +
  scale_color_gradient(low="lightgrey", high = "red", name="-log10(p-adjust)") 
#write.csv(res.df.gst, file = "./Results/reactome.de.gstcell.solid.specific.csv")



###Making concordance plots
de.leiden <- de.leiden.unfiltered
de.leiden <- de.leiden[!is.na(de.leiden$padj),]
de.leiden <- de.leiden[de.leiden$padj < 0.05,]
de.leiden <- de.leiden[!is.na(de.leiden$log2FoldChange),]
de.leiden <- de.leiden[abs(de.leiden$log2FoldChange) > 1,]
de.leiden <- de.leiden[order(de.leiden[,1], -de.leiden[,3]), ]
de.leiden <- de.leiden %>% distinct(ENSEMBL, .keep_all = TRUE)


#Joining dataframes from different timepoints
de.gstcell.solid <- full_join(de.leiden, reference.solid, by = "ENSEMBL")

#Adding value of 0 to NAs in L2FC column
de.gstcell.solid$log2FoldChange.x[is.na(de.gstcell.solid$log2FoldChange.x)] <- 0
de.gstcell.solid$log2FoldChange.y[is.na(de.gstcell.solid$log2FoldChange.y)] <- 0

#Adding value of 1 to NAs in padj column
de.gstcell.solid$padj.x[is.na(de.gstcell.solid$padj.x)] <- 1
de.gstcell.solid$padj.y[is.na(de.gstcell.solid$padj.y)] <- 1


####filtering
de.gstcell.solid <- de.gstcell.solid[abs(de.gstcell.solid$log2FoldChange.x) > 2 | 
                                     abs(de.gstcell.solid$log2FoldChange.y) > 2,]


de.gstcell.solid <- de.gstcell.solid[de.gstcell.solid$padj.x < 0.05 | 
                                     de.gstcell.solid$padj.y < 0.05, ]

#Creating an extra column with the overlap category of the genes
de.gstcell.solid$Overlap <- paste(de.gstcell.solid$padj.x < 0.05 &
                                   abs(de.gstcell.solid$log2FoldChange.x) > 1,
                                 de.gstcell.solid$padj.y < 0.05 & 
                                   abs(de.gstcell.solid$log2FoldChange.y) > 1, sep = "")
de.gstcell.solid$Overlap <- gsub("TRUETRUE", "Both", de.gstcell.solid$Overlap)

de.gstcell.solid$Overlap <- gsub("TRUENA", "gsTcells", de.gstcell.solid$Overlap)
de.gstcell.solid$Overlap <- gsub("TRUEFALSE", "gsTcells", de.gstcell.solid$Overlap)
de.gstcell.solid$Overlap <- gsub("NATRUE", "solid", de.gstcell.solid$Overlap)
de.gstcell.solid$Overlap <- gsub("FALSETRUE", "solid", de.gstcell.solid$Overlap)

de.gstcell.solid$Overlap <- gsub("NANA", NA, de.gstcell.solid$Overlap)
de.gstcell.solid$Overlap <- gsub("NAFALSE", NA, de.gstcell.solid$Overlap)
de.gstcell.solid$Overlap <- gsub("FALSENA", NA, de.gstcell.solid$Overlap)
de.gstcell.solid$Overlap <- gsub("FALSEFALSE", NA, de.gstcell.solid$Overlap)


#Cleaning
de.gstcell.solid <- de.gstcell.solid[!is.na(de.gstcell.solid$Overlap), ]

de.gstcell.solid$Overlap <- as.factor(de.gstcell.solid$Overlap)

de.gstcell.solid$Overlap <- factor(de.gstcell.solid$Overlap, levels = c("solid", "gsTcells", "Both"))

de.gstcell.solid <- de.gstcell.solid %>% distinct(ENSEMBL, .keep_all = TRUE)

#plotting concordant GST cells and Biopsy
annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("Q3", "Q4", "Q2","Q1"),
  hjustvar = c(-1, -1, 2, 2) ,
  vjustvar = c(-1, 2, -1, 2)) #<- adjust

p.solid <- ggplot(data = de.gstcell.solid, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
  geom_point(aes(colour = Overlap), alpha = 0.5, size = 2) + 
  #geom_point(data = subset(de.gstcell.solid, Overlap == "solid"), aes(colour = Overlap), alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #geom_text(data = subset(de.gstcell.solid, Gene == "IL21"), aes(label = Gene),hjust=0,vjust=0) +
  ylab("L2FC in solid") +
  xlab("L2FC in gsTcells") +
  geom_text(data=annotations,aes(x = xpos, y = ypos,hjust=hjustvar,vjust=vjustvar,label= annotateText), size = 4) +
  xlim(-8, 18) +
  ylim(-8, 18) +
  coord_fixed() +
  scale_color_manual(values=c("#8CED3C", "#4B4BCD", "#30BCBC")) +
  ggtitle("")

p.solid <- theme.plain(p.solid) + theme(legend.position = "none")# + guides(colour=guide_legend(nrow = 2, byrow = TRUE))




####Enrichments Solid


de.gstcell.solid <- dplyr::left_join(de.gstcell.solid, ensembl, by = "ENSEMBL")

###Adding quadrant data
de.gstcell.solid$Quadrant <- NA
de.gstcell.solid[de.gstcell.solid$Overlap == "Both" & 
                  de.gstcell.solid$log2FoldChange.x > 1 & 
                  de.gstcell.solid$log2FoldChange.y > 1, ]$Quadrant <- "Q1"
de.gstcell.solid[de.gstcell.solid$Overlap == "Both" & 
                  de.gstcell.solid$log2FoldChange.x > 1 & 
                  de.gstcell.solid$log2FoldChange.y < -1, ]$Quadrant <- "Q2"
de.gstcell.solid[de.gstcell.solid$Overlap == "Both" & 
                  de.gstcell.solid$log2FoldChange.x < -1 & 
                  de.gstcell.solid$log2FoldChange.y < -1, ]$Quadrant <- "Q3"
de.gstcell.solid[de.gstcell.solid$Overlap == "Both" & 
                  de.gstcell.solid$log2FoldChange.x < -1 & 
                  de.gstcell.solid$log2FoldChange.y > 1, ]$Quadrant <- "Q4"


de.gstcell.solid$gene <- NULL
de.gstcell.solid <- de.gstcell.solid %>% distinct(ENSEMBL, .keep_all = TRUE)


#write.csv(de.gstcell.solid, file = "./Results/de.gstcell.solid.v3.csv")

de.gstcell.solid.2 <- de.gstcell.solid[,c("ENSEMBL", "ENTREZID", "SYMBOL.y", "log2FoldChange.x", "padj.x", "Quadrant")]
names(de.gstcell.solid.2) <- c("ENSEMBL", "ENTREZID", "SYMBOL", "log2FoldChange", "padj", "Quadrant")
de.gstcell.solid.2 <- de.gstcell.solid.2[abs(de.gstcell.solid.2$log2FoldChange) > 1 & de.gstcell.solid.2$padj < 0.05,]

de.gstcell.solid.2 <- de.gstcell.solid.2 %>% distinct(ENTREZID, .keep_all = TRUE)

short.de.gstcell.solid <- de.gstcell.solid.2[,c("ENTREZID", "log2FoldChange", "Quadrant")]

list.de.gstcell.solid <- list(Q1 = subset(short.de.gstcell.solid, Quadrant == "Q1"), 
                             Q2 = subset(short.de.gstcell.solid, Quadrant == "Q2"),
                             Q3 = subset(short.de.gstcell.solid, Quadrant == "Q3"),
                             Q4 = subset(short.de.gstcell.solid, Quadrant == "Q4"))
list.de.gstcell.solid <- lapply(list.de.gstcell.solid, function(Q){
  Q$Quadrant <- NULL
  return(Q)
})


#list with vector names as genes and L2FC
list.de.gstcell.solid.reactomeformatted <- lapply(list.de.gstcell.solid , formatting.gene.list.second)
#list with only gene names
list.de.gstcell.solid.reactomeformatted.2 <- lapply(list.de.gstcell.solid.reactomeformatted, gene.list.only.names)

#REACTOME enrichment
res.df.solid <- compareCluster(list.de.gstcell.solid.reactomeformatted.2, fun="enrichPathway", readable=T)

#PLOTS
# Convert adjusted p to -log10 scale
res.df.solid@compareClusterResult$p.adjust <- -log10(res.df.solid@compareClusterResult$p.adjust)

#write.csv(res.df.solid, file = "./Results/reactome.de.solid.gstcellv2.csv")

p.react.solid <- theme.plain(dotplot(res.df.solid, font.size = 8, showCategory = 10), base_size = 8) +
  scale_color_gradient(low="lightgrey", high = "red", name="-log10(p-adjust)") 



#####################

reference.biopsy <- read.table(file = "/Users/R_projects/gsTcells/Biopsy_data/DEseq_DonatellaBiopsies_21072017.txt", sep = "\t", header = T, stringsAsFactors = F)

reference.Tcell.stim <- read.csv(file = "/Users/R_projects/gsTcells/Pattern_genes/CD4_STIM.Vs.CD4_NAIVE.csv", header = T, stringsAsFactors = F)

de.leiden <- de.leiden.unfiltered 

#formatting columng to ENSEMBL for further work
names(reference.biopsy)[1] <- "ENSEMBL"
names(reference.Tcell.stim)[7] <- "ENSEMBL"

names(reference.Tcell.stim)[6] <- "padj"
names(reference.Tcell.stim)[5] <- "log2FoldChange"


#Joining dataframes from different timepoints
de.gstcell.dice <- full_join(de.leiden, reference.Tcell.stim, by = "ENSEMBL")
de.gstcell.biopsy <- full_join(de.leiden, reference.biopsy, by = "ENSEMBL")

#Adding value of 0 to NAs in L2FC column
de.gstcell.dice$log2FoldChange.x[is.na(de.gstcell.dice$log2FoldChange.x)] <- 0
de.gstcell.dice$log2FoldChange.y[is.na(de.gstcell.dice$log2FoldChange.y)] <- 0

de.gstcell.biopsy$log2FoldChange.x[is.na(de.gstcell.biopsy$log2FoldChange.x)] <- 0
de.gstcell.biopsy$log2FoldChange.y[is.na(de.gstcell.biopsy$log2FoldChange.y)] <- 0
#Adding value of 1 to NAs in padj column
de.gstcell.dice$padj.x[is.na(de.gstcell.dice$padj.x)] <- 1
de.gstcell.dice$padj.y[is.na(de.gstcell.dice$padj.y)] <- 1

de.gstcell.biopsy$padj.x[is.na(de.gstcell.biopsy$padj.x)] <- 1
de.gstcell.biopsy$padj.y[is.na(de.gstcell.biopsy$padj.y)] <- 1

####filtering
de.gstcell.dice <- de.gstcell.dice[abs(de.gstcell.dice$log2FoldChange.x) > 1 | 
                                     abs(de.gstcell.dice$log2FoldChange.y) > 1,]
de.gstcell.biopsy <- de.gstcell.biopsy[abs(de.gstcell.biopsy$log2FoldChange.x) > 1 | 
                                         abs(de.gstcell.biopsy$log2FoldChange.y) > 1,]

de.gstcell.dice <- de.gstcell.dice[de.gstcell.dice$padj.x < 0.05 | 
                                     de.gstcell.dice$padj.y < 0.05, ]
de.gstcell.biopsy <- de.gstcell.biopsy[de.gstcell.biopsy$padj.x < 0.05 | 
                                         de.gstcell.biopsy$padj.y < 0.05, ]

#Creating an extra column with the overlap category of the genes
de.gstcell.dice$Overlap <- paste(de.gstcell.dice$padj.x < 0.05 &
                                   abs(de.gstcell.dice$log2FoldChange.x) > 1,
                                 de.gstcell.dice$padj.y < 0.05 & 
                                   abs(de.gstcell.dice$log2FoldChange.y) > 1, sep = "")
de.gstcell.dice$Overlap <- gsub("TRUETRUE", "Both", de.gstcell.dice$Overlap)

de.gstcell.dice$Overlap <- gsub("TRUENA", "gsTcells", de.gstcell.dice$Overlap)
de.gstcell.dice$Overlap <- gsub("TRUEFALSE", "gsTcells", de.gstcell.dice$Overlap)
de.gstcell.dice$Overlap <- gsub("NATRUE", "DICE", de.gstcell.dice$Overlap)
de.gstcell.dice$Overlap <- gsub("FALSETRUE", "DICE", de.gstcell.dice$Overlap)

de.gstcell.dice$Overlap <- gsub("NANA", NA, de.gstcell.dice$Overlap)
de.gstcell.dice$Overlap <- gsub("NAFALSE", NA, de.gstcell.dice$Overlap)
de.gstcell.dice$Overlap <- gsub("FALSENA", NA, de.gstcell.dice$Overlap)
de.gstcell.dice$Overlap <- gsub("FALSEFALSE", NA, de.gstcell.dice$Overlap)


de.gstcell.biopsy$Overlap <- paste(de.gstcell.biopsy$padj.x < 0.05 &
                                     abs(de.gstcell.biopsy$log2FoldChange.x) > 1,
                                   de.gstcell.biopsy$padj.y < 0.05 & 
                                     abs(de.gstcell.biopsy$log2FoldChange.y) > 1, sep = "")
de.gstcell.biopsy$Overlap <- gsub("TRUETRUE", "Both", de.gstcell.biopsy$Overlap)

de.gstcell.biopsy$Overlap <- gsub("TRUENA", "gsTcells", de.gstcell.biopsy$Overlap)
de.gstcell.biopsy$Overlap <- gsub("TRUEFALSE", "gsTcells", de.gstcell.biopsy$Overlap)
de.gstcell.biopsy$Overlap <- gsub("NATRUE", "biopsy", de.gstcell.biopsy$Overlap)
de.gstcell.biopsy$Overlap <- gsub("FALSETRUE", "biopsy", de.gstcell.biopsy$Overlap)

de.gstcell.biopsy$Overlap <- gsub("NANA", NA, de.gstcell.biopsy$Overlap)
de.gstcell.biopsy$Overlap <- gsub("NAFALSE", NA, de.gstcell.biopsy$Overlap)
de.gstcell.biopsy$Overlap <- gsub("FALSENA", NA, de.gstcell.biopsy$Overlap)
de.gstcell.biopsy$Overlap <- gsub("FALSEFALSE", NA, de.gstcell.biopsy$Overlap)

#Cleaning
de.gstcell.dice <- de.gstcell.dice[!is.na(de.gstcell.dice$Overlap), ]
de.gstcell.biopsy <- de.gstcell.biopsy[!is.na(de.gstcell.biopsy$Overlap), ]


de.gstcell.dice$Overlap <- as.factor(de.gstcell.dice$Overlap)
de.gstcell.biopsy$Overlap <- as.factor(de.gstcell.biopsy$Overlap)

de.gstcell.dice$Overlap <- factor(de.gstcell.dice$Overlap, levels = c("DICE", "gsTcells", "Both"))
de.gstcell.biopsy$Overlap <- factor(de.gstcell.biopsy$Overlap, levels = c("biopsy", "gsTcells", "Both"))

de.gstcell.dice <- de.gstcell.dice %>% distinct(ENSEMBL, .keep_all = TRUE)
de.gstcell.biopsy <- de.gstcell.biopsy %>% distinct(ENSEMBL, .keep_all = TRUE)

#plotting concordant GST cells and Biopsy
annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("Q3", "Q4", "Q2","Q1"),
  hjustvar = c(-1, -1, 2, 2) ,
  vjustvar = c(-1, 2, -1, 2)) #<- adjust

p.DICE <- ggplot(data = de.gstcell.dice, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
  geom_point(aes(colour = Overlap), alpha = 0.5, size = 2) + 
  #geom_point(data = subset(de.gstcell.dice, Overlap == "DICE"), aes(colour = Overlap), alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #geom_text(data = subset(de.gstcell.dice, Gene == "IL21"), aes(label = Gene),hjust=0,vjust=0) +
  ylab("L2FC in DICE") +
  xlab("L2FC in gsTcells") +
  geom_text(data=annotations,aes(x = xpos, y = ypos,hjust=hjustvar,vjust=vjustvar,label= annotateText), size = 4) +
  xlim(-8, 18) +
  ylim(-8, 18) +
  coord_fixed() +
  scale_color_manual(values=c("#8CED3C", "#4B4BCD", "#30BCBC")) +
  ggtitle("")
p.DICE <- theme.plain(p.DICE) + theme(legend.position = "none")# + guides(colour=guide_legend(nrow = 2, byrow = TRUE))

p.Biopsy <- ggplot(data = de.gstcell.biopsy, aes(x = log2FoldChange.x, y = log2FoldChange.y)) + 
  geom_point(aes(colour = Overlap), alpha = 0.5, size = 2) + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  ylab("L2FC in Biopsy") +
  xlab("L2FC in gsTcells") +
  geom_text(data=annotations,aes(x = xpos, y = ypos,hjust=hjustvar,vjust=vjustvar,label= annotateText), size = 4) +
  xlim(-8, 18) +
  ylim(-8, 18) +
  coord_fixed() +
  scale_color_manual(values=c("#8CED3C", "#4B4BCD", "#30BCBC")) +
  ggtitle("")
p.Biopsy <- theme.plain(p.Biopsy) + theme(legend.position = "none")# + guides(colour=guide_legend(nrow = 2, byrow = TRUE))

####Enrichments

de.gstcell.dice <- dplyr::left_join(de.gstcell.dice, ensembl, by = "ENSEMBL")
de.gstcell.biopsy <- dplyr::left_join(de.gstcell.biopsy, ensembl, by = "ENSEMBL")

###Adding quadrant data
de.gstcell.dice$Quadrant <- NA
de.gstcell.dice[de.gstcell.dice$Overlap == "Both" & 
                  de.gstcell.dice$log2FoldChange.x > 1 & 
                  de.gstcell.dice$log2FoldChange.y > 1, ]$Quadrant <- "Q1"
de.gstcell.dice[de.gstcell.dice$Overlap == "Both" & 
                  de.gstcell.dice$log2FoldChange.x > 1 & 
                  de.gstcell.dice$log2FoldChange.y < -1, ]$Quadrant <- "Q2"
de.gstcell.dice[de.gstcell.dice$Overlap == "Both" & 
                  de.gstcell.dice$log2FoldChange.x < -1 & 
                  de.gstcell.dice$log2FoldChange.y < -1, ]$Quadrant <- "Q3"
de.gstcell.dice[de.gstcell.dice$Overlap == "Both" & 
                  de.gstcell.dice$log2FoldChange.x < -1 & 
                  de.gstcell.dice$log2FoldChange.y > 1, ]$Quadrant <- "Q4"

de.gstcell.biopsy$Quadrant <- NA
de.gstcell.biopsy[de.gstcell.biopsy$Overlap == "Both" & 
                    de.gstcell.biopsy$log2FoldChange.x > 1 & 
                    de.gstcell.biopsy$log2FoldChange.y > 1, ]$Quadrant <- "Q1"
de.gstcell.biopsy[de.gstcell.biopsy$Overlap == "Both" & 
                    de.gstcell.biopsy$log2FoldChange.x > 1 & 
                    de.gstcell.biopsy$log2FoldChange.y < -1, ]$Quadrant <- "Q2"
de.gstcell.biopsy[de.gstcell.biopsy$Overlap == "Both" & 
                    de.gstcell.biopsy$log2FoldChange.x < -1 & 
                    de.gstcell.biopsy$log2FoldChange.y < -1, ]$Quadrant <- "Q3"
de.gstcell.biopsy[de.gstcell.biopsy$Overlap == "Both" & 
                    de.gstcell.biopsy$log2FoldChange.x < -1 & 
                    de.gstcell.biopsy$log2FoldChange.y > 1, ]$Quadrant <- "Q4"

de.gstcell.dice$gene <- NULL
de.gstcell.dice$Gene <- NULL
de.gstcell.dice <- de.gstcell.dice %>% distinct(ENSEMBL, .keep_all = TRUE)
de.gstcell.biopsy <- de.gstcell.biopsy %>% distinct(ENSEMBL, .keep_all = TRUE)
#write.csv(de.gstcell.dice, file = "./Results/de.gstcell.dice.v3.csv")
#write.csv(de.gstcell.biopsy, file = "./Results/de.gstcell.biopsyv2.csv")

de.gstcell.dice.2 <- de.gstcell.dice[,c("ENSEMBL", "ENTREZID", "SYMBOL", "log2FoldChange.x", "padj.x", "Quadrant")]
names(de.gstcell.dice.2) <- c("ENSEMBL", "ENTREZID", "SYMBOL", "log2FoldChange", "padj", "Quadrant")
de.gstcell.dice.2 <- de.gstcell.dice.2[abs(de.gstcell.dice.2$log2FoldChange) > 1 & de.gstcell.dice.2$padj < 0.05,]

de.gstcell.dice.2 <- de.gstcell.dice.2 %>% distinct(ENTREZID, .keep_all = TRUE)

short.de.gstcell.dice <- de.gstcell.dice.2[,c("ENTREZID", "log2FoldChange", "Quadrant")]

list.de.gstcell.dice <- list(Q1 = subset(short.de.gstcell.dice, Quadrant == "Q1"), 
                             Q2 = subset(short.de.gstcell.dice, Quadrant == "Q2"),
                             Q3 = subset(short.de.gstcell.dice, Quadrant == "Q3"),
                             Q4 = subset(short.de.gstcell.dice, Quadrant == "Q4"))
list.de.gstcell.dice <- lapply(list.de.gstcell.dice, function(Q){
  Q$Quadrant <- NULL
  return(Q)
})



#list with vector names as genes and L2FC
list.de.gstcell.dice.reactomeformatted <- lapply(list.de.gstcell.dice , formatting.gene.list.second)
#list with only gene names
list.de.gstcell.dice.reactomeformatted.2 <- lapply(list.de.gstcell.dice.reactomeformatted, gene.list.only.names)



#Biopsy reactome
de.gstcell.biopsy.2 <- de.gstcell.biopsy[,c("ENSEMBL", "ENTREZID", "SYMBOL", "log2FoldChange.x", "padj.x", "Quadrant")]
names(de.gstcell.biopsy.2) <- c("ENSEMBL", "ENTREZID", "SYMBOL", "log2FoldChange", "padj", "Quadrant")
de.gstcell.biopsy.2 <- de.gstcell.biopsy.2[abs(de.gstcell.biopsy.2$log2FoldChange) > 1 & de.gstcell.biopsy.2$padj < 0.05,]

de.gstcell.biopsy.2 <- de.gstcell.biopsy.2 %>% distinct(ENTREZID, .keep_all = TRUE)

short.de.gstcell.biopsy <- de.gstcell.biopsy.2[,c("ENTREZID", "log2FoldChange", "Quadrant")]

list.de.gstcell.biopsy <- list(Q1 = subset(short.de.gstcell.biopsy, Quadrant == "Q1"), 
                               Q2 = subset(short.de.gstcell.biopsy, Quadrant == "Q2"),
                               Q3 = subset(short.de.gstcell.biopsy, Quadrant == "Q3"),
                               Q4 = subset(short.de.gstcell.biopsy, Quadrant == "Q4"))
list.de.gstcell.biopsy <- lapply(list.de.gstcell.biopsy, function(Q){
  Q$Quadrant <- NULL
  return(Q)
})
#list with vector names as genes and L2FC
list.de.gstcell.biopsy.reactomeformatted <- lapply(list.de.gstcell.biopsy , formatting.gene.list.second)
#list with only gene names
list.de.gstcell.biopsy.reactomeformatted.2 <- lapply(list.de.gstcell.biopsy.reactomeformatted, gene.list.only.names)


#REACTOME enrichment
res.df.dice <- compareCluster(list.de.gstcell.dice.reactomeformatted.2, fun="enrichPathway", readable=T)

res.df.biopsy <- compareCluster(list.de.gstcell.biopsy.reactomeformatted.2, fun="enrichPathway", readable=T)

#PLOTS
# Convert adjusted p to -log10 scale
res.df.dice@compareClusterResult$p.adjust <- -log10(res.df.dice@compareClusterResult$p.adjust)
# Convert adjusted p to -log10 scale
res.df.biopsy@compareClusterResult$p.adjust <- -log10(res.df.biopsy@compareClusterResult$p.adjust)
#write.csv(res.df.dice, file = "./Results/reactome.de.dice.gstcellv2.csv")
#write.csv(res.df.biopsy, file = "./Results/reactome.de.biopsy.gstcellv2.csv")

p.react.dice <- theme.plain(dotplot(res.df.dice, font.size = 8, showCategory = 10), base_size = 8) +
  scale_color_gradient(low="lightgrey", high = "red", name="-log10(p-adjust)") 

p.react.biopsy <- theme.plain(dotplot(res.df.biopsy, font.size = 8, showCategory = 10), base_size = 8) +
  scale_color_gradient(low="lightgrey", high = "red", name="-log10(p-adjust)") 



###############



de.gstcell.dice.gst.2 <- de.gstcell.dice[,c("ENSEMBL", "ENTREZID", "SYMBOL", "log2FoldChange.x", "padj.x", "Overlap", "Direction")]
names(de.gstcell.dice.gst.2) <- c("ENSEMBL", "ENTREZID", "SYMBOL", "log2FoldChange", "padj", "Overlap", "Direction")
de.gstcell.dice.gst.2 <- de.gstcell.dice.gst.2[abs(de.gstcell.dice.gst.2$log2FoldChange) > 2 & de.gstcell.dice.gst.2$padj < 0.05,]

de.gstcell.dice.gst.2 <- de.gstcell.dice.gst.2 %>% distinct(ENTREZID, .keep_all = TRUE)
de.gstcell.dice.gst.2 <- de.gstcell.dice.gst.2[de.gstcell.dice.gst.2$Overlap == "gsTcells",]

short.de.gstcell.dice.gst <- de.gstcell.dice.gst.2[,c("ENTREZID", "log2FoldChange", "Direction")]

list.de.gstcell.dice.gst <- list(Upregulated = subset(short.de.gstcell.dice.gst, Direction == "Upregulated"), 
                                 Downregulated = subset(short.de.gstcell.dice.gst, Direction == "Downregulated"))
list.de.gstcell.dice.gst <- lapply(list.de.gstcell.dice.gst, function(Q){
  Q$Direction <- NULL
  return(Q)
})

#list with vector names as genes and L2FC
list.de.gstcell.dice.gst.reactomeformatted <- lapply(list.de.gstcell.dice.gst , formatting.gene.list.second)
#list with only gene names
list.de.gstcell.dice.gst.reactomeformatted.2 <- lapply(list.de.gstcell.dice.gst.reactomeformatted, gene.list.only.names)

res.df.dice.gst <- compareCluster(list.de.gstcell.dice.gst.reactomeformatted.2, fun="enrichPathway", readable=T)

res.df.dice.gst@compareClusterResult$p.adjust <- -log10(res.df.dice.gst@compareClusterResult$p.adjust)


#write.csv(res.df.dice.gst, file = "./Results/reactome.de.dice.gstcell.gst.csv")

### dice specific genes reactome

#important, using columns with .y and changing overlap for "DICE"
de.gstcell.dice.dice.2 <- de.gstcell.dice[,c("ENSEMBL", "ENTREZID", "SYMBOL", "log2FoldChange.y", "padj.y", "Overlap", "Direction")]
names(de.gstcell.dice.dice.2) <- c("ENSEMBL", "ENTREZID", "SYMBOL", "log2FoldChange", "padj", "Overlap", "Direction")

de.gstcell.dice.dice.2$Direction <- "Upregulated"
de.gstcell.dice.dice.2$Direction[de.gstcell.dice.dice.2$log2FoldChange < 0] <- "Downregulated"

de.gstcell.dice.dice.2 <- de.gstcell.dice.dice.2[abs(de.gstcell.dice.dice.2$log2FoldChange) > 2 & de.gstcell.dice.dice.2$padj < 0.05,]

de.gstcell.dice.dice.2 <- de.gstcell.dice.dice.2 %>% distinct(ENTREZID, .keep_all = TRUE)
de.gstcell.dice.dice.2 <- de.gstcell.dice.dice.2[de.gstcell.dice.dice.2$Overlap == "DICE",]

short.de.gstcell.dice.dice <- de.gstcell.dice.dice.2[,c("ENTREZID", "log2FoldChange", "Direction")]

list.de.gstcell.dice.dice <- list(Upregulated = subset(short.de.gstcell.dice.dice, Direction == "Upregulated"), 
                                  Downregulated = subset(short.de.gstcell.dice.dice, Direction == "Downregulated"))
list.de.gstcell.dice.dice <- lapply(list.de.gstcell.dice.dice, function(Q){
  Q$Direction <- NULL
  return(Q)
})

#list with vector names as genes and L2FC
list.de.gstcell.dice.dice.reactomeformatted <- lapply(list.de.gstcell.dice.dice , formatting.gene.list.second)
#list with only gene names
list.de.gstcell.dice.dice.reactomeformatted.2 <- lapply(list.de.gstcell.dice.dice.reactomeformatted, gene.list.only.names)

res.df.dice.dice <- compareCluster(list.de.gstcell.dice.dice.reactomeformatted.2, fun="enrichPathway", readable=T)

res.df.dice.dice@compareClusterResult$p.adjust <- -log10(res.df.dice.dice@compareClusterResult$p.adjust)

#write.csv(res.df.dice.dice, file = "./Results/reactome.de.dice.gstcell.dice.csv")

####### Specific genes reactome enrichment gstcells compared with biopsy
### gstcells
de.gstcell.biopsy.gst.2 <- de.gstcell.biopsy[,c("ENSEMBL", "ENTREZID", "SYMBOL", "log2FoldChange.x", "padj.x", "Overlap", "Direction")]
names(de.gstcell.biopsy.gst.2) <- c("ENSEMBL", "ENTREZID", "SYMBOL", "log2FoldChange", "padj", "Overlap", "Direction")
de.gstcell.biopsy.gst.2 <- de.gstcell.biopsy.gst.2[abs(de.gstcell.biopsy.gst.2$log2FoldChange) > 1 & de.gstcell.biopsy.gst.2$padj < 0.05,]

de.gstcell.biopsy.gst.2 <- de.gstcell.biopsy.gst.2 %>% distinct(ENTREZID, .keep_all = TRUE)
de.gstcell.biopsy.gst.2 <- de.gstcell.biopsy.gst.2[de.gstcell.biopsy.gst.2$Overlap == "gsTcells",]

short.de.gstcell.biopsy.gst <- de.gstcell.biopsy.gst.2[,c("ENTREZID", "log2FoldChange", "Direction")]

list.de.gstcell.biopsy.gst <- list(Upregulated = subset(short.de.gstcell.biopsy.gst, Direction == "Upregulated"), 
                                   Downregulated = subset(short.de.gstcell.biopsy.gst, Direction == "Downregulated"))
list.de.gstcell.biopsy.gst <- lapply(list.de.gstcell.biopsy.gst, function(Q){
  Q$Direction <- NULL
  return(Q)
})

#list with vector names as genes and L2FC
list.de.gstcell.biopsy.gst.reactomeformatted <- lapply(list.de.gstcell.biopsy.gst , formatting.gene.list.second)
#list with only gene names
list.de.gstcell.biopsy.gst.reactomeformatted.2 <- lapply(list.de.gstcell.biopsy.gst.reactomeformatted, gene.list.only.names)

res.df.biopsy.gst <- compareCluster(list.de.gstcell.biopsy.gst.reactomeformatted.2, fun="enrichPathway", readable=T)

res.df.biopsy.gst@compareClusterResult$p.adjust <- -log10(res.df.biopsy.gst@compareClusterResult$p.adjust)


#write.csv(res.df.biopsy.gst, file = "./Results/reactome.de.biopsy.gstcell.gst.csv")

### biopsy specific genes reactome

#important, using columns with .y and changing overlap for "biopsy"
de.gstcell.biopsy.biopsy.2 <- de.gstcell.biopsy[,c("ENSEMBL", "ENTREZID", "SYMBOL", "log2FoldChange.y", "padj.y", "Overlap", "Direction")]
names(de.gstcell.biopsy.biopsy.2) <- c("ENSEMBL", "ENTREZID", "SYMBOL", "log2FoldChange", "padj", "Overlap", "Direction")

de.gstcell.biopsy.biopsy.2$Direction <- "Upregulated"
de.gstcell.biopsy.biopsy.2$Direction[de.gstcell.biopsy.biopsy.2$log2FoldChange < 0] <- "Downregulated"

de.gstcell.biopsy.biopsy.2 <- de.gstcell.biopsy.biopsy.2[abs(de.gstcell.biopsy.biopsy.2$log2FoldChange) > 1 & de.gstcell.biopsy.biopsy.2$padj < 0.05,]

de.gstcell.biopsy.biopsy.2 <- de.gstcell.biopsy.biopsy.2 %>% distinct(ENTREZID, .keep_all = TRUE)
de.gstcell.biopsy.biopsy.2 <- de.gstcell.biopsy.biopsy.2[de.gstcell.biopsy.biopsy.2$Overlap == "biopsy",]

short.de.gstcell.biopsy.biopsy <- de.gstcell.biopsy.biopsy.2[,c("ENTREZID", "log2FoldChange", "Direction")]

list.de.gstcell.biopsy.biopsy <- list(Upregulated = subset(short.de.gstcell.biopsy.biopsy, Direction == "Upregulated"), 
                                      Downregulated = subset(short.de.gstcell.biopsy.biopsy, Direction == "Downregulated"))
list.de.gstcell.biopsy.biopsy <- lapply(list.de.gstcell.biopsy.biopsy, function(Q){
  Q$Direction <- NULL
  return(Q)
})

#list with vector names as genes and L2FC
list.de.gstcell.biopsy.biopsy.reactomeformatted <- lapply(list.de.gstcell.biopsy.biopsy , formatting.gene.list.second)
#list with only gene names
list.de.gstcell.biopsy.biopsy.reactomeformatted.2 <- lapply(list.de.gstcell.biopsy.biopsy.reactomeformatted, gene.list.only.names)

res.df.biopsy.biopsy <- compareCluster(list.de.gstcell.biopsy.biopsy.reactomeformatted.2, fun="enrichPathway", readable=T)

res.df.biopsy.biopsy@compareClusterResult$p.adjust <- -log10(res.df.biopsy.biopsy@compareClusterResult$p.adjust)

#write.csv(res.df.biopsy.biopsy, file = "./Results/reactome.de.biopsy.gstcell.biopsy.csv")



###Enrichment for specific genes solid

### dice specific genes reactome

#important, using columns with .y and changing overlap for "solid"
de.gstcell.solid.solid.2 <- de.gstcell.solid[,c("ENSEMBL", "ENTREZID", "SYMBOL.y", "log2FoldChange.y", "padj.y", "Overlap", "Direction")]
names(de.gstcell.solid.solid.2) <- c("ENSEMBL", "ENTREZID", "SYMBOL", "log2FoldChange", "padj", "Overlap", "Direction")

de.gstcell.solid.solid.2$Direction <- "Upregulated"
de.gstcell.solid.solid.2$Direction[de.gstcell.solid.solid.2$log2FoldChange < 0] <- "Downregulated"

de.gstcell.solid.solid.2 <- de.gstcell.solid.solid.2[abs(de.gstcell.solid.solid.2$log2FoldChange) > 2 & de.gstcell.solid.solid.2$padj < 0.05,]

de.gstcell.solid.solid.2 <- de.gstcell.solid.solid.2 %>% distinct(ENTREZID, .keep_all = TRUE)
de.gstcell.solid.solid.2 <- de.gstcell.solid.solid.2[de.gstcell.solid.solid.2$Overlap == "solid",]

short.de.gstcell.solid.solid <- de.gstcell.solid.solid.2[,c("ENTREZID", "log2FoldChange", "Direction")]

list.de.gstcell.solid.solid <- list(Upregulated = subset(short.de.gstcell.solid.solid, Direction == "Upregulated"), 
                                  Downregulated = subset(short.de.gstcell.solid.solid, Direction == "Downregulated"))
list.de.gstcell.solid.solid <- lapply(list.de.gstcell.solid.solid, function(Q){
  Q$Direction <- NULL
  return(Q)
})

#list with vector names as genes and L2FC
list.de.gstcell.solid.solid.reactomeformatted <- lapply(list.de.gstcell.solid.solid , formatting.gene.list.second)
#list with only gene names
list.de.gstcell.solid.solid.reactomeformatted.2 <- lapply(list.de.gstcell.solid.solid.reactomeformatted, gene.list.only.names)

res.df.solid.solid <- compareCluster(list.de.gstcell.solid.solid.reactomeformatted.2, fun="enrichPathway", readable=T)

res.df.solid.solid@compareClusterResult$p.adjust <- -log10(res.df.solid.solid@compareClusterResult$p.adjust)

#write.csv(res.df.solid.solid, file = "./Results/reactome.de.solid.gstcell.solid.csv")


#Upset plot
#upset.p


#dotplots of comparisons with quadrants
#p.Biopsy
#p.DICE
#p.solid


#all reactome plots of quadrants genes
#Enrichment of quadrants for dice comparison
#p.react.dice

#Enrichment of quadrants for biopsy comparison
#p.react.biopsy

#Enrichment of quadrants for solid comparison
#p.react.solid
#p.react.gst.shared


#Enrichments of specific sets of DE genes
#gstcells vs DICE: gstcells
p.react.dice.gst <- theme.plain(dotplot(res.df.dice.gst, 
                                        font.size = 8, 
                                        showCategory = 10), 
                                base_size = 8) +
  scale_color_gradient(low="lightgrey", high = "red", name="-log10(p-adjust)") 

#gstcells vs DICE: dice
p.react.dice.dice <- theme.plain(dotplot(res.df.dice.dice, 
                                         font.size = 8, 
                                         showCategory = 10), 
                                 base_size = 8) +
  scale_color_gradient(low="lightgrey", high = "red", name="-log10(p-adjust)") 

#gstcells vs biopsy: gstcells
p.react.biopsy.gst <- theme.plain(dotplot(res.df.biopsy.gst, 
                                          font.size = 8, 
                                          showCategory = 10), 
                                  base_size = 8) +
  scale_color_gradient(low="lightgrey", high = "red", name="-log10(p-adjust)") 

#gstcells vs biopsy: biopsy
p.react.biopsy.biopsy <- theme.plain(dotplot(res.df.biopsy.biopsy, 
                                             font.size = 8, 
                                             showCategory = 10), 
                                     base_size = 8) +
  scale_color_gradient(low="lightgrey", high = "red", name="-log10(p-adjust)") 
#gstcells vs solid vs DICE: gstcells
#p.react.gst.specific 

#gstcells vs solid: solid
p.react.solid.solid <- theme.plain(dotplot(res.df.solid.solid, 
                                         font.size = 8, 
                                         showCategory = 10), 
                                 base_size = 8) +
  scale_color_gradient(low="lightgrey", high = "red", name="-log10(p-adjust)") 



#

#Fig3, Upset, doptplots, enrihcment of specific DE gstcells
pdf(file="Results/Fig3_upset_concordance_enrichmentv7.pdf", width=9, height=12,  family = "ArialMT", paper="a4")
ggarrange(as.ggplot(upset.p, scale = 1), 
          ggarrange(p.DICE, p.solid, p.Biopsy, ncol = 3, nrow = 1, widths = c(1,1,1)), 
          ggarrange(p.react.dice.gst, p.react.gst.specific, ncol = 2, widths = c(1,1)),
          p.react.biopsy.gst, 
          ncol = 1, nrow = 4, heights = c(1.5,1,1,1.5)) 
decorate_annotation("Intersection\nsize", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))})
dev.off()
embedFonts(file = "Results/Fig3_upset_concordance_enrichmentv5.pdf", 
           outfile = "Results/Fig3_upset_concordance_enrichmentv5.pdf")


#SupplFig Enrichments of quadrants 
pdf(file="Results/SupFig3_quadrant_enrichmentv5.pdf", width=9, height=12,  family = "ArialMT", paper="a4")
ggarrange(p.react.dice, 
          p.react.biopsy, 
          p.react.gst.shared,
          ncol = 1, nrow = 3, heights = c(1,2,1), widths = c(0.8, 1, 0.5), labels = c("A", "B", "C")) 
dev.off()

#SupplFig Enrichments of specific DE external
pdf(file="Results/SupFig3_upset_external_specific_enrichmentv3.pdf", width=9, height=12,  family = "ArialMT", paper="a4")
ggarrange(p.react.dice.dice, 
          p.react.biopsy.biopsy, 
          ncol = 1, nrow = 4, heights = c(1.5,2), labels = c("A","B","C")) 
dev.off()

#gst  dice  both
#1469 6836  1926 
#Q1     Q2    Q3    Q4 
#1155   73    687   11 

#biopsy   gsTcells    Both 
#1291     3251        144 
#Q1 Q2 Q3 Q4 
#55 29 32 28 

#Quadrants solid
#Q1 Q2 Q3 Q4 
#47 27 36 13 