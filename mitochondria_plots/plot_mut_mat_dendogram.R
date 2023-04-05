#Rscipt for single cell ECT SNV analysis on mitochondria 
#allele frequency heatmap
#render dendrograms (hierarchical clustering)




library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyverse)

#working dir
setwd("")
#read in Master_mutation_matrix
mut_mat <- read.table("./mutation_matrix.txt", sep = "\t", row.names = 1, header = T, dec=".", stringsAsFactors=FALSE,
                      na.strings=c(".") )
#read in metafile
meta <- read.table("./meta.txt", sep = "\t", header = 1)
#rename TR14 to TR-14 for consistency
meta <- meta %>% 
  mutate(cell_line = str_replace(cell_line, "TR14", "TR-14"))

#get working subset of cells
sample_names <- as.vector(unlist(meta['SNV_samplename']))
Master_filter <- intersect(sample_names,colnames(mut_mat))
meta_all_cells <- meta[meta$SNV_samplename %in% Master_filter, ]



#remove blacklisted rows
blacklisted_pos <- c("310T>C", "3107N>C")
mut_mat[!(row.names(mut_mat) %in% blacklisted_pos), Master_filter]
Mut_mat_all_cells <- mut_mat[!(row.names(mut_mat) %in% blacklisted_pos), Master_filter]

#remove rows with less than 80% information
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}

Mut_mat_all_cells <- delete.na(Mut_mat_all_cells, 35)


#Heatmap_colouring
HM_colour = colorRamp2(c(0, 0.1,1), c("white", "#F29604" ,"#E01D22"))
#Annotation
ha_col <- HeatmapAnnotation("Cell line" = as.character(meta_all_cells$"cell_line"),col = list("Cell line" = c("CHP-212" = "#F04D53", "TR-14" = "#55BADF")))

pdf("./mito_SNV_all_mutations.pdf", width = 15, height = 12)
ComplexHeatmap::Heatmap(Mut_mat_all_cells , use_raster = F, name = "Allele Frequency",top_annotation=ha_col ,col = HM_colour,show_column_names =F,cluster_rows = T,cluster_columns = T,show_row_names = F, clustering_method_columns = "complete")
dev.off()



###################
#Plot only homoplasmic
###################
Mut_mat_all_cells_only_homoplasmic <- Mut_mat_all_cells[rowSums(Mut_mat_all_cells, na.rm=T)>10,]
HM_homoplasmic_colour = colorRamp2(c(0,1), c("white", "#E01D22"))
pdf("./mito_SNV_homoplasmic.pdf", width = 15, height = 12)
ComplexHeatmap::Heatmap(Mut_mat_all_cells_only_homoplasmic , use_raster = F, name = "Allele Frequency",top_annotation=ha_col ,col = HM_homoplasmic_colour,show_column_names =F,cluster_rows = T,cluster_columns = T,show_row_names = T, clustering_method_columns = "complete")
dev.off()

homoplasmic <- row.names(Mut_mat_all_cells_only_homoplasmic)

###################
#Plot only Heteroplsmic CHP212
###################
filter_variants_CHP212 <- read.table("./Filter_CHP212.txt", sep = "\t")
list_CHP212 <- unlist(filter_CHP212)
#subset CHP212
CHP212_meta <- meta_all_cells[meta_all_cells$cell_line == "CHP-212",]
CHP212_cells <- as.vector(CHP212_meta['SNV_samplename'])
CHP212_cells <- unlist(CHP212_cells)
Mut_mat_CHP212 <- Mut_mat_all_cells[list_CHP212,CHP212_cells]
ha_row <- rowAnnotation("Cell line" = as.character(CHP212_meta$"cell_line"),col = list("Cell line" = c("CHP-212" = "#F04D53"))) #row annotation


Mut_mat_CHP212_only_heteroplasmic <- Mut_mat_CHP212[!(row.names(Mut_mat_CHP212) %in% homoplasmic),]
remove <- c("NA")
Mut_mat_CHP212_only_heteroplasmic <- na.omit(Mut_mat_CHP212_only_heteroplasmic)
Mut_mat_CHP212_only_heteroplasmic <- Mut_mat_CHP212_only_heteroplasmic^(1/2) #sqrt transform

af <- Mut_mat_CHP212_only_heteroplasmic 
#filter for variants with max AF > 5% and mean AF < 0.3 and remove cells with mean AF < 0.001 to remove cells with no information 
af3 <- af[rowMaxs(data.matrix(af)) > 0.05 & rowMeans(data.matrix(Mut_mat_CHP212_only_heteroplasmic)) < 0.3, colMeans(data.matrix((Mut_mat_CHP212_only_heteroplasmic))) > 0.001 ]
freq <- apply(af3, 2, table)
freq <- freq[lengths(freq)>2] #filter out cells with only one variant
af3 <- af3[,names(freq)] 
af3 <- af3^(1/2)
CHP212_meta2 <- CHP212_meta # new metatable
rownames(CHP212_meta2) <- CHP212_meta2$SNV_samplename
CHP212_meta3 <- CHP212_meta2[names(freq),]
ha_row <- rowAnnotation("Cell line" = as.character(CHP212_meta3$"cell_line"),col = list("Cell line" = c("CHP-212" = "#F04D53"))) #row annotation




HM_heteroplasmic_colour = colorRamp2(c(0,0.5), c("white", "#E01D22"))
pdf("./revised_mito_SNV_CHP212_heteroplasmic.pdf", width = 15, height = 12)
ComplexHeatmap::Heatmap(t(af3) , use_raster = F, name = "sqrt(Allele Frequency)",row_dend_width = unit(3.5, "cm"),col = HM_heteroplasmic_colour,show_column_names =T,cluster_rows = T,cluster_columns = T,show_row_names = F, clustering_method_columns = "complete", clustering_method_rows = "complete", 
                        row_title = "single cells",left_annotation = ha_row )
dev.off()


Heatmap_obj_chp212 <- ComplexHeatmap::Heatmap(t(af3) , use_raster = F, name = "sqrt(Allele Frequency)",row_dend_width = unit(3.5, "cm"),col = HM_heteroplasmic_colour,show_column_names =T,cluster_rows = T,cluster_columns = T,show_row_names = F, clustering_method_columns = "complete", clustering_method_rows = "complete", 
                        row_title = "single cells",left_annotation = ha_row )
dend <- row_dend(Heatmap_obj_chp212)

library(dendextend)
pdf("./revised_mito_SNV_CHP212_dendrogram.pdf", width = 15, height = 12)
dend %>% 
  set("leaves_pch", 19)  %>% 
  set("leaves_cex", 1) %>% 
  set("leaves_col", "#F04D53") %>% 
  plot(labels=FALSE)
dev.off()


###################
#Plot only Heteroplsmic TR14
###################
filter_variants_TR14 <- read.table("./Filter_TR14.txt", sep = "\t")
list_TR14 <- unlist(filter_CHP212)
#subset TR-14
TR14_meta <- meta_all_cells[meta_all_cells$cell_line == "TR-14",]
TR14_cells <- as.vector(TR14_meta['SNV_samplename'])
TR14_cells <- unlist(TR14_cells)
Mut_mat_TR14 <- Mut_mat_all_cells[list_TR14,TR14_cells]
ha_row <- rowAnnotation("Cell line" = as.character(TR14_meta$"cell_line"),col = list("Cell line" = c("TR-14" = "#55BADF"))) #row annotation

###


af_TR <- Mut_mat_TR14_only_heteroplasmic 
#filter for variants with max AF > 0% and mean AF < 0.3 and remove cells with mean AF < 0.001 to remove cells with no information 


Mut_mat_TR14_only_heteroplasmic <- Mut_mat_TR14[!(row.names(Mut_mat_TR14) %in% homoplasmic),]
Mut_mat_TR14_only_heteroplasmic <- na.omit(Mut_mat_TR14_only_heteroplasmic)
af3_TR <- af_TR[rowMaxs(data.matrix(af_TR)) > 0.01 & rowMeans(data.matrix(Mut_mat_TR14_only_heteroplasmic)) < 0.3, colMeans(data.matrix((Mut_mat_TR14_only_heteroplasmic))) > 0.001 ]
af3_TR <- af3_TR^(1/2)



###

TR_meta2 <- TR14_meta # new metatable
rownames(TR14_meta) <- TR_meta2$SNV_samplename
TR14_meta3 <- TR14_meta[colnames(af3_TR),]

ha_row <- rowAnnotation("Cell line" = as.character(TR14_meta3$"cell_line"),col = list("Cell line" = c("TR-14" = "#55BADF"))) #row annotation

HM_heteroplasmic_colour = colorRamp2(c(0,0.5), c("white", "#E01D22"))
pdf("./revised_mito_SNV_TR_14.pdf", width = 15, height = 12)

ComplexHeatmap::Heatmap(t(af3_TR) , use_raster = F, name = "sqrt(Allele Frequency)",row_dend_width = unit(3.5, "cm"),col = HM_heteroplasmic_colour,show_column_names =T,cluster_rows = T,cluster_columns = T,show_row_names = F, clustering_method_columns = "complete", clustering_method_rows = "complete", 
                        row_title = "single cells", left_annotation = ha_row)
dev.off()

###
##dendrogram
Heatmap_obj_TR <- ComplexHeatmap::Heatmap(t(af3_TR) , use_raster = F, name = "sqrt(Allele Frequency)",row_dend_width = unit(3.5, "cm"),col = HM_heteroplasmic_colour,show_column_names =T,cluster_rows = T,cluster_columns = T,show_row_names = F, clustering_method_columns = "complete", clustering_method_rows = "complete", 
                                          row_title = "single cells", left_annotation = ha_row)
dend <- row_dend(Heatmap_obj_TR)
library(dendextend)
pdf("./revised_mito_SNV_TR_14_dendrogram.pdf", width = 15, height = 12)
dend %>% 
  set("leaves_pch", 19)  %>% 
  set("leaves_cex", 1) %>% 
  set("leaves_col", "#55BADF") %>% 
  plot(labels=FALSE)
dev.off()