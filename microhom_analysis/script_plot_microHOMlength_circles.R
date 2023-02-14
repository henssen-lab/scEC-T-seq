####### BOXPLOTS + PVALUE MICROHOMOLOGY IN CIRCLES LENGTH DISTRIBUTION TR14 & CHP212 ####### 

#libraries
library(ggplot2)

#Open files
tr14 <- read.table("/Users/elias/Documents/ecDNA_projects.CHARITE/Singlecell_ecDNA_Rocio_project.CHARITE/Circle_coords/tr14_exo5/ALL_cells_tr14_exo5.BLAST_RESULTS_HOMOLOGY_min4bp_100bparoundbkp.txt", sep="\t", stringsAsFactors=FALSE, col.names=c("circle_id", "HOM_length", "evalue"))
chp212 <- read.table("/Users/elias/Documents/ecDNA_projects.CHARITE/Singlecell_ecDNA_Rocio_project.CHARITE/Circle_coords/chp212_exo5/ALL_cells_chp212_exo5.BLAST_RESULTS_HOMOLOGY_min4bp_100bparoundbkp.txt", sep="\t", stringsAsFactors=FALSE, col.names=c("circle_id", "HOM_length", "evalue"))

#Join all
tr14$cell_id <- "tr14"
chp212$cell_id <- "chp212"

ALL_circs <- rbind(tr14, chp212)


#### Statistics

#Difference between cells (tr14 vs chp212)?
kruskal.test(HOM_length ~ cell_id, data = ALL_circs)
#Wilcoxon test
pairwise.wilcox.test(jitter(ALL_circs$HOM_length),ALL_circs$cell_id,p.adjust.method ="bonferroni")


#### PLOT

#P <- ggplot(ALL_circs, aes(x=cell_id, y=HOM_length, color=cell_id, fill=cell_id)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, outlier.colour="grey", outlier.size=0.5) + scale_y_continuous(trans='log10') + theme_classic() + labs(title="Distribution of microhomology length per circle",x="Cell line", y = "Length (bp)")
P <- ggplot(ALL_circs, aes(x=cell_id, y=HOM_length, color=cell_id, fill=cell_id)) + geom_violin(alpha = 1/5, color=NA) + geom_boxplot(width=0.2, fill=NA, outlier.size=0.1) + scale_y_continuous(trans='log10') + theme_classic() + labs(title="Distribution of microhomology length per circle",x="Cell line", y = "Length (bp)")
ggsave("plot_microHOMlength_circles_chp212+tr14_22.05.31.pdf", plot=P, width = 5, height = 5,)


ALL_circs$cell_id <- "both"
P_all <- ggplot(ALL_circs, aes(x=cell_id, y=HOM_length, color=cell_id, fill=cell_id)) + geom_violin(alpha = 1/5, color=NA) + geom_boxplot(width=0.2, fill=NA, outlier.size=0.1) + scale_y_continuous(trans='log10') + theme_classic() + labs(title="Distribution of microhomology length per circle",x="Cell line", y = "Length (bp)")
ggsave("plot_microHOMlength_circles_chp212+tr14_TOGETHER_22.05.31.pdf", plot=P_all, width = 5, height = 5,)

