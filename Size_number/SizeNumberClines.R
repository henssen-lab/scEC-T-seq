# Size and number of circles in cell lines
# Rocio Chamorro

library(ggpubr)
library(ggplot2)
library(stringr)
library(superheat)
library(RColorBrewer)

##########1. Create a table including all circular DNA regions identified in all single cells from TR14 and CHP-212

setwd("/Users/rociochamorro/scECTseq/circleCalls/unfiltered/")

## Save all circle calls from all cells 
files<- list.files(pattern = ".bed") # bed files including the coordinates of circles called in each single cell
sc<- lapply(files, read.table)
samples_name_sc<- str_remove(files, ".enriched.merged.b2bRefined.Counts.ThreshFinal.bed")
samples_name_sc<- str_remove(samples_name_sc, "Regions_")

samples_name<- NULL
samples_name_v<- NULL
for (i in 1:length(samples_name_sc)) {
  samples_name<- paste(unlist(strsplit(samples_name_sc[i], "_"))[1],unlist(strsplit(samples_name_sc[i], "_"))[2], sep="_")
  samples_name_v<- append(samples_name_v, samples_name)
}

samples_name_sc<- samples_name_v 
names(sc)<- samples_name_sc

circles<-NULL
all_circles<- NULL
cell_ids<- NULL
cell_id<- NULL

for (i in 1:length(sc)) {
  circles<- sc[[i]]
  all_circles<- rbind(all_circles, circles)
  cell_id<- rep(names(sc[i]), length(circles[,1]))
  cell_ids<- append(cell_ids, cell_id)
}

all_circles<- cbind(all_circles, cell_ids)
all_circles_threshold<-all_circles[which(all_circles$V4>=2),] # We keep only circles with more than 2 circle-supporting reads
all_circles_threshold<-all_circles_threshold[-which(all_circles_threshold$V1=="plasmid2_endo"),]  # Remove plasmid_endo and plasmid linear
all_circles_threshold<-all_circles_threshold[-which(all_circles_threshold$V1=="plasmid3_linear"),]

# Keep only 5days exonuclease digestion condition and include cell line information

metadata<- readRDS("/Users/rociochamorro/scECTseq/RDS/metadata_cLines_Tcell_Tumor.RDS")
metadata<- metadata[which(metadata$protocol == "+ Exo 5days" & metadata$qc == "Yes"),]  
all_circles_threshold_filtered<- all_circles_threshold[which(all_circles_threshold$cell_ids %in% metadata$unique_id),]

cline<- NULL
cline_v<- NULL
for (i in 1:length(all_circles_threshold_filtered$cell_ids)) {
  cline<- metadata[which(metadata$unique_id == all_circles_threshold_filtered[i,]$cell_id),]$subject_id
  cline_v<- append(cline_v, cline)
}

all_circles_threshold_filtered$cell_line<- cline_v
colnames(all_circles_threshold_filtered)<- c("chr","start","end","circle_supportingReads","totalReads","unique_id","cell_line")

# Create table 

#write.csv(all_circles_threshold_filtered, "/Users/rociochamorro/scECTseq/RDS/CircularRegions_clines.csv")
#saveRDS(all_circles_threshold_filtered, "/Users/rociochamorro/scECTseq/RDS/CircularRegions_clines.RDS")

# -------------------------------------------------------------------------------------------------------------------
# Supplementary Fig 4a
# -------------------------------------------------------------------------------------------------------------------

# Calculate circles' length
all_circles_threshold_filtered$circle_length<- all_circles_threshold_filtered$end - all_circles_threshold_filtered$start

# Plot length of all circles per cell line
ggplot(all_circles_threshold_filtered, aes(x=cell_line, y=circle_length, fill=cell_line)) +
  geom_violin() +
  scale_y_continuous(trans='log10') +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 45,hjust = 1)) +
  ylab("Circular DNA length") + 
  xlab("")

# -------------------------------------------------------------------------------------------------------------------
# Figure 2a
# -------------------------------------------------------------------------------------------------------------------

## Create dataframe with number of circles per bin of length
#Each bin consists of 500 bp and it ranges from 1bp to 100,000 bp - from = 1, to = 100000, by = 500

breaks<- seq(from = 1, to = 100000, by = 500) #from 1bp to 1,2MB binsize=10kb
df<- NULL
freq<- NULL
df_length<-NULL
sample<- NULL

for (i in 1:length(unique(all_circles_threshold_filtered$unique_id))) {
  sample<- as.character(unique(all_circles_threshold_filtered$unique_id)[i])
  length<- all_circles_threshold_filtered[which(as.character(all_circles_threshold_filtered$unique_id) == sample),]$circle_length
  df<- t(as.data.frame(table(cut(length, breaks = breaks))))
  freq<- as.numeric(as.character(df[2,]))
  df_length<- rbind(df_length, freq)
}

colnames(df_length)<- df[1,]
row.names(df_length)<- as.character(unique(all_circles_threshold_filtered$unique_id))

## Order rows to separate cell lines

coul <- c("white",colorRampPalette(brewer.pal(5, "Blues"))(25))
row_order<- row.names(df_length)[order(as.numeric(as.character(rowSums(df_length))))]
toMatch<- c("plate5","plate6")
row_order_chp<- row_order[grep(paste(toMatch,collapse="|"), row_order) ]
row_order_tr14<- row_order[grep("plate3", row_order)]

position<- NULL
positions<- NULL
for (i in 1:length(row.names(df_length))) {
  position<- which(row.names(df_length) == row_order_chp[i])
  positions<- append(positions, position)
}

positions_chp<- positions

position<- NULL
positions<- NULL
for (i in 1:length(row.names(df_length))) {
  position<- which(row.names(df_length) == row_order_tr14[i])
  positions<- append(positions, position)
}

positions_tr14<- positions

positions<- c(positions_tr14, positions_chp)

## Plot heatmap 

set.seed(2016113)
superheat(df_length,
          # retain original order of rows/cols
          pretty.order.rows = TRUE,
          pretty.order.cols = FALSE,
          order.rows = positions,
          heat.lim = c(0, 45),
          heat.pal= c("white","#00BFC4","blue"), 
          heat.pal.values = c(0,1),
          yr = as.numeric(as.character(rowSums(df_length))),
          yr.plot.type = "bar",
          yt = as.numeric(as.character(colMeans(df_length))),
          yt.plot.type = "bar")

