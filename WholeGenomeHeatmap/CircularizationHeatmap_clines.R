# Whole-genome circularization heatmaps cell lines
# Rocio Chamorro

library(rtracklayer)
library(stringr)
library(pheatmap)


######### Divide whole genome in bins and calculate the percentage of each bin covered by circle regions 

##1. Create GRegions object that contains the whole genome binned in the same size 

#Interval cutting

chr_length<- read.table("/Users/rociochamorro/Documents/Rocio/AG Henssen/Projects/SCeccDNA & RNA seq/scSeq analyses/Manuscript_analyses/Whole-genome heatmap/R script/hg19_filtered.bed") # bed file which includes the size in bp of each chromosome
chr_length[,1]<- factor(as.character(chr_length[,1]), levels = paste("chr",c(1:22,"X","Y","M"),sep=""))
colnames(chr_length)<- c("chr.name","chr.start","chr.end")
chr_length$divide<- round(as.numeric(as.character(chr_length$chr.end))/3000000) #We picked intervals of 3Mb (bin size)

data<-NULL
start<NULL
end<- NULL
matrix.chromsome <- NULL
matrix.chromsome$chr <- NULL
m.chr<- NULL

for (i in 1:length(chr_length$chr.name)) {
  
  data<-round(seq(chr_length$chr.start[i],chr_length$chr.end[i], length.out =chr_length$divide[i]))
  start=data[1:(chr_length$divide[i]-1)]+1
  end=data[2:chr_length$divide[i]]
  
  matrix.chromsome <- data.frame(start=start,end=end)
  matrix.chromsome$chr <- chr_length$chr.name[i]
  
  m.chr<- rbind(m.chr,matrix.chromsome)
  
}

m.chr<- rbind(m.chr, c(1,16569,"chrM")) #we manually added the length of chrM

#Convert into GRanges

gr.chr<- GRanges(seqnames = m.chr[,3], ranges =IRanges(start = as.numeric(m.chr[,1]), end = as.numeric(m.chr[,2])))

## 2. Find overlap between gr.chr and the circle regions in each single cell and calculate the %region_covered

# We select only exo 5 days condition and quality filtered cells 
# CHP-212
setwd("/Users/rociochamorro/scECTseq/circleCalls/filtered/chp212_exo5/")
files_chp<- dir() # bed files_chp containing the circle regions identified per cell
files_chp<- files_chp[grep(".bed",files_chp)]
m<- data.frame(row.names = length(files_chp))
test_bed<- NULL

for (i in 1:length(files_chp)) {
  
  test_bed<- read.table(files_chp[i])
  test_bed<- test_bed[test_bed$V4>=2,] # Only circles that have more than 2 circle-supporting reads are considered
  colnames(test_bed)<- c("chr","start","end","circle_reads","total_reads")
  test_bed<- GRanges(test_bed)
  
  hits<- NULL
  ov <- NULL
  ov_length<- NULL
  ov_lengths<- NULL
  
  for (j in 1:length(gr.chr)) {
    
    hits<- findOverlaps(gr.chr[j],test_bed) 
    ov <- pintersect(test_bed[subjectHits(hits)], gr.chr[j])
    ov_length<- sum(width(ov)) /sum(width(gr.chr[j])) # fraction of the gr.region that is covered by circular DNA regions
    ov_lengths<- append(ov_lengths, ov_length)
    
  }
  
  m<- cbind(m, ov_lengths)
}

chp_m<- m 
samples_name_chp<- str_remove(files_chp, ".enriched.merged.b2bRefined.Counts.ThreshFinal.bed")
samples_name_chp<- str_remove(samples_name_chp, "Regions_")

samples_name<- NULL
samples_name_v<- NULL
for (i in 1:length(samples_name_chp)) {
  samples_name<- paste(unlist(strsplit(samples_name_chp[i], "_"))[1],unlist(strsplit(samples_name_chp[i], "_"))[2], sep="_")
  samples_name_v<- append(samples_name_v, samples_name)
}

colnames(chp_m)<- samples_name_v # add unique_ids as colnames
rownames(chp_m)<- paste(as.character(seqnames(gr.chr)),as.character(start(gr.chr)),as.character(end(gr.chr)), sep="_") # add coordinates as column names

#write.csv(chp_m,"/Users/rociochamorro/scECTseq/RDS/Fig2b_3mb_chp212.csv",quote = FALSE, col.names = TRUE, row.names = TRUE)


############ TR14 exonuclase 5 days

setwd("/Users/rociochamorro/scECTseq/circleCalls/filtered/tr14_exo5/")
files_tr14<- dir()
files_tr14<- files_tr14[grep(".bed",files_tr14)]
m<- data.frame(row.names = length(files_tr14))
test_bed<- NULL

for (i in 1:length(files_tr14)) {
  
  test_bed<- read.table(files_tr14[i])
  test_bed<- test_bed[test_bed$V4>=2,]
  colnames(test_bed)<- c("chr","start","end","circle_reads","total_reads")
  test_bed<- GRanges(test_bed)
  
  hits<- NULL
  ov <- NULL
  ov_length<- NULL
  ov_lengths<- NULL
  
  for (j in 1:length(gr.chr)) {
    
    hits<- findOverlaps(gr.chr[j],test_bed)
    ov <- pintersect(test_bed[subjectHits(hits)], gr.chr[j])
    ov_length<- sum(width(ov)) /sum(width(gr.chr[j]))
    
    ov_lengths<- append(ov_lengths, ov_length)
    
  }
  
  m<- cbind(m, ov_lengths)
  
}

tr14_m<- m 
samples_name_tr14<- str_remove(files_tr14, ".enriched.merged.b2bRefined.Counts.ThreshFinal.bed")
samples_name_tr14<- str_remove(samples_name_tr14, "Regions_")

samples_name<- NULL
samples_name_v<- NULL
for (i in 1:length(samples_name_tr14)) {
  samples_name<- paste(unlist(strsplit(samples_name_tr14[i], "_"))[1],unlist(strsplit(samples_name_tr14[i], "_"))[2], sep="_")
  samples_name_v<- append(samples_name_v, samples_name)
}

colnames(tr14_m)<- samples_name_v # add unique_ids as colnames
rownames(tr14_m)<- paste(as.character(seqnames(gr.chr)),as.character(start(gr.chr)),as.character(end(gr.chr)), sep="_") # add coordinates as column names

#write.csv(tr14_m,"/Users/rociochamorro/scECTseq/RDS/Fig2b_3mb_tr14.csv",quote = FALSE)


# -------------------------------------------------------------------------------------------------------------------
# Figure 2b - CHP-212
# -------------------------------------------------------------------------------------------------------------------

chp_m<- read.table("/Users/rociochamorro/scECTseq/RDS/Fig2b_3mb_chp212.csv", sep=",", header = TRUE) #we load the created df to avoid running the code again everytime
rownames(chp_m)<- chp_m[,1]
chp_m<-chp_m[,-1]
chp_m<- as.matrix(t(chp_m))

#create chr annotation for annotation's bar
split<-strsplit(as.character(colnames(chp_m)),"_")
chr_number<- NULL
chr_vector<- NULL
for (i in 1:length(split)) {
  chr_number<- split[[i]][1]
  chr_vector<- append(chr_vector,chr_number)
}

my_sample_col <- data.frame(chr_vector)
row.names(my_sample_col) <- colnames(chp_m)

cols <- makeColorRampPalette(c("white", "red",    # distances 0 to 3 colored from white to red
                               "blue"))


# Plot heatmap
chp_wo_mt<- chp_m[,-1009] # we remove mtDNA to visualize all other chromosome better
pheatmap(chp_wo_mt, color=colorRampPalette(c("white","red","darkred"))(200),cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE, annotation_col = my_sample_col)

# -------------------------------------------------------------------------------------------------------------------
# Figure 2b - TR14
# -------------------------------------------------------------------------------------------------------------------

tr14_m<- read.table("/Users/rociochamorro/scECTseq/RDS/Fig2b_3mb_tr14.csv", sep=",", header = TRUE) #we load the created df to avoid running the code again everytime
rownames(tr14_m)<- tr14_m[,1]
tr14_m<-tr14_m[,-1]
tr14_m<- as.matrix(t(tr14_m))

#create a chr annotation for annotation's bar
split<-strsplit(as.character(colnames(tr14_m)),"_")
chr_number<- NULL
chr_vector<- NULL
for (i in 1:length(split)) {
  chr_number<- split[[i]][1]
  chr_vector<- append(chr_vector,chr_number)
}

my_sample_col <- data.frame(chr_vector)
row.names(my_sample_col) <- colnames(tr14_m)

cols <- makeColorRampPalette(c("white", "red",    # distances 0 to 3 colored from white to red
                               "blue"))

# Plot heatmap 
tr14_wo_mt<- tr14_m[,-1009] # we remove mtDNA to visualize all other chromosome better
pheatmap(tr14_wo_mt, color=colorRampPalette(c("white","red","darkred"))(200),cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE, annotation_col = my_sample_col)



