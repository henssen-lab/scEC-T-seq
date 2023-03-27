# Line plots ecDNA regions in cell lines
# Rocio Chamorro

# Load required packages
library(rtracklayer)
library(tidyverse)
library(GenomicRanges)



# -------------------------------------------------------------------------------------------------------------------
# Figure 3a
# -------------------------------------------------------------------------------------------------------------------

# Load read coverage bigwig files from CHP-212 cells
# Only exo5 condidtion and excluding low quality cells that have been filtered out for further analyses

setwd("/Users/rociochamorro/Desktop/Bioinformatics/BIH/work/projects/bigwig/chp")
cells<- list.files(pattern=".bw")

# Filter 50bp intervals from the bw file to keep only regions flanking predicted ecDNA in CHP-212 bulk

bw<- NULL
bw_seq = NULL
bw_start = NULL
bw_end = NULL
bw_score = NULL
bw_bed = NULL
int_region1<- NULL
int_region2<- NULL
int_regions1<- NULL
int_regions2<- NULL
int_regions<- NULL

for (i in 1:length(cells)) {
  bw<- import(cells[i])
  bw_seq = as.character(seqnames(bw))
  bw_start = start(bw)
  bw_end = end(bw)
  bw_score = score(bw)
  bw_bed = data.frame(bw_seq, bw_start, bw_end, bw_score)
  colnames(bw_bed) = c("chr", "start", "end", "score")
  bw_bed$cell = paste("bw",i, sep = "")
  
  bw_bed %>% 
    filter(chr == "chr2") %>% 
    filter(start > 11500000 & end < 13300000) ->int_region1 # Select 1.8 Mb region flanking ecDNA region 1 in chr2
  int_regions1<- rbind(int_regions1,int_region1)
  
  bw_bed %>% 
    filter(chr == "chr2") %>% 
    filter(start > 15000000 & end < 16300000) ->int_region2 # Select 1.8 Mb region flanking ecDNA region 2 in chr2
  int_regions2<- rbind(int_regions2,int_region2)
  
}
int_regions1$region<- "r1"
int_regions2$region<- "r2"

int_regions<- rbind(int_regions1,int_regions2)

colnames(int_regions) = c("chr", "start", "end", "score", "cell","region")

#saveRDS(int_regions, "/Users/rociochamorro/Documents/Rocio/AG Henssen/Projects/SCeccDNA & RNA seq/scSeq analyses/Manuscript_analyses/ecDNA heatmaps/RDS/all_regions_chp212.RDS")

# Load the previously created df to avoid repeating the previous code 
mycn_chp<- readRDS("/Users/rociochamorro/Documents/Rocio/AG Henssen/Projects/SCeccDNA & RNA seq/scSeq analyses/Manuscript_analyses/ecDNA heatmaps/RDS/all_regions_chp212.RDS")

pdf("/Users/rociochamorro/Documents/Rocio/AG Henssen/Projects/SCeccDNA & RNA seq/scSeq analyses/Manuscript_analyses/ecDNA heatmaps/plots/all_regions_chp212_mycn.pdf", width=10, height=5)
mycn_chp %>% 
  mutate(mean = (start + end)/ 2) %>% 
  filter(score > 0) %>% 
  ggplot(aes(x=cell,y=mean))+
  geom_linerange(aes(ymin =start, ymax = end, colour = "#F8766D"))+
  coord_flip()+
  facet_wrap( ~ region, scales = "free_x",ncol = 2) +
  theme_classic()
dev.off()


# -------------------------------------------------------------------------------------------------------------------
# Figure 3b
# -------------------------------------------------------------------------------------------------------------------

# Load read coverage bigwig files from TR14 cells
# Only exo5 condidtion and excluding low quality cells that have been filtered out for further analyses

setwd("/Users/rociochamorro/Desktop/Bioinformatics/BIH/work/projects/bigwig/kk04/tr14_exo")
cells<- list.files(pattern=".bw")

# Filter 50bp intervals from the bw file to keep only regions flanking predicted ecDNA amplicons in TR14 bulk
# In this case, we will plot filtered chr2 and chr12 regions of TR14 ecDNAs separately. 

###############chr2
bw<- NULL
bw_seq = NULL
bw_start = NULL
bw_end = NULL
bw_score = NULL
bw_bed = NULL
int_region1<- NULL
int_region2<- NULL
int_regions1<- NULL
int_regions2<- NULL
int_region3<- NULL
int_regions3<- NULL
int_regions<- NULL

for (i in 1:length(cells)) {
  bw<- import(cells[i])
  bw_seq = as.character(seqnames(bw))
  bw_start = start(bw)
  bw_end = end(bw)
  bw_score = score(bw)
  bw_bed = data.frame(bw_seq, bw_start, bw_end, bw_score)
  colnames(bw_bed) = c("chr", "start", "end", "score")
  bw_bed$cell = paste("bw",i, sep = "")
  
  bw_bed %>% 
    filter(chr == "chr2") %>% 
    filter(start > 2000000 & end < 4000000) ->int_region1
  int_regions1<- rbind(int_regions1,int_region1)
  
  bw_bed %>% 
    filter(chr == "chr2") %>% 
    filter(start > 13600000 & end < 15600000) ->int_region2
  int_regions2<- rbind(int_regions2,int_region2)
  
  bw_bed %>% 
    filter(chr == "chr2") %>% 
    filter(start > 15800000 & end < 17800000) ->int_region3
  int_regions3<- rbind(int_regions3,int_region3)
  
}
int_regions1$region<- "r1"
int_regions2$region<- "r2"
int_regions3$region<- "r3"


int_regions<- rbind(int_regions1,int_regions2)
int_regions<- rbind(int_regions, int_regions3)

colnames(int_regions) = c("chr", "start", "end", "score", "cell","region")

#saveRDS(int_regions, "/Users/rociochamorro/Documents/Rocio/AG Henssen/Projects/SCeccDNA & RNA seq/scSeq analyses/Manuscript_analyses/ecDNA heatmaps/RDS/all_regions_tr14_chr2.RDS")

tr14_chr2<- readRDS("/Users/rociochamorro/Documents/Rocio/AG Henssen/Projects/SCeccDNA & RNA seq/scSeq analyses/Manuscript_analyses/ecDNA heatmaps/RDS/all_regions_tr14_chr2.RDS")

pdf("/Users/rociochamorro/Documents/Rocio/AG Henssen/Projects/SCeccDNA & RNA seq/scSeq analyses/Manuscript_analyses/ecDNA heatmaps/plots/all_regions_tr14_chr2.pdf", width=10, height=3)
tr14_chr2 %>% 
  mutate(mean = (start + end)/ 2) %>% 
  filter(score > 0) %>% 
  ggplot(aes(x=cell,y=mean, color=mean))+
  scale_color_manual(values="#00BFC4")+
  geom_linerange(aes(ymin =start, ymax = end, color="#00BFC4"))+
  coord_flip()+
  facet_wrap( ~ region, scales = "free_x",ncol = 4) +
  theme_classic()
dev.off()

######chr12

bw<- NULL
bw_seq = NULL
bw_start = NULL
bw_end = NULL
bw_score = NULL
bw_bed = NULL
int_region1<- NULL
int_region2<- NULL
int_regions1<- NULL
int_regions2<- NULL
int_region3<- NULL
int_regions3<- NULL
int_regions<- NULL

for (i in 1:length(cells)) {
  bw<- import(cells[i])
  bw_seq = as.character(seqnames(bw))
  bw_start = start(bw)
  bw_end = end(bw)
  bw_score = score(bw)
  bw_bed = data.frame(bw_seq, bw_start, bw_end, bw_score)
  colnames(bw_bed) = c("chr", "start", "end", "score")
  bw_bed$cell = paste("bw",i, sep = "")
  
  bw_bed %>% 
    filter(chr == "chr12") %>% 
    filter(start > 57000000 & end < 59000000) ->int_region1
  int_regions1<- rbind(int_regions1,int_region1)
  
  bw_bed %>% 
    filter(chr == "chr12") %>% 
    filter(start > 68800000 & end < 70800000) ->int_region2
  int_regions2<- rbind(int_regions2,int_region2)
  
  bw_bed %>% 
    filter(chr == "chr12") %>% 
    filter(start > 110500000 & end < 112500000) ->int_region3
  int_regions3<- rbind(int_regions3,int_region3)
  
  
}
int_regions1$region<- "r1"
int_regions2$region<- "r2"
int_regions3$region<- "r3"

int_regions<- rbind(int_regions1,int_regions2)
int_regions<- rbind(int_regions, int_regions3)


colnames(int_regions) = c("chr", "start", "end", "score", "cell","region")

#saveRDS(int_regions, "/Users/rociochamorro/Documents/Rocio/AG Henssen/Projects/SCeccDNA & RNA seq/scSeq analyses/Manuscript_analyses/ecDNA heatmaps/RDS/all_regions_tr14_chr12.RDS")

tr14_chr12<- readRDS("/Users/rociochamorro/Documents/Rocio/AG Henssen/Projects/SCeccDNA & RNA seq/scSeq analyses/Manuscript_analyses/ecDNA heatmaps/RDS/all_regions_tr14_chr12.RDS")

pdf("/Users/rociochamorro/Documents/Rocio/AG Henssen/Projects/SCeccDNA & RNA seq/scSeq analyses/Manuscript_analyses/ecDNA heatmaps/plots/all_regions_tr14_chr12.pdf", width=10, height=3)
tr14_chr12 %>% 
  mutate(mean = (start + end)/ 2) %>% 
  filter(score > 0) %>% 
  ggplot(aes(x=cell,y=mean, color=mean))+
  scale_color_manual(values="#00BFC4")+
  geom_linerange(aes(ymin =start, ymax = end, color="#00BFC4"))+
  coord_flip()+
  facet_wrap( ~ region, scales = "free_x",ncol = 4) +
  theme_classic() 
dev.off()


