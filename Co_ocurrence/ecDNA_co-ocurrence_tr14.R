# Rocio Chamorro
# ecDNA co-ocurrence analysis in TR14 cells

#Load required packages
library(ggplot2)
library(biomaRt)
library(stringr)
library(UpSetR)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

##################################### Find overlaps between circular DNA regions called in single cells and ecDNA regions in TR14

#Load dfs including the genomic coordinates from the three reconstructed ecDNA in TR14

tr14_cdk4<-read.table("/Users/rociochamorro/scECTseq/ecDNA reconstruction coordinates/tr14_cdk4_coordinates.txt",header = TRUE)
tr14_mycn<-read.table("/Users/rociochamorro/scECTseq/ecDNA reconstruction coordinates/tr14_mycn_coordinates.txt",header = TRUE)
tr14_mdm2<-read.table("/Users/rociochamorro/scECTseq/ecDNA reconstruction coordinates/tr14_mdm2_coordinates.txt",header = TRUE)
tr14_amp<- rbind(tr14_cdk4,tr14_mycn,tr14_mdm2)
colnames(tr14_amp)<- c("chr","Start","End","fragment_id","main_gene")

#Load circular DNA regions called in single TR14 cells - after qc filtering and exo5 days condition only

setwd("/Users/rociochamorro/scECTseq/circleCalls/filtered/tr14_exo5/"))

samples<- list.files(pattern = ".bed")
circ_tr14<- lapply(samples, read.delim, header = F)
samples_name<- str_remove(samples, ".enriched.merged.b2bRefined.Counts.ThreshFinal.bed")
samples_name<- str_remove(samples_name, "Regions_")

name<- NULL
names<- NULL
p_name<- NULL
for (i in 1:length(samples_name)) {
  name<- unlist(strsplit(samples_name[i], "_"))
  p_name<- paste(name[1],name[2], sep = "_")
  names<- append(names,p_name)
}

names(circ_tr14)<- names

#Find overlaps among tr14 ecDNA fragments and circular DNA regions called in single cells

overlap.files<- NULL
overlap.file<- NULL
test_bed<- NULL

for( i in 1: length(circ_tr14)) {
  test_bed = circ_tr14[[i]]
  colnames(test_bed) = c("chr", "Start", "End", "value1", "value2")
  test_bed$Start <- as.numeric(as.character(test_bed$Start))
  test_bed$End <- as.numeric(as.character(test_bed$End))
  
  #Find which circular regions overlap with ecDNA regions
  chr.vector = c(paste0("chr",1:22), "chrX","chrY")
  overlap.file = lapply(chr.vector, function(x){
    chrom = x
    tr14_amp %>% 
      filter(chr == chrom) %>% 
      dplyr::select(Start,End,chr,fragment_id,main_gene) %>% 
      as.data.frame()->df1
    colnames(df1)<- c("start1","end1","chr","fragment_id","main_gene")
    
    test_bed %>% 
      filter(chr==chrom) %>% 
      dplyr::select(Start,End) %>% 
      as.data.frame()->df2
    colnames(df2)<- c("start2","end2")
    
    ir1 = with(df1, IRanges(start1, end1))
    ir2 = with(df2, IRanges(start2, end2))
    df1 %>% 
      mutate(overlap = countOverlaps(ir1, ir2))->df1
    return(df1)
  })
  
  overlap.file<- do.call(rbind,overlap.file)
  overlap.file<- overlap.file[which(overlap.file$overlap >=1),] # 
  overlap.file<- list(overlap.file)
  
  overlap.files<- c(overlap.files, overlap.file)
  
}

ov.tr14<- overlap.files
names(ov.tr14)<- names(circ_tr14)

######################## Create recurrence table - binomial 0,1 table, in which 0 means the ecDNA region is absent or 1 if present

r_tr14<- NULL
sc<- NULL

for (i in 1:length(ov.tr14)) {
  sc<- ov.tr14[[i]]
  
  circ<- NULL
  circles_sc<- NULL
  for (j in 1:length(tr14_amp$fragment_id)) {
    circ<- length(which(sc$fragment_id == tr14_amp$fragment_id[j]))
    circles_sc<- append(circles_sc, circ)
  }
  r_tr14<- cbind(r_tr14, circles_sc)
}

colnames(r_tr14)<- names(ov.tr14)
r_tr14<- data.frame(r_tr14)
r_tr14$recurrence<- rowSums(r_tr14)
r_tr14<- cbind(tr14_amp, r_tr14)

r_tr14<- r_tr14[order(-r_tr14$recurrence),]
r_tr14_exo5<- r_tr14

colnames(r_tr14_exo5)<- c("chr","start","end", "fragment_id","main_gene",names, "recurrence")

#write.csv(r_tr14_exo5, "/Users/rociochamorro/scECTseq/RDS/ecDNArecurrence_tr14.csv")

#### Create table for upset plot 
tr14_fragRecurrence<- r_tr14_exo5
tr14_fragRecurrence<- as.data.frame(t(tr14_fragRecurrence))
colnames(tr14_fragRecurrence)<- tr14_fragRecurrence[4,]
tr14_fragRecurrence<- tr14_fragRecurrence[-c(1,2,3,4,5),] #remove coordinates and fragment ids
tr14_fragRecurrence<- tr14_fragRecurrence[-length(tr14_fragRecurrence[,1]),] #remove recurrence row

#### Create table for upset plot with recurrence by ecDNA amplicon 

# As MYCN and CDK4 ecDNAs are multi-fragmented, if any of their fragments overlap with circular regions in one cell, we considered the ecDNA to be present

# CDK4 and MYCN ecDNAs share some segments, these are excluded from the recurrence analysis.
tr14_fragRecurrence_unique<- tr14_fragRecurrence[,-grep("CDK4_1", colnames(tr14_fragRecurrence))] 
tr14_fragRecurrence_unique<- tr14_fragRecurrence_unique[,-grep("MYCN_4", colnames(tr14_fragRecurrence_unique))] 

#CDK4 ecDNA

suma<- NULL
sample<- NULL
samples<- NULL
cdk4_amp<- NULL
cdk4_amps<- NULL
for(i in 1:length(tr14_fragRecurrence_unique[,1])) {
  suma<- sum(as.numeric(tr14_fragRecurrence_unique[i,grep("CDK4",colnames(tr14_fragRecurrence_unique))]))
  sample<- as.character(rownames(tr14_fragRecurrence_unique)[i])
  samples<- append(samples, sample)
  
  if (suma >=1) {
    cdk4_amp<- 1
  }
  else{
    cdk4_amp<- 0
  }
  cdk4_amps<- append(cdk4_amps, cdk4_amp)
}

#MYCN ecDNA

suma<- NULL
sample<- NULL
samples<- NULL
mycn_amp<- NULL
mycn_amps<- NULL

for(i in 1:length(tr14_fragRecurrence_unique[,1])) {
  suma<- sum(as.numeric(tr14_fragRecurrence_unique[i,grep("MYCN",colnames(tr14_fragRecurrence_unique))]))
  sample<- as.character(rownames(tr14_fragRecurrence_unique)[i])
  samples<- append(samples, sample)
  
  if (suma >=1) {
    mycn_amp<- 1
  }
  else{
    mycn_amp<- 0
  }
  mycn_amps<- append(mycn_amps, mycn_amp)
}

#MDM2 (simple ecDNA with one fragment)
mdm2_amps<- as.numeric(tr14_fragRecurrence_unique$MDM2_1)

tr14_recurrence_byamp<- as.data.frame(cbind(cdk4_amps,mycn_amps,mdm2_amps))
colnames(tr14_recurrence_byamp)<- c("CDK4","MYCN","MDM2")
rownames(tr14_recurrence_byamp)<- samples

#write.csv(tr14_recurrence_byamp, "/Users/rociochamorro/scECTseq/RDS/Fig4b_sourceData.csv")

########################################################################################################################

# Fig. 4b - Upset plot

########################################################################################################################

upset(tr14_recurrence_byamp, sets = colnames(tr14_recurrence_byamp), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")


