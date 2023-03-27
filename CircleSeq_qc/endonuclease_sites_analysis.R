## Endo 
## PmeI cutting site "GTTTAAAC"

library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(GenomicRanges)
library(stringr)

#Identify circle regions in single cells which contain the PmeI cutting site "GTTTAAAC" 

setwd("/Users/rociochamorro/scECTseq/circleCalls/unfiltered/")

samples<- list.files(pattern = ".bed") #bed files including the coordinates of circles called in each single cell
circ<- lapply(samples, read.table)
samples_name<- samples

file<-NULL
name<- NULL
seq<-NULL
endo_reads<- NULL
total_reads<- NULL
names<- NULL
endo<- NULL
total<- NULL
endo_c<- NULL
total_c<- NULL
endo_circles<- NULL
total_circles<- NULL


for (i in 1:length(circ)) {
  file<- circ[[i]]
  name<- samples_name[[i]]
  colnames(file)<- c("chrom","start","end","circle_reads","total_reads")
  file<- file[which(file$circle_reads>=2),] # Filtering step --> Include only circles which have at least 2 circleEdge supporting reads 
  file[grep("chrM",file$chrom),3]<- 16569
  file[grep("chrM",file$chrom),2]<- 1
  file[which(file$start == 0),2]<- 1
  file<- file[grep("chr",file$chrom),]
  seq<- Views(Hsapiens, GRanges(file)) #Extract the sequence of each putative circle 
  
  s<- NULL
  cut<- NULL
  cuts<- NULL
  
  for (j in 1:length(file[,1])) {
    s<- as.character(seq[[j]])
    cut<- sum(unlist(str_locate_all(s, "GTTTAAAC"))) != 0 # Identify which of the putative circles contain the cutting sequence
    cuts<- append(cuts,cut)
  }
  
  file$endo_site<- cuts
  endo_reads<- sum(file[which(file$endo_site == TRUE),5]) #count the reads mapping to circles containing the endo site
  endo_circles<- length(file[which(file$endo_site == TRUE),5]) #count how many circles contain cutting sites
  total_circles<- length(file[,1])
  total_reads<- sum(file$total_reads)
  
  names<- append(names,name)
  endo<- append(endo, endo_reads)
  total<- append(total, total_reads)
  endo_c<- append(endo_c, endo_circles)
  total_c<- append(total_c,total_circles)
  
}

table_endo<- as.data.frame(cbind(names, endo))
table_endo<- cbind(table_endo, endo_c)
table_endo<- cbind(table_endo, total)
table_endo<- cbind(table_endo, total_c)
table_endo$p<- as.numeric(as.character(table_endo$endo))*100/as.numeric(as.character(table_endo$total)) #fraction of reads mapping to circles containing endo sites

# We modify the column "names" to include the unique identifiers (plate_well) for each single cell

sample_name<- NULL
sample_name_v<- NULL
for (i in 1:length(table_endo$names)) {
sample_name<- paste(unlist(strsplit(table_endo$names[i], "_"))[2],unlist(strsplit(table_endo$names[i], "_"))[3], sep="_")
sample_name_v<- append(sample_name_v, sample_name)
}

table_endo$names<- sample_name_v


## add 0 reads files --> some bed files are empty, meaning we did not identify any circles, we still need to include them in the final table with all columns = 0
setwd("/Users/rociochamorro/scECTseq/circleCalls/unfiltered/0_reads/")
samples_0<- list.files(pattern = ".bed")

sample0_name<- NULL
sample0_name_v<- NULL
for (i in 1:length(table_endo$names)) {
  sample0_name<- paste(unlist(strsplit(samples_0[i], "_"))[2],unlist(strsplit(samples_0[i], "_"))[3], sep="_")
  sample0_name_v<- append(sample0_name_v, sample0_name)
}

samples0_name<- sample0_name_v
samples0_table<- cbind(samples0_name,0,0,0,0,0) #no circles are identified in these files, therefore, all values should be 0
colnames(samples0_table)<- colnames(table_endo)

table_endo<- rbind(table_endo,samples0_table) #combine both tables
colnames(table_endo)<- c("unique_id","nreads_mapping_endoCircles", "ncircles_with_endoSite", "nreads_total","ncircles_total", "fraction_endoReads")

#write.table(table_endo, "/Users/rociochamorro/Desktop/endo_all_final.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

endo<- read.table("/Users/rociochamorro/Documents/Rocio/AG Henssen/Projects/SCeccDNA & RNA seq/scSeq analyses/Manuscript_analyses/QC circle seq/tables final/endo_all_final.txt", header=TRUE)

plate<- NULL
well<- NULL
endo_name<- NULL
endo_name_v<- NULL
for (i in 1:length(endo$names)) {
  plate<- unlist(strsplit(endo$names[i], "_"))[1]
  well<- unlist(strsplit(endo$names[i], "_"))[2]
  if( plate == "KK03") {
    plate<- "plate2"
  }
  if (plate == "KK04"){
    plate<- "plate3"
  }
  else{
    plate<- plate
  }
  endo_name<- paste(plate,well, sep="_")
  endo_name_v<- append(endo_name_v, endo_name)
}

endo$names<- endo_name_v 
colnames(endo)<- c("unique_id","nreads_mapping_endoCircles", "ncircles_with_endoSite", "nreads_total","ncircles_total", "fraction_endoReads")

#write.table(endo, "/Users/rociochamorro/scECTseq/RDS/EndoCircles.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

endo<- endo[which(endo$unique_id %in% qc_table_cLines$`Well ID`),]
#write.table(endo, "/Users/rociochamorro/scECTseq/RDS/EndoCircles.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
#saveRDS(endo, "/Users/rociochamorro/scECTseq/RDS/EndoCircles.RDS")
