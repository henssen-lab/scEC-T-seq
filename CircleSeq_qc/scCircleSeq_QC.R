# Rocio Chamorro
# scCircleSeq - QC

#Load packages
library(ggplot2)
library(ggpubr)

#Load table including the summary of QC results from all single cells/nuclei 
qc_table<- readRDS("/Users/rociochamorro/scECTseq/RDS/CircleSeq_QC_cLines_Tcell_Tumor.RDS")

# -------------------------------------------------------------------------------------------------------------------
# Supplementary Fig 1c
# -------------------------------------------------------------------------------------------------------------------

qc_table_cLines<- qc_table[which(qc_table$Origin == "Cell line" ),]
qc_table_cLines$Protocol<- factor(qc_table_cLines$Protocol, 
                                  levels = c("no cell",
                                             "+ Endo",
                                             "No digestion",
                                             "+ Exo 1day",
                                             "+ Exo 5days"))
qc_table_cLines$`Cell Line/Patient ID`<- factor(qc_table_cLines$`Cell Line/Patient ID`, 
                                   levels = c("no cell","CHP-212","TR14"))

qc_table_cLines$`mtDNA depth coverage`<- as.numeric(as.character(qc_table_cLines$`mtDNA depth coverage`)) + 0.1 #to be able to convert into log


my_comparisons <- list( c("No digestion","+ Exo 1day"),
                        c("+ Exo 1day", "+ Exo 5days"), 
                        c("+ Exo 5days", "+ Endo") )

#Plotting chrM per-base depth
ggplot(qc_table_cLines, aes(x=Protocol, y = `mtDNA depth coverage`, fill =`Cell Line/Patient ID`)) + 
  geom_violin(position=position_dodge(0.75),width=1)  +
  scale_y_continuous(trans='log10') +
  theme_classic(base_size = 14)+
  geom_jitter(size=0.25,shape=16, position=position_jitterdodge(0.25)) +
  scale_fill_manual(values=c("#999999", "#f04d53", "#55badf")) +
  xlab("") +
  ylab("chrM per-base coverage depth") +
  geom_hline(yintercept  = 10, linetype='dotted') + 
  stat_compare_means(method = "t.test", comparisons=my_comparisons)

# -------------------------------------------------------------------------------------------------------------------
# Fig 1d
# -------------------------------------------------------------------------------------------------------------------

#Filter low quality cells from no digestion and exonuclease-treated conditions
qc_table_cLines_filtered<- qc_table_cLines[-which(qc_table_cLines$`QC pass` == "No" & qc_table_cLines$Protocol == "+ Exo 5days") ,]
qc_table_cLines_filtered<- qc_table_cLines_filtered[-which(qc_table_cLines_filtered$`QC pass` == "No" & qc_table_cLines_filtered$Protocol == "+ Exo 1day" ),]
qc_table_cLines_filtered<- qc_table_cLines_filtered[-which(qc_table_cLines_filtered$`QC pass` == "No" & qc_table_cLines_filtered$Protocol == "No digestion"),]
colnames(qc_table_cLines_filtered)[4]<- "Cell Line"

#harmonize for plotting
qc_table_cLines_filtered$Protocol<- factor(qc_table_cLines_filtered$Protocol, levels = c("No digestion","+ Exo 1day","+ Exo 5days","+ Endo", "no cell"))
qc_table_cLines_filtered$`Fraction of reads mapping to mtDNA`<- as.numeric(as.character(qc_table_cLines_filtered$`Fraction of reads mapping to mtDNA`))

my_comparisons <- list( c("No digestion","+ Exo 1day"),c("+ Exo 1day", "+ Exo 5days"), c("+ Exo 5days", "+ Endo") )

#plotting
ggboxplot(qc_table_cLines_filtered[-which(qc_table_cLines_filtered$`Number of cells sorted` == 0),] , #exclude no cell control
          x="Protocol",
          y = "Fraction of reads mapping to mtDNA", color ="Cell Line", 
          palette= c("#f04d53","#55badf"),
          add ="jitter", ylim = c(0,1.75), 
          add.params = list(size=0.7, jitter=0.4), 
          xlab = "", 
          ylab = "% reads mapping to chrM " ) + stat_compare_means(method = "t.test", comparisons=my_comparisons, label.y = c(1,1.3,1.3)) + scale_y_continuous(breaks = seq(0, 1.75, by=0.25))

# -------------------------------------------------------------------------------------------------------------------
# Fig 1e
# -------------------------------------------------------------------------------------------------------------------

#harmonize for plotting
qc_table_cLines_filtered$`Fraction of reads mapping to circular DNA`<- as.numeric(as.character(qc_table_cLines_filtered$`Fraction of reads mapping to circular DNA`))

my_comparisons <- list( c("No digestion","+ Exo 1day"),c("+ Exo 1day", "+ Exo 5days") )

#plotting
ggboxplot(qc_table_cLines_filtered[-which(qc_table_cLines_filtered$`Number of cells sorted` == 0 | qc_table_cLines_filtered$Protocol == "+ Endo" ),], #exclude no cell control and endo
          x="Protocol", y = "Fraction of reads mapping to circular DNA", color ="Cell Line", 
          palette= c("#f04d53","#55badf"),
          ylim=c(0,100), add="jitter" ,
          ylab = "%reads mapping to circular regions", 
          xlab="", 
          add.params = list(size=0.7, jitter=0.4)) +  stat_compare_means(method = "t.test", comparisons=my_comparisons, label.y = c(95,95,95))

# -------------------------------------------------------------------------------------------------------------------
# Fig 1f
# -------------------------------------------------------------------------------------------------------------------

# Load the table including the number of reads mapping to circles containing endo sites - the code to generate this table can be found in "endonuclease_sites_analysis.R"

table_endo<- read.table("/Users/rociochamorro/scECTseq/RDS/EndoCircles.txt", header=TRUE)

#add metadata info to the table about cell line and protocol
position<- NULL
cline<- NULL
protocol<- NULL
cline_v<- NULL
protocol_v<- NULL

for (i in 1:length(table_endo$unique_id)) {
  position<- which(qc_table_cLines$`Well ID` == table_endo$unique_id[i])
  cline<- as.character(qc_table_cLines[position,]$`Cell Line`)
  protocol<- as.character(qc_table_cLines[position,]$Protocol)
  
  cline_v<- append(cline_v,cline)
  protocol_v<- append(protocol_v, protocol)
  }

table_endo$cell_line<- cline_v
table_endo$protocol<- protocol_v


table_endo_filtered<- table_endo[which(table_endo$unique_id %in% qc_table_cLines_filtered$`Well ID`),] # filter out low quality cells from no digestion and exo conditions

# Plotting

my_comparisons <- list( c("No digestion","+ Endo"),c("+ Exo 1day","+ Endo"), c("+ Exo 5days","+ Endo") )
table_endo_filtered$fraction_endoReads<- as.numeric(as.character(table_endo_filtered$fraction_endoReads))
table_endo_filtered$protocol<- factor(table_endo_filtered$protocol, levels = c("no cell","No digestion", "+ Exo 1day","+ Exo 5days","+ Endo"))


ggboxplot(table_endo_filtered[-which(table_endo_filtered$protocol == "no cell"),], 
          x="protocol", 
          y = "fraction_endoReads", 
          color ="cell_line", 
          palette= c("#f04d53","#55badf"),
          ylim=c(0,100), 
          add="jitter",
          ylab = "%reads circular regions endo-cutting (GTTTAAAC)", 
          xlab="", add.params = list(size=0.7, jitter=0.4)) +  stat_compare_means(method = "t.test", comparisons=my_comparisons, label.y = c(75,80,90)) 

#write.csv(table_endo, "/Users/rociochamorro/scECTseq/RDS/Fig1f_endoAnalysis.csv", row.names = FALSE, quote = FALSE)

