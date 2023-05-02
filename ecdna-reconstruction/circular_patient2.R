library(BioCircos)
library(pafr)
library(ggplot2)
library(gridExtra)

# library(Repitools)
library(StructuralVariantAnnotation)
library(gTrack)
library(gGnome)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(dplyr)
library(parallel)
library(plyranges)
library(diffloop)


standardchrs <- c("chr1", "chr2" , "chr3",  "chr4" , "chr5"  ,"chr6" , "chr7" , "chr8" , "chr9" , "chr10" ,"chr11", "chr12", "chr13", "chr14" ,"chr15", "chr16", "chr17", "chr18", "chr19", "chr20"
                  ,"chr21" ,"chr22", "chrX" )

make_bigwig_gTrack = function(fname, window, bins, name = "NoName", y0=0, y1=1000){
  this_track = read_bigwig(fname)
  this_track = diffloop::addchr(this_track)
  
  seqlevels(this_track , pruning.mode="coarse") = standardchrs
  seqlevels(window , pruning.mode="coarse") = standardchrs
  seqlevels(bins , pruning.mode="coarse") = standardchrs

  this_track = binnedAverage(bins, coverage(this_track, weight=this_track$score), "score")
  this_gt = gTrack(this_track, y.field="score", bar=TRUE, name = name, y0=y0, y1=y1) #, 
  return(this_gt)
}

draw_amplicon = function(file, covfile, window, bins, cov, gencode){

  sv_curated = jJ(file, geno = FALSE, chr.convert = FALSE)
  gg_curated = gG(juncs = sv_curated)
  gt.cov = make_bigwig_gTrack(covfile, window, bins, name = "coverage", y1=cov)
  gg_curated$set(name = 'reconstruction')
  
  print(gg_curated$nodes)
  
  p<-plot(c(gencode, gt.cov, gg_curated$gt), window, links = gg_curated$grl)
  return(p)
}

draw_amplicon_manual = function(nodesfile, edgesfile, covfile, window, bins, cov, gencode){
  
  reconstruct_coordinates <- read.table(nodesfile, header=FALSE, sep="\t")
  colnames(reconstruct_coordinates) <- c("chr", "start","end", "strand",	"target",	"coverage",	"structure","fragment")
  nodes<-makeGRangesFromDataFrame(reconstruct_coordinates,
                                  keep.extra.columns=FALSE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field=c("seqnames", "seqname",
                                                   "chromosome", "chrom",
                                                   "chr", "chromosome_name",
                                                   "seqid"),
                                  start.field="start",
                                  end.field=c("end", "stop"),
                                  strand.field="strand",
                                  starts.in.df.are.0based=FALSE)
  
  edges<-read.table(edgesfile, header=FALSE, sep="\t")
  colnames(edges)<- c("n1","n2","n1.side","n2.side")
  gg_curated <- gG(nodes = nodes, edges = edges)
  
  gt.cov = make_bigwig_gTrack(covfile, window, bins, name = "coverage", y1=cov)
  gg_curated$set(name = 'reconstruction')
  
  p<-plot(c(gencode, gt.cov, gg_curated$gt), window, links = gg_curated$grl, '3')
  return(p)
}

draw_circle = function(file, sample){
  
  df <- read.csv2(file, header=FALSE, sep="\t", row.names = NULL, skip = 1)
  colnames(df) <- c("chr", "start","stop","strand","target","coverage","structure", "fragment")
  df$size <- df$stop - df$start
  
  genomeChr = df$fragment
  lengthChr = df$size
  names(lengthChr) <- genomeChr
  chrLabels <- df$chr
  
  tracklist = BioCircosTextTrack('structure',sample, size = "2em", opacity = 0.5,  x = -0.3, y = 0)
  
  p<-BioCircos(tracklist, genome = lengthChr,
            genomeFillColor = c("tomato2", "darkblue"),
            genomeTicksDisplay = FALSE,
            genomeLabelTextSize = "8pt",
            chrPad = 0.05)
  return(p)
}

draw_dotplot = function(file,sample){
  
  print(file)
  ali <- read_paf(file)
  p<-dotplot(ali, label_seqs=TRUE, alignment_colour="blue", xlab=sample, ylab="reference") + theme_bw()
  return(p)
}