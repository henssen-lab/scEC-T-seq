######### RUN BLAST TO ANALYZE MICROHOMOLOGY AROUND CIRCLE JUNCTIONS ######### 

##Libraries
library (devtools)
library (tidyverse)
library(tidyr)
library (BSgenome.Hsapiens.UCSC.hg19)

#Initialize DFs
FINAL_RESULTS_BLAST = NULL

#Arguments
args <- commandArgs(trailingOnly=TRUE)

input_file <- args[1]


##### Read file

#Open input file
circ_file <- read.table(input_file, sep="\t",stringsAsFactors=FALSE)
names(circ_file) <- c("chr", "start", "end", ".", "..")

#Get a simpler table selecting circles > 100bp (we extend 50bp the sequence around the junction)
circ_file$length <- abs(circ_file$end - circ_file$start)
circ_file_ALL <- unique(circ_file[circ_file$length>=100,c(1,2,3)]) #chr start end
circ_file_ALL <- unique(circ_file[circ_file$start>=50,c(1,2,3)]) #chr start end

#Create flag for each mutation chr_start_end
circ_file_ALL$flag <- paste(circ_file_ALL$chr,"_",circ_file_ALL$start,"_",circ_file_ALL$end,sep="")


###### Generate the DF with all sequences

#Loop through all the circles
for (c in 1:dim(circ_file_ALL)[1]){
        #We extend the sequence +/-50bp around the bkp. 100bp inside the circle and 100bp in the linear genome.
        circ_file_ALL$seq_start[c] <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19, circ_file_ALL[c,1], as.numeric(circ_file_ALL[c,2])-50, as.numeric(circ_file_ALL[c,2])+50))
        circ_file_ALL$seq_end[c] <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19, circ_file_ALL[c,1], as.numeric(circ_file_ALL[c,3])-50, as.numeric(circ_file_ALL[c,3])+50))
}

###### Create the database for BLAST with the end sequences

input_DB <- circ_file_ALL[,c(4,6)]
#Change tab to fasta
sink(paste(input_file,".DB_BLAST_endseq.fa",sep=""))

for (i in 1:dim(input_DB)[1]){
        name <- paste(">",input_DB[i,1])
        sequence <- paste(input_DB[i,2])
        cat(name,sep="\n")
        cat(sequence,sep="\n")
}
  
#this is sink all the console output to the file 
sink()


###### Run DB generator BLAST

cmd_db <- paste("/usr/local/ncbi/blast/bin/makeblastdb -in", paste(input_file,".DB_BLAST_endseq.fa",sep=""), "-dbtype nucl -parse_seqids")
system(cmd_db,intern = FALSE)


###### Run BLAST

for (C in 1:dim(circ_file_ALL)[1]){
        #Create seqid file
        seqid <- circ_file_ALL$flag[C]
        write.table(seqid, "index_seqid.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
        #Create query
        query_input_start <- circ_file_ALL[C,c(4,5)]
        #Change tab to fasta and create files for query
        sink(paste(input_file,".QUERY_BLAST_startseq.fa",sep=""))

        for (i in 1:dim(query_input_start)[1]){
                name <- paste(">",query_input_start[i,1])
                sequence <- paste(query_input_start[i,2])
                cat(name,sep="\n")
                cat(sequence,sep="\n")
        }
        sink()
        #RUN BLAST
        #Create arguments
        arg_query <- paste(input_file, ".QUERY_BLAST_startseq.fa",sep="")
        arg_DB <- paste(input_file, ".DB_BLAST_endseq.fa",sep="")
        #Launch BLAST
        cmd_blast <- paste("/usr/local/ncbi/blast/bin/blastn -task megablast -db", arg_DB, "-query", arg_query, "-seqidlist index_seqid.txt -word_size=4 -evalue=1 -outfmt \"6 qseqid length evalue\" -subject_besthit -reward=1 -penalty=-2")
        result_blast <- data.frame(system(cmd_blast,intern = TRUE))
        if (dim(result_blast)[1]>0){
                result_blast_COLS <- separate(result_blast,col=1,into=c("circle_id", "HOM_length", "evalue"),sep="\t")[1,] #I keep the result with lower evalue
                FINAL_RESULTS_BLAST <- rbind(FINAL_RESULTS_BLAST, result_blast_COLS)
        }
}

###### Write results BLAST
write.table(FINAL_RESULTS_BLAST,paste(input_file,".BLAST_RESULTS_HOMOLOGY_min4bp_100bparoundbkp.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)

##### Remove queries and databases
system(paste("rm ",input_file,".DB_BLAST_endseq.fa*",sep=""))
system(paste("rm ",input_file,".QUERY_BLAST_startseq.fa",sep=""))
system("rm index_seqid.txt")

