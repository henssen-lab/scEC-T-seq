#!/bin/bash
#SBATCH --nodes 1
#SBATCH --mem 200G
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH --time 80:00:00

#$CONDA_PREFIX/var/lib/arriba/download_references.sh hs37d5+GENCODE19
# conda activate arriba

ARRIBA_FILES=$CONDA_PREFIX/var/lib/arriba

root_arriba=~/CircleSeq/singleCellData/chimeric/scripts/arriba
dir=${root_arriba}/STAR_index_hs37d5_GENCODE19
anno=${root_arriba}/GENCODE19.gtf
ref=${root_arriba}/hs37d5.fa

root_data=tr14_plate3
data=tr14_plate3/tr14_plate3_metadata.txt

# 1. Merge fastq generate pseudobulk
set -e
set -x

tail -n+2 $data | awk '{print $18"_R1_001_val_1.fq.gz"}' | xargs gunzip -c > ${root_data}/R1.fq
gzip ${root_data}/R1.fq

tail -n+2 $data | awk '{print $18"_R2_001_val_2.fq.gz"}' | xargs gunzip -c > ${root_data}/R2.fq
gzip ${root_data}/R2.fq

# 2. Run Mapping
cd ${root_data}
STAR --runThreadN 24 --genomeDir $dir  --genomeLoad NoSharedMemory \
--readFilesIn ${root_data}/R1.fq.gz ${root_data}/R2.fq.gz --readFilesCommand zcat \
--outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 \
--alignSplicedMateMapLminOverLmate 0.2 --alignSJstitchMismatchNmax 5 -1 5 5 \
--chimOutType WithinBAM SoftClip \
--chimSegmentMin 10 \
--chimJunctionOverhangMin 10 \
--chimScoreDropMax 100 \
--chimScoreJunctionNonGTAG 0 \
--chimScoreSeparation 1 \
--chimSegmentReadGapMax 3 \
--chimMultimapNmax 50

samtools sort -@ 4 -o Aligned.out.sorted.bam Aligned.out.bam
samtools index Aligned.out.sorted.bam

# 3. Fusion genes calling
arriba -x Aligned.out.sorted.bam -o fusions.tsv -O fusions.discarded.tsv \
-a $ref -g $anno \
-F 150 -U 700 \
-i "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M AC_* NC_*" \
-b ${ARRIBA_FILES}/blacklist_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz \
-k ${ARRIBA_FILES}/known_fusions_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz \
-t ${ARRIBA_FILES}/known_fusions_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz \
-p ${ARRIBA_FILES}/protein_domains_hg19_hs37d5_GRCh37_v2.1.0.gff3

# 4. Filter
tail -n+2 fusions.tsv | awk '$13 >= 50 && $14 >= 50 && 2*($10+$12+$11)/($13+$14) > 0.3 {print $0}' > fusions_confident_tr14.tsv
cat fusions_confident.tsv | awk '{split($5, a, ":"); split($6, b, ":"); print a[1]"\t"a[2]-1"\t"a[2]"\t"b[1]"\t"b[2]-1"\t"b[2]}' > fusions_confident_tr14.bedpe
