#!/bin/bash
#SBATCH --nodes 1
#SBATCH --mem 200G
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH --time 80:00:00

dir="STAR_index_hs37d5_GENCODE19/"
ref="hs37d5.fa"
anno="GENCODE19.gtf"

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $dir --genomeFastaFiles $ref --sjdbGTFfile $anno --sjdbOverhang 74

