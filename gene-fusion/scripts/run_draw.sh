#!/bin/bash
#SBATCH --nodes 1
#SBATCH --mem 10G
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --time 8:00:00

# conda activate arriba
ARRIBA_FILES=$CONDA_PREFIX/var/lib/arriba

# arriba script found under $ARRIBA_FILES
draw_fusions.R --annotation=GENCODE19.gtf \
--fusions=chp212_plate5a6/fusions_candidates.tsv \
--output=chp212_plate5a6/CHP212_fusions_candidates.pdf \
--cytobands=$ARRIBA_FILES/cytobands_hg19_hs37d5_GRCh37_v2.1.0.tsv \
--alignments=chp212_plate5a6/Aligned.out.sorted.bam \
--proteinDomains=$ARRIBA_FILES/protein_domains_hg19_hs37d5_GRCh37_v2.1.0.gff3

draw_fusions.R --annotation=GENCODE19.gtf \
--fusions=tr14_plate3/fusions_candidates.tsv \
--output=tr14_plate3/TR14_fusions_candidates.pdf \
--cytobands=$ARRIBA_FILES/cytobands_hg19_hs37d5_GRCh37_v2.1.0.tsv \
--alignments=tr14_plate3/Aligned.out.sorted.bam \
--proteinDomains=$ARRIBA_FILES/protein_domains_hg19_hs37d5_GRCh37_v2.1.0.gff3

