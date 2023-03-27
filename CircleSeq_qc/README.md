## Quality Control - Single cell Circle-seq data 

scCircle-seq data quality control for CHP212 and TR14 cell lines (Fig. 1 and Suppl. Fig 1). Values for plotting are available as supplementary data (Suppl. Table 1) in Chamorro et al, 2023.

#### Identification of circular DNA regions containing endo-cutting sites

```bash
Rscript endonuclease_sites_analysis.R

```

#### Plot chrM per-base depth, fraction of reads mapping to mtDNA, fraction of reads mapping to all circular DNA regions and fraction of reads mapping to circular DNA regions containing endo-cutting sites.

```bash
Rscript scCircleSeq_QC.R

```
