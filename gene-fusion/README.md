## Gene fusion analysis code

We describe below the gene fusion analysis from pseudobulk scEC&T-seq Illumina short-read data for TR14 and CHP212 neuroblastoma celllines.

### Installation

To create the conda environment, do:

```bash
conda env create -f env.yaml 
export CONDA_PREFIX=<where your conda>
```

### Run analysis

Download reference genome and annotation:

```bash
bash scripts/download.sh
```

Run STAR index:
```bash
bash scripts/run_star_index.sh
```

Run STAR alignment and gene fusion calling using Arriba:

```bash
bash scripts/run_chp212.sh
bash scripts/run_tr14.sh
```

Filter fusions in the proximity of the amplicon boundaries:

```bash
bash scripts/run_chp212_filter_fusions.sh
bash scripts/run_tr14_filter_fusions.sh
```

Draw fusions:

```bash
bash run_draw.sh
```
