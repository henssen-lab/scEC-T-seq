## Structural variants analysis

Analysis of the SVs in single-cell Circle-seq and single-cell WGS data.

### Installation

Clone repository containing the sv calling for single-cell batch.
This calls internally `lumpy` and `svaba`.

```bash
# create conda environment with snakemake>=7.0.2
conda env create -f env.yaml

# clone sv calling pipeline
git clone https://github.com/henssen-lab/sv-tools.git
cd sv-tools/pipelines/scdna-sv-workflow

```
### Run analysis

```bash
conda activate snakemake
ref=/fast/projects/NB_CircleSeq/work/CircleSeq_data/CircleSeq/reference/bwa_plasmid/hg19_plasmid.fa

snakemake --cores 16 \
        --jobs 16 \
        --use-conda \
        --rerun-incomplete --keep-going \
        --config reference=$ref root=</../fastq/DNA/plate5>
```

