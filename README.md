# 2023-bac-phage-rnaseq

This repository contains code to analyze an RNA-seq sample that contains *Escherichia coli* strain B (ATCC Strain 11303), *Bacillus subtilis* strain 168 (ATCC Strain 27370) and their respective infecting phages, T4 phage and SPO1 phage.
It contains notebooks that generate the figures in the pub [DOI 10.57844/arcadia-j03a-mz33](https://research.arcadiascience.com/pub/negative-data-phage-rna-isolation).

## Getting started

This repository uses snakemake to run the pipeline and conda to manage software environments and installations.
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).
After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.
We installed Miniconda3 version `py39_4.12.0` and mamba version `0.15.3`.

```
mamba env create -n env --file environment.yml
conda activate env
```

To start the pipeline, run:
```
snakemake --use-conda -j 1
```

## Pipeline description

The pipeline downloads reference files (genome and annotation files) for the genomes listed above.
It combines the references for mapping so that all genomes are simultaneously mapped against, which should reduce spurious mapping between genomes.
The pipeline uses [BWA mem](https://bio-bwa.sourceforge.net/bwa.shtml) to align reads against the reference genomes, which produces a BAM file that can be viewed in a genome viewer (e.g. IGV).
BWA mem is not a splice-aware aligner, but given that bacteria and phage do not typically splice their transcripts, this should be fine.
The BAM file is post-processed to determine mapping depth and number of reads mapping to each reference.
The pipeline also downloads the GFF and GTF annotation files that accompany these genomes.
