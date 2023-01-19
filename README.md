# 2023-bac-phage-rnaseq

This repository contains code to analyze an RNA-seq sample that contains *Escherichia coli* strain B (ATCC Strain 11303), *Bacillus subtilis* strain 168 (ATCC Strain 27370) and their respective infecting phages, T4 phage and SPO1 phage.

## Getting started

This repository uses snakemake to run the pipeline and conda to manage software environments and installations.
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).
After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.
We installaed Miniconda3 version `py39_4.12.0` and mamba version `0.15.3`.

```
mamba env create -n env --file environment.yml
conda activate env
```

To start the pipeline, run:
```
snakemake --use-conda -j 1
```

## Pipeline description

The pipeline downloads reference files (genome, CDS from genomic, and annotation files) for the genomes listed above.
It combines the references for mapping so that all genomes/transcripts are simultaneously mapped against, which should reduce spurious mapping between genomes.
It uses [salmon](https://salmon.readthedocs.io/en/latest/salmon.html) to "quantify" the number of reads that pseudomap against transcripts (`*CDS_from_genomic.fna.gz`).
This produces a count file, `outputs/salmon_quant/J1_quant/quant.sf`, which contains the name of each transcript and count values.
The pipeline also uses [BWA mem](https://bio-bwa.sourceforge.net/bwa.shtml) to align reads against the reference genomes, which produces a BAM file that can be viewed in a genome viewer (e.g. IGV).
BWA mem is not a splice-aware aligner, but given that bacteria and phage do not typically splice their transcripts, this should be fine.

The pipeline also downloads the GFF and GTF annotation files that accompany these genomes.
These can be joined to the salmon `quant.sf` file or layered as an annotation in genome viewers.
