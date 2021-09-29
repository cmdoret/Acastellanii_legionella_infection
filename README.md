# Genomic changes during infection of amoeba by _L. pneumophila_

## Background

This repository contains the analysis of _Acanthamoeba castellanii_ infection by Legionella pneumophila.
We investigate how the host genome is remodelled during infection by an intracellular bacterium. To investigate these changes, we use Hi-C and RNAseq to measure both 3D changes in chromatin and gene expression changes. We use two biological replicates of uninfected _A. castellanii_ (strain C3) and two infected replicates at 5h post infection.

## Dependencies

The pipeline is written using snakemake and has the following dependencies:

* python >= 3.7
* conda >= 4.8
* snakemake >= 5.5

Each rule is encapsulated in a conda environment where its dependencies are automatically pulled from the internet.
Fastq files containing the Hi-C and RNA-seq reads are also downloaded automatically from SRA.

## Installation

In order to rerun the whole processing pipeline, you need to:

* download the input archive and extract it into `data/input`

The processing pipeline is written using snakemake and relies on conda to manage dependencies. You need to have a working conda installation on your machine and install snakemake via pip or conda.

## Usage

You can then run the pipeline with:

```sh
snakemake -j6 --use-conda
```
And the pipeline should fetch required packages as it runs.

## Configuration

Some metadata files are provided with the pipeline to help understand the design and modify parameters. The following files may be of interest:

* `samples.tsv`: Samples used in analyses and associated informations
* `units.tsv`: sequencing libraries used in the pipeline, file paths for the reads and metadata
* `config.yaml`: path to key files and general parameters to control results of the pipeline.
* `cluster_slurm.json`: Cluster resource requirements in the event that the pipeline is run on a HPC with the SLURM scheduler. In that case, the following command should be used to run the pipeline instead:
  + snakemake --rerun-incomplete --use-conda --cluster-config cluster_slurm.json --cluster "sbatch -n {cluster.ntasks} -c {cluster.ncpus} --mem {cluster.mem} --qos {cluster.queue}" --jobs 30

## Pipeline

The pipeline is subdivided into submodules relating to the processing and downstream analysis of Hi-C and RNAseq data. It starts from fastq files to generate Hi-C matrices and differential expression results. It also computes statistics and does pattern detection on Hi-C contact map to generate figures and tables which will be used by tailored analyses in jupyter notebooks.

Here is a visual summary of pipeline steps and their dependencies:

![](docs/img/rulegraph.svg)

For a more detailed visual summary showing input/output files, see the [filegraph](docs/img/filegraph.svg)

## Analyses

Analyses are described in jupyter notebooks located in the `docs/notebooks` folder. These notebooks are numbered to reflect the logical order in which analyses should be done. They should be executed in that order as some will generate files for the next notebook.

* Statistical exploration of chromatin loop changes [Notebook](docs/notebooks/01_diff_contacts_annot.ipynb)
* Visual exploration of global contact changes during infection [Notebook](docs/notebooks/02_diff_contacts_viz.ipynb)
* Analysis of interchromosomal contacts changes [Notebook](docs/notebooks/03_interchrom_contacts.ipynb)
* Detection and overview of chromatin insulation domains [Notebook](docs/notebooks/04_domains_analyses.ipynb)
* Analysis of the relationship between expression and contacts changes during infection [Notebook](docs/notebooks/05_infection_contacts_diff_expr.ipynb)
* Analysis of gene coexpression versus contact changes using lifted-over expression data from [Li et al. 2020](https://www.frontiersin.org/articles/10.3389/fcimb.2020.00428/full) [Notebook](docs/notebooks/06_Li2020_liftover_c3_infection_contacts_coexpr.ipynb)

