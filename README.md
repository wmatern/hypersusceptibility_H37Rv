# Tn-seq analysis of Mycobacterium Tuberculosis antibiotic hypersusceptibility
This workflow, implemented in a Jupyter notebook, seeks to reproduce the results presented in our upcoming paper: "Functional whole genome screen of nutrient-starved Mycobacterium tuberculosis identifies genes involved in antibiotic tolerance". Following the steps below should allow you to get identical results to those reported in our paper.

## Repo Authors

* William Matern (@wmatern)

## Dependencies

You will need to install the following packages (and their dependencies) in order to run this workflow:
* conda (tested with version: 23.1.0)

## Raw Data
The necessary Tn-seq data can be found in NCBI SRA under BioProject: PRJNA946182. You will need to download all of the associated fastq files (SRR23906849-SRR23906902), 108 files, with filenames ending in \*\_1.fastq and \*\_2.fastq) into the appropriate input/fastq folder before running the Jupyter notebooks. There are two analysis folders: 7H9_analysis (rich medium) and PBS_analysis (nutrient starvation). View input/Samples\_Labels.csv to see the list of expected fastq files for each workflow.  Additionally, the genome information contained in the Genbank file (.gb, .gbff) along with corresponding FASTA files (.fa, .fna) can be downloaded at https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2/.

The expected genome file names in order to run this workflow are:

    input/Genome/GCF_000195955.2_ASM19595v2_genomic.fna
    input/Genome/GCF_000195955.2_ASM19595v2_genomic.gbff

## Usage
There are two workflows used to analyze the two different experiments: 7H9_analysis (rich medium) and PBS_analysis (nutrient starvation). You'll need to run the appropriate Jupyter notebook inside the respective folder to reproduce our results.

### Step 1: Download workflow
If you simply want to run these two workflows, download and extract the [latest release](https://github.com/).

### Step 2: Install conda environment
conda env create -f environment.yml

### Step 3: Execute workflow
Activate the conda environment and run each cell of appropriate Jupyter notebook in order. This will reproduce the plots and supplementary tables.

## Report Issue
If you have any questions or issues reproducing our results please send me an email at maternwill@gmail.com.
