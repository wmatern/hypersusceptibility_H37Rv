# Tn-seq analysis of Mycobacterium Tuberculosis antibiotic hypersusceptibility
This workflow, implemented in a Jupyter notebook, seeks to reproduce the results presented in our upcoming paper: "Functional whole genome screen of nutrient-starved Mycobacterium tuberculosis identifies genes involved in antibiotic tolerance". Following the steps below should allow you to get identical results to those reported in our paper.

## Repo Authors

* William Matern (@wmatern)

## Dependencies

You will need to install the following packages (and their dependencies) in order to run this workflow:
* conda (tested with version: 23.1.0)

## Raw Data
The necessary Tn-seq data can be found in NCBI SRA under BioProject: PRJNA946182. You will need to add all of the associated fastq files (SRR23906849-SRR23906902), 108 files, with filenames ending in \*\_1.fastq and \*\_2.fastq) into the input/fastq folder before running the Jupyter notebook. Additionally, the genome information contained in the Genbank file (.gb, .gbff) along with corresponding FASTA files (.fa, .fna) can be downloaded at https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2/.

The expected file names in order to run this workflow are:

    input/Genome/GCF_000195955.2_ASM19595v2_genomic.fna
    input/Genome/GCF_000195955.2_ASM19595v2_genomic.gbff

You can view input/Samples\_Labels.csv to see the list of sample information.

## Usage

### Step 1: Download workflow
If you simply want to run this workflow, download and extract the [latest release](https://github.com/).

### Step 2: Execute workflow
Run each cell of Jupyter notebook in order.

## Report Issue
If you have any questions or issues reproducing our results please send me an email at maternwill@gmail.com.
