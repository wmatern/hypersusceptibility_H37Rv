# Tn-seq analysis of Mycobacterium Avium antibiotic hypersusceptibility
This workflow, implemented in a Jupyter notebook and two scripts (make_barplots.py and make_venn.R), seeks to reproduce the results presented in our upcoming paper: "Genetic determinants of intrinsic antibiotic tolerance in Mycobacterium avium". Following the steps below should allow you to get identical results to those reported in our paper.

## Authors

* William Matern (@wmatern)

## Dependencies

You will need to install the following packages (and their dependencies) in order to run this workflow:
* conda (tested with version: 4.8.3)
* R and eulerr

## Raw Data
The necessary Tn-seq data can be found in NCBI SRA under BioProject: XXXXX. You will need to add all of the associated fastq files (SRRXXX-SRRXXX, XXX files, with files ending in \*\_1.fastq and \*\_2.fastq) into the input/ folder before running the Jupyter notebook. Additionally, the genome information contained in the Genbank file (.gb, .gbff) along with corresponding FASTA files (.fa, .fna) can be downloaded at https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2/. You will need to rename the genome information files to use the prefix XXXX. 

The expected file names in order to run this workflow are:

    input/Genome/XXXX.fa
    input/Genome/XXXX.gb

You can view input/Samples\_Labels.csv to see the list of sample information.

## Usage

### Step 1: Download workflow
If you simply want to run this workflow, download and extract the [latest release](https://github.com/).

### Step 2: Execute workflow
Run each cell of Jupyter notebook in order.

### Step 3: Run make_barplots.py to create barplots
`python3 make_barplots.py`

### Step 4: Change working directory and run make_venn.R
XXXX
`vim make_venn.R`

Set the first line of make_venn.R to the appropriate location.

`Rscript make_venn.R`

## Report Issue
If you have any questions or issues reproducing our results please send me an email at maternwill@gmail.com.
