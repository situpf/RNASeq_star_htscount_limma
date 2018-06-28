# RNASeq_star_htscount_limma
Pipeline for running RNASeq data. In this directory you'll find all necessary scripts for analyzing RNASeq data and extracting some important information as coverage plots and a summary of counts. The starting point for doing so are the Fastq files. 

Scripts in this directory: 
* x1_align.sh
* x2_get_counts.sh
* x3_rnaseq_limma.R
* x4_counts_and_fpkm.sh
* get_counts_and_fpkm.R
* create-table.py

# WORKFLOW

## 1. Aligning to the reference genome
First thing to do is to align the reads in FASTQ files to a reference genome. Here we have the script called **x1_align.sh** which will generate an index of the reference genome in the first place, and then it will align our input reads giving us a BAM file (among other output files) for each fastq file (or pair of fastq files if the data is paired end) aligned. 

For running this script we have to call it from the terminal as follows: 
```bash
sbatch x1_align.sh </path/with/fastq> <file.gtf> <reference.fasta> <Overhang> </output/directory> [-p]
```
The script needs 5 mandatory arguments and an extra one specifying the type of data we're going to analyze. 

1. </path/with/fastq> : This argu
