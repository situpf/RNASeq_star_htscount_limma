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

## 1. Aligning to the reference genome - STAR
First thing to do is to align the reads in FASTQ files to a reference genome. Here we have the script called **x1_align.sh** which will generate an index of the reference genome in the first place, and then it will align our input reads giving us a BAM file (among other output files) for each fastq file (or pair of fastq files if the data is paired end) aligned. 

For running this script we have to call it from the terminal as follows: 
```bash
sbatch x1_align.sh </path/with/fastq> <file.gtf> <reference.fasta> <Overhang> </output/directory> [-p]
```
The script needs 5 mandatory arguments and an extra one specifying the type of data we're going to analyze. 

1. </path/with/fastq> : In here you have to specify the path to the directory whith the FASTQ files. 
2. <file.gtf> : GTF file for creating the index and aligning the reads. Please enter the full path. 
3. <reference.fasta> : FASTA file with the reference genome. Please enter the full path. 
4. Overhang : Specify a VALUE for the --sjdbOverhang parameter for aligning with STAR.
5. </output/directory> : Specify the directory where you want to store all the output files. 

Example of running the script with human paired data: 
```bash 
sbatch /woliveros/RNASeq_star_htscount_limma/x1_align.sh /woliveros/scratch/Test_RNASeq /woliveros/scratch/Test_RNASeq/Schizosaccharomyces_pombe.ASM294v2.39.gtf /woliveros/scratch/Test_RNASeq/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa 49 /woliveros/scratch/Test_RNASeq/results2
```
After running this script you'll find 2 new directories created inside the specified output directory with some files inside: 
* **star**:
  * sample1Aligned.sortedByCoord.out.bam
  * sample1Log.final.out
  * sample1Log.out
  * sample1Log.progress.out
  * sample1SJ.out.tab
* **star_index**:
  * All files generated when creating the indes for the reference genome

## 2. HTSeq-Count 
