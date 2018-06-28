# RNASeq_star_htscount_limma


## **TABLE OF CONTENTS**

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->
- [RNASeq_star_htscount_limma](#rnaseq_star_htscount_limma)
  - [Table of contents](#table-of-contents)
  - [Description](#description)
  - [Workflow](#workflow)
   - [Aligning to the reference genome](#aligning-to-the-reference-genome-star)
   - [HTSeq-Count](#htseq-count)
  

  
<!-- /TOC -->

## Description

Pipeline for running RNASeq data. In this directory you'll find all necessary scripts for analyzing RNASeq data and extracting some important information as coverage plots and a summary of counts. The starting point for doing so are the Fastq files. 

Scripts in this directory: 
* x1_align.sh
* x2_get_counts.sh
* x3_rnaseq_limma.R
* x4_counts_and_fpkm.sh
* get_counts_and_fpkm.R
* create-table.py

## WORKFLOW

### Aligning to the reference genome STAR
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
6. \[-p] : This last optional argument is to tell the script if you're working with paired data, if it's not specified it wil assume that you're working with single-end data.  

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

### HTSeq-Count 
Now that we have aligned the reads to a reference genome and we have a BAM file with the results, it's time to use HTSeq to get the counts of the reads to each feature on a gtf file. Finally thanks to the python script available in this directory it will create a .txt file with a tabulated table with a summary of the STAR and HTSeq results. 

For running this script we have to call it from the terminal as follows: 
```bash 
sbatch x2_get_counts.sh </path/with/results> <file.gtf> </yes/no/reverse>
```
This bash script needs 3 mandatory arguments: 

1. </path/with/results> : Specify the directory where you have the results from the previous step (same argument as number 5 of the previous script). 
2. <file.gtf> : GTF file with the features you want to get the counts. 
3. </yes/no/reverse> : Specify if the data is stranded, no stranded or reverse stranded. 

