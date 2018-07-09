#!/bin/bash
#SBATCH --job-name=x2
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G 

export PATH:$PATH="/path/to/RNASeq_star_htscount_limma"

    if [  $# -le 2 ] || [ $1 == "--help" ] || [ $1 == "-h" ]
	then 
		echo -e "\nUsage: ${0} </path/with/fastq> <file.gtf> [yes|no|reverse] \n" 
		echo 
		echo "OPTIONS: "
		echo
		echo "</path/with/results/>:			        Enter the whole path for results with star directory. "
		echo "<file.gtf>:			                Specify the GTF file for generating the Index. "
		echo "[yes|no|reverse]: 				Enter whether the data is from a strand-specific assay."
		echo
		exit 1
	fi 


mkdir $1/htscount/
module load parallel/20151222
module load HTSeq/0.9.1-foss-2016b-Python-2.7.12
###  reverse gene counts with htseq-counts
function counts() {

    base=`basename "$1" | sed 's/Aligned.sortedByCoord.out.bam//g' `
    echo "processing "$base
    ### counts from bam sorted by position and not strand-specific assay
    htseq-count -f bam -r pos -s $3 $1 $2 > $4/htscount/$base\_htseq.csv
}
export -f counts

ls $1/star/*Aligned.sortedByCoord.out.bam | parallel --progress  -k counts {} $2 $3 $1


module load Python/3.5.2-foss-2016b 

create-table.py -c $1/htscount/ -s $1/star/ -o $1 -v
