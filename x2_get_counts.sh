#!/bin/bash
#SBATCH --job-name=x2
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G 

#SBATCH --mail-type=END
#SBATCH --mail-user=winona.oliveros01@estudiant.upf.edu

#module load R/3.2.3
#Rscript x2_rnaseq_get-counts.R

#$1 = directory with sample files (fastq)
#$2 = GTF file ; input the whole path
#$3 = Directory with star output files [/path/to/results/]star -> last part /star/ not necessary

    if [  $# -le 2 ] || [ $1 == "--help" ] || [ $1 == "-h" ]
	then 
		echo -e "\nUsage: x2_get_counts.sh </path/with/fastq> <file.gtf> </star/output/directory> \n" 
		echo 
		echo "OPTIONS: "
		echo
		echo "</path/with/fastq>:			        Enter the whole path with the sample .fastq files. "
		echo "<file.gtf>:			                Specify the GTF file for generating the Index. "
        	echo "<star/output/directory/>:         		Enter the directory with the star files. Example: /homes/users/username/results/star is the directory -> specify /homes/users/username/results/ "
		echo "[yes|no|reverse]: 				Enter whether the data is from a strand-specific assay."
		echo
		exit 1
	fi 


mkdir $3htscount/
module load parallel/20151222
module load HTSeq/0.9.1-foss-2016b-Python-2.7.12
###  reverse gene counts with htseq-counts
function counts() {
    base=`basename "$1" | sed 's/.fastq.gz//g'`
    echo "processing "$base
    ### counts from bam sorted by position and not strand-specific assay
    htseq-count -f bam -r pos -s $4 $3star/$base\Aligned.sortedByCoord.out.bam $2 > $3htscount/$base\_htseq.csv
}
export -f counts

ls $1*.fastq.gz | parallel --progress  -k counts {} $2 $3 $4
