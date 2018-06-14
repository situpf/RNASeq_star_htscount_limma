#!/bin/bash
#SBATCH --job-name=x2
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G 

#SBATCH --mail-type=END
#SBATCH --mail-user=marc.tormo@upf.edu

#module load R/3.2.3
#Rscript x2_rnaseq_get-counts.R

#$1 = directory with sample files (fastq)
#$2 = GTF file ; input the whole path


mkdir htscount/
module load parallel/20151222
module load HTSeq/0.9.1-foss-2016b-Python-2.7.12
###  reverse gene counts with htseq-counts
function counts() {
    base=`basename "$1" | sed 's/.fastq.gz//g'`
    echo "processing "$base
    ### counts from bam sorted by position and not strand-specific assay
    htseq-count -f bam -r pos -s reverse star/$base\Aligned.sortedByCoord.out.bam $2 > htscount/$base\_htseq.csv
}
export -f counts

ls $1*.fastq.gz | parallel --progress  -k counts {}
