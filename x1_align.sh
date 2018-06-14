#!/bin/bash
# # to submit sbatch, sinfo, scancel, squeue 
# Le ponemos nombre: 
#SBATCH --job-name=star
#SBATCH --cpus-per-task=4 
#SBATCH --mem-per-cpu=16G 

# #SBATCH -o slurm.%j.out 
# #SBATCH -e slurm.%j.err 
#SBATCH --mail-type=END
#SBATCH --mail-user=winona.oliveros01@estudiant.upf.edu

#$1 = path with sample files (fastq)
#$2 = GTF file (whole path)
#$3 = REF file (FASTA - whole path)
#$4 = Specify the --sjdbOverhang for generating the index with STAR (ENTER A VALUE)
#$5 = Directory for output files

    if [  $# -le 4 ] || [ $1 == "--help" ] || [ $1 == "-h" ]
	then 
		echo -e "\nUsage: x1_align.sh </path/with/fastq> <file.gtf> <reference.fasta> <Overhang> </output/directory> \n" 
		echo 
		echo "OPTIONS: "
		echo
		echo "</path/with/fastq>:			Enter the whole path with the sample .fastq files. "
		echo "<file.gtf>:			        Specify the GTF file for generating the Index. "
		echo "<reference.fasta>:		    	Specify the FASTA file with the reference genome (only one allowed)."
       		echo "<Overhang>:                   		Specify a VALUE for the --sjdbOverhang parameter for aligning with STAR."
        	echo "</output/directory/>:         		Enter the desired directory where will be stored output files. "
		echo
		exit 1
	fi 

#how many CPUs to use on the current machine?
NUMCPUS=4

### list of samples
#FASTQLOC="results"
reads=`ls $1*.fastq.gz`


GTF=$2
REF=$3


##if these programs are not in any PATH directories, please edit accordingly:
module load STAR/2.5.2b-foss-2016b


### GENERATE INDEX
if [ ! -d star_index ]; then mkdir -p $5star_index;fi
STAR --runThreadN $NUMCPUS --runMode genomeGenerate --genomeDir $5star_index --genomeFastaFiles $REF --sjdbGTFfile $GTF --sjdbOverhang $4

### ALIGNMENT

if [ ! -d star ]; then mkdir -p $5star;fi

for i in $reads; do
    ### remove extension
    sample=`echo $i | sed 's/\.fastq\.gz//g'`

    STAR --runThreadN $NUMCPUS --genomeDir $5star_index --readFilesIn $i --readFilesCommand zcat --outFileNamePrefix $5star/$sample --outSAMtype BAM SortedByCoordinate
done


#### get stats with multiqc
module load MultiQC
multiqc -o $5star -f $5star


