#!/bin/bash
# # to submit sbatch, sinfo, scancel, squeue 
# Le ponemos nombre: 
#SBATCH --job-name=star
#SBATCH --cpus-per-task=4 
#SBATCH --mem-per-cpu=16G 

# #SBATCH -o slurm.%j.out 
# #SBATCH -e slurm.%j.err 
#SBATCH --mail-type=END
#SBATCH --mail-user=marc.tormo@upf.edu


#how many CPUs to use on the current machine?
NUMCPUS=4

### list of samples
#FASTQLOC="results"
reads=`ls *.fastq.gz`


GTF="gtf/Schizosaccharomyces_pombe.ASM294v2.35.chr.gtf"
REF="ref/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa"


##if these programs are not in any PATH directories, please edit accordingly:
module load STAR/2.5.2b-foss-2016b


### GENERATE INDEX
if [ ! -d star_index ]; then mkdir -p star_index;fi
STAR --runThreadN $NUMCPUS --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles $REF --sjdbGTFfile $GTF --sjdbOverhang 49

### ALIGNMENT

if [ ! -d star ]; then mkdir -p star;fi

for i in $reads; do
    ### remove extension
    sample=`echo $i | sed 's/\.fastq\.gz//g'`

    STAR --runThreadN $NUMCPUS --genomeDir star_index --readFilesIn $i --readFilesCommand zcat --outFileNamePrefix star/$sample --outSAMtype BAM SortedByCoordinate
done


#### get stats with multiqc
module load MultiQC
multiqc -o star -f star


