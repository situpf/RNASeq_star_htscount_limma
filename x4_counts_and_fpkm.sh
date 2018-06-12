#!/bin/bash
#SBATCH --job-name=x3
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G 

#SBATCH --mail-type=END
#SBATCH --mail-user=marc.tormo@upf.edu

module load R/3.3.3-foss-2016b-X11-20160819
Rscript get_counts_and_fpkm.R

