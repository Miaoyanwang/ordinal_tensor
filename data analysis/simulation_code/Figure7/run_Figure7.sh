#!/bin/bash
#SBATCH --mail-user=miaoyan@stat.wisc.edu
#SBATCH -p long
#SBATCH -t 5-00:00:00
#SBATCH -n 4
#SBATCH --mem-per-cpu=10000M
module load R/R-3.5.3
export R_LIBS=/workspace/miaoyan/x86_64-pc-linux-gnu-library/3.4/
R CMD BATCH Figure8.R output_new.Rout
