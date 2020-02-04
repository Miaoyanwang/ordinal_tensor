#!/bin/bash
#SBATCH --mail-user=miaoyan@stat.wisc.edu
#SBATCH -p long
#SBATCH -t 5-00:00:00
#SBATCH -n 4
#SBATCH --mem-per-cpu=10000M
export MATLABPATH=/workspace/miaoyan/TensorCompletion_1bit_noisy-master/
export MATLABPATH=/workspace/miaoyan/TensorCompletion_1bit_noisy-master/tensorlab/
matlab -nodisplay -nosplash < One_bit_simulation_Fig8.m 




