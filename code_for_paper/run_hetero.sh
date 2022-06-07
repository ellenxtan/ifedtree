#!/bin/bash

#SBATCH --job-name tau-CT-ni500
#SBATCH --mail-user=XXXX
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12g
#SBATCH --time=48:00:00
#SBATCH --array=1-8
#SBATCH --account=XXXX
#SBATCH --partition=XXXX

#  Put your job commands after this line
module load gcc/8.2.0 r/4.1.0
Rscript --vanilla sbatch_hetero.R

# sbatch run_hetero.sh  #run batch script
