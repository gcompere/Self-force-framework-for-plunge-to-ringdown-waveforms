#!/bin/bash
# Submission script for Lyra
#SBATCH --job-name=Plunge_C_HPC
#SBATCH --time=99:00:00 # hh:mm:ss
# Best if small job: #SBATCH --cpus-per-task=6
# Best if large job: #SBATCH --ntasks=60 # maximum 120
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=5000 # megabytes

module load Mathematica/12.3.1

wolframscript -script Plunge_C_HPC.wls
