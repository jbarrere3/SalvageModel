#!/bin/bash
#SBATCH --job-name=SalvageModel
#SBATCH --account=project_2005142
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --partition=small

module load r-env

srun singularity_wrapper exec Rscript --no-save cmd.R