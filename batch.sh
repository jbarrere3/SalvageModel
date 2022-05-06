#!/bin/bash
#SBATCH --job-name=FinnishNFI
#SBATCH --account=project_2005142
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=small

module load r-env-singularity/4.1.1

srun singularity_wrapper exec Rscript --no-save cmd.R