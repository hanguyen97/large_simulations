#!/bin/bash
#SBATCH --job-name=attglasso
#SBATCH --output=logs/attglasso_%A_%a.out
#SBATCH --error=logs/attglasso_%A_%a.err
#SBATCH --array=1-20%2   # Only 2 tasks will run at the same time
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tn325@cornell.edu 


Rscript ~/ATTglasso/03_compare_arr.R $SLURM_ARRAY_TASK_ID
