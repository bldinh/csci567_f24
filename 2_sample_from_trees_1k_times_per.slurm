#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --array=0-99
#SBATCH --time=2-00:00:00
#SBATCH --account=rohs_108
#SBATCH --partition=main
#SBATCH --output=logs/2/extract_data%a.out
#SBATCH --error=logs/2/extract_data%a.err
#SBATCH --job-name=2_extract_data
#SBATCH --mem-per-cpu=10gb
#SBATCH --mail-user=bldinh@usc.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module load python/3.7.6
python3 _2_extract_features.py --array $SLURM_ARRAY_TASK_ID --replicates 2000

