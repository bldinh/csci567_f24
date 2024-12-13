#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --array=1-100
#SBATCH --time=2-00:00:00
#SBATCH --partition=main
#SBATCH --output=logs/1_sim_hardsweep_arr%a.out
#SBATCH --error=logs/1_sim_hardsweep_arr%a.err
#SBATCH --job-name=1_simhardsweep
#SBATCH --mem-per-cpu=10gb

module load python/3.7.6
python3 _1_hardsweep_slim.py --array $SLURM_ARRAY_TASK_ID --jobs 1000


