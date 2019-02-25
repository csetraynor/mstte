#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=4571mb
#SBATCH --time=46:00:00

module load mro/3.5.1-foss-2017a

# Use R CMD BATCH to run your script
R CMD BATCH ~/rfactory/mstte-data/baseline_sim.R
