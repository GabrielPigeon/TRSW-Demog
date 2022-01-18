#!/bin/bash
#SBATCH --account=def-pelleti2
#SBATCH --time=6-23:55           # time (DD-HH:MM)
#SBATCH --ntasks=3               # 3 core(CPU)
#SBATCH --job-name=313_m0b # smnsible name for the job
#SBATCH --mem-per-cpu=42G                 # Default memory per CPU is 3GB.
#SBATCH --mail-user=gabriel.pigeon@usherbrooke.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# If you would like to use more please adjust this.

## Below you can put your scripts
# If you want to load module
# module load nixpkgs/16.09
# module load gcc/7.3.0
module load r/4.1.2

# module list # List loaded modules

# Other commands can be included below
cd ~/projects/def-pelleti2/pigeonga/Trsw_CMR/

Rscript --verbose batchJob/4_runModel.R $SLURM_JOB_NAME


# R CMD BATCH batchJob/4_runModel.R
# salloc --time=1:0:0 --ntasks=1 --account=def-pelleti2
