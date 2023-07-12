#!/bin/bash
#SBATCH --account=rrg-zovoilis
#SBATCH --time=144:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

# non pyenv tools
module load StdEnv/2020 gcc/9.3.0
module load bwa/0.7.17
module load samtools/1.13
module load apptainer/1.1.8
# load conda capabilities
source ~/.bashrc
# activate env
conda activate mbdmbd
# run
snakemake --cores 48 --profile profiles/slurm/ --use-singularity
