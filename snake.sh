#!/bin/bash
#SBATCH -J test
#SBATCH --mem=64G
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -c 4

module load system/singularity-2.5.1
singularity run --bind data:/scif/data --bind snakevirome/:/snakevirome --bind /work/aexbrayat/samples/:/samples  snakevir.img  run snakemake -p
