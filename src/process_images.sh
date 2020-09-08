#!/bin/bash

#SBATCH --job-name=binarize_images    # Job name
#SBATCH --nodes=1                     # Use one node
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --mem-per-cpu=1gb             # Memory per processor
#SBATCH --time=00:10:00               # Time limit hrs:min:sec
#SBATCH --output=log/%j_%A-%a.out     # Standard output and error log
#SBATCH --array=1-5                   # Array range

# iterate through list of input images
file=$(ls data/raw/*.tif | sed -n ${SLURM_ARRAY_TASK_ID}p)

outdir="data/processed"

./src/binarize_images.py -i ${file} -o ${outdir}
