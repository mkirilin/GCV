#!/bin/bash --login
#SBATCH --job-name=gcv_blur
#SBATCH --output=gcv_blur.out
#SBATCH --error=gcv_blur.err
#SBATCH --time=24:00:00
#SBATCH --mem=400GB
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=kirilin@math.tu-berlin.de
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --exclusive

module load matlab/2024a
matlab -nodisplay -nosplash -r "\
  addpath('/work/kirilin/tools/AIRToolsII'); \
  AIRToolsII_setup(); \
  addpath('/work/kirilin/tools/IRtools'); \
  IRtools_setup(); \
  run('Tests.m'); \
  exit;
"
sacct --format=JobName,JobID,Elapsed,MaxRSS,MaxVMSize,AllocTRES