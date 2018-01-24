#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xxxx@ccf.org
#SBATCH --job-name=star
#SBATCH --nodes=1
#SBATCH -o STAR_jan1_run.o%j

module load STAR/2.5.2b

python3.6 run_STAR_aligner.py
