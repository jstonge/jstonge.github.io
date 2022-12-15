#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=8
#SBATCH --output=sourcesink2_output/res_8.out 
julia models/sourcesink2.jl --db 'sourcesink-db' -O 280 -L 320 -o sourcesink2_output