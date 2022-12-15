#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=7
#SBATCH --output=sourcesink2_output/res_7.out 
julia models/sourcesink2.jl --db 'sourcesink-db' -O 240 -L 280 -o sourcesink2_output