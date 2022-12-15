#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=4
#SBATCH --output=sourcesink2_output/res_4.out 
julia models/sourcesink2.jl --db 'sourcesink-db' -O 120 -L 160 -o sourcesink2_output