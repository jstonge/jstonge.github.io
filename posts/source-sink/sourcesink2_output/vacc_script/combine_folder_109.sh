#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=109
julia models/sourcesink2.jl --beta 0.07:0.05:0.22 -a 0.6 -g 0.9 -r 0.1 -b 0.12 -c 0.55 -m 0.0001 -o sourcesink2_output