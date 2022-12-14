#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=52
julia models/sourcesink2.jl --beta 0.07:0.05:0.22 -a 0.5 -g 1.0 -r 0.25 -b 0.12 -c 2.05 -m 0.0001 -o sourcesink2_output