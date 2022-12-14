#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=42
julia models/sourcesink2.jl --beta 0.07:0.05:0.22 -a 0.5 -g 1.0 -r 0.1 -b 0.17 -c 1.05 -m 0.0001 -o sourcesink2_output