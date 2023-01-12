#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=5
#SBATCH --output=sourcesink1_output/res_5.out 
julia models/sourcesink1.jl --db source-sink.db -O 120 -L 30 -o sourcesink1_output