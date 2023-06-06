#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=1
#SBATCH --output=sourcesink2_output/res_1.out 
julia models/sourcesink2.jl --db source-sink.db -O 0 -L 25 -o sourcesink2_output