#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=5
#SBATCH --output=sourcesink2_output/res_5.out 
julia models/sourcesink2.jl --db source-sink.db -O 100 -L 25 -o sourcesink2_output