#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=2
#SBATCH --output=sourcesink2_output/res_2.out 
julia models/sourcesink2.jl --db source-sink.db -O 30 -L 60 -o sourcesink2_output