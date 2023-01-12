#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=14
#SBATCH --output=sourcesink2_output/res_14.out 
julia models/sourcesink2.jl --db source-sink.db -O 390 -L 30 -o sourcesink2_output