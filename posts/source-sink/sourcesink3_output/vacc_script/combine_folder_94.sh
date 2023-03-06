#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=94
#SBATCH --output=sourcesink3_output/res_94.out 
julia models/sourcesink3.jl --db source-sink.db -O 2790 -L 30 -o sourcesink3_output