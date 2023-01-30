#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=118
#SBATCH --output=sourcesink3_output/res_118.out 
julia models/sourcesink3.jl --db source-sink.db -O 7020 -L 60 -o sourcesink3_output