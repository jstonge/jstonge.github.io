#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=47
#SBATCH --output=sourcesink3_output/res_47.out 
julia models/sourcesink3.jl --db source-sink.db -O 1380 -L 30 -o sourcesink3_output