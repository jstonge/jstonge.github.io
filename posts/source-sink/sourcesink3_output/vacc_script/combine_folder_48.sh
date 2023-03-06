#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=48
#SBATCH --output=sourcesink3_output/res_48.out 
julia models/sourcesink3.jl --db source-sink.db -O 1410 -L 30 -o sourcesink3_output