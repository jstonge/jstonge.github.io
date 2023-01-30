#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=39
#SBATCH --output=sourcesink3_output/res_39.out 
julia models/sourcesink3.jl --db source-sink.db -O 2280 -L 60 -o sourcesink3_output