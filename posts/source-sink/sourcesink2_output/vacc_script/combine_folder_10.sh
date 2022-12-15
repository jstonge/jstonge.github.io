#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=10
#SBATCH --output=sourcesink2_output/res_10.out 
julia models/sourcesink2.jl --db 'source-sink.db' -O 360 -L 400 -o sourcesink2_output