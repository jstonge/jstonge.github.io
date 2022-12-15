#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=6
#SBATCH --output=sourcesink2_output/res_6.out 
julia models/sourcesink2.jl --db 'source-sink.db' -O 200 -L 240 -o sourcesink2_output