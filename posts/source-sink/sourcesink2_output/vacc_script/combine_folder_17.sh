#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=17
#SBATCH --output=sourcesink2_output/res_17.out 
julia models/sourcesink2.jl --db 'source-sink.db' -O 640 -L 680 -o sourcesink2_output