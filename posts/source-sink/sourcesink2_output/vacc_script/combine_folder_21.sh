#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=21
#SBATCH --output=sourcesink2_output/res_21.out 
julia models/sourcesink2.jl --db 'source-sink.db' -O 800 -L 840 -o sourcesink2_output