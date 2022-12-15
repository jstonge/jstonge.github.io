#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=24
#SBATCH --output=sourcesink2_output/res_24.out 
julia models/sourcesink2.jl --db 'source-sink.db' -O 690 -L 720 -o sourcesink2_output