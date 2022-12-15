#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=23
#SBATCH --output=sourcesink2_output/res_23.out 
julia models/sourcesink2.jl --db 'source-sink.db' -O 660 -L 690 -o sourcesink2_output