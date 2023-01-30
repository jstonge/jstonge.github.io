#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=81
#SBATCH --output=sourcesink3_output/res_81.out 
julia models/sourcesink3.jl --db source-sink.db -O 4800 -L 60 -o sourcesink3_output