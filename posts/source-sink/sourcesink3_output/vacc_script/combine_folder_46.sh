#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=46
#SBATCH --output=sourcesink3_output/res_46.out 
julia models/sourcesink3.jl --db source-sink.db -O 2700 -L 60 -o sourcesink3_output