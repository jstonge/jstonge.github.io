#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=12
#SBATCH --output=sourcesink2_output/res_12.out 
julia --project=../.. InstitutionalDynamics.jl/src/sourcesink2.jl --db source-sink.db -O 330 -L 30 -o sourcesink2_output