#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=107
#SBATCH --output=sourcesink3_output/res_107.out 
julia --project=../.. InstitutionalDynamics.jl/src/sourcesink3.jl --db source-sink.db -O 3180 -L 30 -o sourcesink3_output