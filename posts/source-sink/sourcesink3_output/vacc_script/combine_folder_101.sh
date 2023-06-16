#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=101
#SBATCH --output=sourcesink3_output/res_101.out 
julia --project=../.. InstitutionalDynamics.jl/src/sourcesink3.jl --db source-sink.db -O 3000 -L 30 -o sourcesink3_output