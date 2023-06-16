#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=02:59:59
#SBATCH --job-name=16
#SBATCH --output=sourcesink3_output/res_16.out 
julia --project=../.. InstitutionalDynamics.jl/src/sourcesink3.jl --db source-sink.db -O 450 -L 30 -o sourcesink3_output