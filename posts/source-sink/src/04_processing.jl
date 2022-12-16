# This script creates a database that contains all the results that we want
# to visualize.

using Pkg; Pkg.activate("../..");
using SQLite, DataFrames, CSV, ProgressMeter

db2 = SQLite.DB("sourcesink2-res.db")
OUTPUT_DIR = "/home/jstonge/OneDrive/teenyverse/sourcesink2/sourcesink2_output"
fnames = filter(x -> endswith(x, "txt"),  readdir(OUTPUT_DIR, join=true))

dfs = []
counter = 1
@showprogress for fname in fnames
  fname = fnames[160]  
  sol = CSV.read(fname, DataFrame; header=["timestep", "L", "value"])

  fname = replace(fname, "$(OUTPUT_DIR)/" => "")
  modelname = split(fname, "_")[1]
  p_str = replace(join(split(fname, "_")[2:end], "-"), ".txt" => "")

  gd = groupby(sol, [:timestep, :L])
  n = nrow(gd[1])
  
  # fname
  # solfromdb = CSV.read("sourcesink2_0.07_0.6_1.0_0.25_0.17_1.05_0.0001.txt", DataFrame; header=["timestep", "L", "value"])
  # sol6 = CSV.read("sourcesink2_0.07_0.6_1.0_0.25_0.17_1.05_0.0001_6digits.txt", DataFrame; header=["timestep", "L", "value"])
  # sol10 = CSV.read("sourcesink2_0.07_0.6_1.0_0.25_0.17_1.05_0.0001_10digits.txt", DataFrame; header=["timestep", "L", "value"])
  # solnoRounding = CSV.read("sourcesink2_0.07_0.6_1.0_0.25_0.17_1.05_0.0001_noRounding.txt", DataFrame; header=["timestep", "L", "value"])
  
  # subset(sol, :timestep => x -> x .== 27, :L => x -> x .== 4)
  # subset(solfromdb, :timestep => x -> x .== 27, :L => x -> x .== 4)
  # subset(sol6, :timestep => x -> x .== 27, :L => x -> x .== 4)
  # subset(sol10, :timestep => x -> x .== 27, :L => x -> x .== 4)
  # subset(solnoRounding, :timestep => x -> x .== 27, :L => x -> x .== 4)

  df_agg = combine(gd, :value => x -> iszero(sum(x)) ? 0.0 : sum((collect(0:(n-1)) / n) .* x) / sum(x)) 
  rename!(df_agg, Dict(:value_function => "value")) 
  unique!(df_agg, :value)
  df_agg[!, :name] .= p_str
  push!(dfs, df_agg)
  counter += 1
end

counter

all_dfs = vcat(dfs...)

all_dfs |> SQLite.load!(db2,"$(p_str)") 


# SQLite.execute(db2, """
# DROP TABLE `0.07-0.5-0.9-0.1-0.12-0.55-0.0001`
# """)
