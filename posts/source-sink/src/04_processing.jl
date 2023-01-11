# This script creates a database that contains all the results that we want
# to visualize.

using Pkg; Pkg.activate("../..");
using SQLite, DataFrames, CSV, ProgressMeter, ArgParse

function parse_commandline()
  s = ArgParseSettings()

  @add_arg_table! s begin
      "-d"
      default = "."
      help = "Name of the db containing params. It should only include one model."
      "-o"
      help = "Name of the output DB."
    end

  return parse_args(s)
end

function create_res_db()
  SQLite.execute(db, """
  CREATE TABLE sourcesink2 (
    name STRING,
    timestep INT,
    L INT,
    value REAL,
    PRIMARY KEY (name, timestep, L)
  )
  """)
end

"""
Combine all sols in `sourcesink2_output/`, using only unique value, into database.
"""
function main()
  args = parse_commandline()
  
  # DATA_DIR = "sourcesink2_output"
  DATA_DIR = args["d"]
  fnames = filter(x -> endswith(x, "txt"),  readdir(DATA_DIR, join=true))

  @assert length(fnames) > 0 println("There are no data files at the given directory")

  db = SQLite.DB(args["o"])
  # modelname = split(fnames[1], "_")[1]
  
  try
    create_res_db()
  catch
    println("DB already exists")
  finally
    already_done = DBInterface.execute(db, """SELECT DISTINCT name FROM sourcesink2""") |> DataFrame
  end
  

  @showprogress for fname in fnames
    sol = CSV.read(fname, DataFrame; header=["timestep", "L", "value"])
    fname = split(fname, "/")[end]
    p_str = replace(join(split(fname, "_")[2:end], "_"), ".txt" => "")
    
    if p_str ∉ already_done.name
      
      gd = groupby(sol, [:timestep, :L])
      
      # Each combination of timestep and level must have the same number of rows. 
      # If not, we are unable to combine them. To ensure this is the case,
      # we take the min of each combination. 
      #!TODO Verify that this will not throw away important information down the line
      n = nrow(gd[1])
      # n = minimum([nrow(g) for g in gd])
      
      df_agg = combine(gd, :value => x -> iszero(sum(x)) ? 0.0 : sum((collect(0:(n-1)) / n) .* x) / sum(x)) 
      rename!(df_agg, Dict(:value_function => "value")) 
      unique!(df_agg, :value)
  
      # write to db
      for row in eachrow(df_agg)
        params = (p_str, row["timestep"], row["L"], row["value"])
        SQLite.execute(db, """INSERT INTO sourcesink2 VALUES (?, ?, ?, ?)""", params)
      end
      
    end

  end
end

main()



# SQLite.execute(db, """
# DROP TABLE sourcesink2
# """)

# OUTPUT_DIR = "/home/jstonge/OneDrive/teenyverse/sourcesink2/sourcesink2_output"

# dfs = []
# counter = 1
# @showprogress for fname in fnames
#   fname = fnames[160]  
#   sol = CSV.read(fname, DataFrame; header=["timestep", "L", "value"])

#   fname = replace(fname, "$(OUTPUT_DIR)/" => "")
#   modelname = split(fname, "_")[1]
#   p_str = replace(join(split(fname, "_")[2:end], "-"), ".txt" => "")

#   gd = groupby(sol, [:timestep, :L])
#   n = nrow(gd[1])
  
#   # fname
#   # solfromdb = CSV.read("sourcesink2_0.07_0.6_1.0_0.25_0.17_1.05_0.0001.txt", DataFrame; header=["timestep", "L", "value"])
#   # sol6 = CSV.read("sourcesink2_0.07_0.6_1.0_0.25_0.17_1.05_0.0001_6digits.txt", DataFrame; header=["timestep", "L", "value"])
#   # sol10 = CSV.read("sourcesink2_0.07_0.6_1.0_0.25_0.17_1.05_0.0001_10digits.txt", DataFrame; header=["timestep", "L", "value"])
#   # solnoRounding = CSV.read("sourcesink2_0.07_0.6_1.0_0.25_0.17_1.05_0.0001_noRounding.txt", DataFrame; header=["timestep", "L", "value"])
  
#   # subset(sol, :timestep => x -> x .== 27, :L => x -> x .== 4)
#   # subset(solfromdb, :timestep => x -> x .== 27, :L => x -> x .== 4)
#   # subset(sol6, :timestep => x -> x .== 27, :L => x -> x .== 4)
#   # subset(sol10, :timestep => x -> x .== 27, :L => x -> x .== 4)
#   # subset(solnoRounding, :timestep => x -> x .== 27, :L => x -> x .== 4)

#   df_agg = combine(gd, :value => x -> iszero(sum(x)) ? 0.0 : sum((collect(0:(n-1)) / n) .* x) / sum(x)) 
#   rename!(df_agg, Dict(:value_function => "value")) 
#   unique!(df_agg, :value)
#   df_agg[!, :name] .= p_str
#   push!(dfs, df_agg)
#   counter += 1
# end

# counter

# all_dfs = vcat(dfs...)

# all_dfs |> SQLite.load!(db2,"$(p_str)") 


# # SQLite.execute(db2, """
# # DROP TABLE `0.07-0.5-0.9-0.1-0.12-0.55-0.0001`
# # """)
