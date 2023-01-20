using Pkg; Pkg.activate("../..");
using SQLite, DataFrames, CSV, ProgressMeter, ArgParse

function parse_commandline()
  s = ArgParseSettings()

  @add_arg_table! s begin
      "-d"
      default = "."
      help = "Directory containing the txt files."
    end

  return parse_args(s)
end

function create_res_db(db, modelname)
  SQLite.execute(db, """
  CREATE TABLE IF NOT EXISTS $(modelname) (
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
  
  # DATA_DIR = "sourcesink1_output"
  DATA_DIR = args["d"]
  fnames = filter(x -> endswith(x, "txt"),  readdir(DATA_DIR, join=true))

  @assert length(fnames) > 0 println("There are no data files at the given directory")

  db = SQLite.DB("source-sink-res.db")
  modelname = split(fnames[1], "_")[1]
  
  create_res_db(db, modelname)
  already_done = DBInterface.execute(db, """SELECT DISTINCT name FROM sourcesink""") |> DataFrame

  dfs = []
  @showprogress for fname in fnames
    sol = CSV.read(fname, DataFrame; header=["timestep", "L", "value"])
    fname_parts = split(fname, "/")
    fname = fname_parts[end]
    p_str = replace(join(split(fname, "_")[2:end], "_"), ".txt" => "")
    set_already_done = Set(already_done.name)
    if !(p_str in set_already_done)
      
      gd = groupby(sol, [:timestep, :L])
      n = nrow(gd[1])
      
      df_agg = combine(gd, :value => x -> iszero(sum(x)) ? 0.0 : sum((collect(0:(n-1)) / n) .* x) / sum(x)) 
      rename!(df_agg, Dict(:value_function => "value")) 
      unique!(df_agg, :value)

      df_agg[!, :name] .= p_str
      push!(dfs, df_agg)      
    end
  end

    all_dfs = vcat(dfs...) 
    all_dfs |> SQLite.load!(db, "$(modelname)") 

end

main()


# SQLite.execute(db, """
# DROP TABLE sourcesink2
# """)