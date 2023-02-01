using Pkg; Pkg.activate("../..");
using DataFrames, CSV, ProgressMeter, ArgParse, Parquet

function parse_commandline()
  s = ArgParseSettings()

  @add_arg_table! s begin
      "-d"
      default = "."
      help = "Directory containing the txt files."
    end

  return parse_args(s)
end

processing1(x, n) = sum((collect(0:(n-1)) / n) .* x) / sum(x)

"""
Combine all sols and proportion in `sourcesink_output/`, using only unique value, into a parquet file.
"""
function main()
  args = parse_commandline()
  
  # DATA_DIR = "sourcesink3_output"
  DATA_DIR = args["d"]
  fnames = filter(x -> endswith(x, "txt"),  readdir(DATA_DIR, join=true))

  @assert length(fnames) > 0 println("There are no data files at the given directory")

  modelname = split(fnames[1], "_")[1]
  
  dfs = []
  @showprogress for fname in fnames   
    # fname = fnames[1]
    sol = CSV.read(fname, DataFrame; header=["timestep", "L", "value"], types=Dict(:timestep => Int, :L => Int, :value => Float32))
    fname_parts = split(fname, "/")
    fname = fname_parts[end]
    p_str = replace(join(split(fname, "_")[2:end], "_"), ".txt" => "")
   
    gd = groupby(sol, [:timestep, :L])
    n = nrow(gd[1])

    df_agg = combine(gd, :value => (x -> round(sum(x), digits=8)) => :value_prop, 
                         :value => (x -> iszero(sum(x)) ? 0.0 : round(processing1(x,n), digits=8)) => :value)

    unique!(df_agg, :value)
    
    df_agg[!, :name] .= p_str
    push!(dfs, df_agg)      
  end

  all_dfs = vcat(dfs...) 
  write_parquet("$(modelname).parquet", all_dfs)

end

main()


# SQLite.execute(db, """
# DROP TABLE sourcesink2
# """)
