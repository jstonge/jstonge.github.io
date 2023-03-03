using Pkg; Pkg.activate("../..");
using DataFrames, CSV, ProgressMeter, ArgParse
using Parquet2: writefile

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
  
  # DATA_DIR = "sourcesink2_output"
  DATA_DIR = args["d"]
  fnames = filter(x -> endswith(x, "txt"),  readdir(DATA_DIR, join=true))[4500:5500]

  @assert length(fnames) > 0 println("There are no data files at the given directory")

  modelname = split(fnames[1], "_")[1]

  out_f = "$(modelname).parquet"
  lookup_f = "$(modelname)_lookup.parquet"

  dfs = []
  p = ProgressMeter.Progress(length(fnames))
  Threads.@threads for i=eachindex(fnames)
    fname = fnames[i]

    p_str = split(fname, "/")[end]
    p_str = replace(join(split(p_str, "_")[2:end], "_"), ".txt" => "")
  
    sol = CSV.read(fname, DataFrame; header=["timestep", "L", "value"], 
                   types=Dict(:timestep => Int, :L => Int, :value => Float32))
    
    gd = groupby(sol, [:timestep, :L])
    n = nrow(gd[1])
  
    # process functions
    df_agg = combine(gd, :value => (x -> round(sum(x), digits=8)) => :value_prop, 
                         :value => (x -> iszero(sum(x)) ? 0.0 : round(processing1(x,n), digits=8)) => :value)
  
    # Keep only unique value to save space
    unique!(df_agg, :value)
    
    df_agg[!, :name] .= p_str

    # Take timestep maxmin so all levels have same length
    gd2 = groupby(df_agg, [:L])
    minmax_timestep = minimum(combine(gd2,:timestep=>maximum=>:timestep)[!, :timestep])
    df_agg = filter(:timestep => x -> x <= minmax_timestep, df_agg)

    push!(dfs, df_agg)
    ProgressMeter.next!(p)
  end
  
  all_dfs = vcat(dfs...)

  # Write lookup for rowids -> all_dfs.name to disk
  names_params = unique(all_dfs.name)
  row_ids = Int32.(1:length(names_params))
  
  lookup_name = Dict()
  [get!(lookup_name, n, row_id) for (n, row_id) in zip(names_params, row_ids)];

  writefile(lookup_f, (param_str=names_params, row_id=row_ids))
  
  # Write output to disk
  all_dfs.name = [lookup_name[n] for n in all_dfs.name]
  all_dfs.L = Int8.(all_dfs.L)
  all_dfs.timestep = Int16.(all_dfs.timestep)

  writefile(out_f, all_dfs)
end

main()