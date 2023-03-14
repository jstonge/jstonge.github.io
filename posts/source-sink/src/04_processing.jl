using Pkg; Pkg.activate("../..");
using DataFrames, CSV, ProgressMeter, ArgParse, Parquet
using Pipe: @pipe
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

processing1(x, n) = round(sum((collect(0:(n-1)) / n) .* x) / sum(x), digits=7)

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

  out_f = "$(modelname).parquet"
  lookup_f = "$(modelname)_lookup.parquet"

  dfs = []
  p = ProgressMeter.Progress(length(fnames))
  Threads.@threads for i=eachindex(fnames)
    fname = fnames[i]

    p_str = @pipe split(fname, "/")[end] |> 
      split(_, "_")[2:end] |> 
      join(_, "_") |>
      replace(_, ".txt" => "")
  
    sol = CSV.read(fname, DataFrame; header=["timestep", "L", "value"], 
                   types=Dict(:timestep => Int, :L => Int, :value => Float32))

    gd = groupby(sol, [:timestep, :L])
    n = nrow(gd[1])
  
    # process functions
    df_agg = combine(gd, :value => (x -> round(sum(x), digits=7)) => :value_prop, 
                         :value => (x -> iszero(sum(x)) ? 0.0 : processing1(x,n)) => :value)
  
    # Take timestep maxmin so all levels have same length
    minmax_timestep = @pipe df_agg |>
      unique(_, :value) |>
      groupby(_, :L) |>
      combine(_,:timestep => maximum=> :timestep)[!, :timestep] |>
      minimum(_)

    df_agg = filter(:timestep => x -> x < minmax_timestep, df_agg)

    df_agg[!, :name] .= p_str

    push!(dfs, df_agg)
    ProgressMeter.next!(p)
  end
  
  all_dfs = vcat(dfs...)

  # Write lookup for rowids -> all_dfs.name to disk
  names_params = unique(all_dfs.name)
  row_ids = Int32.(1:length(names_params))
   
  writefile(lookup_f, (row_id=row_ids, param_str=names_params))
  
  # lookup
  lookup_name = Dict()
  [get!(lookup_name, n, row_id) for (n, row_id) in zip(names_params, row_ids)];

  # Write output to disk
  all_dfs.row_id = [lookup_name[n] for n in all_dfs.name]
  select!(all_dfs, Not(:name))
  all_dfs.L = Int8.(all_dfs.L)
  all_dfs.timestep = Int32.(all_dfs.timestep)
  all_dfs.value = round.(all_dfs.value, digits=4)
  all_dfs.value_prop = round.(all_dfs.value_prop, digits=4)

  # writefile(out_f, all_dfs)
  writefile(out_f, all_dfs)
end

main()
