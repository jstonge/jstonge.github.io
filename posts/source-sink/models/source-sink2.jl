using Pkg; Pkg.activate("../../");
using ArgParse, Distributions, StatsBase, OrdinaryDiffEq, JSON, RecursiveArrayTools

function ArgParse.parse_item(::Type{StepRangeLen}, x::AbstractString)
  start, step, stop = split(x, ":")
  start = parse(Float64, start)
  stop = parse(Float64, stop)
  step = parse(Float64, step)
  return range(start,stop,step=step)

end

function initialize_u0(;n::Int=20, L::Int=6, M::Int=20, p::Float64=0.01)
    #!TODO
end

function source_sink2!(du, u, p, t)
    #!TODO
end

function parse_commandline()
  s = ArgParseSettings()

  @add_arg_table! s begin
      "--param1"
      arg_type = StepRangeLen
      default = 0.07:0.05:0.22 # initial value grid
      help = "Explain param1"
      "--param2"
      arg_type = StepRangeLen
      default = 0.9:0.1:1.1
      help = "Explain param1"
    end

  return parse_args(s)
end


# ------------------------------ run individual ------------------------------ #


function main()
  
  args = parse_commandline()

  for β=args["beta"], γ=args["gamma"]
    
    # Init
    #!TODO

    # Load old run
    # old_run = isfile("../data.json") ? JSON.parsefile("../data.json") : Dict()
    
    # Check if already done
    if haskey(old_run, join(p, "_")) && length(old_run[join(p, "_")]) >= 20
      println("Already done.")
    else
      # tspan = (1.0, 4000) # can be modified
      
      # Solve problem
      # prob = ODEProblem(source_sink2!, u₀, tspan, p)
      # sol = solve(prob, DP5(), saveat = 1., reltol=1e-8, abstol=1e-8)
      
      # Get summary statistics
      #!TODO

      # Wrangle summary stats into list of dicts
      #!TODO
      
      # Update old runs with new run
      if isfile("../data2.json")
        get!(old_run, join(p, "_"), vcat(new_run...))    # If already exists, update the old run Dict()
        new_run = old_run
      else
        new_run = Dict(join(p, "_") => vcat(new_run...)) # new Dict()
      end
      
      # Write to disk
      open("../data2.json", "w") do f 
        write(f, JSON.json(new_run))
      end

  end

  end
end

main()