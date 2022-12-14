using Pkg; Pkg.activate("../../");
using ArgParse, Distributions, StatsBase, OrdinaryDiffEq, JSON, RecursiveArrayTools

function ArgParse.parse_item(::Type{StepRangeLen}, x::AbstractString)
  start, step, stop = split(x, ":")
  start = parse(Float64, start)
  stop =  parse(Float64, stop)
  step =  parse(Float64, step)
  return range(start,stop,step=step)

end

function initialize_u0(;n::Int=20, L::Int=6, M::Int=20, p::Float64=0.01)
  G = zeros(L, n+1)
  for _ in 1:M
    ℓ = rand(1:L) # pick a level
    i = sum(collect(rand(Binomial(1, p), n))[1]) # how many total adopters?
    G[ℓ, i+1] += 1 # everytime combination G[ℓ,i], count +1
  end

  G = G ./ M # normalized by tot number of groups
  
  # ArrayPartition are nice because we can still access the level such as x[ℓ][i]
  return ArrayPartition(Tuple([G[ℓ,:] for ℓ=1:L]))
end

function source_sink!(du, u, p, t)
    G, L, n = u, length(u.x), length(first(u.x))
    β, γ, ρ, b, c, μ = p
    Z, pop, R = zeros(L), zeros(L), 0.

    # Calculate mean-field coupling and observed fitness landscape
    for ℓ in 1:L
      n_adopt = collect(0:(n-1))
      Z[ℓ]    = sum(exp.(b*n_adopt .- c*(ℓ-1)) .* G.x[ℓ]) 
      pop[ℓ]  = sum(G.x[ℓ])
      R      += sum(ρ * n_adopt .* G.x[ℓ]) # R = ρ ∑_{i',ℓ} i'G_{i',ℓ}                                [global difusion]
      pop[ℓ] > 0.0 && ( Z[ℓ] /= pop[ℓ] ) # Z_{ℓ} = ∑_{i'} exp(b*i' - cℓ) G_{i',ℓ} / ∑_{i'} G_{i',ℓ} [perceived fitness]
    end

    for ℓ = 1:L, i = 1:n
      n_adopt, gr_size = i-1, n-1

      # Diffusion events
      du.x[ℓ][i] = -γ*n_adopt*G.x[ℓ][i] - (ℓ-1)*β*(n_adopt+R)*(gr_size-n_adopt)*G.x[ℓ][i]

      n_adopt > 0 && ( du.x[ℓ][i] += β*(ℓ-1)*(n_adopt-1+R)*(gr_size-n_adopt+1)*G.x[ℓ][i-1])
      n_adopt < gr_size && ( du.x[ℓ][i] +=  γ*(n_adopt+1)*G.x[ℓ][i+1] )

      # Group selection process
      ℓ > 1 && ( du.x[ℓ][i] += ρ*G.x[ℓ-1][i]*(Z[ℓ] / Z[ℓ-1] + μ) - ρ*G.x[ℓ][i]*(Z[ℓ-1] / Z[ℓ]+μ) )
      ℓ < L && ( du.x[ℓ][i] += ρ*G.x[ℓ+1][i]*(Z[ℓ] / Z[ℓ+1] + μ) - ρ*G.x[ℓ][i]*(Z[ℓ+1] / Z[ℓ]+μ) )
    end
end

function parse_commandline()
  s = ArgParseSettings()

  @add_arg_table! s begin
      "--beta"
      arg_type = StepRangeLen
      default = 0.07:0.05:0.22
      help = "Spreading rate from non-adopter to adopter"
      "--gamma"
      arg_type = StepRangeLen
      default = 0.9:0.1:1.1
      help = "Recovery rate, i.e. rate at which adopters loose behavioral trait"
      "--rho"
      arg_type = StepRangeLen
      default = 0.1:0.15:0.40
      help = "Global behavioral diffusion (allows the behaviour to spread between groups)"
      "-b"
      arg_type = StepRangeLen
      default = 0.12:0.05:0.22
      help = "Group benefits"
      "-c"
      arg_type = StepRangeLen
      default = .55:0.5:2.05
      help = "Institutional cost"
    end

  return parse_args(s)
end


# ------------------------------ run individual ------------------------------ #


function main()
  
  args = parse_commandline()
  for β=args["beta"], γ=args["gamma"], ρ=args["rho"], b=args["b"], c=args["c"]
    
    # Init
    μ = 1e-4 #  constant rate of transition regardless of fitness
    p = [ρ, β, γ, b, c, μ]
    
    println("Doing $p")
    
    n, M = 20, 1000
    u₀ = initialize_u0(n=n, L=6, M=M, p=0.01)
    
    # Load old run
    old_run = isfile("../data.json") ? JSON.parsefile("../data.json") : Dict()
    
    # Check if already done
    if haskey(old_run, join(p, "_")) && length(old_run[join(p, "_")]) >= 20
      println("Already done.")
    else
      tspan = (1.0, 4000)
      
      # Solve problem
      prob = ODEProblem(source_sink!, u₀, tspan, p)
      sol = solve(prob, DP5(), saveat = 1., reltol=1e-8, abstol=1e-8)
      
      # Get summary statistics
      L = length(sol[1].x)
      n = length(sol[1].x[1])
      I = zeros(L, length(sol.t))
      for t in 1:length(sol.t)
        for ℓ in 1:L
          G_nil = sol[t].x[ℓ]
          I[ℓ, t] = sum((collect(0:(n-1)) / n) .* G_nil) / sum(G_nil)
        end
      end
      
      # Wrangle summary stats into list of dicts where each dict is a level
      new_run = []
      for ℓ in 1:L
        tmp_dat = [Dict("L" => ℓ, "value" => I[ℓ, i], "timesteps" => i) for i in 1:3000]
        x = round.([i["value"] for i in tmp_dat], digits=5) # To lighten the output file we only keep unique `y` values. 
        idx = unique(z -> x[z], 1:length(x))
        push!(new_run, tmp_dat[idx] )
      end
      
      # Update old runs with new run
      if isfile("../data.json")
        get!(old_run, join(p, "_"), vcat(new_run...))    # If already exists, update the old run Dict()
        new_run = old_run
      else
        new_run = Dict(join(p, "_") => vcat(new_run...)) # new Dict()
      end
      
      # Write to disk
      open("../data.json", "w") do f 
        write(f, JSON.json(new_run))
      end

  end

  end
end

main()