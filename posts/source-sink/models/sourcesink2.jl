using Pkg; Pkg.activate("../../");
using ArgParse, Distributions, StatsBase, OrdinaryDiffEq, RecursiveArrayTools, JLD2

function ArgParse.parse_item(::Type{StepRangeLen}, x::AbstractString) where { StepRangeLen }
  if tryparse(Float64, x) !== nothing
    start = parse(Float64, x)
    stop = parse(Float64, x)
    step = 1.
  else
    start, step, stop = split(x, ":")
    start = parse(Float64, start)
    stop =  parse(Float64, stop)
    step =  parse(Float64, step)
  end
  return range(start,stop,step=step)
end

function parse_commandline()
  s = ArgParseSettings()

  @add_arg_table! s begin
      "--beta"
      arg_type = StepRangeLen
      default = 0.07:1.:0.07
      help = "Spreading rate from non-adopter to adopter beta"
      "-a"
      arg_type = StepRangeLen
      default = 0.05:1.:0.05
      help = "Negative benefits alpha"
      "-g"
      arg_type = StepRangeLen
      default = 1.:1.:1.
      help = "Recovery rate gamma, i.e. rate at which adopters loose behavioral trait"
      "-r"
      arg_type = StepRangeLen
      default = 0.1:1.:0.1
      help = "Global behavioral diffusion rho (allows the behaviour to spread between groups)"
      "-b"
      arg_type = StepRangeLen
      default = 0.18:1.:0.18
      help = "Group benefits b"
      "-c"
      arg_type = StepRangeLen
      default = 1.05:1.:1.05
      help = "Institutional cost c"
      "-m"
      arg_type = StepRangeLen
      default = 1e-4:1.:1e-4
      help = "Noise u"
      "-o"
      default = "."
      help = "Output file for results"
    end

  return parse_args(s)
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

function source_sink2!(du, u, p, t)
    G, L, n = u, length(u.x), length(first(u.x))
    β, α, γ, ρ, b, c, μ = p
    Z, pop, R = zeros(L), zeros(L), 0.

    # Calculate mean-field coupling and observed fitness landscape
    for ℓ in 1:L
        n_adopt = collect(0:(n-1))
        Z[ℓ]    = sum(exp.(b*n_adopt .- c*(ℓ-1)) .* G.x[ℓ]) 
        pop[ℓ]  = sum(G.x[ℓ])
        R      += sum(ρ * n_adopt .* G.x[ℓ]) 
        pop[ℓ] > 0.0 && ( Z[ℓ] /= pop[ℓ] ) 
      end
      
      for ℓ = 1:L, i = 1:n
        n_adopt, gr_size = i-1, n-1
        # Diffusion events
        du.x[ℓ][i] = -γ*n_adopt*G.x[ℓ][i] - β*(ℓ^-α)*(n_adopt+R)*(gr_size-n_adopt)*G.x[ℓ][i]
        n_adopt > 0 && ( du.x[ℓ][i] += β*(ℓ^-α)*(n_adopt-1+R)*(gr_size-n_adopt+1)*G.x[ℓ][i-1])
        n_adopt < gr_size && ( du.x[ℓ][i] +=  γ*(n_adopt+1)*G.x[ℓ][i+1] )
        # Group selection process
        ℓ > 1 && ( du.x[ℓ][i] += ρ*G.x[ℓ-1][i]*(Z[ℓ] / Z[ℓ-1] + μ) - ρ*G.x[ℓ][i]*(Z[ℓ-1] / Z[ℓ]+μ) )
        ℓ < L && ( du.x[ℓ][i] += ρ*G.x[ℓ+1][i]*(Z[ℓ] / Z[ℓ+1] + μ) - ρ*G.x[ℓ][i]*(Z[ℓ+1] / Z[ℓ]+μ) )
      end
end

function main()
  args = parse_commandline()
  for β=args["beta"], α=args["a"], γ=args["g"], ρ=args["r"], b=args["b"], c=args["c"], μ=args["m"]
    # Init
    # β, α, γ, ρ, b, c = 0.07, 0.5, 1, 0.1, 0.18, 1.05
    p = [β, α, γ, ρ, b, c, μ]
    
    n, M = 20, 1000
    u₀ = initialize_u0(n=n, L=6, M=M, p=0.01)
    tspan = (1.0, 4000)
    
    # Solve problem
    prob = ODEProblem(source_sink2!, u₀, tspan, p)
    sol = solve(prob, DP5(), saveat = 1., reltol=1e-8, abstol=1e-8)

    @save "$(args["o"])/sourcesink2_$(join(p, "_")).jld2" sol
  end
end

main()