using Pkg; Pkg.activate("../../");
using ArgParse, Distributions, StatsBase, OrdinaryDiffEq, RecursiveArrayTools, DataFrames, SQLite

"""
    parse_commandline()

Function for the commandline interface.
"""
function parse_commandline()
  s = ArgParseSettings()

  #!TODO: Update the params and defaults as needed.
  @add_arg_table! s begin
      "--db"
      help = "Use Database to query parameters"
      "-L"
      arg_type = Int
      default = 5
      help = "LIMIT of rows"
      "-O"
      arg_type = Int
      default = 0
      help = "The OFFSET clause after LIMIT specifies how many rows to skip at the beginning of the result set."
      "--beta"
      arg_type = Float64
      default = 0.07
      help = "Spreading rate from non-adopter to adopter beta"
      "-g"
      arg_type = Float64
      default = 1.
      help = "Spreading rate from adopters to non adopters"
      "-r"
      arg_type = Float64
      default = 0.1
      help = "Global behavioral diffusion rho (allows the behaviour to spread between groups)"
      "-b"
      arg_type = Float64
      default = 0.18
      help = "Group benefits b"
      "-c"
      arg_type = Float64
      default = 1.05
      help = "Institutional cost c"
      "-m"
      arg_type = Float64
      default = 1e-4
      help = "Relative rate between individual and group evolution, μ"
      "-o"
      default = "."
      help = "Output file for results"
    end

  return parse_args(s)
end

"""
  write_sol2txt(path, sol)

Function to write solution to textfile. Nothing to do here.
"""
function write_sol2txt(path, sol)
  L = length(sol.u[1].x)
  open(path, "a") do io
    for t=1:length(sol.u), ℓ=1:L
      for val in sol.u[t].x[ℓ]
          write(io, "$(t) $(ℓ) $(round(val, digits=6))\n")
      end
    end
  end
end

"""
  initialize_u0(;n, L, M, p)

Function to initialize the model.
"""
function initialize_u0(;n::Int=20, L::Int=6, M::Int=20, p::Float64=0.01)::ArrayPartition
  G = zeros(L, n+1)
  for _ in 1:M
    ℓ = rand(1:L)
    i = sum(rand(Binomial(1, p), n))
    G[ℓ, i+1] += 1
  end

  G = G ./ M

  return ArrayPartition(Tuple([G[ℓ,:] for ℓ=1:L]))
end

f(x; a1=1) = 1 / (1 + exp(-a1*x)) 
h(x; a2=1) = (1 - exp(-a2*x)) / (1 - exp(-a2))

function source_sink3!(du, u, p, t)
  G, L, n = u, length(u.x), length(u.x[1])
  β, γ, ρ, b, c, μ = p
  Z, pop, R = zeros(L), zeros(L), 0.

  # Calculate mean-field coupling and observed fitness landscape
  for ℓ in 1:L
      n_adopt = collect(0:(n-1))
      Z[ℓ]    = sum(exp.(b*n_adopt .- c*(ℓ-1)) .* G.x[ℓ])
      pop[ℓ]  = sum(G.x[ℓ])
      R      += sum(n_adopt .* G.x[ℓ]) # Global diffusion
      pop[ℓ] > 0.0 && ( Z[ℓ] /= pop[ℓ] )
  end

    for ℓ = 1:L, i = 1:n
      n_adopt, gr_size = i-1, n-1
      # Imitation 
      du.x[ℓ][i] = -γ*n_adopt*(gr_size-n_adopt + ρ*(gr_size-R))*G.x[ℓ][i] - β*(gr_size-n_adopt)*(n_adopt+ρ*R)*G.x[ℓ][i]
      # Cost-benefits
      du.x[ℓ][i] += -n_adopt*f(1-h(ℓ))*G.x[ℓ][i] - (gr_size-n_adopt)*f(h(ℓ)-1)*G.x[ℓ][i]
      # n_adopt > 0 && ( du.x[ℓ][i] += β*(ℓ^-α)*(n_adopt-1+R)*(gr_size-n_adopt+1)*G.x[ℓ][i-1])
      # n_adopt < gr_size && ( du.x[ℓ][i] +=  γ*(n_adopt+1)*G.x[ℓ][i+1] )
      # # Group selection process
      # ℓ > 1 && ( du.x[ℓ][i] += ρ*G.x[ℓ-1][i]*(Z[ℓ] / Z[ℓ-1] + μ) - ρ*G.x[ℓ][i]*(Z[ℓ-1] / Z[ℓ]+μ) )
      # ℓ < L && ( du.x[ℓ][i] += ρ*G.x[ℓ+1][i]*(Z[ℓ] / Z[ℓ+1] + μ) - ρ*G.x[ℓ][i]*(Z[ℓ+1] / Z[ℓ]+μ) )
    end
end

#!TODO: Modify the function where needed. Fct name should correspond to the same number of your model.
function run_source_sink(p)
  n, M = 20, 1000
  u₀ = initialize_u0(n=n, L=6, M=M, p=0.01)
  tspan = (1.0, 4000)

  # Solve problem
  prob = ODEProblem(source_sink3!, u₀, tspan, p)
  return solve(prob, DP5(), saveat = 1., reltol=1e-8, abstol=1e-8)
end


function main()
  # β, α, γ, ρ, b, c = 0.07, 0.5, 1, 0.1, 0.18, 1.05
  args = parse_commandline()

  #TODO: Add your model name here.
  modelname = ""

  # Init
  if isnothing(args["db"])
    #TODO: Update param list
    β = args["beta"]
    α = args["a"]
    γ = args["g"]
    ρ = args["r"]
    b = args["b"]
    c = args["c"]
    μ = args["m"]

    p = [β, α, γ, ρ, b, c, μ]
    #TODO: Update function name
    sol = run_source_sink(p)
    write_sol2txt("$(args["o"])$(modelname)_$(join(p, "_")).txt", sol)
  else
    db = SQLite.DB(args["db"])
    con = DBInterface.execute(db, """SELECT * from $(modelname) LIMIT $(args["L"]) OFFSET $(args["O"])""") |> DataFrame

    for row in eachrow(con)

      β = row["beta"]
      α = row["alpha"]
      γ = row["gamma"]
      ρ = row["rho"]
      b = row["b"]
      c = row["cost"]
      μ = row["mu"]

      p = [β, α, γ, ρ, b, c, μ]
      sol = run_source_sink(p)
      write_sol2txt("$(args["o"])/$(modelname)_$(join(p, "_")).txt", sol)
    end
  end
end

main()


# prototyping -------------------------------------------------------------------------------

# When you are developping your model

β, γ, ρ, b, c, μ = 0.07, 1., 0.1, 0.18, 1.05, 0.0001
p  = [β, γ, ρ, b, c, μ]

n, M = 20, 1000
u₀ = initialize_u0(n=n, L=6, M=M, p=0.01)
tspan = (0., 4000.)

prob = ODEProblem(source_sink3!, u₀, tspan, p)
sol = solve(prob, DP5(), saveat=1, reltol=1e-8, abstol=1e-8)

L = length(sol.u[1].x)
inst_level = Dict()
for ℓ=1:L
  values = []
  for t=1:200
    n = length(sol.u[t].x[ℓ])
    x = sol.u[t].x[ℓ]
    out = sum((collect(0:(n-1)) / n) .* x) / sum(x)
    push!(values, out)
  end
  inst_level[ℓ] = values
end

p =scatter(1:200, inst_level[1], xaxis = :log, legendtitle="grsize", 
        legend=:outertopright, label="0")
for i=2:6
  scatter!(1:200, inst_level[i], label="$(i-1)")
end

p
ylims!(0,0.3)