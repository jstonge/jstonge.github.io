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
      help = "Imitation rate from non-adopter to adopter β"
      "-g"
      arg_type = Float64
      default = 1.
      help = "Imitation rate from adopter to non-adopter γ"
      "-r"
      arg_type = Float64
      default = 0.1
      help = "Global behavioral imitation ρ (allows a behaviour to spread between groups)"
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
      help = "Endogenous rate of institutional change μ"
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
          write(io, "$(t) $(ℓ) $(round(val, 10))\n")
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

f(x; a1=2.2) = 1 / (1 + exp(-a1*x)) # cost-benefits for individuals
h(x; a2=0.3) = (1 - exp(-a2*x)) / (1 - exp(-a2)) # dependency of synergy on institutional level

function source_sink3!(du, u, p, t)
  G, L, n = u, length(u.x), length(u.x[1])
  β, γ, ρ, b, c, μ = p
  Z, pop, R = zeros(L), zeros(L), 0.

  # Calculate mean-field coupling and observed fitness landscape
  # In the following, the functions g (cost-benefits for groups) and g̃ (fitness function) are taken equal to function f. The three have similar properties.
    for ℓ in 1:L
      n_adopt = collect(0:(n-1))
      Z[ℓ]    = sum(f.(b*n_adopt .- c*(ℓ-1)) .* G.x[ℓ])
      pop[ℓ]  = sum(G.x[ℓ])
      R      += sum(n_adopt .* G.x[ℓ]) # Global diffusion
      pop[ℓ] > 0.0 && ( Z[ℓ] /= pop[ℓ] )
  end

    for ℓ = 1:L, i = 1:n
      n_adopt, gr_size = i-1, n-1
      # Inndividual selection process
      du.x[ℓ][i] = -n_adopt*f(1-h(ℓ))*G.x[ℓ][i] - (gr_size-n_adopt)*f(h(ℓ)-1)*G.x[ℓ][i]
      du.x[ℓ][i] += -n_adopt*(gr_size-n_adopt)*(β+γ)*G.x[ℓ][i] - ρ*(gr_size-n_adopt)*β*R*G.x[ℓ][i] - ρ*n_adopt*γ*(gr_size-R)*G.x[ℓ][i]
      n_adopt > 0 && ( du.x[ℓ][i] += (gr_size-n_adopt+1)*f(h(ℓ)-1)*G.x[ℓ][i-1] + β*(n_adopt-1+ρ*R)*(gr_size-n_adopt+1)*G.x[ℓ][i-1] )
      n_adopt < gr_size && ( du.x[ℓ][i] += (n_adopt+1)*f(1-h(ℓ))*G.x[ℓ][i+1] + γ*(gr_size-n_adopt-1+ρ*(gr_size-R))*(n_adopt+1)*G.x[ℓ][i+1] )
      # Group selection process
      ℓ > 1 && ( du.x[ℓ][i] += f(b*n_adopt-c*(ℓ-1))*(μ+ρ*Z[ℓ]/Z[ℓ-1])*G.x[ℓ-1][i] - (μ*f(c*(ℓ-1)-b*n_adopt)+ρ*f(b*n_adopt-c*(ℓ-2))*Z[ℓ-1]/Z[ℓ])*G.x[ℓ][i] )
      ℓ < L && ( du.x[ℓ][i] += (μ*f(c*ℓ-b*n_adopt)+ρ*f(b*n_adopt-c*(ℓ-1))*Z[ℓ]/Z[ℓ+1])*G.x[ℓ+1][i] - f(b*n_adopt-c*ℓ)*(μ+ρ*Z[ℓ+1]/Z[ℓ])*G.x[ℓ][i] )
    end
end

function run_source_sink3(p)
  n, M = 20, 1000
  u₀ = initialize_u0(n=n, L=6, M=M, p=0.01)
  tspan = (1.0, 4000)

  # Solve problem
  prob = ODEProblem(source_sink3!, u₀, tspan, p)
  return solve(prob, DP5(), saveat = 1., reltol=1e-8, abstol=1e-8)
end

function main()
  # β, γ, ρ, b, c, μ = 0.07, 0.5, 1, 0.1, 0.18, 1.05, 0.2
  args = parse_commandline()

  modelname = "sourcesink3"

  # Init
  if isnothing(args["db"])
    β = args["beta"]
    γ = args["g"]
    ρ = args["r"]
    b = args["b"]
    c = args["c"]
    μ = args["m"]

    #!TODO: don't forget to change run_source_sink
    p = [β, γ, ρ, b, c, μ]
    println(p)
    sol = run_source_sink3(p)
    write_sol2txt("$(args["o"])$(modelname)_$(join(p, "_")).txt", sol)
  else
    db = SQLite.DB(args["db"])
    con = DBInterface.execute(db, """SELECT * from $(modelname) LIMIT $(args["L"]) OFFSET $(args["O"])""") |> DataFrame

    for row in eachrow(con)

      β = row["beta"]
      γ = row["gamma"]
      ρ = row["rho"]
      b = row["b"]
      c = row["cost"]
      μ = row["mu"]

      #!TODO: don't forget to change run_source_sink
      p = [β, γ, ρ, b, c, μ]
      sol = run_source_sink3(p)
      write_sol2txt("$(args["o"])/$(modelname)_$(join(p, "_")).txt", sol)
    end
  end
end

main()


# prototyping -------------------------------------------------------------------------------

# β, γ, ρ, b, c, μ = 0.1, 0.1, 0.2, 0, 2, 0.1
# p  = [β, γ, ρ, b, c, μ]

# n, M = 20, 1000
# u₀ = initialize_u0(n=n, L=6, M=M, p=0.01)
# tspan = (0., 4000.)

# prob = ODEProblem(source_sink3!, u₀, tspan, p)
# sol = solve(prob, DP5(), saveat=1, reltol=1e-8, abstol=1e-8)


# function from_file()
#   sol = CSV.read(".sourcesink3_0.1_0.1_0.2_0.0_2.0_0.1.txt", DataFrame; header=["timestep", "L", "value"])
#   L = 6
#   inst_level = Dict()
#   lower_limit = 1
#   upper_limit = 21
#   for t=1:4000
#     for ℓ=1:L
#       myrange = UnitRange(lower_limit:upper_limit)
#       n = length(sol.value[myrange])
#       x = sol.value[myrange]
#       out = sum((collect(0:(n-1)) / n) .* x) / sum(x)
#       if haskey(inst_level, ℓ)
#         inst_level[ℓ] = [inst_level[ℓ]; out]
#       else
#         inst_level[ℓ] = out
#       end

#       lower_limit += 21
#       upper_limit += 21

#     end
#   end
#   return inst_level
# end

# inst_level = from_file()

# when sol is available
# L = length(sol.u[1].x)
# for ℓ=1:L
#   values = []
#   for t=1:4000
#     n = length(sol.u[t].x[ℓ])
#     x = sol.u[t].x[ℓ]
#     out = sum((collect(0:(n-1)) / n) .* x) / sum(x)
#     push!(values, out)
#   end
#   inst_level[ℓ] = values
# end

# scatter(1:4000, [inst_level[i] for i in 1:L], xaxis = :log, legendtitle="grsize", 
#       legend=:outertopright, labels=["0" "1" "2" "3" "4" "5"], palette = palette(:Reds)[2:7],
#       markerstrokewidth=0, markersize = 3.)