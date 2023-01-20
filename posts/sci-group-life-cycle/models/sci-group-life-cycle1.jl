using Pkg; Pkg.activate("../../");
using ArgParse, Distributions, StatsBase, OrdinaryDiffEq, RecursiveArrayTools, DataFrames, SQLite

function parse_commandline()
  s = ArgParseSettings()

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
      "-m"
      arg_type = Float64
      default = 1e-4
      help = "inflow new students-non coders"
      "--nn"
      arg_type = Float64
      default = 1.
      help = "death rate non-coders"
      "--np"
      arg_type = Float64
      default = 0.1
      help = "death rate coders"
      "-b"
      arg_type = Float64
      default = 0.18
      help = "benefits non coders"
      "-a"
      arg_type = Float64
      default = 0.5
      help = "benefits coders"
      "-o"
      default = "."
      help = "Output file for results"
    end

  return parse_args(s)
end

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

function initialize_u0(;N::Int=20)
  N_plus_1 = N + 1 # add column for zeroth case
  G = zeros(N_plus_1, N_plus_1)
  for i=1:N_plus_1, j=1:N_plus_1
    G[i,j] = 1/(N_plus_1*N_plus_1)
  end
  return ArrayPartition(Tuple([G[n,:] for n=1:N_plus_1]))
end

c(n, i) = 0.95 * exp(-i / n)                # cost function
τ(n, i, α, β) = n*exp(-α + β*(1 - c(n, i))) # group benefits

function life_cycle_research_groups!(du, u, p, t)

  G, N, P = u, length(u.x), length(first(u₀.x)) # Note that there can be no coders but not non-coders
  μ, νₙ, νₚ, α, β = p
  for n=1:N, i=1:P
    println("n:$(n), i:$(i), G.x[n][i]:$(G.x[n][i])")
    coder, non_coder = i-1, n-1   # we distinguish indices from actual values.
    
    du.x[n][i] = 0

    non_coder > 0 && ( du.x[n][i] += μ*(G.x[n-1][i]) )                # 1st term
    
    # for everybody
    println("2: $(νₙ*non_coder*G.x[n][i])")
    du.x[n][i] -= νₙ*non_coder*G.x[n][i]
    println("3: $(νₚ*coder*G.x[n][i])")
    du.x[n][i] -= νₚ*coder*G.x[n][i]

    # upper boxes don't exist 
    if i < P
      # non_coder > 0 && println("4: $(τ(non_coder, coder, α, β)*G.x[n][i] )")
      # We don't want to pass non_coder = 0 to τ()
      non_coder > 0 && ( du.x[n][i] -= τ(non_coder, coder, α, β)*G.x[n][i] )               # 4th term
      # println("5: $(νₚ*(coder+1)*G.x[n][i+1])")
      du.x[n][i] += νₚ*(coder+1)*G.x[n][i+1]  # 5th term
    end
    
    # the bottom boxes don't exist
    if n < N
      # println("6: $(μ*G.x[n][i])")
      du.x[n][i] -= μ*G.x[n][i]                                       # 1st term
      du.x[n][i] += τ(non_coder+1, coder, α, β)*(c(non_coder+1, coder))*G.x[n+1][i]     # 6th term
      du.x[n][i] += νₙ*(non_coder+1)*G.x[n+1][i]                                            # 2nd term
      coder > 0 && ( du.x[n][i] += τ(non_coder+1, coder-1, α, β)*(1-c(non_coder+1, coder-1))*G.x[n+1][i-1] ) # 3rd term 
    end
  end
end

function run_life_cycle_research_groups(p)
  u₀ = initialize_u0(N=3)
  tspan = (1.0, 4000)
  
  # Solve problem
  prob = ODEProblem(life_cycle_research_groups!, u₀, tspan, p)
  return solve(prob, DP5(), saveat = 1., reltol=1e-8, abstol=1e-8)
end


function main()
    # β, γ, ρ, b, c, μ = 0.07, 1., 0.1, 0.18, 1.05, 0.0001
    args = parse_commandline()
    
    if isnothing(args["db"])
      μ  = args["mu"] 
      νₙ = args["nn"] 
      νₚ = args["np"] 
      α  = args["a"]  
      β  = args["b"]  
      
      p  = [μ, νₙ, νₚ, α, β]
      sol = run_life_cycle_research_groups(p)
      write_sol2txt("$(args["o"])/sci-group-life-cycle1_$(join(p, "_")).txt", sol) 

    else
      
      db = SQLite.DB(args["db"])
      con = DBInterface.execute(db, """SELECT * from sci-group-life-cycle1 LIMIT $(args["O"]), $(args["L"])""") |> DataFrame
    
      for row in eachrow(con)
        μ  = args["mu"]   # inflow new students-non coders
        νₙ = args["nn"]   # death rate non-coders
        νₚ = args["np"]    # death rate coders
        α  = args["a"]   # benefits non coders
        β  = args["b"]     # benefits coders
        
        p  = [μ, νₙ, νₚ, α, β]
        sol = run_life_cycle_research_groups(p)    
        write_sol2txt("$(args["o"])/sci-group-life-cycle1_$(join(p, "_")).txt", sol) 
      end
    end  
end

main()


# prototyping -------------------------------------------------------------------------------

μ  = 0.001   # inflow new students-non coders
νₙ = 0.01    # death rate non-coders
νₚ = 0.05    # death rate coders
α  = 0.01    # benefits non coders
β  = 0.1     # benefits coders
p  = [μ, νₙ, νₚ, α, β]

u₀ = initialize_u0(N=3)
tspan = (0., 160.)

prob = ODEProblem(life_cycle_research_groups!, u₀, tspan, p)
sol = solve(prob, Tsit5(), saveat=1, reltol=1e-8, abstol=1e-8)

N = length(sol[1].x)
P = length(sol[1].x[1]) # i programmers ~ i adopters
maxsize = (N-1)+(P-1)
I = zeros(maxsize, length(sol.t))
weighted_avg = zeros(maxsize, length(sol.t))
size_dis = zeros(maxsize, length(sol.t))
dist_gsize = zeros(maxsize, length(sol.t))
sommeGNI = zeros(maxsize, length(sol.t))

for t=1:length(sol.t)
  for n=1:N
    for i=1:P
      
      coder, noncoder = i-1, n-1 
      G_nil = sol[t].x[n] # sol of group with n non-coders at time t
      gsize = coder+noncoder
      
      gsize > 0 && ( sommeGNI[gsize, t] += gsize*G_nil[i] )
      gsize > 0 && ( dist_gsize[gsize, t] = G_nil[i] )  
      gsize > 0 && ( weighted_avg[gsize, t] += (coder/gsize)*G_nil[i] )
      gsize > 0 && ( size_dis[gsize, t] += G_nil[i] )

    end
  end
  
  for gsize=1:maxsize
      I[gsize,t] = weighted_avg[gsize, t]/size_dis[gsize, t]
  end
end

# Plotting -------------------------------------------------------------

for i=1:160
  dist_gsize[:,i] = dist_gsize[:,i] / sum(dist_gsize[:,i])
end

scatter(dist_gsize[1,1:100], legendtitle="grsize", legend=:outertopright, label="1")
for i=2:6
  scatter!(dist_gsize[i,1:100], label="$(i)")
end
vline!([10], label="")
xlabel!("time")
ylabel!("Fraction of gsize (%)")


# heatmap(sommeGNI)
# xlabel!("time")
# ylabel!("grsize")
# title!("Couleur: (coder+noncoder)*G_nil[i] ")

# Writing to file -------------------------------------------------------------


old_run = isfile("data.json") ? JSON.parsefile("data.json") : Dict()

new_run = []
for n in 1:N
  tmp_dat = [Dict("N" => n, "value" => I[n, i], "timesteps" => i) for i in 1:length(sol.t)]
  x = round.([i["value"] for i in tmp_dat], digits=5) # To lighten the output file we only keep unique `y` values.
  push!(new_run, tmp_dat )
end

if isfile("data.json")
  get!(old_run, join(p, "_"), vcat(new_run...))
  new_run = old_run
else
  new_run = Dict(join(p, "_") => vcat(new_run...))
end

new_run = Dict("$(N)_$(join(p, "_"))" => vcat(new_run...))

open("data.json", "w") do f
  write(f, JSON.json(new_run))
end

