using Pkg; Pkg.activate("../..")
using OrdinaryDiffEq, Plots, Distributions, RecursiveArrayTools, JSON

function initialize_u0(;N::Int=20)
  N_plus_1 = N + 1 # add column for zeroth case
  G = zeros(N_plus_1, N_plus_1)
  for i=1:N_plus_1, j=1:N_plus_1
    G[i,j] = 1/(N_plus_1*N_plus_1)
  end
  return ArrayPartition(Tuple([G[n,:] for n=1:N_plus_1]))
end

μ  = 0.001   # inflow new students-non coders
νₙ = 0.01    # death rate non-coders
νₚ = 0.05    # death rate coders
α  = 0.01    # benefits non coders
β  = 0.1     # benefits coders
p  = [μ, νₙ, νₚ, α, β]

u₀ = initialize_u0(N=3)
tspan = (0., 160.)

c(n, i) = 0.95 * exp(-i / n)              # cost function
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
    # println("2: $(νₙ*non_coder*G.x[n][i])")
    du.x[n][i] -= νₙ*non_coder*G.x[n][i]
    # println("3: $(νₚ*coder*G.x[n][i])")
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

prob = ODEProblem(life_cycle_research_groups!, u₀, tspan, p)
sol = solve(prob, Tsit5(), saveat=1, reltol=1e-8, abstol=1e-8)

N = length(sol[1].x)
P = length(sol[1].x[1]) # i programmers ~ i adopters
maxsize = (N-1)+(P-1)
I = zeros(maxsize, length(sol.t))
weighted_avg = zeros(maxsize, length(sol.t))
size_dis = zeros(maxsize, length(sol.t))


dist_gsize = zeros(maxsize, 160)
sommeGNI = zeros(maxsize, length(sol.t))
for t in 1:160
  for n in 1:N
    
    G_nil = sol[t].x[n] # sol of group with n non-coders at time t

    for i in 1:P

      coder, noncoder = i-1, n-1 
      
      gsize = coder+noncoder
      
      gsize > 0 && ( sommeGNI[gsize, t] += (coder+noncoder)*G_nil[i] )
      
      gsize > 0 && ( dist_gsize[gsize, t] = G_nil[i] )   # sol of group with n non-coders and i coders at time t, e.g.
                                                         # proportion of groups with n non-coders and i coders.

      gsize > 0 && ( weighted_avg[gsize, t] += (coder/gsize)*G_nil[i] )
      # gsize > 0 && ( println( "t=$(t), n=$(N),i=$(P) == $(G_nil[i])" ) )
      gsize > 0 && ( size_dis[gsize, t] += G_nil[i] )

    end
  end
  
  for gsize in 1:maxsize
      # println("$(gsize),  $(weighted_avg[gsize, t]), $(size_dis[gsize, t])")
      I[gsize,t] = weighted_avg[gsize, t]/size_dis[gsize, t]
  end
end

# round.(dist_gsize[:, 10].*100, digits=3) # at time 10, 31.9% of groups had grsize = 0
                                           # and 1.3% had groups of size 6

# for i=1:160
#   dist_gsize[:,i] = dist_gsize[:,i] / sum(dist_gsize[:,i])
# end

# scatter(dist_gsize[1,1:100], legendtitle="grsize", legend=:outertopright, label="1")
# for i=2:6
#   scatter!(dist_gsize[i,1:100], label="$(i)")
# end
# vline!([10], label="")
# xlabel!("time")
# ylabel!("Fraction of gsize (%)")


heatmap(sommeGNI)
xlabel!("time")
ylabel!("grsize")
title!("Couleur: (coder+noncoder)*G_nil[i] ")

# if isfile("data.json")
#     old_run = JSON.parsefile("data.json")
# end

new_run = []
for n in 1:N
  tmp_dat = [Dict("N" => n, "value" => I[n, i], "timesteps" => i) for i in 1:1000]
  x = round.([i["value"] for i in tmp_dat], digits=5) # To lighten the output file we only keep unique `y` values.
  # idx = unique(z -> x[z], 1:length(x))
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

keys(new_run)