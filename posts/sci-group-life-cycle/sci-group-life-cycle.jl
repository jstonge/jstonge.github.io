using OrdinaryDiffEq, Plots, Distributions, RecursiveArrayTools, JSON

function initialize_u0(;N::Int=20, M::Int=20, p::Float64=0.01)

  function make_draw()
    G = zeros(N+1, N+1)
    
    for _ in 1:M
      n = rand(1:N+1)
      i = sum(collect(rand(Binomial(1, p), N+1))[1]) # how many total programmers?
      G[n, i+1] += 1 # everytime combination G[ℓ,i], count +1
    end
    return G
  end

  G = make_draw()

  # making sure we start with zero programmers but we have varying group size
  while sum(G[:,2:N+1]) != 0
    G = make_draw()
  end

  G = G ./ M # normalized by tot number of groups

  hm = heatmap(G, xlabel="# non-coders", ylabel="# coders", title="Frac groups with p coders and n non-coder", c = :heat)

  # ArrayPartition are nice because we can access the level such as G.x[ℓ][i]
  return ArrayPartition(Tuple([G[n,:] for n=1:N+1])), hm
end

μ  = 0.001   # inflow new students-non coders
νₙ = 0.01    # death rate non-coders
νₚ = 0.05    # death rate coders
α  = 0.01    # benefits non coders
β  = 0.1     # benefits coders
p  = [μ, νₙ, νₚ, α, β]

n,M = 2, 1000
u₀, hm = initialize_u0(N=n, M=M, p=0.01)
tspan = (0., 1000.)

u₀ = zeros(n+1, n+1)
for i=1:3, j=1:3
  u₀[i,j] = 1/9
end
u₀ = ArrayPartition(Tuple([u₀[n,:] for n=1:3]))

length(first(u₀.x))
length(u₀.x)

# hm # show p_n init mat

c(n, i) = 0.95 * exp(-i / n)              # cost function
τ(n, i, α, β) = n*exp(-α + β*(1 - c(n, i))) # group benefits

function life_cycle_research_groups!(du, u, p, t)

  G, N, P = u, length(u.x), length(first(u₀.x)) # Note that there can be no coders but not non-coders
  μ, νₙ, νₚ, α, β = p
  
  for n=1:N, i=1:P
    
    # println("n:$(n), i:$(i), G.x[n][i]:$(G.x[n][i])")
    
    coder, non_coder = i-1, n-1   # we distinguish indices from actual values.
    
    du.x[n][i] = 0

    non_coder > 0 && println("1: $(μ*(G.x[n-1][i]))")
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

τ(1, 1, α, β)*2.9548258834291983

# sol[1].x[1]

N = length(sol[1].x)
P = length(sol[1].x[1]) # i programmers ~ i adopters
maxsize = (N-1)*(P-1)
I = zeros(maxsize, length(sol.t))
weighted_avg = zeros(maxsize, length(sol.t))
size_dis = zeros(maxsize, length(sol.t))

# for t in 1:length(sol.t)
for t in 1:1000
  for n in 1:N
    G_nil = sol[t].x[n]
    for i in 1:P
      coder = i-1
      noncoder = n-1
      gsize = coder+noncoder
      gsize > 0 && ( weighted_avg[gsize, t] += (coder/gsize)*G_nil[i] )
      gsize > 0 && ( size_dis[gsize, t] += G_nil[i] )
    end
  end
  for gsize in 1:maxsize
      I[gsize,t] = weighted_avg[gsize, t]/size_dis[gsize, t]
  end
end

I

heatmap(I)

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

# if isfile("data.json")
#   get!(old_run, join(p, "_"), vcat(new_run...))
#   new_run = old_run
# else
#   new_run = Dict(join(p, "_") => vcat(new_run...))
# end

new_run = Dict(join(p, "_") => vcat(new_run...))

open("data.json", "w") do f
  write(f, JSON.json(new_run))
end
