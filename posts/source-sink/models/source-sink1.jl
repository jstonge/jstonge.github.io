using Distributions, StatsBase, OrdinaryDiffEq, Plots, JSON, RecursiveArrayTools

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

"""
source_sink

Params:
=======
- β: spreading rate from non-adopter to adopter
- γ: recovery at rate gamma (loose behavioral trait)
- ρ: global behavioral diffusion (allows the behaviour to spread between groups)
- μ: constant rate of transition regardless of fitness
- b: benefits
- c: institutional cost
"""
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


function main(p)
  n, M = 20, 1000
  u₀ = initialize_u0(n=n, L=6, M=M, p=0.01)

  # Load old run
  if isfile(".../data.json")
    old_run = JSON.parsefile("../data.json")
  end
  
  if haskey(old_run, join(p, "_")) == false
    
    tspan = (1.0, 4000)

    prob = ODEProblem(source_sink!, u₀, tspan, p)
    sol = solve(prob, DP5(), saveat = 1., reltol=1e-8, abstol=1e-8)
  
    L = length(sol[1].x)
    n = length(sol[1].x[1])
    I = zeros(L, length(sol.t))
    
    for t in 1:length(sol.t)
      for ℓ in 1:L
        G_nil = sol[t].x[ℓ]
        I[ℓ, t] = sum((collect(0:(n-1)) / n) .* G_nil) / sum(G_nil)
      end
    end
    
    function update_runs()     
      new_run = []
      for ℓ in 1:L
        tmp_dat = [Dict("L" => ℓ, "value" => I[ℓ, i], "timesteps" => i) for i in 1:3000]
        x = round.([i["value"] for i in tmp_dat], digits=5) # To lighten the output file we only keep unique `y` values. 
        idx = unique(z -> x[z], 1:length(x))
        push!(new_run, tmp_dat[idx] )
      end
  
      if isfile("../data.json")
        get!(old_run, join(p, "_"), vcat(new_run...))    # If already exists, update the old run Dict()
        new_run = old_run
      else
        new_run = Dict(join(p, "_") => vcat(new_run...)) # new Dict()
      end
      
      open("../data.json", "w") do f 
        write(f, JSON.json(new_run))
      end
    end
    
    update_runs()
  
  end

end

# ------------------------------ run individual ------------------------------ #

β, γ, ρ, b, c = 0.07, 1, 0.1, 0.18, 1.05 # base case from the paper
γ = 1.1
μ = 1e-4
params = [β, γ, ρ, b, c, μ]

main(params)

# --------------------------------- run batch -------------------------------- #

function get_tot_it()
  tot_it = 0
  for β=0.07:0.05:0.22, c = .55:0.5:2.05, b=0.12:0.05:0.22, ρ=0.1:0.15:0.40
    tot_it += 1
  end
  return tot_it
end

tot_it = get_tot_it()

counter = 0 
for β=0.07:0.05:0.22, c = .55:0.5:2.05, b=0.12:0.05:0.22, ρ=0.1:0.15:0.40
  # γ = 0.9:0.1:1.0
  γ = 1.1
  μ = 1e-4
  params = [β, γ, ρ, b, c, μ]

  main(params)
 
  counter += 1
  println((counter  / tot_it) * 100)
end

