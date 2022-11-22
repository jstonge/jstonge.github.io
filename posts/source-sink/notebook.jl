### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 973344d6-6039-11ed-0d52-fdda616eb8ee
# ╠═╡ show_logs = false
using Pkg; Pkg.activate(".")

# ╔═╡ 1f62ee0b-2148-4f50-9868-13efa26eaa11
begin
	using RecursiveArrayTools: ArrayPartition # to flatten our arrays
	using JSON, Plots, OrdinaryDiffEq, Distributions, StatsBase, PlutoUI
end

# ╔═╡ 5e8f2365-2347-429f-8a70-a5a4cb5d174d
md"""
## Julia model
"""

# ╔═╡ 20cd6f5e-43ba-4bdc-88ee-91af4c7554b7
function source_sink!(du, u, p, t)
    G, L, n = u, length(u.x), length(first(u.x))
    β, γ, ρ, b, c, μ = p
    Z, pop, R = zeros(L), zeros(L), 0.

    # Calculate mean-field coupling and observed fitness landscape
    for ℓ in 1:L
      n_adopt = collect(0:(n-1))
      Z[ℓ]    = sum(exp.(b*n_adopt .- c*(ℓ-1)) .* G.x[ℓ])
      pop[ℓ]  = sum(G.x[ℓ])
      R      += sum(ρ*n_adopt .* G.x[ℓ])
      pop[ℓ] > 0.0 && ( Z[ℓ] /= pop[ℓ] )
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

# ╔═╡ aee67c51-5240-4230-81c0-3207d4f82905
md"""
## Output 

β: $(@bind β Slider(0.07:0.03:0.21, show_value=true, default=0.07))

γ: $(@bind γ Slider(0.9:0.1:1.1, show_value=true, default=1.))

c: $(@bind c Slider(1.05:0.5:2.05, show_value=true, default=1.05))

b: $(@bind b Slider(0.12:0.05:0.22, show_value=true, default=1.05))

ρ: $(@bind ρ Slider(0.1:0.15:0.40, show_value=true, default=1.05))
"""

# ╔═╡ 89eae328-1148-424b-a431-fffc4f9d371b
md"""
## Appendix
"""

# ╔═╡ 7647ba81-1519-4ebd-9ef9-6951f0157cd0
function initialize_u0(;n::Int=20, L::Int=6, M::Int=20, p::Float64=0.01)
  G = zeros(L, n+1)

  for _ in 1:M
    ℓ = rand(1:L) # pick a level
    i = sum(collect(rand(Binomial(1, p), n))[1]) # how many total adopters?
    G[ℓ, i+1] += 1 # everytime combination G[ℓ,i], count +1
  end

  G = G ./ M # normalized by tot number of groups

  # !TODO: find better way to flatten matrix.
  @assert L == 6 "Number of lvl must equal 6 for now"
  G = ArrayPartition(G[1,:], G[2,:], G[3,:], G[4,:], G[5,:], G[6,:])

  return G
end

# ╔═╡ a7dff3fb-4667-4404-9a6c-ede7768f14ee
begin
	n, M = 20, 1000
	u₀ = initialize_u0(n=n, L=6, M=M, p=0.01)
	μ = 1e-4
	p = [β, γ, ρ, b, c, μ]
	tspan = (1.0, 4000)

	prob = ODEProblem(source_sink!, u₀, tspan, p)
 	sol = solve(prob, DP5(), saveat=1., reltol=1e-8, abstol=1e-8)

	L = length(sol[1].x)
	n_sol = length(sol[1].x[1])
	I = zeros(L, length(sol.t))

	for t in 1:length(sol.t)
	    for ℓ in 1:L
	      G_nil = sol[t].x[ℓ]
	      I[ℓ, t] = sum((collect(0:(n_sol-1)) / n_sol) .* G_nil) / sum(G_nil)
	    end
	end

	p_sol = scatter(I[1,1:3000], xaxis=:log, ylims = (0,1.))
	[scatter!(I[ℓ,1:3000], xaxis=:log, ylims = (0,1.)) for ℓ=2:6]
	
	p_sol
end

# ╔═╡ Cell order:
# ╟─973344d6-6039-11ed-0d52-fdda616eb8ee
# ╠═1f62ee0b-2148-4f50-9868-13efa26eaa11
# ╟─5e8f2365-2347-429f-8a70-a5a4cb5d174d
# ╠═20cd6f5e-43ba-4bdc-88ee-91af4c7554b7
# ╟─aee67c51-5240-4230-81c0-3207d4f82905
# ╟─a7dff3fb-4667-4404-9a6c-ede7768f14ee
# ╟─89eae328-1148-424b-a431-fffc4f9d371b
# ╟─7647ba81-1519-4ebd-9ef9-6951f0157cd0
