using Pkg; Pkg.activate("../..");
using OrdinaryDiffEq, Plots

function plot_lv_solution(dims=2)
  ys = sum(sol.u[2], dims=dims)
  ys = dims==1 ? transpose(ys) : ys
  xs = 1:length(ys)
  p=plot(xs, ys, marker = 2, label="t=2")
  for t=5:3:22
    ys = sum(sol.u[t], dims=dims)
    ys = dims==1 ? transpose(ys) : ys
    plot!(xs, ys, marker = 2, label="t=$(t)")
  end
  type_ent = dims == 1 ? "non-programmers" : "programmers"
  xlabel!("Number of $(type_ent)")
  ylabel!("Occupation number")
  p    
end

function init(;N::Int=20)
  N_plus_1 = 20 + 1 # add column for zeroth case
  G = zeros(N_plus_1, N_plus_1)
  for i=1:N_plus_1, j=1:N_plus_1
    G[i,j] = 1/(N_plus_1*N_plus_1)
  end
  return G
end

function life_cycle_research_groups!(du, u, p, t)
    """
    Similar to Lotka-Volterra model

    Params
    ======
     - μ: const inflow of non-programmers
     - νₚ: death rate programmers
     - νₙ: death rate non-programmers
     - β: conversion rate from non-programmers to programmers

    Note
    ====
     - νₚ > νₙ : programmers get out of academia at a faster rate.
     - Non-programmer reproduction depends on a sigmoid
     - Academic deaht of non-programmers can happen on its own.
     - death rate of programmers is the same
     - No cost-benefits of learning to code, just a conversion rate.
    """
    N1, N2 = size(u)
    K = size(u, 1)
    α, νₙ, νₚ, β = p
  
    for n₁=1:N1, n₂=1:N2
      non_coder, coder = n₁-1, n₂-1   # we distinguish indices from actual values.
      du[n₁,n₂] = 0
  
      non_coder > 0 && ( du[n₁,n₂] += α*non_coder*(K-non_coder)*u[n₁-1,n₂]/K )   # non-prog birth input
      du[n₁,n₂] -= α*non_coder*(K-non_coder)*u[n₁,n₂]/K                          # non-prog birth output
      
      n₁ < N1 && ( du[n₁,n₂] += νₙ*(non_coder+1)*u[n₁+1,n₂] ) # non-prog death input
      du[n₁,n₂] -= νₚ*coder*u[n₁,n₂]                          # non-prog death output
     
      n₂ < N2 && ( du[n₁,n₂] += νₚ*(coder+1)*u[n₁,n₂+1] ) # prog death input
      du[n₁,n₂] -= νₙ*non_coder*u[n₁,n₂]                  # prog death output

      (n₁ < N1 && coder > 0) && ( du[n₁,n₂] += β*(non_coder+1)*(coder-1)*u[n₁+1,n₂-1] ) # nonProg2prog input
      n₂ < N2 && ( du[n₁,n₂] -= β*non_coder*coder*u[n₁,n₂] ) # nonProg2prog output
    end
end

α  = 0.9   # inflow new students-non coders
νₙ = 0.01   # death rate non-coders
νₚ = 0.01    # death rate coders
β  = 0.03     # conversion rates non-prog to prog
p  = [α, νₙ, νₚ, β]

tspan = (0., 100)
nb_of_states = 25

u₀ = zeros(nb_of_states,nb_of_states)
u₀[12,8] = 1

# u₀ = init(N=20)

prob = ODEProblem(life_cycle_research_groups!, u₀, tspan, p)
sol = solve(prob, Tsit5(), saveat=1)

p1 = plot_lv_solution(1)
title!(param_str)
p2 = plot_lv_solution(2)
p3 = heatmap(log10.(sol.u[2]))
plot(p1, p3, p2, layout=(2,2))

param_str = join(["$(x)$(y)" for (x,y) in zip(p, ["α", "νₙ", "νₚ", "β"])], "_")
savefig("figs/prototypes/m2_$(param_str)_init2.png")
