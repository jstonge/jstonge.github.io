using Pkg; Pkg.activate("../../");
using ArgParse, Distributions, StatsBase, OrdinaryDiffEq, RecursiveArrayTools, DataFrames, SQLite, Plots, LaTeXStrings

include("helpers.jl")

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
        "--beta"
        arg_type = Float64
        default = 0.07
        help = "Spreading rate from non-adopter to adopter beta"
        "-a"
        arg_type = Float64
        default = 0.5
        help = "Negative benefits alpha"
        "-g"
        arg_type = Float64
        default = 1.
        help = "Recovery rate gamma, i.e. rate at which adopters loose behavioral trait"
        "-r"
        arg_type = Float64
        default = 0.1
        help = "rate between groups to spread the contagion"
        "-e"
        arg_type = Float64
        default = 0.1
        help = "rate between groups to spread the institution level."
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
        help = "Noise u"
        "-o"
        default = "."
        help = "Output file for results"
      end
  
    return parse_args(s)
end

function initialize_u0(;n::Int=20, L::Int=6, M::Int=20, p::Float64=0.01, lvl_1_inf::Bool=false)
    G = zeros(L, n+1)
    for _ in 1:M
      ℓ = rand(1:L) # pick a level
      i = sum(collect(rand(Binomial(1, p), n))) # how many total adopters?
      G[ℓ, i+1] += 1 # everytime combination G[ℓ,i], count +1
    end
  
    G = G ./ M # normalized by tot number of groups
    
    if lvl_1_inf # ("almost" only lowest level populated at t = 0)
      u0 = ArrayPartition(Tuple([G[ℓ,:] for ℓ=1:L]))
      u0.x[1][1] = 0.99993
      u0.x[1][2] = 0.00001
      [u0.x[l] .= 0 for l in 2:L]
      [u0.x[l][1] = 0.00001 for l in 2:L]
      return u0
    else 
      # ArrayPartition are nice because we can still access the level such as x[ℓ][i]
      return ArrayPartition(Tuple([G[ℓ,:] for ℓ=1:L]))
    end
end

function source_sink2!(du, u, p, t)
    G, L, n = u, length(u.x), length(first(u.x))
    β, α, γ, ρ, η, b, c, μ = p
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
        ℓ > 1 && ( du.x[ℓ][i] += η*G.x[ℓ-1][i]*(Z[ℓ] / Z[ℓ-1] + μ) - η*G.x[ℓ][i]*(Z[ℓ-1] / Z[ℓ]+μ) )
        ℓ < L && ( du.x[ℓ][i] += η*G.x[ℓ+1][i]*(Z[ℓ] / Z[ℓ+1] + μ) - η*G.x[ℓ][i]*(Z[ℓ+1] / Z[ℓ]+μ) )
      end
end

function run_source_sink2(p; lvl_1_inf::Bool=false, susc_pop::Bool=true)
  n, M = 20, 1000
  L = 6
  if susc_pop
    u₀ = initialize_u0(n=n, L=L, M=M, p=0.001, lvl_1_inf=lvl_1_inf)
  else
    u₀ = initialize_u0(n=n, L=L, M=M, p=0.3, lvl_1_inf=lvl_1_inf)
  end

  tspan = (1.0, 10000)
  
  # Solve problem
  prob = ODEProblem(source_sink2!, u₀, tspan, p)
  return solve(prob, DP5(), saveat = 1., reltol=1e-8, abstol=1e-8)
end

function main()
  # β, α, γ, ρ, b, c = 0.07, 0.5, 1, 0.1, 0.18, 1.05
  args = parse_commandline()
  
  # Init
  if isnothing(args["db"])
    
    β = args["beta"]
    α = args["a"]
    γ = args["g"]
    ρ = args["r"]
    η = args["e"]
    b = args["b"]
    c = args["c"]
    μ = args["m"]
    
    p = [β, α, γ, ρ, η, b, c, μ]
    sol = run_source_sink2(p)
    write_sol2txt("$(args["o"])/sourcesink2_$(join(p, "_")).txt", sol) 
  else
    db = SQLite.DB(args["db"])
    con = DBInterface.execute(db, """SELECT * from sourcesink2 LIMIT $(args["L"]) OFFSET $(args["O"])""") |> DataFrame

    for row in eachrow(con)
      
      β = row["beta"]
      α = row["alpha"]
      γ = row["gamma"]
      ρ = row["rho"]
      η = row["eta"]
      b = row["b"]
      c = row["cost"]
      μ = row["mu"]
      
      p = [β, α, γ, ρ, η, b, c, μ]
      sol = run_source_sink2(p)
      write_sol2txt("$(args["o"])/sourcesink2_$(join(p, "_")).txt", sol)
    end
  end  
end

main()

# prototyping -------------------------------------------------------------------------------

default(legendfont = ("Computer modern", 12),
        tickfont = ("Computer modern", 12),
        guidefontsize = 12, markerstrokewidth=0., markersize = 4.,
        linewidth=1, framestyle=:axis,
        titlefontsize=10)


# abstract figure

# β, α, γ, ρ, η, b, c, μ = 0.5, 1., 1., 0.2, 0.1, -0.3, 1., 0.0001

# p = [β, α, γ, ρ, η, b, c, μ]
# sol = run_source_sink2(p)

# inst_level, inst_level_prop = parse_sol(sol)  # params: β, γ, ρ, η, b, c, μ, δ


# regimes

# tmax = 10000

# function plot_regimes(ηs, t_max, lvl_1_inf)
#   # 0.13, 2., 1., 0.05, -1., 1., 0.0001
#   p = [0.17, 2., 1., 0.05, ηs[1], -1., 1., 0.0001]  # β, α, γ, ρ, η, b, c, μ
#   sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf)
#   res, res_prop = parse_sol(sol)
#   L = length(res)
#   global_freq = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
#   pl = scatter(1:t_max, global_freq, width = 4, xscale=:log, xlabel = L"\textrm{time}",
#             ylabel = L"\textrm{prevalence}", legend=:right, label = L"\eta/\rho = %$(round(ηs[1]/0.05, digits = 2))",
#             palette = palette(:Greys)[4:(4+size(ηs,1))], grid =:none,
#             xlims = (2,t_max))
  
#   for i in 2:size(ηs,1)
#     p[5] = ηs[i]
#     sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf)
#     res, res_prop = parse_sol(sol)
#     global_freq = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
#     scatter!(2:t_max, global_freq, width = 4, label = L"\eta/\rho = %$(round(ηs[i]/0.05, digits = 2))")
#   end
#   title!(join([L"%$(p)=%$(v)" for (p,v) in zip(["β","α","γ","ρ","η","b","c","μ"], p)], ", "))
#   return pl
# end

# ηs = [5.,0.05,0.0005]
# p1 = plot_regimes(ηs, 5000, true)
# p2 = plot_regimes(ηs, 5000, false)


t_max = 5500
# i.c. 1
lvl_1_inf = true
susc_pop = true
ηs = [5.,0.05,0.0005]
p = [0.17, 2., 1., 0.05, ηs[1], -1., 1., 0.0001]  # β, α, γ, ρ, η, b, c, μ
sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf, susc_pop=susc_pop)
res, res_prop = parse_sol(sol)
L = length(res)
global_freq = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
pl = plot(1:t_max, global_freq, width = 4, xscale=:log, xlabel = L"\textrm{time}",
            ylabel = L"\textrm{prevalence}", legend=:right, label = L"\ \eta/\rho = %$(round(ηs[1]/0.05, digits = 2)),\ \textrm{i.c.}\ 1",
            palette = palette(:Reds)[[3,6,8]], grid =:none);
for i in 2:size(ηs,1)
  p[5] = ηs[i]
  sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf, susc_pop=susc_pop)
  res, res_prop = parse_sol(sol)
  global_freq = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
  plot!(1:t_max, global_freq, width = 4, label = L"\ \eta/\rho = %$(round(ηs[i]/0.05, digits = 2)),\ \textrm{i.c.}\ 1");
end
# i.c. 2
lvl_1_inf = false
susc_pop = true
p = [0.17, 2., 1., 0.05, ηs[3], -1., 1., 0.0001]  # β, α, γ, ρ, η, b, c, μ
sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf, susc_pop=susc_pop)
res, res_prop = parse_sol(sol)
global_freq = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
plot!(1:t_max, global_freq, width = 4,
      label = L"\ \eta/\rho = %$(round(ηs[3]/0.05, digits = 2)),\ \textrm{i.c.}\ 2",
      palette = palette(:Blues)[[8]]);
# i.c. 3
lvl_1_inf = false
susc_pop = false
p = [0.17, 2., 1., 0.05, ηs[2], -1., 1., 0.0001]  # β, α, γ, ρ, η, b, c, μ
sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf, susc_pop=susc_pop)
res, res_prop = parse_sol(sol)
global_freq = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
plot!(1:t_max, global_freq, width = 4,
      label = L"\ \eta/\rho = %$(round(ηs[2]/0.05, digits = 2)),\ \textrm{i.c.}\ 3",
      palette = palette(:Greens)[[8]], grid =:none)

savefig("NetSci_abstract_fig___.pdf")