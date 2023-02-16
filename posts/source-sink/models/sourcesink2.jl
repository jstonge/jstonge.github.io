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
        "--xi"
        arg_type = Float64
        default = 1.
        help = "Simple-complex contagion parameter"
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

function initialize_u0(;n::Int=20, L::Int=6, M::Int=20, p::Float64=0.001,
  lvl_1_inf::Bool=false)
  G = zeros(L, n+1)

  if lvl_1_inf # ≈ only (99.995%) lowest level populated at t = 0
    for _ in 1:M
      i = sum(collect(rand(Binomial(1, 2*p), n))) # how many total adopters?
      G[1, i+1] += 1 # everytime combination [ℓ,i], count +1
    end
    for _ in 1:M
      ℓ = rand(2:L) # pick a level
      i = sum(collect(rand(Binomial(1, 0.01*p), n))) # how many total adopters?
      G[ℓ, i+1] += 1 # everytime combination [ℓ,i], count +1
    end
    G = G ./ (2*M) # normalized by tot number of groups
    return ArrayPartition(Tuple([G[ℓ,:] for ℓ=1:L]))
  else 
      for _ in 1:M
        ℓ = rand(1:L) # pick a level
        i = sum(collect(rand(Binomial(1, p), n))) # how many total adopters?
        G[ℓ, i+1] += 1 # everytime combination [ℓ,i], count +1
      end
    G = G ./ M # normalized by tot number of groups
    return ArrayPartition(Tuple([G[ℓ,:] for ℓ=1:L]))
  end
end

# function to switch between simple (ξ = 1) and complex contagion (ξ ≠ 1)
g(x; ξ=1.) = x^ξ # for x ∈ ℕ ⇒ ξ = 1: linear growth; 0 < ξ < 1: sublinear growth; ξ > 1: superlinear growth

function source_sink2!(du, u, p, t)
    G, L, n = u, length(u.x), length(first(u.x))
    β, ξ, α, γ, ρ, η, b, c, μ = p
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
        du.x[ℓ][i] = -γ*n_adopt*G.x[ℓ][i] - β*(ℓ^-α)*g(n_adopt+R, ξ=ξ)*(gr_size-n_adopt)*G.x[ℓ][i]
        n_adopt > 0 && ( du.x[ℓ][i] += β*(ℓ^-α)*g(n_adopt-1+R, ξ=ξ)*(gr_size-n_adopt+1)*G.x[ℓ][i-1])
        n_adopt < gr_size && ( du.x[ℓ][i] +=  γ*(n_adopt+1)*G.x[ℓ][i+1] )
        # Group selection process
        ℓ > 1 && ( du.x[ℓ][i] += η*G.x[ℓ-1][i]*(Z[ℓ] / Z[ℓ-1] + μ) - η*G.x[ℓ][i]*(Z[ℓ-1] / Z[ℓ] + μ) )
        ℓ < L && ( du.x[ℓ][i] += η*G.x[ℓ+1][i]*(Z[ℓ] / Z[ℓ+1] + μ) - η*G.x[ℓ][i]*(Z[ℓ+1] / Z[ℓ] + μ) )
      end
end

function run_source_sink2(p; perc_inf::Float64=0.001, lvl_1_inf::Bool=false)
  n, M = 20, 1000
  L = 6
  u₀ = initialize_u0(n=n, L=L, M=M, p=perc_inf, lvl_1_inf=lvl_1_inf)

  tspan = (1.0, 10000)
  
  # Solve problem
  prob = ODEProblem(source_sink2!, u₀, tspan, p)
  return solve(prob, DP5(), saveat = 1., reltol=1e-8, abstol=1e-8)
end

function main()
  args = parse_commandline()
  
  # Init
  if isnothing(args["db"])
    
    β = args["beta"]
    ξ = args["xi"]
    α = args["a"]
    γ = args["g"]
    ρ = args["r"]
    η = args["e"]
    b = args["b"]
    c = args["c"]
    μ = args["m"]
    
    p = [β, ξ, α, γ, ρ, η, b, c, μ]
    sol = run_source_sink2(p)
    write_sol2txt("$(args["o"])/sourcesink2_$(join(p, "_")).txt", sol) 
  else
    db = SQLite.DB(args["db"])
    con = DBInterface.execute(db, """SELECT * from sourcesink2 LIMIT $(args["L"]) OFFSET $(args["O"])""") |> DataFrame

    for row in eachrow(con)
      
      β = row["beta"]
      ξ = row["xi"]
      α = row["alpha"]
      γ = row["gamma"]
      ρ = row["rho"]
      η = row["eta"]
      b = row["b"]
      c = row["cost"]
      μ = row["mu"]
      
      p = [β, ξ, α, γ, ρ, η, b, c, μ]
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
        titlefontsize=10, grid=:none,
        bottom_margin = 1mm, left_margin = 1mm, right_margin = 0mm)
gr(size=(650,400))

# # abstract figure
# run_source_sink2(p, lvl_1_inf=lvl_1_inf)
# t_max = 3500
# # i.c. 1
# lvl_1_inf = true
# ηs = [5.,0.05,0.0005]
# p = [0.17, 1., 2., 1., 0.05, ηs[1], -1., 1., 0.0001]  # β, ξ, α, γ, ρ, η, b, c, μ
# sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf)
# res, res_prop = parse_sol(sol)
# L = length(res)
# global_freq1 = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
# pl = plot(1:t_max, global_freq1[1:t_max], width = 4, xscale=:log, xlabel = L"\textrm{time}",
#             ylabel = L"\textrm{prevalence}", legend=:right, label = L"\ \eta = %$(ηs[1]),\ \beta = %$(p[1]),\ \textrm{i.c.}\ 1",
#             color = "seagreen", grid =:none);
# p[6] = ηs[2]
# sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf)
# res, res_prop = parse_sol(sol)
# global_freq2 = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
# plot!(1:t_max, global_freq2[1:t_max], width = 4, label = L"\ \eta = %$(ηs[2]),\ \beta = %$(p[1]),\ \textrm{i.c.}\ 1",
#         color = "deepskyblue");
# p[6] = ηs[3]
# sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf)
# res, res_prop = parse_sol(sol)
# global_freq3 = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
# plot!(1:t_max, global_freq3[1:t_max], width = 4, label = L"\ \eta = %$(ηs[3]),\ \beta = %$(p[1]),\ \textrm{i.c.}\ 1",
#       color = "darkorchid2");
# # i.c. 2
# lvl_1_inf = false
# p[6] = ηs[3]
# sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf)
# res, res_prop = parse_sol(sol)
# global_freq4 = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
# plot!(1:t_max, global_freq4[1:t_max], width = 4,
#       label = L"\ \eta = %$(ηs[3]),\ \beta = %$(p[1]),\ \textrm{i.c.}\ 2",
#       color = "purple");
# # i.c. 3
# lvl_1_inf = false
# perc_inf = 0.3
# p[6] = ηs[2]
# sol = run_source_sink2(p, perc_inf = perc_inf, lvl_1_inf=lvl_1_inf)
# res, res_prop = parse_sol(sol)
# global_freq5 = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
# plot!(1:t_max, global_freq5[1:t_max], width = 4,
#       label = L"\ \eta = %$(ηs[2]),\ \beta = %$(p[1]),\ \textrm{i.c.}\ 3",
#       color = "blue4");
# p[1] = 0.06
# sol = run_source_sink2(p, perc_inf = perc_inf, lvl_1_inf=lvl_1_inf)
# res, res_prop = parse_sol(sol)
# global_freq6 = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
# plot!(1:t_max, global_freq6[1:t_max], width = 4,
#       label = L"\ \eta = %$(ηs[2]),\ \beta = %$(p[1]),\ \textrm{i.c.}\ 3",
#       color = "orange3")

# savefig("NetSci_abstract_fig.pdf")


params_name = "β", "ξ", "α", "γ", "ρ", "η", "b", "c", "μ"

lvl_1_inf = false
perc_inf = 0.002
p = [0.08, 1., 1., 1., 0.01, 0.01, -1., 1., 0.0001]  # β, ξ, α, γ, ρ, η, b, c, μ
sol = run_source_sink2(p, perc_inf = perc_inf, lvl_1_inf=lvl_1_inf)
res, res_prop = parse_sol(sol)
t_max = 6500
global_freq = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
plot([res[l][1:t_max] for l=1:L], xscale=:log, ylabel = L"\textrm{prevalence}", labels = " " .* string.([1:L;]'),
      width = 3., legendtitle = L"\textrm{level}", palette = palette(:Reds)[3:9], legend=:left,
      xticks = 10 .^ [0,1,2,3,4]);
plot!(1:t_max, global_freq[1:t_max], width = 3, color =:black, ls =:dash, label = L"\textrm{global}",
      title = join([params_name[i] * "=" * string.(p)[i] for i in 1:length(params_name)], "  "))
plot([res_prop[l][1:t_max] for l=1:L], xscale=:log, ylabel = L"\textrm{level\ proportion}", labels = " " .* string.([1:L;]'),
      width = 3., legendtitle = L"\textrm{level}", palette = palette(:Blues)[3:9], legend=:left,
      title = join([params_name[i] * "=" * string.(p)[i] for i in 1:length(params_name)], "  "),
      xticks = 10 .^ [0,1,2,3,4])