using Pkg; Pkg.activate("../../");
using ArgParse, Distributions, StatsBase, OrdinaryDiffEq, RecursiveArrayTools, DataFrames, SQLite, Plots
using OrdinaryDiffEq:ODESolution
using Plots.PlotMeasures

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

function initialize_u0(;N=20, M::Int=100, p::Float64=0.01)
  N_plus_one = N+1
  G = zeros(N_plus_one, N_plus_one)
  for _ in 1:M
    ℓ = rand(1:N_plus_one)
    i = sum(rand(Binomial(1, p), N_plus_one))
    G[ℓ, i+1] += 1
  end

  G = G ./ M

  return G
end

function plot_cost(a=3; p=5)
  p1=plot(x -> c(x, p, a=a), 1, 21, label="p=$(p)")
  xlabel!("# non programmers")
  
  p2=plot(x -> c(1, x, a=a), 0, 21, label="n=1")
  plot!(x -> c(5, x, a=a),  0, 21, label="n=5")
  plot!(x -> c(10, x, a=a), 0, 21, label="n=10")
  plot!(x -> c(15, x, a=a), 0, 21, label="n=15")
  xlabel!("# programmers")
  
  p3=plot(p1,p2, layout=(2,1))
  ylabel!("cost")
  title!("a=$(a)")  
  return p3
end

function plot_tryptic_cost(p)
  p1=plot_cost(5; p=p)
  p2=plot_cost(3; p=p)
  p3=plot_cost(1; p=p)
  plot(p1,p2,p3, layout=(1,3), bottom_margin = 10mm)
  plot!(size=(700,500))
end

plot_tryptic_cost(1)

c(n, i; a=3) = n == i == 0 ? 0.95 : 0.95 * exp(-a*i / n)  # cost function

τ(n, i, α, β) = exp(-α + β*(1 - c(n, i))) # group benefits


# plot(x -> τ(1,x,α,β)*c(1,x), 0, 20, label="#nonprog=1")
# plot!(x -> τ(3,x,α,β)*c(3,x), 0, 20, label="#nonprog=3")
# plot!(x -> τ(5,x,α,β)*c(5,x), 0, 20, label="#nonprog=5")
# xlabel!("# prog")
# ylabel!("cost")
# title!("plotting τ*c")


function life_cycle_research_groups!(du, u, p, t)
  N, P = size(u) # Note that there can be no coders but not non-coders
  G = u # G for groups
  μ, νₙ, νₚ, α, β, a = p
  for n=1:N, i=1:P
    coder, noncoder = i-1, n-1 
    du[n,i] = 0

    noncoder > 0 && ( du[n,i] += μ*G[n-1,i] ) # non-prog repro input
    n < N && ( du[n,i] -= G[n,i]*μ )                       # non-prog repro output

    n < N && ( du[n,i] += G[n+1,i]*νₙ*(noncoder+1) ) # non-prog death input
    du[n,i] -= G[n,i]*noncoder*νₙ                    # non-prog death output
    
    i < P && ( du[n,i] += G[n,i+1]*νₚ*(coder+1) ) # prog death input
    du[n,i] -= G[n,i]*coder*νₚ                    # prog death output
    
    (coder > 0 && n < N) && ( du[n,i] += G[n+1,i-1]*(noncoder+1)*τ(noncoder+1, coder-1, α, β)*(1-c(noncoder+1,coder-1, a=a)) ) # non-prog to prog predation input
    i < P && ( du[n,i] -= G[n,i]*noncoder*τ(noncoder, coder, α, β)*(1-c(noncoder,coder, a=a)) ) # non-prog to prog predation output
    
    i < P && ( du[n,i] -= G[n,i]*noncoder*τ(noncoder, coder, α, β)*c(noncoder,coder, a=a) )           # non-prog death output due to cost
    (n < N && i < P) && ( du[n,i] += G[n+1,i]*(noncoder+1)*τ(noncoder+1, coder, α, β)*c(noncoder+1,coder, a=a) ) # non-prog death input due to cost
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

# main()


# prototyping -------------------------------------------------------------------------------

μ  = 0.1   # inflow new students-non coders
νₙ = 0.01    # death rate non-coders
νₚ = 0.05    # death rate coders
α  = 0.01    # benefits non coders
β  = 0.02     # benefits coders
a = 1       # parameter cost function
p  = [μ, νₙ, νₚ, α, β, a]

u₀ = initialize_u0(N=20)

t_max = 4000
tspan = (0., t_max)

prob = ODEProblem(life_cycle_research_groups!, u₀, tspan, p)
sol = solve(prob, Tsit5(), saveat=1, reltol=1e-8, abstol=1e-8)


# ------------------------------ wrangling n1+n2 ----------------------------- #

# for s in group size
# for p in 1:max number of programers OR s
# (nb prog / group size) * (sol for that nb prog and p)
# normalized by (sol for that nb prog and p)
function wrangle()
  tot_out = []
  for t=1:t_max
    out_num = zeros(41)
    out_denum = zeros(41)
    #!TODO: not just max but for all ts.
    for s=1:20
      for p=1:minimum([21,s+1])
        out_num[s+1] += ((p-1) / s) * sol[t][s-p+2,p]
        out_denum[s+1] += sol[t][s-p+2,p]
      end
    end
    push!(tot_out, out_num[2:20] ./ out_denum[2:20])
  end
  return tot_out
end

out = wrangle()

function plot_sol()

  param_str = join(["$(pname)=$(p);" for (pname, p) in zip(["μ", "νₙ", "νₚ", "α", "β", "a", "p"], p)], ' ')

  ps=plot(2:20, out[2], label="t=2", legend=:outerright, top_margin = 20mm)
  for t=collect(5:5:30)
    plot!(2:20, out[t], label="t=$(t)") 
  end
  plot!(2:20, out[t_max], label="t=$(t_max)")
  ps
  title!("Many programmers in large groups while\nwe have few programmers in small groups\n($(param_str))")
  xlabel!("group size")
  ylabel!("proportion programmers")
  # ylims!(0,1)
  plot!(size=(650,400))
  # savefig("mymodel.pdf")
end

plot_sol()


# checks
round(sum(sol[1]), digits= 2)
round(sum(sol[t_max]), digits= 2)

round.(sum(sol[t_max], dims=1), digits=4)
round.(sum(sol[t_max], dims=2), digits=4)
