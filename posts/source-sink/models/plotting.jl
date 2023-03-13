using Pkg; Pkg.activate("../../");
# using Distributions, StatsBase, OrdinaryDiffEq, RecursiveArrayTools, DataFrames, SQLite, Plots, LaTeXStrings
using Plots, Measures

include("helpers.jl")
include("sourcesink2.jl")

default(legendfont = ("Computer modern", 12),
        tickfont = ("Computer modern", 12),
        guidefontsize = 12, markerstrokewidth=0., markersize = 4.,
        linewidth=1, framestyle=:axis,
        titlefontsize=10, grid=:none)

gr(size=(650,400))

# ----------------------------- NETSCI MAIN PLOT ------------------------------ #

t_max = 3500

# i.c. = 1
lvl_1_inf = true
ηs = [5.,0.05,0.0005]
p = [0.17, 1., 2., 1., 0.05, ηs[1], -1., 1., 0.0001]  # β, ξ, α, γ, ρ, η, b, c, μ
sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf)
res, res_prop = parse_sol(sol)

L = length(res)

global_freq1 = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
pl = plot(1:t_max, global_freq1[1:t_max], width = 4, xscale=:log, xlabel = L"\textrm{time}",
            ylabel = L"\textrm{prevalence}", legend=:top, label = L"\ \eta = %$(ηs[1]),\ \beta = %$(p[1]),\ \textrm{i.c.}\ 1",
            color = "seagreen", grid =:none);

## slower inst
p[6] = ηs[2]
sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf)
res, res_prop = parse_sol(sol)
global_freq2 = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
plot!(1:t_max, global_freq2[1:t_max], width = 4, label = L"\ \eta = %$(ηs[2]),\ \beta = %$(p[1]),\ \textrm{i.c.}\ 1",
        color = "deepskyblue");

## slowest inst but still i.c = 1
p[6] = ηs[3]
sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf)
res, res_prop = parse_sol(sol)
global_freq3 = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
plot!(1:t_max, global_freq3[1:t_max], width = 4, label = L"\ \eta = %$(ηs[3]),\ \beta = %$(p[1]),\ \textrm{i.c.}\ 1",
      color = "darkorchid2");

pl

# i.c. 2
lvl_1_inf = false
p[6] = ηs[3]
sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf)
res, res_prop = parse_sol(sol)
global_freq4 = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
plot!(1:t_max, global_freq4[1:t_max], width = 4,
      label = L"\ \eta = %$(ηs[3]),\ \beta = %$(p[1]),\ \textrm{i.c.}\ 2",
      color = "purple");

# i.c. 3
lvl_1_inf = false
perc_inf = 0.3
p[6] = ηs[2]
sol = run_source_sink2(p, perc_inf = perc_inf, lvl_1_inf=lvl_1_inf)
res, res_prop = parse_sol(sol)
global_freq5 = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
plot!(1:t_max, global_freq5[1:t_max], width = 4,
      label = L"\ \eta = %$(ηs[2]),\ \beta = %$(p[1]),\ \textrm{i.c.}\ 3",
      color = "blue4");
# p[1] = 0.06
# sol = run_source_sink2(p, perc_inf = perc_inf, lvl_1_inf=lvl_1_inf)
# res, res_prop = parse_sol(sol)
# global_freq6 = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
# plot!(1:t_max, global_freq6[1:t_max], width = 4,
#       label = L"\ \eta = %$(ηs[2]),\ \beta = %$(p[1]),\ \textrm{i.c.}\ 3",
#       color = "orange3")


savefig("figs/NetSci_abstract_fig_no_orange.pdf")


# -------------------------------- Extra fig 1 ------------------------------- #

t_max = 3500

ηs = [5.,0.05,0.0005]
p = [0.17, 1., 2., 1., 0.05, ηs[1], -1., 1., 0.0001]  # β, ξ, α, γ, ρ, η, b, c, μ

# i.c. = 1
## η = 5.0; β = 0.17
### fast imitation, strong but costly institutions can fail to emerge
### epidemic becomes highly endemic. Not great.
lvl_1_inf = true
sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf)
res, res_prop = parse_sol(sol)

L = length(res)
tmax = length(res[L[1]])

pl1 = plot_scatter(res, res_prop, size=(500, 350))
plot!([], [], legend=false)
annotate!(repeat([t_max/0.15], inner=L), [last(res[i]) for i=eachindex(res)], [L"$\ell=$"*"$(i)" for i=eachindex(res)], 10)
xlabel!("timestep")
ylabel!(L"$I(\ell$)")

pl2 = plot_scatter(res, res_prop, plot_prop=true, size=(500, 350))
plot!([], [], legend=false)
annotate!(repeat([t_max/0.15], inner=L), [last(res_prop[i]) for i=eachindex(res)], [L"$\ell=$"*"$(i)" for i=eachindex(res)], 10)
xlabel!("timestep")
ylabel!(L"prop($\ell$)")

pl_tot = plot(pl1, pl2, layout = @layout([a b]), size=(1000, 350))
plot!(margin=10mm)
savefig("figs/fast_imit_ic1_green.pdf")



## η = 0.05; β = 0.17
p[6] = ηs[2]
lvl_1_inf = true
sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf)
res, res_prop = parse_sol(sol)

L = length(res)
tmax = length(res[L[1]])

pl3 = plot_scatter(res, res_prop, size=(500, 350))
plot!([], [], legend=false)
annotate!(repeat([t_max/0.15], inner=L), [i == 4 ? last(res[i])-0.03 : last(res[i]) for i=eachindex(res)], [L"$\ell=$"*"$(i)" for i=eachindex(res)], 10)
xlabel!("timestep")
ylabel!(L"$I(\ell$)")

pl4 = plot_scatter(res, res_prop, plot_prop=true, size=(500, 350))
plot!([], [], legend=false)
annotate!(repeat([t_max/0.15], inner=L), [i == 4 ? last(res_prop[i])-0.01 : last(res_prop[i]) for i=eachindex(res)], [L"$\ell=$"*"$(i)" for i=eachindex(res)], 10)
xlabel!("timestep")
ylabel!(L"prop($\ell$)")

pl_tot2 = plot(pl3, pl4, layout = @layout([a b]), size=(1000, 350))
plot!(margin=10mm)
savefig("figs/slower_inst_ic1_lightblue.pdf")

# i.c. = 1
## η = 0.0005; β = 0.17
lvl_1_inf = true
p[6] = ηs[3]
sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf)
res, res_prop = parse_sol(sol)

L = length(res)
tmax = length(res[L[1]])

pl5 = plot_scatter(res, res_prop, size=(500, 350))
plot!([], [], legend=false)
annotate!(repeat([t_max/0.15], inner=L), [i == 4 ? last(res[i])-0.03 : last(res[i]) for i=eachindex(res)], [L"$\ell=$"*"$(i)" for i=eachindex(res)], 10)
xlabel!("timestep")
ylabel!(L"$I(\ell$)")

pl6 = plot_scatter(res, res_prop, plot_prop=true, size=(500, 350))
plot!([], [], legend=false)
annotate!(repeat([t_max/0.15], inner=L), [i == 4 ? last(res_prop[i])-0.01 : last(res_prop[i]) for i=eachindex(res)], [L"$\ell=$"*"$(i)" for i=eachindex(res)], 10)
xlabel!("timestep")
ylabel!(L"prop($\ell$)")

pl_tot3 = plot(pl5, pl6, layout = @layout([a b]), size=(1000, 350))
plot!(margin=10mm)
savefig("figs/slowest_inst_ic1_darkorchid.pdf")


# i.c. = 2
## η = 0.0005; β = 0.17
lvl_1_inf = false
p[6] = ηs[3]
sol = run_source_sink2(p, lvl_1_inf=lvl_1_inf)
res, res_prop = parse_sol(sol)

L = length(res)
tmax = length(res[L[1]])

pl6 = plot_scatter(res, res_prop, size=(500, 350))
plot!([], [], legend=false)
annotate!(repeat([t_max/0.15], inner=L), [i == 4 ? last(res[i])-0.03 : last(res[i]) for i=eachindex(res)], [L"$\ell=$"*"$(i)" for i=eachindex(res)], 10)
xlabel!("timestep")
ylabel!(L"$I(\ell$)")

pl7 = plot_scatter(res, res_prop, plot_prop=true, size=(500, 350))
plot!([], [], legend=false)
annotate!(repeat([t_max/0.15], inner=L), [i == 4 ? last(res_prop[i])-0.01 : last(res_prop[i]) for i=eachindex(res)], [L"$\ell=$"*"$(i)" for i=eachindex(res)], 10)
xlabel!("timestep")
ylabel!(L"prop($\ell$)")

pl_tot4 = plot(pl5, pl6, layout = @layout([a b]), size=(1000, 350))
plot!(margin=10mm)
savefig("figs/slowest_inst_ic2_purple.pdf")


# i.c. 3 (blue4 / darkblue)
lvl_1_inf = false ## random level not just first level infected
perc_inf = 0.3 ## we start with an outbreak
p[6] = ηs[2]
sol = run_source_sink2(p, perc_inf = perc_inf, lvl_1_inf=lvl_1_inf)
res, res_prop = parse_sol(sol)

L = length(res)
tmax = length(res[L[1]])

pl8 = plot_scatter(res, res_prop, size=(500, 350))
plot!([], [], legend=false)
annotate!(repeat([t_max/0.15], inner=L), [i == 4 ? last(res[i])-0.03 : last(res[i]) for i=eachindex(res)], [L"$\ell=$"*"$(i)" for i=eachindex(res)], 10)
xlabel!("timestep")
ylabel!(L"$I(\ell$)")

pl9 = plot_scatter(res, res_prop, plot_prop=true, size=(500, 350))
plot!([], [], legend=false)
annotate!(repeat([t_max/0.15], inner=L), [i == 4 ? last(res_prop[i])-0.01 : last(res_prop[i]) for i=eachindex(res)], [L"$\ell=$"*"$(i)" for i=eachindex(res)], 10)
xlabel!("timestep")
ylabel!(L"prop($\ell$)")

pl_tot5 = plot(pl5, pl6, layout = @layout([a b]), size=(1000, 350))
plot!(margin=10mm)
savefig("figs/slower_inst_ic3_darkblue.pdf")
