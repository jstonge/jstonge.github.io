using Pkg; Pkg.activate("../../");
using JLD2, JSON

fname = "sourcesink1_0.07_1.0_0.1_0.18_1.05_0.0001.jld2"
modelname = split(fname, "_")[1]
p = replace(join(split(fname, "_")[2:end], "_"), ".jld" => "")
JLD2.@load fname sol

L = length(sol.u[1].x)
n = length(sol.u[1].x[1])
I = zeros(L, length(sol.t))

for t in 1:length(sol.t)
  for ℓ in 1:L
    G_nil = sol.u[t].x[ℓ]
    I[ℓ, t] = sum((collect(0:(n-1)) / n) .* G_nil) / sum(G_nil)
  end
end

# Wrangle summary stats into list of dicts where each dict is a level
new_run = []
for ℓ in 1:L
  tmp_dat = [Dict("L" => ℓ, "value" => I[ℓ, i], "timesteps" => i) for i in 1:3000]
  x = round.([i["value"] for i in tmp_dat], digits=5) # To lighten the output file we only keep unique `y` values. 
  idx = unique(z -> x[z], 1:length(x))
  push!(new_run, tmp_dat[idx] )
end

# Update old runs with new run
if isfile("data2.json")
  get!(old_run, p, vcat(new_run...))    # If already exists, update the old run Dict()
  new_run = old_run
else
  new_run = Dict(p => vcat(new_run...)) # new Dict()
end

# Write to disk
open("data2.json", "w") do f 
  write(f, JSON.json(new_run))
end