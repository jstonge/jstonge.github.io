using CSV

"""
write_sol2txt(path, sol)

Function to write solution to textfile. Nothing to do here.
"""
function write_sol2txt(path, sol)
    L = length(sol.u[1].x)
    open(path, "a") do io
        for t=1:length(sol.u), ℓ=1:L
            for val in sol.u[t].x[ℓ]
                write(io, "$(t) $(ℓ) $(round(val, digits=14))\n")
            end
        end
    end
end


processing_sol1(x, n) = sum((collect(0:(n-1)) / n) .* x) / sum(x) 

function parse_sol(s::String)
    
    sol = CSV.read(s, DataFrame; header=["timestep", "L", "value"])
    L = 6
    inst_level = Dict()
    lower_limit, upper_limit = 1, 21
    for t=1:4000
        for ℓ=1:L
            myrange = UnitRange(lower_limit:upper_limit)
            n = length(sol.value[myrange])
            x = sol.value[myrange]
            out = processing_sol1(x, n)
        if haskey(inst_level, ℓ)
            inst_level[ℓ] = [inst_level[ℓ]; out]
        else
            inst_level[ℓ] = out
        end

        lower_limit += 21
        upper_limit += 21

        end
    end
    return inst_level
end

function parse_sol(s::ODESolution)
    L = length(s.u[1].x)
    for ℓ=1:L
      values = []
      for t=1:4000
        n = length(s.u[t].x[ℓ])
        x = s.u[t].x[ℓ]
        out = processing_sol1(x,n)
        push!(values, out)
      end
      inst_level[ℓ] = values
    end
    return inst_level
end

function plot_scatter(res::Dict)
    L = length(res)
    scatter(1:length(res[1]), [res[i] for i in 1:L], xaxis = :log, legendtitle="grsize", 
            legend=:outertopright, labels=collect(1:L)', palette = palette(:Reds)[2:(L-1)],
            markerstrokewidth=0, markersize = 3.)
end

# using SQLite

# function plot_phase_diagram(res_db)
    
#     L = length(res)
#     ys = [last(inst_level[i]) for i=1:L]

#     ps = [heatmap(1:length(res[1]), ys[lvl]) for lvl=1:6]
# end