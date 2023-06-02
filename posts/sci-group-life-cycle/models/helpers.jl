using CSV
using OrdinaryDiffEq:ODESolution

processing_sol1(x, n) = sum((collect(0:(n-1)) / n) .* x) / sum(x) 

function parse_sol(s::ODESolution)
    t1 = 1
    @assert s.u[t1] isa ArrayPartition
    N2 = length(first(s.u[t1].x)) # N2 can be anything; institution level, # predators, # programmers
    tmax = length(s)-1
    n2_indices = Dict()
    n2_indices_prop = Dict()
      
    for n₂=1:N2
      values = []
      values_prop = []
      for t=1:tmax
        n2_probs = s.u[t].x[n₂] # probs to find system in states n₂, regardless of n₁
        out = processing_sol1(n2_probs,N2)
        push!(values, out)
        out = sum(n2_probs)
        push!(values_prop, out)
      end
      n2_indices[n₂] = values
      n2_indices_prop[n₂] = values_prop
    end
    return n2_indices, n2_indices_prop
  end