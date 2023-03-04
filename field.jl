#####################################################################
# File 'field.jl' -- hydrodynamic variables fields
# FSLA.JL project
# Copyright (c) Vadim V. Shepelev, Ph.D.
#####################################################################
using StaticArrays

include("problem.jl")


mutable struct Field_primitive
    """ w = (rho, u, p)^T
    """
    N::Int64
    x::Vector{Float64}
    w::Vector{SVector{3,Float64}}
    w_new::Vector{SVector{3,Float64}}
    t::Float64
    dm::Float64
    dt::Float64

    function Field_primitive(pr::Dict)
        N = pr["N"]
        xmin, xmax = pr["task"].x_span 
        x = collect(range(xmin, xmax, N+1))
        x_new = collect(range(xmin, xmax, N+1))
        rho = pr["task"].w[1]
        dm = (xmax - xmin) / N / rho
        new(N, x, [zeros(SVector{3}) for _ in 1:N+1], [zeros(SVector{3}) for _ in 1:N+1], 0.0, dm, 0.0)
    end
end


    
