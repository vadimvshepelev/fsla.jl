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
    xmin::Float64
    xmax::Float64
    w::Vector{SVector{3,Float64}}
    w_new::Vector{SVector{3,Float64}}
    t::Float64
    tmin::Float64
    tmax::Float64
    dm::Float64
    dt::Float64

    function Field_primitive(pr::Dict)
        N = pr["N"]
        xmin, xmax = pr["task"].x_span
        tmin, tmax = pr["task"].t_span 
        x = collect(range(xmin, xmax, N+1))
        x_new = collect(range(xmin, xmax, N+1))
        rho = pr["task"].w[1]
        dm = (xmax - xmin) / N / rho
        w = [copy(pr["task"].w) for _ in 1:N+1]
        w_new = copy(w)
        new(N, x, xmin, xmax, w, w_new, tmin, tmin, tmax, dm, 0.0)
    end
end


####### Функция setics для постановки начальных условий #######
function setics!(pr::Dict, fld::Field_primitive)
# Не нужна, все делает конструктор Field_primitive    
end
