#####################################################################
# File 'field.jl' -- hydrodynamic variables fields
# FSLA.JL project
# Copyright (c) Vadim V. Shepelev, Ph.D.
#####################################################################


mutable struct Field_primitive
    """ w = (rho, u, p)^T
    """
    N::Int64
    x::Vector{Float64}
    x_new::Vector{Float64}
    xmin::Float64
    xmax::Float64
    w::Vector{Vector{Float64}}
    w_new::Vector{Vector{Float64}}
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
        w = [deepcopy(pr["task"].w) for _ in 1:N+1]
        w_new = deepcopy(w)
        new(N, x, x_new, xmin, xmax, w, w_new, tmin, tmin, tmax, dm, 0.0)
    end
end