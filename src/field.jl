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
    w_prev::Vector{Vector{Float64}}
    g::Vector{Float64}
    t::Float64
    tmin::Float64
    tmax::Float64
    dm::Float64
    dt::Float64
    CFL::Float64

    function Field_primitive(pr::Dict)
        N = pr["N"]
        xmin, xmax = pr["task"].x_span
        tmin, tmax = pr["task"].t_span 
        dx::Float64 = (xmax - xmin)/N
        x = collect(range(xmin-dx, xmax+dx, N+3))        
        x_new = deepcopy(x)
        rho = pr["task"].w[1]
        dm = (xmax - xmin) / N / rho
        w = [deepcopy(pr["task"].w) for _ in 1:N+3]
        w_new = deepcopy(w)
        w_prev = deepcopy(w)
        g = zeros(N+3)
        CFL = pr["CFL"]
        new(N, x, x_new, xmin, xmax, w, w_new, w_prev, g, tmin, tmin, tmax, dm, 0.0, CFL)
    end
end