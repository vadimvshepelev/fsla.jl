#####################################################################
# File 'problem.jl' -- hydrodynamic problems formulation
# FSLA.JL project
# Copyright (c) Vadim V. Shepelev, Ph.D.
#####################################################################

abstract type Problem end;


###### 1D Riemann Problem in primitive variables formulation
struct RiemannProblem <: Problem
    x
    wl::Vector{Float64}
    wr::Vector{Float64}
    xbnd::Float64

    function RiemannProblem(; x::Tuple{Float64, Float64}, wl::Vector{Float64}, wr::Vector{Float64}, xbnd::Float64)
        new(x, wl, wr, xbnd) 
    end
end

function Base.show(io::IO, pr::RiemannProblem)
    print(io, "RiemannProblem: [$(pr.x[1]), $(pr.x[2])], wl=$(pr.wl), wr=$(pr.wr), x0=$(pr.xbnd)")
end




