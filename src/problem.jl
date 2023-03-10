#####################################################################
# File 'problem.jl' -- hydrodynamic problems formulation
# FSLA.JL project
# Copyright (c) Vadim V. Shepelev, Ph.D.
#####################################################################

abstract type Problem end

###### 1D Vacuum Problem in primitive variables formulation
struct VacuumProblem <: Problem
    x_span::Tuple{Float64, Float64}
    t_span::Tuple{Float64, Float64}
    w::Vector{Float64}

    function VacuumProblem(; 
        x_span::Tuple{Float64, Float64}, 
        t_span::Tuple{Float64, Float64},
        w::Vector{Float64}
        )
        new(x_span, t_span, w) 
    end
end

function Base.show(io::IO, pr::VacuumProblem)
    print(io, "VacuumProblem: [$(pr.x[1]), $(pr.x[2])], wl=$(pr.wl), wr=$(pr.wr), x0=$(pr.xbnd)")
end


###### 1D Riemann Problem in primitive variables formulation
struct RiemannProblem <: Problem
    x_span::Tuple{Float64, Float64}
    t_span::Tuple{Float64, Float64}
    wl::Vector{Float64}
    wr::Vector{Float64}
    xbnd::Float64

    function RiemannProblem(; 
        x_span::Tuple{Float64, Float64}, 
        t_span::Tuple{Float64, Float64},
        wl::Vector{Float64}, 
        wr::Vector{Float64}, 
        xbnd::Float64)

        new(x_span, t_span, wl, wr, xbnd) 
    end
end

function Base.show(io::IO, pr::RiemannProblem)
   print(io, "RiemannProblem: [$(pr.x_span[1]), $(pr.x_span[2])], [$(pr.t_span[1]), $(pr.t_span[2])], wl=$(pr.wl), wr=$(pr.wr), x0=$(pr.xbnd)")
end




