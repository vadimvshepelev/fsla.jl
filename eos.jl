#####################################################################
# File 'eos.jl' -- equations of state
# FSLA.JL project
# Copyright (c) Vadim V. Shepelev, Ph.D.
#####################################################################

abstract type EOS end


####### Ideal gas EOS
struct eos_ideal<: EOS  
    gamma::Float64

    function eos_ideal(_gamma::Float64)
        new(_gamma)
    end
end


function getp(eos::eos_ideal, rho::Float64, e::Float64)::Float64 
    (eos.gamma-1)*rho*e
end

function gete(eos::eos_ideal, rho::Float64, p::Float64)::Float64
    p/(eos.gamma-1)/rho
end

function getc(eos::eos_ideal, rho::Float64, p::Float64)::Float64
    sqrt(eos.gamma*p/rho)
end


####### Mie-Gruneisen LiF EOS