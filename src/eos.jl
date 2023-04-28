module EquationsOfState

export EOS, eos_ideal, eos_LiF, getp, gete, getc

#####################################################################
# File 'eos.jl' -- equations of state
# FSLA.JL project
# Copyright (c) Vadim V. Shepelev, Ph.D.
#####################################################################

abstract type EOS end


# Ideal gas EOS
struct eos_ideal <: EOS  
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


# Mie-Gruneisen LiF EOS
struct eos_LiF <: EOS  
    rho0::Float64
    G::Float64
    c::Float64
    s::Float64    
    
    # rho0(2650.), G(.71), c(5150.), s(1.35) 
    function eos_LiF()
        new(2650., .71, 5150., 1.35)
    end
end

function getp(eos::eos_LiF, rho::Float64, e::Float64)::Float64  # [Pa]
	p = p_cold(eos, rho) + rho*eos.G*(e - e_cold(eos, rho)) 
end

function gete(eos::eos_LiF, rho::Float64, p::Float64)::Float64  # [J/kg]
	e::Float64 = e_cold(eos, rho) + (p - p_cold(eos, rho)) / (rho * eos.G)
end

function getc(eos::eos_LiF, rho::Float64, p::Float64)::Float64    # [m/s]
	dp_drho::Float64 = - eos.rho0 / (rho * rho) * dp_cold_dx(eos, rho) + 
        eos.G*gete(eos, rho, p) - eos.G*e_cold(eos, rho) + eos.G*eos.rho0/rho*de_cold_dx(eos, rho)
	dp_de = eos.G * rho
	c = sqrt(dp_drho + p*dp_de / (rho * rho))
end

function p_cold(eos::eos_LiF, rho::Float64)::Float64
	x::Float64 = eos.rho0 / rho 
	pc::Float64 = eos.rho0 * eos.c * eos.c * (1.0 - x)/((1.0 - eos.s*(1.0 - x))*(1.0 - eos.s*(1.0 - x))) 
end

function e_cold(eos::eos_LiF, rho::Float64)::Float64
	x::Float64 = eos.rho0 / rho 
	ec::Float64 = p_cold(eos, rho) / eos.rho0 * (1.0 - x) / 2.0; 
end

function dp_cold_dx(eos::eos_LiF, rho::Float64)::Float64 
	x::Float64 = eos.rho0 / rho
	dpc_dx::Float64 = -eos.rho0*eos.c*eos.c*(1.0+eos.s*(1.0-x)) / ((1.0-eos.s*(1.0-x))*(1.0-eos.s*(1.0-x))*(1.0-eos.s*(1.0-x)))
end

function de_cold_dx(eos::eos_LiF, rho::Float64)::Float64 
	x::Float64 = eos.rho0 / rho
	dec_dx::Float64 = 1.0 / (2.0 * eos.rho0) * ((1.0 - x) * dp_cold_dx(eos, rho) - p_cold(eos, rho))
end


end # module EquationsOfState