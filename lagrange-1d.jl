# include("problem.jl")
include("eos.jl")
include("field.jl")
include("method.jl")


rp = VacuumProblem(x_span=(0., 1.), t_span=(0.0, 0.2), w=[2650., 0., 35.6])

problem = Dict(
    "task" => rp, 
    "method" => "samarskii", 
    "N" => 100, 
    "CFL" => 0.5)

fld = Field_primitive(problem)
eos = eos_ideal(1.4)

println("fsLA.jl hydrodynamic 1D Lagrangian code, (c)2023 Vadim V. Shepelev, Ph.D.")

iter = 1

while fld.t < fld.tmax
    fld.dt = calcdt(problem, eos, fld)
    if fld.t + fld.dt > fld.tmax 
        fld.dt = pr.t_max - fld.t 
    end
    # println("step=$iter: $(calc(problem, eos, fld)) iterations, t=$(fld.t), dt=$(fld.dt)")
end



