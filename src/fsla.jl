include("eos.jl")
include("field.jl")
include("method.jl")
include("problem.jl")

# rp = VacuumProblem(x_span=(0., 1.), t_span=(0.0, 0.2), w=[2650., 0., 35.6])
rp = VacuumProblem(x_span=(0., 1.), t_span=(0.0, 0.2), w=[1.0, 0.0, 1.0])

problem = Dict(
    "task" => rp, 
    "method" => "samarskii", 
    "N" => 100, 
    "CFL" => 0.5)

fld = Field_primitive(problem)
eos = eos_ideal(1.4)

println("fsLA.jl hydrodynamic 1D Lagrangian code, (c) Vadim V. Shepelev, Ph.D., 2023")
println()
println("Starting simulation...")

iter::Int64 = 1
while fld.t < fld.tmax
    fld.dt = calcdt(problem, eos, fld)
    if fld.t + fld.dt > fld.tmax 
        fld.dt = pr.t_max - fld.t 
    end
    t0 = time()
    it_num = calc(problem, eos, fld)
    t1 = time()
    delta_t = round(t1-t0, sigdigits=3)
    println("step=$iter: $(it_num) iters, t=$(round(fld.t, digits=3)), dt=$(round(fld.dt, sigdigits=2)), time=$(delta_t)s")
    global iter += 1
    fld.t += fld.dt
end



