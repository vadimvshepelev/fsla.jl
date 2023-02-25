# include("problem.jl")
include("eos.jl")
include("field.jl")


rp = VacuumProblem(x_span=(0., 1.), t_span=(0.0, 0.2), w=[2650., 0., 35.6])

problem = Dict(
    "task" => rp, 
    "method" => "samarskii", 
    "N" => 100, 
    "CFL" => 0.5)

fld = field_primitive(problem)

println("Hello, world!!!")
