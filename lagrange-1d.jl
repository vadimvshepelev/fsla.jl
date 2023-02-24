include("problem.jl")
include("eos.jl")
include("field.jl")


rp = RiemannProblem(x_span=(0., 1.), t_span=(0.0, 0.2) wl=[0., 0., 0.], wr=[2650., 0., 35.6], xbnd=0.)

problem = Dict(
    "task" => rp, 
    "method" => "samarskii", 
    "N" => 100, 
    "CFL" => 0.5)

fld = field_primitive(problem)

println("Hello, world!!!")
