include("problem.jl")
include("eos.jl")


rp = RiemannProblem(x=(0., 1.), wl=[0., 0., 0.], wr=[2650., 0., 35.6], xbnd=0.)
println(rp)