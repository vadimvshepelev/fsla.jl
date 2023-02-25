include("problem.jl")
include("eos.jl")
include("field.jl")

abstract type Method end

######## Samarskii explicit-implicit finite difference method in Lagrangian mass coordinates

function calc{:samarskii}(pr::VacuumProblem, eos<:EOS, fld::Field_primitive)  



end


