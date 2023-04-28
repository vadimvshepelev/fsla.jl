######## Samarskii explicit-implicit finite difference method in Lagrangian mass coordinates #######

function calc(pr::Dict, eos::eos_ideal, fld::Field_primitive)::Int8  
	# double h = ms[0].dm;
	# int itCounter = 0;
	# double *g = new double[nSize];
    N::Int64 = pr["N"]
    it_num::Int8, max_it::Int8 = 0, 30
    dm, dt = fld.dm, fld.dt
	visc_const::Float64 = 0.001 # 20.0 
	x::Vector{Float64}, x_new::Vector{Float64} = fld.x, fld.x_new
	w::Vector{Vector{Float64}}, w_new::Vector{Vector{Float64}} = fld.w, fld.w_new
    g::Vector{Float64} = cat([(visc_const * w[i][1] * (w[i+1][2] - w[i][2])) for i in 1:N], [0.0], dims=1)
    w1 = deepcopy(w)
	flag::Bool = (w1 === w)
    while it_num < 31
        # w1 = deepcopy(w_new)
		for i in 1:N, j in 1:3 
			w1[i][j] = w_new[i][j]
		end
        for i = 1:N       
			if i==1  
				p_next_plus  = w1[i][3]
				p_next_minus = 0.0
				p_plus       = w[i][3]
				p_minus      = 0.0
			elseif i==N 
				p_next_plus  = 0.0
				p_next_minus = w1[i-1][3]
				p_plus       = 0.0
				p_minus      = w[i-1][3]
            else
				p_next_plus  = w1[i][3] + g[i]
				p_next_minus = w1[i-1][3] + g[i-1]
				p_plus       = w[i][3] + g[i]
				p_minus      = w[i-1][3] + g[i-1]
			end
			w_new[i][2] = w[i][2] - 0.5*dt/dm*(p_next_plus - p_next_minus + p_plus - p_minus)
			x_new[i] = x[i] + 0.5*dt*(w_new[i][2] + w[i][2])
		end  		
		for i=1:N 
			w_new[i][1] = 1.0 / (1.0/w[i][1] + 0.5*dt/dm *(w_new[i+1][2] + w[i+1][2] - w_new[i][2] - w[i][2]))        
			rho::Float64, p::Float64 = w[i][1], w[i][3]
			e::Float64 = gete(eos, rho, p)
			e_new::Float64 = e - 0.25*dt/dm*(w1[i][3] + w[i][3] + g[i])*(w_new[i+1][2] + w[i+1][2] - w_new[i][2] - w[i][2]);
			w_new[i][3] = getp(eos, w_new[i][1], e_new)
		end
		it_num += 1
        all(isapprox.(w_new, w1; atol=0.001)) && break
 	end
	it_num == 31 && println("No convergence in $(max_it) iterations")  # throw(DomainError(-1, "No convergence in $(max_it) iterations"))
	# w = deepcopy(w_new)
	for i in 1:N, j in 1:3 
		w[i][j] = w_new[i][j]
	end


	i = findall(x -> x .> 0.001, w_new-w1)
		println(i)


	it_num
end

######## Calculation of time step according to CFL condition #######
function calcdt(pr::Dict, eos::T, fld::Field_primitive)::Float64 where {T<:EOS}
	x::Vector{Float64}, w::Vector{Vector{Float64}}, N::Int64, CFL::Float64 = fld.x, fld.w, fld.N, fld.CFL
	c_local::Float64 = getc(eos, w[1][1], w[1][3])
	imin::Int64, imax::Int64 = fld.imin, fld.imax
	dt::Float64 = CFL * (x[imin+1] - x[imin]) / (abs(w[imin][2]) + c_local)
	for i::Int64 = imin:imax+1
		c_local, dt = getc(eos, w[i][1], w[i][3]), min(dt, CFL * (x[i+1] - x[i]) / (abs(w[i][2]) + c_local))
	end
	dt
end
