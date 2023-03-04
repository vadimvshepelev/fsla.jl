include("problem.jl")
include("eos.jl")
include("field.jl")

abstract type Method end

######## Samarskii explicit-implicit finite difference method in Lagrangian mass coordinates

function calc{:samarskii}(pr::VacuumProblem, eos<:EOS, fld::Field_primitive)  
	# double h = ms[0].dm;
	# int itCounter = 0;
	# double *g = new double[nSize];
    N::Int64 = pr["N"]
    it_num::Int8, max_it::Int8 = 0, 30
    dm, dt = fld.dm, fld.dt
	visc_const::Float64 = 20.0 
    g::Vector{Float64} = [visc_const*w[i][0]*(w[i+1][1]-w[i][1]) for i in 1:N] * [0.]
    w1=w.copy()
    while it_num < 31
        w_new=w1.copy
        for i = 1:N+1  begin       
			if i==1  
				p_next_plus  = w1[i][3];
				p_next_minus = 0.0;
				p_plus       = w[i][3];
				p_minus      = 0.0;
			elseif i==N 
				p_next_plus  = 0.0;
				p_next_minus = w1[i-1][3];
				p_plus       = 0.0;
				p_minus      = w[i-1][3];
            else
				p_next_plus  = w1[i][3] + g[i];
				p_next_minus = w1[i-1][3] + g[i-1];
				p_plus       = w[i][3] + g[i];
				p_minus      = w[i-1][3] + g[i-1];
			end
			w_new[i][2] = w[i][2] - 0.5*dt/dm*(p_next_plus - p_next_minus + p_plus - p_minus)
			x_new[i] = x[i] + 0.5*dt*(w_new[i][2] + w[i][2])
            it_num += 1
            approx.(w_new, w) && break
		end
    end
    it_num == 31 && throw(ConvergenceError("No convergence in $(max_it) iterations"))
    for i=1:N+1 
        w_new[i][0] = 1.0 / (1.0/w[i][0] + 0.5*dt/dm *(w_new[i+1][2] + w[i+1][2] - w_new[i][2] - w[i][2]))
        
        w_new[ms_temp[i].ei = ms[i].ei - tau/4.0/h*(ms_prev[i].pi + ms[i].pi + g[i]) *
                        (ms_temp[i+1].v + ms[i+1].v - ms_temp[i].v - ms[i].v);
        ms_temp[i].ee = ms[i].ee - tau/4.0/h*(ms_prev[i].pe + ms[i].pe + g[i]) *
                        (ms_temp[i+1].v + ms[i+1].v - ms_temp[i].v - ms[i].v);
        ms_temp[i].e  = ms_temp[i].ei + ms_temp[i].ee;
    }


    # do {


		for(int i=0; i<nSize; i++) {
			ms_temp[i].ro = 1.0/(1.0/ms[i].ro + tau/2.0/h*
							(ms_temp[i+1].v + ms[i+1].v - ms_temp[i].v-ms[i].v));
			ms_temp[i].ei = ms[i].ei - tau/4.0/h*(ms_prev[i].pi + ms[i].pi + g[i]) *
							(ms_temp[i+1].v + ms[i+1].v - ms_temp[i].v - ms[i].v);
			ms_temp[i].ee = ms[i].ee - tau/4.0/h*(ms_prev[i].pe + ms[i].pe + g[i]) *
							(ms_temp[i+1].v + ms[i+1].v - ms_temp[i].v - ms[i].v);
			ms_temp[i].e  = ms_temp[i].ei + ms_temp[i].ee;
		}
		for(int i=0; i<nSize; i++)	{
			ms_temp[i].ti = eos.getti(ms_temp[i].ro, ms_temp[i].ei);
			if(ms_temp[i].ti == -1.) {
				cout << "Ooops! No solution of non-linear equation of hydrodynamics in node " << i << endl
					 << "ro = " << ms[i].ro << endl
					 << "ti = " << ms[i].ti << endl
					 << "te = " << ms[i].te << endl
					 << "e = "  << ms[i].e  << endl
					 << "ei = " << ms[i].ei << endl
					 << "ee = " << ms[i].ee << endl
					 << "p = "  << ms[i].p  << endl
					 << "pi = " << ms[i].pi << endl
					 << "pe = " << ms[i].pi << endl;
				saveSolution("error.dat", t);
				cin.get();
			}
			ms_temp[i].te = eos.gette(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].ee);
		}
		for(int i=0; i<nSize; i++) 	{
			ms_temp[i].pi = eos.getpi(ms_temp[i].ro, ms_temp[i].ti);
			ms_temp[i].pe = eos.getpe(ms_temp[i].ro, ms_temp[i].ti, ms_temp[i].te);
			ms_temp[i].p  = ms_temp[i].pi + ms_temp[i].pe;
		}
	}
	while (compTi(ms_temp, ms_prev) > .01);
	for(int i=0; i<nSize; i++) {
		 ms[i].x = ms_temp[i].x;
		 ms[i].v = ms_temp[i].v;
		ms[i].ro = ms_temp[i].ro;
		 ms[i].e = ms_temp[i].e;
		ms[i].ee = ms_temp[i].ee;
		ms[i].ei = ms_temp[i].ei;
		ms[i].p  = ms_temp[i].p;
		ms[i].pi = ms_temp[i].pi;
		ms[i].pe = ms_temp[i].pe;
		ms[i].te = ms_temp[i].te;
		ms[i].ti = ms_temp[i].ti;
		ms[i].C     = eos.getC(ms[i].ro, ms[i].ti, ms[i].te);
		ms[i].Alphaei	    = eos.getAlpha(ms[i].ro, ms[i].ti, ms[i].te);
		ms[i].kappa = eos.getkappa(ms[i].ro, ms[i].ti, ms[i].te);
		ms[i].ce    = eos.getce(ms[i].ro, ms[i].te);
	}
	ms[nSize].x = ms_temp[nSize].x;
	ms[nSize].v = ms_temp[nSize].v;
	delete[] g;
	return 1;
}


end


