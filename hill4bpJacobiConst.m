function C = hill4bpJacobiConst(state, params)
%HILL4BPJACOBICONST  Jacobi constant C = -2H for planar Hill 4BP states.

C = -2*hill4bpHamiltonian(state, params);

end
