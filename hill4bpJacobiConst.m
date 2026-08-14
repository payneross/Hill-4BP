function C = hill4bpJacobiConst(state, params)
% Jacobi constant C = -2H for planar Hill 4BP 

C = -2*hill4bpHamiltonian(state, params);

end
