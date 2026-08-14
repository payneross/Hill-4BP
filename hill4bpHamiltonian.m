function H = hill4bpHamiltonian(state, params)
%HILL4BPHAMILTONIAN  Hamiltonian for planar Hill 4BP velocity states.
%
% state = [x1; x2; x1dot; x2dot].

q = state(1:2);
v = state(3:4);
H = 0.5*(v(1)^2 + v(2)^2) - hill4bpEffectivePotential(q, params);

end
