function equilibria = hill4bpEquilibria(params)
%HILL4BPEQUILIBRIA  Planar Hill 4BP equilibria and their energies.

gamma = params.gamma;
nu = params.nu;
alpha = params.alpha;
c = params.c;

r1 = equilibriumRadius(params.lambda2, params);
r2 = equilibriumRadius(params.lambda1, params);

rho1 = r1^gamma;
rho2 = r2^gamma;

equilibria.E1.mcgehee = [r1; 0; 0; r1];
equilibria.E2.mcgehee = [r1; pi; 0; r1];
equilibria.E3.mcgehee = [r2; pi/2; 0; r2];
equilibria.E4.mcgehee = [r2; 3*pi/2; 0; r2];

equilibria.E1.cartesian = [rho1; 0; 0; 0];
equilibria.E2.cartesian = [-rho1; 0; 0; 0];
equilibria.E3.cartesian = [0; rho2; 0; 0];
equilibria.E4.cartesian = [0; -rho2; 0; 0];

equilibria.E1.h = hill4bpHamiltonian(equilibria.E1.cartesian, params);
equilibria.E2.h = equilibria.E1.h;
equilibria.E3.h = hill4bpHamiltonian(equilibria.E3.cartesian, params);
equilibria.E4.h = equilibria.E3.h;

equilibria.collisionPlus = [0; NaN; sqrt(2*c); 0];
equilibria.collisionMinus = [0; NaN; -sqrt(2*c); 0];

end

function r = equilibriumRadius(lambda, params)
nu = params.nu;
alpha = params.alpha;
c = params.c;
gamma = params.gamma;
power = 2 - gamma*(nu + 2);

f = @(r) lambda*r.^2 - nu*r.^power - alpha*c;

lo = eps;
hi = max(1, lambda^(-1/(2 - power))*2);
while f(hi) <= 0
    hi = 2*hi;
end

r = fzero(f, [lo, hi]);
end
