function stateDot = Hill4BP(~, state, params)
%HILL4BP  Singular planar Hill 4BP equations in Cartesian velocity form.
%
% This is useful for reference and diagnostics.  The Poincare script uses
% Hill4BPRegularized instead because the Cartesian equations are singular
% at collisions with the tertiary.

x = state(1);
y = state(2);
xdot = state(3);
ydot = state(4);

rho = sqrt(x^2 + y^2);
if rho == 0
    error('Hill4BP:collision', ...
        'The Cartesian Hill 4BP vector field is singular at the origin.');
end

gravityFactor = 1/rho^3 + 3*params.c/rho^5;

xddot = 2*ydot + (params.lambda2 - gravityFactor)*x;
yddot = -2*xdot + (params.lambda1 - gravityFactor)*y;

stateDot = [xdot; ydot; xddot; yddot];

end
