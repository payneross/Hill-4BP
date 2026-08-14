function Omega = hill4bpEffectivePotential(q, params)
% Planar effective potential Omega(q).
% In velocity coordinates, H = 1/2*|qdot|^2 - Omega(q), so a point is
% energetically allowed on H=h when h + Omega(q) >= 0.

x = q(1);
y = q(2);
rho = sqrt(x^2 + y^2);

if rho == 0
    Omega = Inf;
    return
end

Omega = 0.5*(params.lambda2*x^2 + params.lambda1*y^2) ...
    + 1/rho + params.c/rho^3;

end
