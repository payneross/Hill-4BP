function s = hill4bpCartesianToMcGehee(state, params)
%HILL4BPCARTESIANTOMCGEHEE  Convert [x1;x2;x1dot;x2dot] to [r;theta;v;w].

q1 = state(1);
q2 = state(2);
xdot = state(3);
ydot = state(4);

rho = sqrt(q1^2 + q2^2);
if rho == 0
    error('hill4bpCartesianToMcGehee:collision', ...
        'McGehee coordinates require a nonzero Cartesian position.');
end

theta = atan2(q2, q1);
r = rho^(1/params.gamma);

% Canonical momenta from qdot = [p1 + x2; p2 - x1].
p1 = xdot - q2;
p2 = ydot + q1;

scale = r^(params.gamma*params.beta);
v = scale*(p1*cos(theta) + p2*sin(theta));
w = scale*(-p1*sin(theta) + p2*cos(theta));

s = [r; theta; v; w];

end
