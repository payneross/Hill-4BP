function state = hill4bpMcGeheeToCartesian(s, params)
%HILL4BPMCGEHEETOCARTESIAN  Convert [r;theta;v;w] to [x1;x2;x1dot;x2dot].

r = max(s(1), 0);
theta = s(2);
v = s(3);
w = s(4);

rho = r^params.gamma;
cTheta = cos(theta);
sTheta = sin(theta);

q1 = rho*cTheta;
q2 = rho*sTheta;

if r == 0
    pScale = Inf;
else
    pScale = r^(-params.gamma*params.beta);
end

p1 = pScale*(v*cTheta - w*sTheta);
p2 = pScale*(v*sTheta + w*cTheta);

xdot = p1 + q2;
ydot = p2 - q1;

state = [q1; q2; xdot; ydot];

end
