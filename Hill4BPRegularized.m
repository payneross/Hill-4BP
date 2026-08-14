function ds = Hill4BPRegularized(~, s, params)
%HILL4BPREGULARIZED  McGehee-regularized planar Hill 4BP equations.
%
% State s = [r; theta; v; w], and prime denotes d/dtau with dt = r dtau.

r = s(1);
theta = s(2);
v = s(3);
w = s(4);

if r < 0
    r = 0;
end

beta = params.beta;
gamma = params.gamma;
nu = params.nu;
alpha = params.alpha;
A = params.A;
B = params.B;
c = params.c;

radialPower = 2 - gamma*(nu + 2);

ds = [
    (beta + 1)*v*r;
    w - r;
    beta*v^2 + w^2 - alpha*c ...
        - nu*r^radialPower ...
        - 2*A*r^2*cos(theta)^2 ...
        - 2*B*r^2*sin(theta)^2;
    (beta - 1)*v*w + 2*(A - B)*r^2*sin(theta)*cos(theta)
];

end
