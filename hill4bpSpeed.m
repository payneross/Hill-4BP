function [speed, allowed] = hill4bpSpeed(q, h, params)
% Velocity magnitude on the Hamiltonian level H=h

speedSquared = 2*(h + hill4bpEffectivePotential(q, params));
allowed = speedSquared >= 0 && isfinite(speedSquared);

if allowed
    speed = sqrt(max(speedSquared, 0));
else
    speed = NaN;
end

end
