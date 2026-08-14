clear
format long

% Planar Hill 4BP with oblate tertiary, regularized with McGehee variables.
% Energy convention: H = h and Jacobi C = -2*h.
params = hill4bpParameters();
equilibria = hill4bpEquilibria(params);
poincare = params.poincare;

% Default section parameters come from ../Hill_4BP_parameters.cfg.
% The E1/E2 threshold is the first neck-opening energy.  
h = poincare.hE1Scale * equilibria.E1.h;
jacobiC = -2*h;

fprintf('Hill 4BP Hamiltonian h = %.15g, Jacobi C = %.15g\n', h, jacobiC);
fprintf('E1/E2 h = %.15g, E3/E4 h = %.15g\n', equilibria.E1.h, equilibria.E3.h);

% Increase k for finer initial-condition sampling.
% Increase iterates to collect more intersections per initial condition.
k = poincare.k;
l = poincare.l;
iterates = poincare.iterates;

% Initial section y = 0, ydot > 0, plotted as (x, xdot).
x_begin = poincare.xBegin;
x_end = poincare.xEnd;
xdot_lower = poincare.xdotLower;
xdot_upper = poincare.xdotUpper;

y_section = poincare.ySection;
tauStep = poincare.tauStep;

if k > 1
    xStep = (x_end - x_begin)/(k - 1);
else
    xStep = 0;
end

if l > 1
    xdotStep = (xdot_upper - xdot_lower)/(l - 1);
else
    xdotStep = 0;
end

maxIterates = l*k*iterates;
PoincareMap = NaN(maxIterates, 2);
totalIterates = 0;
skippedInitialConditions = 0;
truncatedInitialConditions = 0;
maxSegmentsPerInitialCondition = 5*iterates + 100;

odeOptions = odeset( ...
    'RelTol', poincare.relTol, ...
    'AbsTol', poincare.absTol, ...
    'Events', @(tau, s) hill4bpSectionEvent(tau, s, params));

for m = 1:l
    for n = 1:k
        fprintf('velocity row %d/%d, x node %d/%d, crossings %d\n', m, l, n, k, totalIterates);

        x0 = x_begin + (n - 1)*xStep;
        xdotDirection = xdot_lower + (m - 1)*xdotStep;

        q0 = [x0; y_section];
        direction = [xdotDirection; 1.0];
        direction = direction/norm(direction);

        [speed, allowed] = hill4bpSpeed(q0, h, params);
        if ~allowed
            skippedInitialConditions = skippedInitialConditions + 1;
            continue
        end

        v0 = speed*direction;
        initialCartesian = [q0; v0];
        initialMcGehee = hill4bpCartesianToMcGehee(initialCartesian, params);

        localIterates = 0;
        segmentCount = 0;
        while localIterates < iterates && segmentCount < maxSegmentsPerInitialCondition
            segmentCount = segmentCount + 1;

            [~, trajectory, ~, eventStates, ~] = ode113( ...
                @(tau, s) Hill4BPRegularized(tau, s, params), ...
                [0 tauStep], initialMcGehee, odeOptions);

            if isempty(trajectory) || any(~isfinite(trajectory(end, :)))
                break
            end

            for eventIndex = 1:size(eventStates, 1)
                eventCartesian = hill4bpMcGeheeToCartesian(eventStates(eventIndex, :).', params);

                if eventCartesian(4) <= 0
                    continue
                end

                localIterates = localIterates + 1;
                totalIterates = totalIterates + 1;

                if totalIterates > size(PoincareMap, 1)
                    PoincareMap = [PoincareMap; NaN(maxIterates, 2)];
                end

                PoincareMap(totalIterates, :) = eventCartesian([1, 3]).';

                if localIterates >= iterates
                    break
                end
            end

            initialMcGehee = trajectory(end, :).';
        end

        if localIterates < iterates
            truncatedInitialConditions = truncatedInitialConditions + 1;
        end
    end
end

PoincareGrid = PoincareMap(~isnan(PoincareMap(:, 1)), :);

% Zero-velocity boundary in the plotted section coordinates.
% On y = 0, H = h gives ydot^2 = 2*(h + Omega(x,0)) - xdot^2.
% The boundary of the allowed region in the (x, xdot) plot is therefore
% xdot = +/-sqrt(2*(h + Omega(x,0))). Outside these curves, ydot is not real.
zvbX = linspace(x_begin, x_end, 4000).';
zvbSpeedSquared = NaN(size(zvbX));
for i = 1:numel(zvbX)
    zvbSpeedSquared(i) = 2*(h + hill4bpEffectivePotential([zvbX(i); 0], params));
end
zvbAllowed = zvbSpeedSquared >= 0 & isfinite(zvbSpeedSquared);
ZeroVelocityBoundary = [ ...
    zvbX(zvbAllowed), sqrt(max(zvbSpeedSquared(zvbAllowed), 0)); ...
    NaN, NaN; ...
    zvbX(zvbAllowed), -sqrt(max(zvbSpeedSquared(zvbAllowed), 0))];

save closeData2 PoincareGrid h jacobiC params equilibria ...
    skippedInitialConditions truncatedInitialConditions ZeroVelocityBoundary

writematrix(PoincareGrid, 'poincare_grid.csv')
writematrix(ZeroVelocityBoundary, 'zero_velocity_boundary.csv')

figure
plot(PoincareGrid(:, 1), PoincareGrid(:, 2), 'b.', 'MarkerSize', 4)
hold on
plot(ZeroVelocityBoundary(:, 1), ZeroVelocityBoundary(:, 2), 'r-', 'LineWidth', 1.25)
hold off
xlabel('x_1')
ylabel('dx_1/dt')
title(sprintf('Planar regularized Hill 4BP Poincare section, h = %.6g', h))
legend('Poincare iterates', 'zero-velocity boundary', 'Location', 'best')
grid on
