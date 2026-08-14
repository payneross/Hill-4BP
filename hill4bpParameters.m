function params = hill4bpParameters()
%HILL4BPPARAMETERS  Shared parameters for the planar oblate Hill 4BP model.

values = readHill4BPParameterFile();

params.G = requiredValue(values, 'G');
params.mu = requiredValue(values, 'mu');
params.u1 = requiredValue(values, 'u1');
params.u2 = requiredValue(values, 'u2');
params.c1 = requiredValue(values, 'c1');
params.c2 = requiredValue(values, 'c2');
params.c3 = requiredValue(values, 'c3');
params.c = -params.c3;

mu = params.mu;
u1 = params.u1;
u2 = params.u2;
Delta = (mu*u1^3 + (1 - mu)*u2^3)^2 - mu*(1 - mu)*u1*u2 ...
    * (-u1^4 - u2^4 + 2*u1^2 + 2*u2^2 + 2*u1^2*u2^2 - 1);
base = 2 - 2*(1 - mu)/u1^5 - 2*mu/u2^5 + 3*(1 - mu)/u1^3 + 3*mu/u2^3;
computedLambda1 = 0.5*(base - 3*sqrt(Delta)/(u1^3*u2^3));
computedLambda2 = 0.5*(base + 3*sqrt(Delta)/(u1^3*u2^3));

params.lambda1 = optionalValue(values, 'lambda1', computedLambda1);
params.lambda2 = optionalValue(values, 'lambda2', computedLambda2);
params.A = (1 - params.lambda2)/2;
params.B = (1 - params.lambda1)/2;

params.nu = requiredValue(values, 'nu');
params.alpha = requiredValue(values, 'alpha');
params.beta = params.alpha/2;
params.gamma = 2/(params.alpha + 2);

params.poincare = struct( ...
    'hE1Scale', requiredValue(values, 'poincare_h_e1_scale'), ...
    'k', integerValue(values, 'poincare_k'), ...
    'l', integerValue(values, 'poincare_l'), ...
    'iterates', integerValue(values, 'poincare_iterates'), ...
    'xBegin', requiredValue(values, 'poincare_x_begin'), ...
    'xEnd', requiredValue(values, 'poincare_x_end'), ...
    'xdotLower', requiredValue(values, 'poincare_xdot_lower'), ...
    'xdotUpper', requiredValue(values, 'poincare_xdot_upper'), ...
    'ySection', requiredValue(values, 'poincare_y_section'), ...
    'tauStep', requiredValue(values, 'poincare_tau_step'), ...
    'relTol', requiredValue(values, 'poincare_rel_tol'), ...
    'absTol', requiredValue(values, 'poincare_abs_tol'));

end

function values = readHill4BPParameterFile()
configPath = hill4bpParameterFile();
lines = regexp(fileread(configPath), '\r\n|\n|\r', 'split');
values = struct();

for i = 1:numel(lines)
    line = regexprep(lines{i}, '#.*$', '');
    line = strtrim(line);
    if isempty(line)
        continue
    end

    tokens = regexp(line, '^([A-Za-z][A-Za-z0-9_]*)\s*=\s*([-+0-9.eEdD]+)$', 'tokens', 'once');
    if isempty(tokens)
        error('hill4bpParameters:invalidConfigLine', ...
            'Invalid line in %s at line %d: %s', configPath, i, lines{i});
    end

    value = str2double(strrep(strrep(tokens{2}, 'D', 'e'), 'd', 'e'));
    if isnan(value)
        error('hill4bpParameters:invalidConfigValue', ...
            'Invalid numeric value in %s at line %d: %s', configPath, i, lines{i});
    end
    values.(tokens{1}) = value;
end
end

function configPath = hill4bpParameterFile()
thisDir = fileparts(mfilename('fullpath'));
candidates = { ...
    fullfile(thisDir, 'Hill_4BP_parameters.cfg'), ...
    fullfile(thisDir, '..', 'Hill_4BP_parameters.cfg'), ...
    fullfile(thisDir, '..', '..', 'Hill_4BP_parameters.cfg')};

for i = 1:numel(candidates)
    candidate = candidates{i};
    if exist(candidate, 'file') == 2
        configPath = candidate;
        return
    end
end

error('hill4bpParameters:missingConfig', ...
    'Could not find Hill_4BP_parameters.cfg near %s.', thisDir);
end

function value = requiredValue(values, name)
if ~isfield(values, name)
    error('hill4bpParameters:missingConfigValue', ...
        'Missing required Hill 4BP parameter "%s".', name);
end
value = values.(name);
end

function value = optionalValue(values, name, defaultValue)
if isfield(values, name)
    value = values.(name);
else
    value = defaultValue;
end
end

function value = integerValue(values, name)
value = requiredValue(values, name);
if value < 1 || abs(value - round(value)) > eps(value)
    error('hill4bpParameters:invalidIntegerConfigValue', ...
        'Hill 4BP parameter "%s" must be a positive integer.', name);
end
value = round(value);
end
