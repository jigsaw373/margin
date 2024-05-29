% System parameters
R = 10; L1 = 1e-3; L2 = 1e-3; C1 = 10e-6; C2 = 1e-6; Vin = 100; D = 0.5;

% System transfer function
nm = [R*Vin*C1*L1/(1-D) -(L1*Vin*D^2)/(1-D)^2 R*Vin];
dn = [R*C1*C2*L1*L2 C1*L1*L2 (R*C2*L2*(1-D)^2+R*C2*L1*D^2+R*C1*L1) (L2*(1-D)^2+L1*D^2) R*(1-D)^2];
G = tf(nm,dn);

% PSO parameters
numParticles = 50;
numIterations = 100;
lb = [0, 0]; % Lower bounds for kp and kr
ub = [10, 10]; % Upper bounds for kp and kr

% Number of PSO runs
numRuns = 15;

% Initialize best results
bestKp = NaN;
bestKr = NaN;
bestCost = Inf;

% Objective function
objectiveFcn = @(params) evaluateMargins(params, G);

% Run PSO multiple times
for run = 1:numRuns
    fprintf('Running PSO iteration %d...\n', run);
    options = optimoptions('particleswarm', 'SwarmSize', numParticles, 'MaxIterations', numIterations);
    [optimalParams, fval] = particleswarm(objectiveFcn, 2, lb, ub, options);
    
    % Check if this run's result is better
    if fval < bestCost
        bestCost = fval;
        bestKp = optimalParams(1);
        bestKr = optimalParams(2);
    end
    
    fprintf('Iteration %d: Optimal kp: %.4f, Optimal kr: %.4f, Cost: %.4f\n', run, optimalParams(1), optimalParams(2), fval);
end

% Display best results
fprintf('Best kp: %.4f\n', bestKp);
fprintf('Best kr: %.4f\n', bestKr);
fprintf('Best Cost: %.4f\n', bestCost);

% Function to evaluate gain and phase margins
function cost = evaluateMargins(params, G)
    kp = params(1);
    kr = params(2);
    Gpr = tf([kp kr kp*((2*pi*50)^2)], [1 0 (2*pi*50)^2]);
    sys = series(Gpr, G);
    [Gm, Pm] = margin(sys);
    desiredGm = 20; % Desired Gain Margin in dB
    desiredPm = 60; % Desired Phase Margin in degrees
    if isnan(Gm) || isnan(Pm) || Gm < 0 || Pm < 0
        cost = 1e6; % Large penalty for non-computable margins
    else
        % Cost function considering desired Gm and Pm
        cost = (abs(Gm - desiredGm))^2 + (abs(Pm - desiredPm))^2;
    end
end

% Plot Bode plot of best system
figure;
Gpr_best = tf([bestKp bestKr bestKp*((2*pi*50)^2)], [1 0 (2*pi*50)^2]);
sys_best = series(Gpr_best, G);
bode(sys_best);
grid on;
margin(sys_best);
grid on;

% Manual verification
kp_manual = 0;
kr_manual = 0.95;
Gpr_manual = tf([kp_manual kr_manual kp_manual*((2*pi*50)^2)], [1 0 (2*pi*50)^2]);
sys_manual = series(Gpr_manual, G);
figure;
bode(sys_manual);
grid on;
margin(sys_manual);
grid on;