% System parameters
R = 10; L1=1e-3; L2 = 1e-3; C1=10e-6; C2 = 1e-6; Vin = 100; D = 0.5;

% System transfer function
nm = [R*Vin*C1*L1/(1-D) -(L1*Vin*D^2)/(1-D)^2 R*Vin];
dn = [R*C1*C2*L1*L2 C1*L1*L2 (R*C2*L2*(1-D)^2+R*C2*L1*D^2+R*C1*L1) (L2*(1-D)^2+L1*D^2) R*(1-D)^2];
G = tf(nm,dn);

% PSO parameters
numParticles = 30;
numIterations = 100;
lb = [0, 0]; % Lower bounds for kp and kr
ub = [10, 10]; % Upper bounds for kp and kr

% Objective function
objectiveFcn = @(params) evaluateMargins(params, G);

% PSO optimization
options = optimoptions('particleswarm', 'SwarmSize', numParticles, 'MaxIterations', numIterations, 'Display', 'iter');
optimalParams = particleswarm(objectiveFcn, 2, lb, ub, options);

% Optimal kp and kr
optimalKp = optimalParams(1);
optimalKr = optimalParams(2);

% Display results
fprintf('Optimal kp: %.4f\n', optimalKp);
fprintf('Optimal kr: %.4f\n', optimalKr);

% Function to evaluate gain and phase margins
function cost = evaluateMargins(params, G)
    kp = params(1);
    kr = params(2);
    Gpr = tf([kp kr kp*((2*pi*50)^2)], [1 0 (2*pi*50)^2]);
    sys = series(Gpr, G);
    [Gm, Pm] = margin(sys);
    desiredGm = 20; % Desired Gain Margin in dB
    desiredPm = 60; % Desired Phase Margin in degrees
    if isnan(Gm) || isnan(Pm)
        cost = 1e6; % Large penalty for non-computable margins
    else
        % Cost function considering desired Gm and Pm, prioritizing Gm first
        cost = (max(0, desiredGm - Gm))^2 + (max(0, desiredPm - Pm))^2;
    end
end

% Plot Bode plot of optimal system
figure;
Gpr_opt = tf([optimalKp optimalKr optimalKp*((2*pi*50)^2)], [1 0 (2*pi*50)^2]);
sys_opt = series(Gpr_opt, G);
bode(sys_opt);
grid on;
margin(sys_opt);
grid on;
