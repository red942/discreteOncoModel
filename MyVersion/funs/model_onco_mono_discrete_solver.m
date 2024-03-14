% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
%
% Difference equations for onco model monotherapy and the variational equations
%

function soltn_vals = model_onco_mono_discrete_solver(time_start, time_end, h, w_0, parameters)

global g_actual_curve;
global gcell_dosing;

l0 = parameters(1);
l1 = parameters(2);
k1 = parameters(3);
k2 = parameters(4);

switch gcell_dosing{1,1,g_actual_curve}
    case 0
        fun_c = @(t) 0;
    case 1 
        fun_c = @(t) fun_c_PKpo_1Komp(t);
    case 2
        fun_c = @(t) fun_c_PKpo_2Komp(t);
end

% Time vector (input vals)
t = linspace(time_start, time_end, (time_end-time_start)/h + 1);

% Initialize vectors to store results
x1 = zeros(size(t));
x2 = zeros(size(t));
x3 = zeros(size(t));
x4 = zeros(size(t));

% Set initial values
x1(1) = w_0;
x2(1) = 0;
x3(1) = 0;
x4(1) = 0;

% Time-stepping loop
for i = 1:length(t)-1
    w = x1(i) + x2(i) + x3(i) + x4(i);
    
    % Update the difference equations with the time
    x1(i+1) = x1(i) + h * ((2*l0*l1*x1(i)^2)/((l1+2*l0*x1(i))*w) - k2*fun_c(i)*x1(i));
    x2(i+1) = x2(i) + h * (k2*fun_c(i)*x1(i) - k1*x2(i));
    x3(i+1) = x3(i) + h * (k1*(x2(i) - x3(i)));
    x4(i+1) = x4(i) + h * (k1*(x3(i) - x4(i)));
end

% Store the results
soltn_vals = [x1', x2', x3', x4'];

