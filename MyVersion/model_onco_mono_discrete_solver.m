% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
%
% Difference equations for onco model monotherapy and the variational equations
%
function soltn_vals = model_onco_mono_discrete_solver(time_range, h, w_0, parameters)

l0 = parameters(1);
l1 = parameters(2);
k1 = parameters(3);
k2 = parameters(4);
%

%switch gcell_dosing{1,1,g_actual_curve}
%    case 0
%        fun_c = 0;
%    case 1 
%        fun_c = @(t) fun_c_PKpo_1Komp(t);
%    case 2
%        fun_c = @(t) fun_c_PKpo_2Komp(t);
%end
fun_c = @(t) 0;

% Time vector
t = linspace(time_range(1), time_range(2), (time_range(2)-time_range(1))/h + 1);

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
    
    % Evaluate c(t) at the current time
    current_c = fun_c(i);
    
    % Update the difference equations with the time-dependent c(t)
    x1(i+1) = x1(i) + h * ((2*l0*l1*x1(i)^2)/((l1+2*l0*x1(i))*w) - k2*current_c*x1(i));
    x2(i+1) = x2(i) + h * (k2*current_c*x1(i) - k1*x2(i));
    x3(i+1) = x3(i) + h * (k1*(x2(i) - x3(i)));
    x4(i+1) = x4(i) + h * (k1*(x3(i) - x4(i)));
end

% Store the results in matrix dxx
soltn_vals = [x1', x2', x3', x4'];

