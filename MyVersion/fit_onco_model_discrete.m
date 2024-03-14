%temp test file
clear all;  
clear global all;

global g_t_data g_y_data g_fid g_param_temp; 
global gv_start_vals;
global gs_c_info;
global gs_wr_2_disc gs_wr_2_disp;
global g_n_curve g_anz_data_all gv_anz_data;
global g_actual_curve;
global g_model;
global gs_LargeScale;
global timescale;
timescale = 1;

addpath('data','funs'); %lets matlab see these folders

% =====   Load data file
% 
% Choose a data set and run the program
%

%struct_data = fun_data_control_RKO();
%struct_data = fun_data_control_PC3();
%struct_data = fun_data_control_MDA();
%struct_data = fun_data_control_A549();

%struct_data = fun_data_control_1();
%struct_data = fun_data_control_2();

%struct_data = fun_data_mono_drug_A1_180();
%struct_data = fun_data_mono_drug_A2_120();
struct_data = fun_data_mono_drug_B_100();
%struct_data = fun_data_mono_drug_C_100();
%struct_data = fun_data_mono_drug_C_150();
%struct_data = fun_data_mono_drug_C_100_150();

%struct_data = fun_data_comb_drug_A1_180_and_C_100();
%struct_data = fun_data_comb_drug_A2_120_and_C_100();
%struct_data = fun_data_comb_drug_B_100_and_C_150();
%
%=================================


% Plot the object data
figure(1); plot(struct_data(1).t_data,struct_data(1).y_data,'red x'); hold on;
for i=2:g_n_curve
    figure(1); plot(struct_data(i).t_data,struct_data(i).y_data,'black s'); hold on;
end
xlabel('t (days)'); ylabel('w(t) (cm^3)');
if (g_model == 0)
    h_leg = legend(gs_c_info);
else
    h_leg = legend('Control Data',gs_c_info);
end
set(h_leg,'Interpreter','none');

% Read data and create column vectors g_t_data and g_y_data         
g_anz_data_all = 0;
for i=1:g_n_curve
    gv_anz_data(i) = struct_data(i).anz_data;     % Number of data of i-th object
    for j=1:gv_anz_data(i)
        g_t_data(j+(i-1)*gv_anz_data(i)) = struct_data(i).t_data(j);    
        g_y_data(j+(i-1)*gv_anz_data(i)) = struct_data(i).y_data (j);
    end
    g_anz_data_all = g_anz_data_all + gv_anz_data(i);
end





%Displays certain things based on if variables are set to on or off

%case off
if strcmp(gs_wr_2_disp,'off')
    if (g_model == 0)
        text = sprintf('Inital values\nl0 = %12.8f, l1 = %12.8f, w0 = %12.8f', ...
                       gv_start_vals(1),gv_start_vals(2),gv_start_vals(5));
        disp(text);
    end
    if (g_model == 1)
        text = sprintf('Inital values\nl0 = %12.8f, l1 = %12.8f, k1 = %12.8f, k2 = %12.8f, w0 = %12.8f', ...
                       gv_start_vals(1),gv_start_vals(2),gv_start_vals(3),gv_start_vals(4),gv_start_vals(5));
        disp(text);
    end
    if (g_model == 2)
        text = sprintf('Inital values\nl0 = %12.8f, l1 = %12.8f, k1 = %12.8f, psi = %12.8f, w0 = %12.8f', ...
                       gv_start_vals(1),gv_start_vals(2),gv_start_vals(3),gv_start_vals(4),gv_start_vals(5));
        disp(text);
    end
    disp('calculating...');
end

%case on
if strcmp (gs_wr_2_disc,'on')
    cd results2;
    datei_name = strcat('fit_process_',gs_c_info,'.txt');
    g_fid = fopen(datei_name,'w');
    cd ..;
end

% Lower and upper bounds for parameter estimation 
lb = [0, 0, 0, 0, 0]; 
ub = [2, 2, 2, 2, 2];  

%[optimized_params,resnorm,residual,exitflag,output]  = ... 
%    lsqcurvefit(@fun_discrete_eqn_vals,gv_start_vals,...
%                g_t_data,g_y_data,lb,ub,options_lsq);

% Define the objective function
objective_function = @fun_discrete_eqn_vals;

% Define optimization options
lsq_ops = [1e-10, 1e-10, 1e7, 1e7];   % TolFun, TolX, MaxIter,MaxFunEvals
options_lsq = optimset('TolFun',lsq_ops(1),...
                       'TolX',lsq_ops(2),'MaxIter',lsq_ops(3),...
                       'MaxFunEvals',lsq_ops(4),...
                       'LargeScale',gs_LargeScale);
tic;

% Perform optimization
optimized_params = lsqcurvefit(objective_function, gv_start_vals, g_t_data, g_y_data, lb, ub, options_lsq);

% Calculate the fitted values using the optimized parameters
fitted_values = fun_discrete_eqn_vals(optimized_params, g_t_data);

% An atttempt at rerunning to get it to not stop until it works
% Define custom optimization stopping criteria
%threshold = 0.05; % 5% threshold
%tolerance = threshold * ones(size(g_y_data)); % Tolerance for each data point
%while any(abs(fitted_values - g_y_data) > tolerance)
%    % If any fitted value is not within the tolerance, re-optimize
%    optimized_params = lsqcurvefit(objective_function, optimized_params, g_t_data, g_y_data, lb, ub, options);%
%    fitted_values = fun_discrete_eqn_vals(optimized_params, g_t_data);
%end

time = toc;

if strcmp (gs_wr_2_disc,'on')
    fclose(g_fid);
end

% Write resulting solver properties 
text = sprintf('\nCalculation time: %8.0f (secs)',time);
disp(text);
disp(g_param_temp)


% Compute solutions of different objects and store them in a struct      
for k=1:g_n_curve
    g_param_temp = [optimized_params(1), optimized_params(2), optimized_params(3), optimized_params(4)];    
    t_end = struct_data(k).t_data(length(struct_data(k).t_data));
    g_actual_curve = k;
    %x0 = zeros(24,1);
    %x0(1) = optimized_params(5); x0(21) = 1;
    x0 = optimized_params(5);
    if (g_model == 0 || g_model == 1)
        %struct_sol(k) = ode45(@model_onco_mono, [0,t_end+1], ...    z1 
        solutions = model_onco_mono_discrete_solver(0, t_end + 1, timescale, x0, g_param_temp);
        w_t = sum(solutions', 1)';
        t_values = linspace(0, t_end, length(w_t));
        struct_sol(k) = struct('x', t_values, 'y', w_t);
    else
        struct_sol(k) = ode45(@model_onco_comb, [0,t_end+1], ... 
                        x0,options);
    end
    
    
    %T = struct_sol(k).x ;
    %Y = (struct_sol(k).y)';
    %figure(1); plot(T,Y(:,1)+Y(:,2)+Y(:,3)+Y(:,4),'black'); hold on;
    plot(t_values, w_t, 'black'); hold on;
end


% Calculate statistics
%stat = fun_stats(optimized_params,struct_sol,struct_data);


% Plot fitting results
%void = fun_print_result(optimized_params,g_model,stat);