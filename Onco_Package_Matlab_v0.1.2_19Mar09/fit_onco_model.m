% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
%
% Fitting program for the onco model 
% (Optimization toolbox is required)
%
%
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
struct_data = fun_data_mono_drug_A2_120();
%struct_data = fun_data_mono_drug_B_100();
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

disp(g_n_curve)
% Solver options
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
lsq_jac = 'off';                     % Jacobian from user or not
lsq_ops = [1e-8, 1e-8, 1e5, 1e5];   % TolFun, TolX, MaxIter,MaxFunEvals
options_lsq = optimset('Jacobian',lsq_jac,'TolFun',lsq_ops(1),...
                       'TolX',lsq_ops(2),'MaxIter',lsq_ops(3),...
                       'MaxFunEvals',lsq_ops(4),...
                       'LargeScale',gs_LargeScale);


% Lower and upper bounds for parameter estimation 
lb = []; 
ub = [];  


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
if strcmp (gs_wr_2_disc,'on')
    cd results;
    datei_name = strcat('fit_process_',gs_c_info,'.txt');
    g_fid = fopen(datei_name,'w');
    cd ..;
end


% Perform fitting process with lsqcurvefit
tic;
[erg,resnorm,residual,exitflag,output]  = ... 
    lsqcurvefit(@fun_F_and_Jac_mono_or_comb,gv_start_vals,...
                g_t_data,g_y_data,lb,ub,options_lsq);                
time = toc;

if strcmp (gs_wr_2_disc,'on')
    fclose(g_fid);
end


% Write resulting solver properties 
text = sprintf('\nCalculation time: %8.0f (secs)',time);
disp(text);
text = sprintf('\nlsqcurvefit results: Iterations: %8.0f, FuncCount: %8.0f, Algortihm: %s',...
output.iterations,output.funcCount,output.algorithm);
disp(text);
text = sprintf('lsqcurvefit settings: Jacobi: %s, TolFun: %e, TolX: %e, MaxIter: %e, MaxFunEval: %e', ...
               lsq_jac,lsq_ops(1),lsq_ops(2),lsq_ops(3),lsq_ops(4));
disp(text);
disp(g_param_temp)

% Compute solutions of different objects and store them in a struct      
for k=1:g_n_curve
    g_param_temp = [erg(1), erg(2), erg(3), erg(4)];    
    t_end = struct_data(k).t_data(length(struct_data(k).t_data));
    g_actual_curve = k;
    x0 = zeros(24,1);
    x0(1) = erg(5); x0(21) = 1;
    if (g_model == 0 || g_model == 1)
        struct_sol(k) = ode45(@model_onco_mono, [0,t_end+1], ... 
                        x0,options);
    else
        struct_sol(k) = ode45(@model_onco_comb, [0,t_end+1], ... 
                        x0,options);
    end
    T = struct_sol(k).x ;
    Y = (struct_sol(k).y)';
    figure(1); plot(T,Y(:,1)+Y(:,2)+Y(:,3)+Y(:,4),'black'); hold on;
end


% Calculate statistics
stat = fun_stats(erg,struct_sol,struct_data);


% Plot fitting results
void = fun_print_result(erg,g_model,stat);