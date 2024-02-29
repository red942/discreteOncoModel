% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
% 
% Task: Returns the function values and the values of the jacobian
%       at the points in time given in 't_data'
%
function [F,J] = fun_F_and_Jac_mono_or_comb(param,t_data)

global g_param_temp g_y_data;
global g_fid;
global gs_wr_2_disc gs_wr_2_disp;
global g_anz_data_all;
global gv_anz_data;
global g_actual_curve
global g_n_curve;
global g_model;
global gv_SofS_CT;

g_param_temp = param;     

options = odeset('RelTol',1e-8,'AbsTol',1e-8);

Y = zeros(24,1);   
F = zeros(1,g_anz_data_all);
J = zeros(g_anz_data_all,5);

act_pos = 0;

% Loop over number of objects
for k=1:g_n_curve
    
    g_actual_curve = k;                    % Actual object
    
    % act_pos+1 is the index of the first point in time of object k 
    % in column vector t_data  
    if (k>1)
        act_pos = act_pos + gv_anz_data(k-1);
    end
    
    % Initial values 
    x0 = zeros(1,24);
    x0(1) = param(5); x0(21) = 1;

    % Calculate solution of monotherapy or combination therapy in [0,t(k,nk)]
    if (g_model == 0 || g_model == 1)   % Monotherapy
        sol = ode45(@model_onco_mono,[0, t_data(gv_anz_data(k) + act_pos)],x0,options);
    else                                % Combination therapy
        sol = ode45(@model_onco_comb,[0, t_data(gv_anz_data(k) + act_pos)],x0,options);
    end

    % Evaluate solution at the specified points in time
    for i=1:gv_anz_data(k)
    
        Y(:,1) = deval(sol,t_data(i+act_pos));
            
        F(i+act_pos)   = Y(1,1)  + Y(2,1)  + Y(3,1)  + Y(4,1);
        J(i+act_pos,1) = Y(5,1)  + Y(6,1)  + Y(7,1)  + Y(8,1);
        J(i+act_pos,2) = Y(9,1)  + Y(10,1) + Y(11,1) + Y(12,1);
        J(i+act_pos,3) = Y(13,1) + Y(14,1) + Y(15,1) + Y(16,1);
        J(i+act_pos,4) = Y(17,1) + Y(18,1) + Y(19,1) + Y(20,1);
        J(i+act_pos,5) = Y(21,1) + Y(22,1) + Y(23,1) + Y(24,1);
    end

end

% Calculate sum of squares
h = 0;
for i=1:g_anz_data_all
    h = h + (F(i) - g_y_data(i))^2;
end

% Calculate threshold concentration for monotherapy
if (g_model == 1)
    C_T = (-param(3) + param(3)*sqrt(1+24*(param(1)/param(3))))/(6*param(4));
else 
    C_T = 0;
end

% Output 
switch g_model
    case 0
        if strcmp(gs_wr_2_disp,'on')
            text = sprintf('SofS = %18.12f | %12.6f %12.6f %12.6f', ...
                           h,param(1),param(2),param(5));               
            disp(text);
        end
        if strcmp(gs_wr_2_disc,'on')
            fprintf(g_fid, '%18.12f | %12.6f %12.6f %12.6f \n\r\n\r', ...
                    h,param(1),param(2),param(5));
        end
        gv_SofS_CT = [h; 0];
    case 1
        if strcmp(gs_wr_2_disp,'on')
            text = sprintf('SofS = %18.12f | %12.6f %12.6f %12.6f %12.6f %12.6f, C_T = %12.6f', ...
                           h,param(1),param(2),param(3),param(4),param(5),C_T);               
            disp(text);
        end
        if strcmp(gs_wr_2_disc,'on')
            fprintf(g_fid, '%18.12f | %12.6f %12.6f %12.6f %12.6f %12.6f, C_T = %12.6f \n\r\n\r', ...
                    h,param(1),param(2),param(3),param(4),param(5),C_T);
        end
        gv_SofS_CT = [h;C_T];
    case 2
        if strcmp(gs_wr_2_disp,'on')
            text = sprintf('SofS = %18.12f | %12.6f %12.6f %12.6f %12.6f %12.6f', ...
                           h,param(1),param(2),param(3),param(4),param(5));               
            disp(text);
        end
        if strcmp(gs_wr_2_disc,'on')
            fprintf(g_fid, '%18.12f | %12.6f %12.6f %12.6f %12.6f %12.6f\n\r\n\r', ...
                    h,param(1),param(2),param(3),param(4),param(5));
        end
        gv_SofS_CT = [h;0];
end