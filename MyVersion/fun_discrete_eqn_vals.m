% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
% 
% Task: Returns the function values and the values of the jacobian
%       at the points in time given in 't_data'
%
function [F] = fun_discrete_eqn_vals(param, t_data)

global g_param_temp g_y_data;
global gs_wr_2_disc;
global g_anz_data_all;
global gv_anz_data;
global g_actual_curve
global g_n_curve;
global g_model;
global gv_SofS_CT;
global timescale;

g_param_temp = param;

F = zeros(1, g_anz_data_all);

act_pos = 0;

for k = 1:g_n_curve
    g_actual_curve = k;

    if (k > 1)
        act_pos = act_pos + gv_anz_data(k - 1);
    end

    %x0 = zeros(1, 24);
    %x0(1) = param(5); x0(21) = 1;
    x0 = param(5);

    % Call the discrete solver
    %assumes x0 = w0 from the paper
    soltn_vals = model_onco_mono_discrete_solver(t_data(1), t_data(end), timescale, x0, param);

    % Evaluate function values
    for i = 1:gv_anz_data(k)
        F(i + act_pos) = sum(soltn_vals(i, :));  % Sum across columns for actual values
        % No need to compute Jacobian separately since it is not used
    end
end

% Calculate sum of squares
h = 0;
for i = 1:g_anz_data_all
    h = h + (F(i) - g_y_data(i))^2;
end

% Calculate threshold concentration for monotherapy
if (g_model == 1)
    C_T = (-param(3) + param(3) * sqrt(1 + 24 * (param(1) / param(3)))) / (6 * param(4));
else
    C_T = 0;
end

% Output 
switch g_model
    case 0
        %if strcmp(gs_wr_2_disp,'on')
        %    text = sprintf('SofS = %18.12f | %12.6f %12.6f %12.6f', ...
        %                   h,param(1),param(2),param(5));               
        %    disp(text);
        %end
        if strcmp(gs_wr_2_disc,'on')
            fprintf(g_fid, '%18.12f | %12.6f %12.6f %12.6f \n\r\n\r', ...
                    h,param(1),param(2),param(5));
        end
        gv_SofS_CT = [h; 0];
    case 1
        %if strcmp(gs_wr_2_disp,'on')
        %    text = sprintf('SofS = %18.12f | %12.6f %12.6f %12.6f %12.6f %12.6f, C_T = %12.6f', ...
        %                   h,param(1),param(2),param(3),param(4),param(5),C_T);               
        %    disp(text);
        %end
        if strcmp(gs_wr_2_disc,'on')
            fprintf(g_fid, '%18.12f | %12.6f %12.6f %12.6f %12.6f %12.6f, C_T = %12.6f \n\r\n\r', ...
                    h,param(1),param(2),param(3),param(4),param(5),C_T);
        end
        gv_SofS_CT = [h;C_T];
    case 2
        %if strcmp(gs_wr_2_disp,'on')
        %    text = sprintf('SofS = %18.12f | %12.6f %12.6f %12.6f %12.6f %12.6f', ...
        %                   h,param(1),param(2),param(3),param(4),param(5));               
        %    disp(text);
        %end
        if strcmp(gs_wr_2_disc,'on')
            fprintf(g_fid, '%18.12f | %12.6f %12.6f %12.6f %12.6f %12.6f\n\r\n\r', ...
                    h,param(1),param(2),param(3),param(4),param(5));
        end
        gv_SofS_CT = [h;0];
end