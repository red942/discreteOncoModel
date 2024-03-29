% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
%
% Difference equations for onco model combination therapy
% (k2a and k2b are fixed) and the variational equations
%
%NOT YET IMPLEMENTD
function dxx = model_onco_comb_discrete(t,x)

global g_param_temp;
global gv_k2a_k2b_fix; 
global gcell_dosing;
global g_actual_curve;

l0  = g_param_temp(1);
l1  = g_param_temp(2);
k1  = g_param_temp(3);
k2a = gv_k2a_k2b_fix(1);     
k2b = gv_k2a_k2b_fix(2);     
psi = g_param_temp(4);

dxx = zeros(24,1);

if (g_actual_curve == 1)    % Control fit ==> c1,c2 = 0
    c1 = 0;
    c2 = 0;
else                        % Combination therapy 
    % First drug
    g_actual_curve = 2;
    switch gcell_dosing{1,1,2}
        case 0
            c1 = 0;
        case 1
            c1 = fun_c_PKpo_1Komp(t);
        case 2
            c1 = fun_c_PKpo_2Komp(t);
    end
    % Second drug
    g_actual_curve = 3;
    switch gcell_dosing{1,1,3}
        case 0
            c2 = 0;
        case 1    
            c2 = fun_c_PKpo_1Komp(t);
        case 2
            c2 = fun_c_PKpo_2Komp(t);
    end
end


% Onco model
%
w = x(1)+x(2)+x(3)+x(4);

dx(1) = (2*l0*l1*x(1)^2)/((l1+2*l0*x(1))*w) - (k2a*c1+psi*k2b*c2)*x(1);
dx(2) = (k2a*c1+psi*k2b*c2)*x(1) - k1*x(2);
dx(3) = k1*(x(2)-x(3));
dx(4) = k1*(x(3)-x(4));


% Variational equations 
%
x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);


df1_dx1 = 4*l0*l1*x1 / ((l1+2*l0*x1)*w) - 4*l0^2*l1*x1^2 / ((l1+2*l0*x1)^2*w) ... 
          -2*l0*l1*x1^2 / ((l1+2*l0*x1)*w^2) - (k2a*c1+psi*k2b*c2);

df1_dx2 = (-2*l0*l1*x1^2) / ((l1+2*l0*x1)*w^2);
df1_dx3 = df1_dx2;
df1_dx4 = df1_dx2;

df_dx = [df1_dx1,             df1_dx2, df1_dx3, df1_dx4;
         (k2a*c1+psi*k2b*c2), -k1    , 0      , 0;
         0                  ,  k1    , -k1    , 0;
         0                  , 0      , k1     , -k1;];

% Create df_deta matrix     
df1_dl0 = 2*l1*x1^2 / ((l1+2*l0*x1)*w) - 4*l0*l1*x1^3 / ((l1+2*l0*x1)^2*w);
df1_dl1 = 2*l0*x1^2 / ((l1+2*l0*x1)*w) - 2*l0*l1*x1^2 / ((l1+2*l0*x1)^2*w);

df_deta = [df1_dl0, df1_dl1, 0       , -k2b*c2*x1  ,    0;
           0      , 0      , -x2     , k2b*c2*x1   ,    0;
           0      , 0      , x2 - x3 , 0           ,    0;
           0      , 0      , x3 - x4 , 0           ,    0];

v_x = [ x(5), x(9),  x(13), x(17), x(21);
        x(6), x(10), x(14), x(18), x(22);
        x(7), x(11), x(15), x(19), x(23);
        x(8), x(12), x(16), x(20), x(24); ];
    
% Calculate variational equation    
du = df_dx*v_x + df_deta;


% Output is a column vector 
%
dxx = [dx(1);dx(2);dx(3);dx(4);du(:)];