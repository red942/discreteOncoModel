% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
%
% Differential equations for onco model monotherapy and the variational equations
%
function dxx = model_onco_mono(t,x)

global g_param_temp;
global g_actual_curve;
global gcell_dosing;

l0 = g_param_temp(1);
l1 = g_param_temp(2);
k1 = g_param_temp(3);
k2 = g_param_temp(4);
dxx = zeros(24,1);

c = 0;
switch gcell_dosing{1,1,g_actual_curve}
    case 0
        c = 0;
    case 1 
        c = fun_c_PKpo_1Komp(t);
    case 2
        c = fun_c_PKpo_2Komp(t);
end


% Onco model
%
w = x(1)+x(2)+x(3)+x(4);

dx(1) = (2*l0*l1*x(1)^2)/((l1+2*l0*x(1))*w) - k2*c*x(1);
dx(2) = k2*c*x(1)-k1*x(2);
dx(3) = k1*(x(2)-x(3));
dx(4) = k1*(x(3)-x(4));


% Variational equations
%
x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);

% Create df_dx matrix
df1_dx1 = 4*l0*l1*x1 / ( (l1+2*l0*x1) * w)  - 4*l0^2*l1*x1^2 / ( (l1+2*l0*x1)^2 * w ) ... 
          - 2*l0*l1*x1^2 / ( (l1+2*l0*x1) * w^2)  - k2*c;
df1_dx2 = -2*l0*l1*x1^2 / ( (l1+2*l0*x1) * w^2) ;

df1_dx3 = df1_dx2;
df1_dx4 = df1_dx2;

df_dx = [df1_dx1, df1_dx2, df1_dx3, df1_dx4;
         k2*c   , -k1    , 0      , 0;
         0      ,  k1    , -k1    , 0;
         0      , 0      , k1     , -k1;];

% Create df_deta matrix     
df1_dl0 = 2*l1*x1^2 / ( (l1+2*l0*x1) * w)  - 4*l0*l1*x1^3 / ( (l1+2*l0*x1)^2 * w );
df1_dl1 = 2*l0*x1^2 / ( (l1+2*l0*x1) * w) - 2*l0*l1*x1^2 / ( (l1+2*l0*x1)^2 * w);

df_deta = [df1_dl0, df1_dl1, 0       , -c*x1    ,    0;
           0      , 0      , -x2     , c*x1     ,    0;
           0      , 0      , x2 - x3 , 0        ,    0;
           0      , 0      , x3 - x4 , 0        ,    0];
     
v_x = [ x(5), x(9),  x(13), x(17), x(21);
        x(6), x(10), x(14), x(18), x(22);
        x(7), x(11), x(15), x(19), x(23);
        x(8), x(12), x(16), x(20), x(24); ];
    
% Calculate variational equation    
du = df_dx*v_x + df_deta;


% Output is a column vector 
%
dxx = [dx(1);dx(2);dx(3);dx(4);du(:)];

