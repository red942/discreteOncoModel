% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
% 
% Task: Multiple dosing with a one compartment model p.o.
% 
% Based on the code of J.Fischer 
% (see http://www.boomer.org/pkin/PK02/PK2002218.html)
%
function c = fun_c_PKpo_1Komp(t)

global g_actual_curve;
global gcell_dosing;

% Get parameters and dosing intervall
dose = gcell_dosing{2,1,g_actual_curve}(1);
ka   = gcell_dosing{2,1,g_actual_curve}(2);
k    = gcell_dosing{2,1,g_actual_curve}(3);
V    = gcell_dosing{2,1,g_actual_curve}(4);
con  = gcell_dosing{3,1,g_actual_curve};


t = t*24;              % Change time unit in hours
j = 1;
ndose = con(1);        % Number of doses

% Count number of given doses until the present time point t
for i=1:ndose           
    j = j + 1;          
    if (t <= con(j))
        i = i - 1;      % Correct index
        break;          % Goto (*) 
    end
end
% (*)

ndose = i;
c = 0;
h = 1;
    
% Superposition Principle: Summing up the concentration
for m=1:ndose
    h = h + 1;          
    td = t - con(h);
    coef = ((ka*dose)/((ka-k)*V))*(exp(-k*td) - exp(-ka*td));
    c    = c + coef;
end
 
    
   