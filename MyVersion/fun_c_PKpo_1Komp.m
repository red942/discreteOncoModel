% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
% 
% Task: Multiple dosing with a one compartment model p.o.
% 
% Based on the code of J.Fischer 
% (see http://www.boomer.org/pkin/PK02/PK2002218.html)
%
function c = fun_c_PKpo_1Komp(t)

g_actual_curve = 1;
global gcell_dosing;

% Get parameters and dosing intervall
dose = 100;
ka   = 103.9673;
k    = 0.1052;
V    = 2.788;
con  = 100;


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
 
    
   