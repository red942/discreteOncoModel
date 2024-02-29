% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
% 
% Task: Return t-values
% 
% Based on
% [DRAPER NR, SMITH H, Applied Regression Analysis, Third Edition, 1998, p.686
%                      Wiley Series in Probability and Statistics]
%
function t_val = fun_calc_t_val_95per(v)

a_t_values = [12.71, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, ...
              2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086, ...
              2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042, ...
              2.021, 2.009, 2.000, 1.990, 1.984, 1.980, 1.960];
        
if (v>= 1 && v< 30)  
    t_val = a_t_values(v);
end

if (v >= 30 && v < 35)
    t_val = a_t_values(30);
end

if (v >= 35 && v < 45)
    t_val = a_t_values(31);
end

if (v >= 45 && v < 55)
    t_val = a_t_values(32);
end

if (v >= 55 && v < 70)
    t_val = a_t_values(33);
end

if (v >= 70 && v < 90)
    t_val = a_t_values(34);
end

if (v >= 90 && v < 110)
    t_val = a_t_values(35);
end

if (v >= 110 && v < 130)
    t_val = a_t_values(36);
end

if (v >= 130)
    t_val = a_t_values(37);
end

    