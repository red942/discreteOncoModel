clear all;  
clear global all;
addpath('data','funs'); %lets matlab see these folders

global g_actual_curve
global gcell_dosing

struct_data = fun_data_mono_drug_B_100();

g_actual_curve = 1;

model_onco_mono_discrete_solver([0, 10], 1, 0.0567, [0.127, 0.293, 7.22, 0.18])
