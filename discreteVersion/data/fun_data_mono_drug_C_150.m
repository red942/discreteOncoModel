% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
%
% Monotherapy: Control 2 + Drug C 150
% 
% Fitting results:
% l0 = 0.19624830, l1 = 0.23757468, k1 = 10.5855802, k2 = 0.01155475, w0 = 0.01046438
% (CV = 36.6397)  (CV = 21.8378)   (CV = 499.1216)  (CV = 24.4490)   (CV = 123.5606)
% SofS = 0.03116588, MSE = 0.00346288
%
% 1 curve: Var = 0.00595832, R^2 = 0.99169530
% 2 curve: Var = 0.00096743, R^2 = 0.99707680  
%
% C_T = 30.8515
%
function struct_data = fun_data_mono_drug_C_150()

global gv_start_vals;
global gs_c_info;
global gs_wr_2_disc gs_wr_2_disp;
global g_n_curve;
global g_model;
global gs_LargeScale;
global gcell_dosing;
global gs_excel;


% Initial Settings 
%
g_model = 1;                      % 0 = Control Data
                                  % 1 = Monotherapy
                                  % 2 = Combination therapy
g_n_curve = 2;                    % Number of objects
gs_wr_2_disc = 'off';             % Write fitting progress to disc
gs_wr_2_disp = 'on';              % Write fitting progress to display
gs_c_info = 'C_150';              % String with drug related information
gs_LargeScale = 'on';             % 'on' = LargeScale, 'off' = MediumScale
gs_excel     = 'on';              % Write results as Excel File
                                  % (works only on Windows with Excel)
                                  
% Initial values                    Control Data:  (l0,l1,0, 0,  w0)
%                                   Monotherapy:   (l0,l1,k1,k2, w0)
%                                   Comb. therapy: (l0,l1,k1,psi,w0)
gv_start_vals = [0.1871, 0.2426, 1, 0.01, 0.0124];
                                  
         
% Fitting Data
%
% Obj_1 : Control Data 2 (Obj_1 is always control data)
% Obj_2 : Drug C 150mg/kg po,  SID1-7, Start day 12
t_data_obj_1 = [9; 12; 14; 16; 19; 21; 23];
y_data_obj_1 = [0.2362101711; 0.4684588444; 0.8166467914; 1.070723326; ...
                1.488022049;  1.750415778;  2.306176399; ];
t_data_obj_2 = [9; 12; 14; 16; 19; 21; 23];
y_data_obj_2 = [0.2448071886; 0.4912571692; 0.6548858702; 0.8480012311; ...
                0.9808028857; 1.337384611; 1.682867548; ];

% Write fitting data in a struct                            
struct_data(1) = struct('anz_data',length(t_data_obj_1),'t_data',t_data_obj_1,'y_data',y_data_obj_1);
struct_data(2) = struct('anz_data',length(t_data_obj_2),'t_data',t_data_obj_2,'y_data',y_data_obj_2);
 

% PK Info
%
% pk_obj_x contains the PK parameter and comp_obj_x defines the Comp.-model
% 0 = No dosing
% 1 = 1-Comp: dose,ka,k,V
% 2 = 2-Comp: Aoral,Boral,alpha,beta,ka           
pk_obj_1 = 0;
pk_obj_2 = [25.2929, 40.3747, 0.17046905, 4.96189217, 84.90413544 ];
            
% Dosing Schedule, [Number of doses, point in time (in hr), point in time (in hr), ... ]
sched_obj_1 = 0;
sched_obj_2 = [7, 288, 312, 336, 360, 384, 408, 432];
comp_obj_1 = 0;
comp_obj_2 = 2;

% Write PK data in a cell
dosing(1) = struct('Comp',comp_obj_1,'pk_obj',pk_obj_1,'sched_obj',sched_obj_1); 
dosing(2) = struct('Comp',comp_obj_2,'pk_obj',pk_obj_2,'sched_obj',sched_obj_2);
gcell_dosing = struct2cell(dosing);
