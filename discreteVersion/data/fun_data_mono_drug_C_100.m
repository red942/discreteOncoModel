% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
% 
% Monotherapy: Control 1 + Drug C 100
%
% Fitting results: 
% l0  = 0.13601299, l1 = 0.28019425, k1 = 13.746919, k2 = 0.01540936, w0 = 0.04924977 
% (CV = 28.4785)    (CV = 17.8520)   (CV = 1662)     (CV = 46.1078)   (CV = 65.6133)
% SofS = 0.16590534, MSE = 0.01106035 
%
% 1 curve: Var = 0.01494051, R^2 = 0.99140434
% 2 curve: Var = 0.00718020, R^2 = 0.98853457
%
% C_T = 16.713884
%
function struct_data = fun_data_mono_drug_C_100()

global gv_start_vals;
global gs_c_info;
global gs_wr_2_disc gs_wr_2_disp;
global g_n_curve;
global g_model;
global gs_LargeScale;
global gcell_dosing;
global gs_excel;


% Initial settings 
%
g_model = 1;                      % 0 = Control Data
                                  % 1 = Monotherapy
                                  % 2 = Combination therapy
g_n_curve = 2;                    % Number of objects
gs_wr_2_disc = 'off';             % Write fitting progress to disc
gs_wr_2_disp = 'on';              % Write fitting progress to display
gs_c_info = 'C_100';              % String with drug related information
gs_LargeScale = 'on';             % 'on' = LargeScale, 'off' = MediumScale
gs_excel = 'on';                  % Write results as Excel file
                                  % (works only on Windows with Excel)

                                  
% Initial values                    Control Data:  (l0,l1,0, 0,  w0)
%                                   Monotherapy:   (l0,l1,k1,k2, w0)
%                                   Comb. therapy: (l0,l1,k1,psi,w0)
gv_start_vals = [0.1356, 0.2744, 1, 0.01, 0.0528];       


% Fitting data
%
% Obj_1 : Control Data 1 (Obj_1 is always control data)
% Obj_2 : Drug C 100mg/kg po,  SID1-12, Start day 11
t_data_obj_1 = [6; 11; 13; 15; 18; 20; 22; 25; 27; 29; ];
y_data_obj_1 = [0.289324513; 0.5201940347; 0.8787273154; 1.034623291; ...
                1.636494591; 1.883689601;  2.273302339;  3.036456225; ...
                3.055183057; 3.789235538; ];
t_data_obj_2 = [6; 11; 13; 15; 18; 20; 22; 25; 27; 29; ];
y_data_obj_2 = [0.2874112509; 0.5317808964; 0.6744040076; 0.7878502488; ...
                0.9540505736; 1.290836649;  1.297906082;  1.79016148; ...
                2.205434541;  2.408760284; ];

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
pk_obj_2 = [16.8619, 26.9165, 0.17046905, 4.96189217, 84.90413544; ];
comp_obj_1 = 0;
comp_obj_2 = 2;

% Dosing schedule, [Number of doses, point in time (in hr), point in time (in hr), ... ]
sched_obj_1 = 0;
sched_obj_2 = [12, 264, 288, 312, 336, 360, 384, 408, 432, 456, 480, 504, 528];

% Write PK data in a cell
dosing(1) = struct('Comp',comp_obj_1,'pk_obj',pk_obj_1,'sched_obj',sched_obj_1); 
dosing(2) = struct('Comp',comp_obj_2,'pk_obj',pk_obj_2,'sched_obj',sched_obj_2);
gcell_dosing = struct2cell(dosing);
    