% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
% 
% Monotherapy: Control 1 + Drug C 100 + Drug C 150
%
% Fitting results: 
% l0 = 0.13637831, l1 = 0.28561818, k1 = 13.8367239, k2 = 0.01468196, w0 = 0.04569168 
%( CV =  22.8305)  (CV =  15.0475)  (CV = 984.6345)  (CV =  29.5018)  (CV = 52.6519)
% SofS = 0.22446708, MSE = 0.00897869 
%
% 1 curve: Var = 0.01441829, R^2 = 0.99139104
% 2 curve: Var = 0.00682655, R^2 = 0.98983040
% 3 curve: Var = 0.00569121, R^2 = 0.98737822  
%
function struct_data = fun_data_mono_drug_C_100_150()

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
g_n_curve = 3;                    % Number of objects                                  
gs_wr_2_disc = 'off';             % Write fitting progress to disc
gs_wr_2_disp = 'on';              % Write fitting progress to display
gs_c_info = 'C_100_150';          % String with drug related information
gs_LargeScale = 'on';             % 'on' = LargeScale, 'off' = MediumScale
gs_excel = 'on';                  % Write results as Excel File
                                  % (works only on Windows with Excel)
                                  
                                  
% Initial values                    Control Data:  (l0,l1,0, 0,  w0)
%                                   Monotherapy:   (l0,l1,k1,k2, w0)
%                                   Comb. therapy: (l0,l1,k1,psi,w0)                                  
gv_start_vals = [0.1356, 0.2744, 1, 0.01, 0.0528];


% Fitting Data
%
% Obj_1 : Control Data (Obj_1 is always control data)
% Obj_2 : Drug C 100mg/kg po, SID1-12, Start day 11
% Obj_3 : Drug C 150mg/kg po, SID1-12, Start day 11
t_data_obj_1 = [6; 11; 13; 15; 18; 20; 22; 25; 27; 29; ];
y_data_obj_1 = [0.289324513; 0.5201940347; 0.8787273154; 1.034623291; ...
                1.636494591; 1.883689601;  2.273302339;  3.036456225; ...
                3.055183057; 3.789235538; ];
t_data_obj_2 = [6; 11; 13; 15; 18; 20; 22; 25; 27; 29; ];
y_data_obj_2 = [0.2874112509; 0.5317808964; 0.6744040076; 0.7878502488; ...
                0.9540505736; 1.290836649;  1.297906082;  1.79016148; ...
                2.205434541;  2.408760284; ];
t_data_obj_3 = [6; 11; 13; 15; 18; 20; 22; 25; 27; 29; ];
y_data_obj_3 = [0.3133178994; 0.5234641689; 0.5870236276; 0.6787753731; ...
                0.819572133;  0.9483855993; 1.068650395;  1.489703; ...
                1.857222075;  2.162818767; ];

% Write fitting data in a struct                    
struct_data(1) = struct('anz_data',length(t_data_obj_1),'t_data',t_data_obj_1,'y_data',y_data_obj_1);
struct_data(2) = struct('anz_data',length(t_data_obj_2),'t_data',t_data_obj_2,'y_data',y_data_obj_2);
struct_data(3) = struct('anz_data',length(t_data_obj_3),'t_data',t_data_obj_3,'y_data',y_data_obj_3);


% PK Info
%
% pk_obj_x contains the PK parameter and comp_obj_x defines the Comp.-model
% 0 = No dosing
% 1 = 1-Comp: dose,ka,k,V
% 2 = 2-Comp: Aoral,Boral,alpha,beta,ka           
pk_obj_1 = 0;
pk_obj_2 = [16.8619,     26.9165,     0.17046905, 4.96189217, 84.90413544];
pk_obj_3 = [16.8619*1.5, 26.9165*1.5, 0.17046905, 4.96189217, 84.90413544];
comp_obj_1 = 0;
comp_obj_2 = 2;
comp_obj_3 = 2;

% Dosing Schedule, [Number of doses, point in time (in hr), point in time (in hr), ... ]
sched_obj_1 = 0;
sched_obj_2 = [12, 264, 288, 312, 336, 360, 384, 408, 432, 456, 480, 504, 528];
sched_obj_3 = [12, 264, 288, 312, 336, 360, 384, 408, 432, 456, 480, 504, 528];

% Write PK data in a cell
dosing(1) = struct('Comp',comp_obj_1,'pk_obj',pk_obj_1,'sched_obj',sched_obj_1); 
dosing(2) = struct('Comp',comp_obj_2,'pk_obj',pk_obj_2,'sched_obj',sched_obj_2);
dosing(3) = struct('Comp',comp_obj_3,'pk_obj',pk_obj_3,'sched_obj',sched_obj_3);
gcell_dosing = struct2cell(dosing);
