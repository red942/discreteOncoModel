% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
% 
% Combination therapy: Control 2 + Drug B 100 + Drug C 150
%
% Fitting results
% l0 = 0.17205175, l1 = 0.26048346, k1 = 0.49127820,  psi = 1.78319345, w0 = 0.01541515 
% (CV = 33.7399)   (CV = 24.3909)   (CV = 16.2409)    (CV = 19.6498)    (CV = 102.0861)
% SofS = 0.03096626, MSE = 0.00344070
%
% 1 curve: Var = 0.00610387, R^2 = 0.99156812
% 2 curve: Var = 0.00077752, R^2 = 0.97801849  
%
function struct_data = fun_data_comb_drug_B_100_and_C_150()

global gv_start_vals;
global gs_c_info;
global gs_wr_2_disc gs_wr_2_disp;
global g_n_curve;
global g_model
global gv_k2a_k2b_fix;
global gs_LargeScale;
global gcell_dosing;
global gs_excel;


% Initial Settings 
%
g_model = 2;                      % 0 = Control Data
                                  % 1 = Monotherapy
                                  % 2 = Combination therapy
g_n_curve = 2;                    % Number of objects
gs_wr_2_disc = 'off';             % Write fitting progress to disc
gs_wr_2_disp = 'on';              % Write fitting progress to display
gs_c_info = 'B_100_+_C_150';      % String with drug related information
gs_LargeScale = 'off';            % 'on' = LargeScale, 'off' = MediumScale
gs_excel = 'on';                  % Write results as Excel File
                                  % (works only on Windows with Excel)
                                  
                                  
% Initial values                    Control Data:  (l0,l1,0, 0,  w0)
%                                   Monotherapy:   (l0,l1,k1,k2, w0)
%                                   Comb. therapy: (l0,l1,k1,psi,w0)
gv_start_vals = [0.19, 0.24, 1, 1, 0.01];

gv_k2a_k2b_fix = [0.0077, 0.0116]; % Fixed potency parameters of Drug 1 and 2


% Fitting Data
%                                  
% Obj_1 : Control Data 2 (Obj_1 is always control data)
% Obj_2 : Drug B 100mg/kg po,  SID1-5,  Start day 12 
%        +Drug C 150mg/kg po,  SID1-7,  Start day 12
t_data_obj_1 = [9; 12; 14; 16; 19; 21; 23];
y_data_obj_1 = [0.2362101711; 0.4684588444; 0.8166467914; 1.070723326; ...
                1.488022049;  1.750415778;  2.306176399; ];
t_data_obj_2 = [9; 12; 14; 16; 19; 21; 23];
y_data_obj_2 = [0.2636992162; 0.4850432728; 0.6450617936; 0.6426685872; ...
                0.6351207896; 0.6160567835; 0.6280185988; ];
                              
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
pk_obj_2_drug_1 = [100, 103.967302, 0.10521723, 2.78825022];
pk_obj_2_drug_2 = [25.2929, 40.3747, 0.17046905, 4.96189217, 84.90413544];
comp_obj_1 = 0;
comp_obj_2_drug_1 = 1;
comp_obj_2_drug_2 = 2;             

% Dosing Schedule, [Number of doses, point in time (in hr), point in time (in hr), ... ]
sched_obj_1 = 0;
sched_obj_2_drug_1 = [5, 288, 312, 336, 360, 384,];
sched_obj_2_drug_2 = [7, 288, 312, 336, 360, 384, 408, 432];
      
% Write PK data in a cell
dosing(1) = struct('Comp',comp_obj_1,'pk_obj',pk_obj_1,'sched_obj',sched_obj_1); 
dosing(2) = struct('Comp',comp_obj_2_drug_1,'pk_obj',pk_obj_2_drug_1,'sched_obj',sched_obj_2_drug_1);
dosing(3) = struct('Comp',comp_obj_2_drug_2,'pk_obj',pk_obj_2_drug_2,'sched_obj',sched_obj_2_drug_2);
gcell_dosing = struct2cell(dosing);             
             









