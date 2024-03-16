% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
% 
% Control 2
%
% Fitting results: 
% l0 = 0.18713875, l1 = 0.24262374, w0 = 0.01242141 
% (CV =  54.3395)  (CV =  31.2105)  (CV = 182.9584)
% SofS = 0.02674929, MSE =   0.00668732 
% Var  = 0.00668732, R^2 =   0.99171780
%
function struct_data = fun_data_control_2()

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
g_model = 0;                      % 0 = Control Data
                                  % 1 = Monotherapy
                                  % 2 = Combination therapy
g_n_curve = 1;                    % Number of objects
gs_wr_2_disc = 'off';             % Write fitting progress to disc
gs_wr_2_disp = 'on';              % Write fitting progress to display
gs_c_info = 'Control_2';          % String with drug related information
gs_LargeScale = 'off';            % 'on' = LargeScale, 'off' = MediumScale
gs_excel     = 'on';              % Write results as Excel file
                                  % (works only on Windows with Excel)
                                  
                                  
% Initial values                    Control Data:  (l0,l1,0, 0,  w0
%                                   Monotherapy:   (l0,l1,k1,k2, w0)
%                                   Comb. therapy: (l0,l1,k1,psi,w0)
gv_start_vals = [0.1, 0.3, 0, 0, 0.01];


% Fitting data
%
% Obj_1 : Control Data 2 (Obj_1 is always control data)
t_data_obj_1 = [9; 12; 14; 16; 19; 21; 23];
y_data_obj_1 = [0.2362101711; 0.4684588444; 0.8166467914; 1.070723326; ...
                1.488022049;  1.750415778;  2.306176399; ];
            
% Write fitting data in a struct   
struct_data(1) = struct('anz_data',length(t_data_obj_1),'t_data',t_data_obj_1,'y_data',y_data_obj_1);


% PK Info
%
% pk_obj_x contains the PK parameter and comp_obj_x defines the Comp.-model
% 0 = No dosing
% 1 = 1-Comp: dose,ka,k,V
% 2 = 2-Comp: Aoral,Boral,alpha,beta,ka
pk_obj_1 = 0;
sched_obj_1 = 0;
comp_obj_1 = 0;

% Write PK data in a cell
dosing(1) = struct('Comp',comp_obj_1,'pk_obj',pk_obj_1,'sched_obj',sched_obj_1); 
gcell_dosing = struct2cell(dosing);

