% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
% 
% Control 1
%
% Fitting results: 
% l0 = 0.13559699, l1 = 0.27442345, w0 = 0.05283816  
% (CV = 40.2811)   (CV = 22.4084)   (CV = 95.6730)
% SofS = 0.10947629, MSE = 0.01563947
% Var  = 0.01563947, R^2 = 0.99135962
%
function struct_data = fun_data_control_1()

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
gs_wr_2_disc = 'off';             % Write fitting process to disc
gs_wr_2_disp = 'on';              % Write fitting process to display
gs_c_info = 'Control_1';          % String with drug related information
gs_LargeScale = 'off';            % 'on' = LargeScale, 'off' = MediumScale 
gs_excel = 'on';                  % Write results as Excel file
                                  % (works only on Windows with Excel)

                                  
% Initial values                    Control Data:  (l0,l1,0, 0,  w0)
%                                   Monotherapy:   (l0,l1,k1,k2, w0)
%                                   Comb. therapy: (l0,l1,k1,psi,w0)
gv_start_vals = [0.1, 0.3, 0, 0, 0.01];    


% Fitting data
%
% Obj_1 : Control Data 1 (Obj_1 is always control data)
t_data_obj_1 = [6; 11; 13; 15; 18; 20; 22; 25; 27; 29; ];
y_data_obj_1 = [0.289324513; 0.5201940347; 0.8787273154; 1.034623291; ...
                1.636494591; 1.883689601;  2.273302339;  3.036456225; ...
                3.055183057; 3.789235538; ];

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

