% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
% 
% Control PC3
%
% Fitting results: 
% l0 = 0.10364973, l1 = 0.27454333, w0 = 0.10973956 
% (CV = 27.7508),  (CV = 43.5648),  (CV = 28.7388)
% SofS = 0.00956218, MSE = 0.00191244, 
% Var  = 0.00191244, R^2 = 0.99379644
%
function struct_data = fun_data_control_PC3()

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
gs_c_info = 'Control_PC3';        % String with drug related information
gs_LargeScale = 'off';            % 'on' = LargeScale, 'off' = MediumScale 
gs_excel     = 'on';              % Write results as Excel File
                                  % (works only on Windows with Excel)
                                  
                                  
% Initial values                    Control Data:  (l0,l1,0, 0,  w0)
%                                   Monotherapy:    (l0,l1,k1,k2, w0)
%                                   Comb. therapy: (l0,l1,k1,psi,w0)    
gv_start_vals = [0.2, 1, 0, 0, 0.015];


% Fitting Data
%
% Obj_1 : Control Data PC3 (Obj_1 is always control data)
t_data_obj_1 = [4; 6; 8; 10; 12; 14; 16; 18;];
y_data_obj_1 = [0.2240800; 0.3091580; 0.4915150; 0.6049460; 0.7189090; ...
               1.0582420; 1.2686990; 1.5395230; ];
           
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


