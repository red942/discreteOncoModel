% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
% 
% Control RKO
%
% Fitting results: 
% l0  = 0.24976992, l1  = 0.22304630, w0  = 0.00143812
% (CV = 32.3767)  , (CV = 18.3026)  , (CV = 167.7960) 
% SofS = 0.01403851, MSE = 0.00233975, 
% Var  = 0.00233975, R^2 = 0.99587053
%
function struct_data = fun_data_control_RKO()

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
gs_c_info = 'Control_RKO';        % String with drug related information
gs_LargeScale = 'off';            % 'on' = LargeScale, 'off' = MediumScale 
gs_excel = 'on';                  % Write results as Excel File
                                  % (works only on Windows with Excel)
                                  
                                  
% Initial values                    Control Data:  (l0,l1,0, 0,  w0)
%                                   Monotherapy:   (l0,l1,k1,k2, w0)
%                                   Comb. therapy: (l0,l1,k1,psi,w0)
gv_start_vals = [0.2, 1, 0, 0, 0.015];


% Fitting Data
%
% Obj_1 : Control Data RKO (Obj_1 is always control data)
t_data_obj_1 = [4; 6; 8; 11; 13; 15; 18; 20; 22];
y_data_obj_1 = [0.0427544; 0.0624530; 0.0800580; 0.1666590; 0.4211890; ...
                0.6233760; 1.0344770; 1.4683570; 1.6991660; ];
            
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


