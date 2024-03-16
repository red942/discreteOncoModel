% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
% 
% Control MDA-MB468
%
% Fitting results: 
% l0 = 0.03360087, l1 = 0.10423437, w0 = 0.05207537
% (CV = 29.8185) , (CV = 152.8675), (CV = 26.1595)
% SofS = 0.00604198, MSE = 0.00075525, 
% Var  = 0.00075525, R^2 = 0.97605735
%
function struct_data = fun_data_control_MDA()

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
gs_c_info = 'Control_MDA';        % String with drug related information
gs_LargeScale = 'off';            % 'on' = LargeScale, 'off' = MediumScale 
gs_excel     = 'on';              % Write results as Excel File
                                  % (works only on Windows with Excel)
                                  
                                  
% Initial values                    Control Data:  (l0,l1,0, 0,  w0)
%                                   Monotherapy:   (l0,l1,k1,k2, w0)
%                                   Comb. therapy: (l0,l1,k1,psi,w0)
gv_start_vals = [0.2, 1, 0, 0, 0.015];


% Fitting Data
%
% Obj_1 : Control Data MDA-MB468 (Obj_1 is always control data)
t_data_obj_1 = [4; 7; 11; 14; 18; 21; 25; 28; 32; 34; 40];
y_data_obj_1 = [0.0327050; 0.0670580; 0.0864680; 0.1503430; 0.2031840; ...
                0.2238060; 0.2361360; 0.2892050; 0.3351280; 0.4018410; 
                0.5667850; ];           
           
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


