% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
% 
% Monotherapy: Control 1 + Drug A2 120
%
% Fitting results: 
% l0 =  0.12067358, l1 = 0.31253305, k1 = 2.96803993, k2 = 0.00732078, w0 = 0.05871917 
% (CV = 30.5039)    (CV = 22.8123)   (CV = 60.6706)   (CV = 12.1140)   (CV = 67.5035)
% SofS = 0.26903953, MSE = 0.01793597
%
% 1 curve: Var = 0.01767965, R^2 = 0.99089278
% 2 curve: Var = 0.01819229, R^2 = 0.91345654
%
% C_T = 27.4085821
%
function struct_data = fun_data_mono_drug_A2_120()

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
gs_c_info = 'A2_120';             % String with drug related information
gs_LargeScale = 'off';            % 'on' = LargeScale, 'off' = MediumScale
gs_excel = 'on';                  % Write results as Excel file
                                  % (works only on Windows with Excel)
                                  
                                  
% Initial values                    Control Data:  (l0,l1,0, 0,  w0)
%                                   Monotherapy:   (l0,l1,k1,k2, w0)
%                                   Comb. therapy: (l0,l1,k1,psi,w0)
gv_start_vals = [0.1356, 0.2744, 1, 0.01, 0.0528];   


% Fitting data
%
% Obj_1 : Control Data 1 (Obj_1 is always control data)
% Obj_2 : Drug A2 120mg/kg po,  SID1-8, Start day 15
t_data_obj_1 = [6; 11; 13; 15; 18; 20; 22; 25; 27; 29; ];
y_data_obj_1 = [0.289324513; 0.5201940347; 0.8787273154; 1.034623291; ...
                1.636494591; 1.883689601;  2.273302339;  3.036456225; ...
                3.055183057; 3.789235538; ];
t_data_obj_2 = [6; 11; 13; 15; 18; 20; 22; 25; 27; 29; ];
y_data_obj_2 = [ 0.3126685895; 0.5021769672;  0.6207679758; 0.890892579; ...
                 0.871297959;  0.9944266331;  1.102733728;  1.146657959; ...
                 1.313193387;  1.608065881; ]; 

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
pk_obj_2 = [120, 5.54270775, 0.15492953, 1.51569635 ];
comp_obj_1 = 0;
comp_obj_2 = 1;

% Dosing schedule, [Number of doses, point in time (in hr), point in time (in hr), ... ]
sched_obj_1 = 0;
sched_obj_2 = [8, 360, 384, 408, 432, 456, 480, 504, 528 ];

% Write PK data in a cell
dosing(1) = struct('Comp',comp_obj_1,'pk_obj',pk_obj_1,'sched_obj',sched_obj_1); 
dosing(2) = struct('Comp',comp_obj_2,'pk_obj',pk_obj_2,'sched_obj',sched_obj_2);
gcell_dosing = struct2cell(dosing);

         

