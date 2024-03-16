% 18.03.09 Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
% 
% Combination therapy: Control 1 + Drug A1 180 + Drug C 100
%
% Fitting results: 
% l0 = 0.15065595, l1 = 0.27124539, k1 = 1.60987451, psi = 0.83705221, w0 = 0.03581512 
% (CV = 31.8044)   (CV = 17.8324)   (CV = 31.4004)   (CV = 28.1088)    (CV = 82.0609)
% SofS = 0.19713797, MSE = 0.01314258
%
% 1 curve: Var = 0.01588778, R^2 = 0.99150777
% 2 curve: Var = 0.01039739, R^2 = 0.94993342
%
function struct_data = fun_data_comb_drug_A1_180_and_C_100()

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
gs_c_info = 'A1_180_+_C_100';     % String with drug related information
gs_LargeScale = 'off';            % 'on' = LargeScale, 'off' = MediumScale    
gs_excel = 'on';                  % Write results as Excel File
                                  % (works only on Windows with Excel)


% Initial Values                    Control Data:  (l0,l1,0, 0,  w0)
%                                   Monotherapy:   (l0,l1,k1,k2, w0)
%                                   Comb. therapy: (l0,l1,k1,psi,w0)
gv_start_vals = [0.14, 0.28, 1, 1, 0.04];

gv_k2a_k2b_fix = [0.1803, 0.0154]; % Fixed potency parameters of Drug 1 and 2


% Fitting Data
%
% Obj_1 : Control Data (Obj_1 is always control data)
% Obj_2 : Drug A1 180mg/kg po,  SID5-8,  Start day 11    
%        +Drug C  100mg/kg po,  SID1-12, Start day 11
t_data_obj_1 = [6; 11; 13; 15; 18; 20; 22; 25; 27; 29];
y_data_obj_1 = [0.289324513; 0.5201940347; 0.8787273154; 1.034623291;...
                1.636494591; 1.883689601;  2.273302339;  3.036456225;...
                3.055183057; 3.789235538];           
t_data_obj_2 = [6; 11; 13; 15; 18; 20; 22; 25; 27; 29; ];
y_data_obj_2 = [ 0.2592458009; 0.5132396492; 0.6465548943; 0.7440010459; ...
                 0.755396527;  0.826979724;  0.8199702174; 1.029438376; ...
                 1.255879676;  1.720928164 ];            

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
pk_obj_2_drug_1 = [115.7367, 11.1723, 3.11424791, 0.66324719,  4.41789241];
pk_obj_2_drug_2 = [16.8619,  26.9165, 0.17046905, 4.96189217, 84.90413544];
comp_obj_1 = 0;
comp_obj_2_drug_1 = 2;
comp_obj_2_drug_2 = 2;

% Dosing Schedule, [Number of doses, point in time (in hr), point in time (in hr), ... ]
sched_obj_1 = 0;
sched_obj_2_drug_1 = [4, 360, 384, 408, 432];
sched_obj_2_drug_2 = [12, 264, 288, 312, 336, 360, 384, 408, 432, 456, 480, 504, 528];
            
% Write PK data in a cell
dosing(1) = struct('Comp',comp_obj_1,'pk_obj',pk_obj_1,'sched_obj',sched_obj_1); 
dosing(2) = struct('Comp',comp_obj_2_drug_1,'pk_obj',pk_obj_2_drug_1,'sched_obj',sched_obj_2_drug_1);
dosing(3) = struct('Comp',comp_obj_2_drug_2,'pk_obj',pk_obj_2_drug_2,'sched_obj',sched_obj_2_drug_2);
gcell_dosing = struct2cell(dosing);
            








