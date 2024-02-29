% 18.03.09, Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
%
% Task: Print results 
%
function void = fun_print_result(erg,g_model,stat)
global g_n_curve;
global gs_excel;
global gs_c_info;
global gv_SofS_CT;

if (g_model == 0) 
    
    Covariance_Matrix = stat.cov;
    Correlation_Matrix = stat.corr;
    
    disp(' ');
    Covariance_Matrix
    disp(' ');
    Correlation_Matrix
    disp(' ');
    
    for i=1:g_n_curve
        text = sprintf('Objective curve %4.0f : Variance = %12.8f,   R^2 = %12.8f (%5.2e)',... 
                       i,stat.var(i),stat.R2(i),stat.R2(i)); 
        disp(text);
    end
    text = sprintf('\nParameter          (CV)            [Confidence Interval]');
    disp(text);
    text = sprintf('l0  = %12.8f (CV = %8.4f) [%8.4f, %8.4f]   (%5.2e)',... 
                   erg(1),stat.varcoeff(1),stat.conf_iv(1,1),stat.conf_iv(1,2),erg(1)); 
                   disp(text);
    text = sprintf('l1  = %12.8f (CV = %8.4f) [%8.4f, %8.4f]   (%5.2e)',... 
                   erg(2),stat.varcoeff(2),stat.conf_iv(2,1),stat.conf_iv(2,2),erg(2)); 
                   disp(text);
    text = sprintf('w0  = %12.8f (CV = %8.4f) [%8.4f, %8.4f]   (%5.2e)',...
                   erg(5),stat.varcoeff(5),stat.conf_iv(3,1),stat.conf_iv(3,2),erg(5)); 
                   disp(text);
    disp(' ');
    text = sprintf('MSE = %12.8f (%5.2e)',stat.mse,stat.mse);   
    disp(text);
                   
    % Create excel sheet               
    if strcmp(gs_excel,'on')
        M = zeros(3,4);
        M(:,1) = [erg(1); erg(2); erg(5)];
        M(:,2) = [stat.varcoeff(1); stat.varcoeff(2); stat.varcoeff(5)];
        M(:,3) = stat.conf_iv(:,1);
        M(:,4) = stat.conf_iv(:,2);
        d = {'','Parameter','CV','Conf.I.V. L','Conf.I.V. R'; ... 
             'l0','','','','';'l1','','','','';'w0','','','','';
             'SofS','','MSE','','';};
        if ispc
            name = strcat('results_',gs_c_info);
            warning off MATLAB:xlswrite:AddSheet
            cd results;
            xlswrite(name,d,'Table1','A1');
            xlswrite(name,M,'Table1','B2');
            xlswrite(name,gv_SofS_CT(1),'Table1','B5');
            xlswrite(name,stat.mse,'Table1','D5');
            cd ..;
        else
            disp('Excel output only supported on Windows systems');
        end
    end
                   
end
if (g_model == 1 || g_model == 2)
    Covariance_Matrix = stat.cov;
    Correlation_Matrix = stat.corr;
    
    disp(' ');
    Covariance_Matrix
    disp(' ');
    Correlation_Matrix
    disp(' ');
    
    
    for i=1:g_n_curve
        text = sprintf('Objective curve %4.0f : Variance = %12.8f,   R^2 = %12.8f  (%5.2e)',...
                       i,stat.var(i),stat.R2(i),stat.R2(i)); 
        disp(text);
    end
     
    text = sprintf('\nParameter          (CV)            [Confidence Interval]'); 
                   disp(text);
    text = sprintf('l0  = %12.8f (CV = %8.4f) [%8.4f, %8.4f]   (%5.2e)',...
                   erg(1),stat.varcoeff(1),stat.conf_iv(1,1),stat.conf_iv(1,2),erg(1)); 
                   disp(text);
    text = sprintf('l1  = %12.8f (CV = %8.4f) [%8.4f, %8.4f]   (%5.2e)',...
                   erg(2),stat.varcoeff(2),stat.conf_iv(2,1),stat.conf_iv(2,2),erg(2)); 
                   disp(text);
    text = sprintf('k1  = %12.8f (CV = %8.4f) [%8.4f, %8.4f]   (%5.2e)',...
                   erg(3),stat.varcoeff(3),stat.conf_iv(3,1),stat.conf_iv(3,2),erg(3)); 
                   disp(text);
    if (g_model == 1)
        text = sprintf('k2  = %12.8f (CV = %8.4f) [%8.4f, %8.4f]   (%5.2e)',...
                       erg(4),stat.varcoeff(4),stat.conf_iv(4,1),stat.conf_iv(4,2),erg(4)); 
                       disp(text);
    end
    if (g_model == 2)
         text = sprintf('psi = %12.8f (CV = %8.4f) [%8.4f, %8.4f]   (%5.2e)',...
                   erg(4),stat.varcoeff(4),stat.conf_iv(4,1),stat.conf_iv(4,2),erg(4)); 
                   disp(text);
    end
    text = sprintf('w0  = %12.8f (CV = %8.4f) [%8.4f, %8.4f]   (%5.2e)',...
                   erg(5),stat.varcoeff(5),stat.conf_iv(5,1),stat.conf_iv(5,2),erg(5)); 
                   disp(text);
    disp(' ');
    text = sprintf('MSE = %12.8f (%5.2e)',stat.mse,stat.mse);   
    disp(text);
                   
    if strcmp(gs_excel,'on') && (g_model == 1)
        M = zeros(5,4);
        M(:,1) = erg(:)';
        M(:,2) = stat.varcoeff(:)';
        M(:,3) = stat.conf_iv(:,1);
        M(:,4) = stat.conf_iv(:,2);
        d = {'','Parameter','CV','Conf.I.V. L','Conf.I.V. R'; ... 
            'l0','','','','';'l1','','','','';'k1','','','',''; ... 
            'k2','','','','';'w0','','','','';
            'SofS','','MSE','','';'CT','','','',''};
        if ispc
            cd results;
            name = strcat('results_',gs_c_info);
            warning off MATLAB:xlswrite:AddSheet
            xlswrite(name,d,'Table1','A1');
            xlswrite(name,M,'Table1','B2');
            xlswrite(name,gv_SofS_CT,'Table1','B7');
            xlswrite(name,stat.mse,'Table1','D7');
            cd ..;
         else
            disp('Excel output only supported on Windows systems');
         end
    end
    if strcmp(gs_excel,'on') && (g_model == 2)
        M = zeros(5,4);
        M(:,1) = erg(:)';
        M(:,2) = stat.varcoeff(:)';
        M(:,3) = stat.conf_iv(:,1);
        M(:,4) = stat.conf_iv(:,2);
        d = {'','Parameter','CV','Conf.I.V. L','Conf.I.V. R'; ... 
             'l0','','','','';'l1','','','','';'k1','','','',''; ... 
             'psi','','','','';'w0','','','','';
             'SofS','','MSE','','';'CT','','','',''};
        if ispc
            cd results;
            name = strcat('results_',gs_c_info);
            warning off MATLAB:xlswrite:AddSheet
            xlswrite(name,d,'Table1','A1');
            xlswrite(name,M,'Table1','B2');
            xlswrite(name,gv_SofS_CT,'Table1','B7');
            xlswrite(name,stat.mse,'Table1','D7');
            cd ..;
        else
            disp('Excel output only supported on Windows systems');
        end
    end
                    
end

void = 0;