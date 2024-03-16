% 18.03.09, Gilbert Koch (University of Konstanz), gilbert.koch@uni-konstanz.de
% http://www.math.uni-konstanz.de/numerik/personen/koch/
%
% Task: Calculates the statistics after a successful fitting process
%
% This function is based on
% [D'ARGENIO D.Z., SCHUMITZKY A., ADAPT User Guide, Chapter 3 ]
%
function struct_erg = fun_stats(erg,struct_sol,struct_data)

global g_n_curve;                    % Number of fitting curves
global gv_anz_data;                  % Vector with number of measurements for each 
                                     % curve/object k
global g_anz_data_all;               % Number of all measurements
global g_model;                      % Model 

if (g_model == 0) 
    l = g_n_curve;                   % Number of objects (=curves)
    p = 3;                           % Number of parameter
    P = zeros(g_anz_data_all,p);     % Allocate memory for the derivatives of the
                                     % objective function with respect to
                                     % the parameter    
    act_pos = 0;
    for k=1:l                        % Loop over number of objects

        m = gv_anz_data(k);          % Number of observations of object k 

        % Calculate position in t- and y-vector
        if (k>1)
            act_pos = act_pos + gv_anz_data(k-1);
        end
  
     
        % Calculate the derivatives of the objective function with respect to
        % the parameter
        sum = 0;                  
        for i=1:m
            y = deval(struct_sol(k),struct_data(k).t_data(i));
            w = y(1)+y(2)+y(3)+y(4);
            
            % Calculate the sum of squares 
            sum = sum + (struct_data(k).y_data(i) - w)^2;

            % Calculate derivatives
            J(1) = y(5)  + y(6)  + y(7)  + y(8);
            J(2) = y(9)  + y(10) + y(11) + y(12);
            %J(3) = y(13) + y(14) + y(15) + y(16);
            %J(4) = y(17) + y(18) + y(19) + y(20);
            J(5) = y(21) + y(22) + y(23) + y(24);
    
            P(i+act_pos,:) = [J(1), J(2), J(5)]; 
        end
        
        % Variance of the residuals with m-(p/l) degrees of freedom
        var(k) = (1/(m-(p/l)))*sum;

        % Calculate predictions
        y_pred = zeros(1,m);
        for i=1:m
            y = deval(struct_sol(k),struct_data(k).t_data(i));
            y_pred(i) = y(1)+y(2)+y(3)+y(4);       
        end
        
        % Calculate R^2 
        sum1 = 0; sum2 = 0; sum3 = 0;
        for i=1:m
            sum1 = sum1 + (y_pred(i)-mean(y_pred))*(struct_data(k).y_data(i)-mean(struct_data(k).y_data));
            sum2 = sum2 + (y_pred(i)-mean(y_pred))^2;
            sum3 = sum3 + (struct_data(k).y_data(i)-mean(struct_data(k).y_data))^2;
        end
        v_R2(k) = sum1^2 /(sum2*sum3);         
        
        % Observed vs. predicted plot
        clear obse_plot pred_plot;
        for i=1:m
            obse_plot(i) = struct_data(k).y_data(i);
            pred_plot(i) = y_pred(i);
        end
        figure(10+k); plot(obse_plot,pred_plot,'black s'); hold on; 
        figure(10+k); plot([0,y_pred(m)],[0,y_pred(m)],'red'); hold on; 
        xlabel('Observed'); ylabel('Predicted'); 
        text = sprintf('Curve %4.0f',k);
        title(text);
        
    end
    
    % Calculate covariance matrix 
    R = zeros(g_anz_data_all,g_anz_data_all);
    for k=1:l
        m = gv_anz_data(k); 
        R(m*(k-1)+1:k*m,m*(k-1)+1:k*m) = var(k) .* eye(m,m);
    end
    try
        A = (P'*P)^(-1);
        cov =  A * (P'*R*P) * A;
    catch
        disp('PT*P singular, calculation of CV(%) not possible'); 
        cov = zeros(p,p);
    end
    
    % Standard deviation of the parameters
    stddev_params = [sqrt(cov(1,1)),sqrt(cov(2,2)),sqrt(cov(3,3))];

    % Correlation matrix
    for i=1:p
        for j=1:p
            corr(i,j) = cov(i,j) / (stddev_params(i)*stddev_params(j));
        end
    end
    
    % Coefficient of variation 
    varcoeff_params = [(stddev_params(1)*100)/erg(1),(stddev_params(2)*100)/erg(2),0,0,(stddev_params(3)*100)/erg(5)]; 
               
    % 95% confidence interval
    t_val = fun_calc_t_val_95per(round(m-(p/l)));      
    conf_iv(1,1) = erg(1) - sqrt(cov(1,1))*t_val; conf_iv(1,2) = erg(1) + sqrt(cov(1,1))*t_val;
    conf_iv(2,1) = erg(2) - sqrt(cov(2,2))*t_val; conf_iv(2,2) = erg(2) + sqrt(cov(2,2))*t_val;
    conf_iv(3,1) = erg(5) - sqrt(cov(3,3))*t_val; conf_iv(3,2) = erg(5) + sqrt(cov(3,3))*t_val;
    
    % MSE     
    mse = sum /(g_anz_data_all-p);     
        
    % Store statistic results in a struct
    struct_erg = struct('var',var,'varcoeff',varcoeff_params,'conf_iv',conf_iv,'cov',cov,'corr',corr,...
                        'R2',v_R2,'mse',mse);                   
end
if (g_model == 1 || g_model == 2)

    l = g_n_curve;                   % Number of objects
    p = 5;                           % Number of parameter
    P = zeros(g_anz_data_all,p);     % Allocate memory for the derivatives of the
                                     % objective function with respect to
                                     % the parameter
    
    act_pos = 0;
    for k=1:l

        m = gv_anz_data(k);          % Number of observations of object k 

        % Calculate position in t- and y-vector
        if (k>1)
            act_pos = act_pos + gv_anz_data(k-1);
        end
  
        sum = 0;
        for i=1:m
            y = deval(struct_sol(k),struct_data(k).t_data(i));
            
            w = y(1)+y(2)+y(3)+y(4);
            
            w_tmp(i) = w;
            
            % Calculate sum of squares
            res(i) = (struct_data(k).y_data(i) - w);
            sum = sum + res(i)^2;
           
            J(1) = y(5)  + y(6)  + y(7)  + y(8);
            J(2) = y(9)  + y(10) + y(11) + y(12);
            J(3) = y(13) + y(14) + y(15) + y(16);
            J(4) = y(17) + y(18) + y(19) + y(20);
            J(5) = y(21) + y(22) + y(23) + y(24);
    
            P(i+act_pos,:) = [J(1), J(2), J(3), J(4), J(5)]; 
        end
               
        % Variance of the residuals
        var(k) = (1/(m-(p/l)))*sum;      
            
        % Observed vs. predicted plot
        clear obse_plot pred_plot;
        for i=1:m
            obse_plot(i) = struct_data(k).y_data(i);
            pred_plot(i) = w_tmp(i);
        end
        figure(10+k); plot(obse_plot,pred_plot,'black s'); hold on; 
        figure(10+k); plot([0,w_tmp(m)],[0,w_tmp(m)],'red'); hold on; 
        xlabel('Observed'); ylabel('Predicted'); 
        text = sprintf('Curve %4.0f',k);
        title(text);
          
        % R^2
        y_pred = zeros(1,m);
        for i=1:m
            y = deval(struct_sol(k),struct_data(k).t_data(i));
            y_pred(i) = y(1)+y(2)+y(3)+y(4);       
        end
        sum1 = 0; sum2 = 0; sum3 = 0;
        for i=1:m
            sum1 = sum1 + (y_pred(i)-mean(y_pred))*(struct_data(k).y_data(i)-mean(struct_data(k).y_data));
            sum2 = sum2 + (y_pred(i)-mean(y_pred))^2;
            sum3 = sum3 + (struct_data(k).y_data(i)-mean(struct_data(k).y_data))^2;
        end       
        v_R2(k) = sum1^2 /(sum2*sum3);      
        
    end

    % MSE
    sum_res = 0;
    for k=1:l
        for i=1:m
            y = deval(struct_sol(k),struct_data(k).t_data(i));        
            w = y(1)+y(2)+y(3)+y(4);

            sum_res = sum_res + (struct_data(k).y_data(i) - w)^2;
        end
    end
    mse = sum_res/(g_anz_data_all-p);
      
    
    % Calculate covariance matrix
    R = zeros(g_anz_data_all,g_anz_data_all);
    for k=1:l
        m = gv_anz_data(k); 
        R(m*(k-1)+1:k*m,m*(k-1)+1:k*m) = var(k) .* eye(m,m);
    end
    try
        A = (P'*P)^(-1);
    catch 
        disp('PT*P singular, calculation of CV(%) not possible'); 
        A = zeros(p,p);
    end
    cov =  A * (P'*R*P) * A;  
        
    % Standard deviation
    stddev_params = [sqrt(cov(1,1)),sqrt(cov(2,2)),sqrt(cov(3,3)),sqrt(cov(4,4)),sqrt(cov(5,5))];

    % Coefficient of variation
    varcoeff_params = [(stddev_params(1)*100)/erg(1),(stddev_params(2)*100)/erg(2),(stddev_params(3)*100)/erg(3),... 
                        stddev_params(4)*100/erg(4),(stddev_params(5)*100)/erg(5)]; 
    
    % Correlation matrix
    corr = zeros(p,p);
    for i=1:p
        for j=1:p
            corr(i,j) = cov(i,j) / (stddev_params(i)*stddev_params(j));
        end
    end
    
    % 95% confidence interval
    t_val = fun_calc_t_val_95per(round(m-(p/l)));  
    conf_iv = zeros(p,2);    
    for i=1:p
        conf_iv(i,1) = erg(i) - sqrt(cov(i,i))*t_val; 
        conf_iv(i,2) = erg(i) + sqrt(cov(i,i))*t_val;
    end
    
    % Store statistic results in a struct
    struct_erg = struct('var',var,'varcoeff',varcoeff_params,'conf_iv',conf_iv,'cov',cov,'corr',corr,...
                        'R2',v_R2,'mse',mse);        

end