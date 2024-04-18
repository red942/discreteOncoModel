function [yVals] = discreteModelValsMono(parameters, data, h)
%DISCRETEMODELVALS solves the difference equations and evaluates at given time values
    
    %assigns c(t) concentration function based on data type (drug one)
    c = @(t) concentrationVals(data, t, 1);
    
    %grabs and names the values
    tVals = data.tVals;
    l0 = parameters(1);
    l1 = parameters(2);
    k1 = parameters(3);
    k2 = parameters(4);
    w0 = parameters(5);
    
    %splits the interval of times based on the value of h
    hPointsNum = (tVals(size(tVals, 1)) - tVals(1))/h;

    %if the value is an integer, add one (it loses the end point otherwise)
    if (hPointsNum == ceil(hPointsNum))
        hPointsNum = hPointsNum + 1;
    else %if its not an integer, take the ceiling
        hPointsNum = ceil(hPointsNum);
    end

    %creates each x_i from the model
    x1 = zeros(1, hPointsNum);
    x2 = zeros(1, hPointsNum);
    x3 = zeros(1, hPointsNum);
    x4 = zeros(1, hPointsNum);
    
    x1(1) = w0;
    
    %Loops from the 1st element to the 2nd to last, assigning the next
    %value in each loop; So the 3rd element's loop assigns the 4th's value 
    for i = 1:hPointsNum-1
        w = x1(i) + x2(i) + x3(i) + x4(i);
        disp("w " + w + "x1 " + x1(i))
        
        x1(i+1) = x1(i) + h * ((2*l0*l1*x1(i)^2)/((l1+2*l0*x1(i))*w) - k2*c(tVals(1) + h*(i-1))*x1(i));
        if (x1(i + 1) == 0)
            disp("x1, uh oh")
        end
        x2(i+1) = x2(i) + h * (k2*c(tVals(1) + h*(i-1))*x1(i) - k1*x2(i));
 
        x3(i+1) = x3(i) + h * (k1*(x2(i) - x3(i)));

        x4(i+1) = x4(i) + h * (k1*(x3(i) - x4(i)));

    end

    w = x1 + x2 + x3 + x4;

    yVals = zeros(1, size(tVals, 1));
   


    for i = 1:size(tVals, 1)
        
        t = tVals(i);
        %finds the index of the nearest value to the left of t (also works
        %if it is the first value)
        hIndex = (t - tVals(1))/h;
        if (hIndex == ceil(hIndex))
            hIndex = hIndex + 1;
        else 
            hIndex = ceil(hIndex);
        end

        if ((hIndex == size(w, 2)))
            yVals(i) = w(hIndex);
            return
        end
        %Does a weighted average of the nearest two evaluated values
        yVals(i) = (w(hIndex)*(h - abs(t - w(hIndex))) + w(hIndex + 1)*(h - abs(t - w(hIndex + 1))))/h;

    end
end

