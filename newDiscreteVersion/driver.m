%DRIVER This is the driver file for the discrete version of the mono vs combination drug therapy code
%   If you wish to run the code, run this file
%   Made by Cameron Miller
%   Originally based on code by Gilbert Koch

%clears old variables and adds the folders to the path
clear;  
clear global;

%gets input for data type
inputNum = input("Which type of data: type 1 for control, 2 for mono, 3 for comb\n");

%gets input for which data within type and retrieves it
switch (inputNum)
    case 1
        inputNum = input("Which data: 1 -> 1, 2 -> 2, 3 -> A549, 4 -> MDA, 5 -> PC3, 6 -> RKO\n");
        data = getData("cont", inputNum);
    case 2
        inputNum = input("Which data: 1 -> A1_180, 2 -> A2_120, 3 -> B_100, 4 -> C_100, 5 -> C_150_1, 6 -> C_150_2\n");
        data = getData("mono", inputNum);
    case 3 
        inputNum = input("Which data: 1 -> A1 + C, 2 -> A2 + C, 3 -> B + C\n");
        data = getData("comb", inputNum);
    otherwise
        error("invalid selection, run program again")
end

%TEST PLOT, TEMPORARY
%plot(data.tControl, data.yControl, 'rx', data.tVals, data.yVals, "bx")


% l0 =  0.12067358, l1 = 0.31253305, k1 = 2.96803993, k2 = 0.00732078, w0 = 0.05871917 
testVals = discreteModelValsMono([0.12067358, 0.31253305, 2.96803993, 0.00732078, 0.05871917], data, 1/24);
%disp(testVals)
plot(data.tControl, data.yControl, 'rx', data.tVals, testVals, "bx")
