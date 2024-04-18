classdef dataWrapper
    %DATAWRAPPER Class to store input data
    
    properties
        paramStartVals %starting guess for parameter values
        potencyParam %paremeters concerning the potency of drug in combination case
        tControl %control data t values (days)
        yControl %control data y values
        tVals %data's t values (days)
        yVals %data's y values
        type %1 is control, 2 is mono, 3 is combo
        pkParameters1 %paremeters used in the pk model for drug 1
        pkType1 %stores which pk values are provided for drug 1
        pkParameters2 %paremeters used in the pk model for drug 2
        pkType2 %stores which pk values are provided for drug 2
        dosingTime1 %times dose was recieved in hours
        dosingTime2 %times dose was recieved in hours
    end
    
    methods

        function obj = dataWrapper(paramStartVals, potencyParam, tControl, ...
                yContorl, tVals, yVals, type, pkParameters1, pkParameters2, ...
                dosingTime1, dosingTime2)
            %DATAWRAPPER Construct an instance of this class

            %universal parameters
            obj.paramStartVals = paramStartVals;
            obj.tControl = tControl;
            obj.yControl = yContorl;
            obj.type = type;

            %mono and comb
            if (type ~= 1)
                obj.tVals = tVals;
                obj.yVals = yVals;
                obj.pkParameters1 = pkParameters1;

                %assigns pkType
               if (size(pkParameters1, 2) == 4)
                    obj.pkType1 = 1;
               else
                   obj.pkType1 = 2;
               end

                obj.dosingTime1 = dosingTime1;
            end
            
            %only comb
            if (type == 3)
                obj.dosingTime2 = dosingTime2;
                obj.potencyParam = potencyParam;
                obj.pkParameters2 = pkParameters2;

               %assigns pkType
               if (size(pkParameters2, 2) == 4)
                    obj.pkType2 = 1;
               else
                   obj.pkType2 = 2;
               end
            end
        end

    end

end

