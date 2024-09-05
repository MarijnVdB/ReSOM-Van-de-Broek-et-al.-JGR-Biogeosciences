% Function to initialize the total C inputs for the soil
function [avg_daily_C_input] = initialize_Cinputs_soil(C_input, annualLitterCInputs, ...
            dates_spinup, dates_treatmentRun, dt, date_start_spinup, ...
            date_startWarming, endDate, seasonalCarbonInputs, variableTimeStep)
        
% If a variable timestep is used, these different timesteps are isolated
if variableTimeStep == 1
    d1 = max(dt); % The largest timestep (used for the spinup)
    d2 = min(dt); % The smallest timestep (used for the treatment run)
end

% First check if a constant input is used, or if data is scaled based on
% the litter input data
n = numel(annualLitterCInputs);

% If no annualLitterCInputs is provided, this means a constant C input for the
% entire run
if n == 0
    
    % If the inputs are constant throughout the year
    if seasonalCarbonInputs == 0
        
        % The number of timesteps is retrieved
        nT = numel([dates_spinup dates_treatmentRun]);
        avg_daily_C_input = ones(nT,1).*(C_input./365);
    
    % If the C inputs are distributed within the year
    elseif seasonalCarbonInputs == 1
        
        % The relative fraction of monthly C inputs
        relMonth = [0.004082665 0.003734601 0.025604134 0.057220096	...
            0.150748113 0.196453025 0.1723213 0.171196473 0.105433194	...
            0.0841401 0.024065782 0.005000516];

        % An id per month of the year
        relMonth_id = 1:1:12;
        
        % The complete simulated time series is constructed
        allDates = date_start_spinup:caldays(1):endDate;
        % The months are isolated
        month = allDates.Month;
        year = allDates.Year;
        % An array for the final C inputs in initialized
        inputsPerDay = NaN(size(allDates,2), 1);
        
        % A loop through all years to calculate the litter inputs
        uniqueYears = unique(allDates.Year);
        for i = 1:numel(uniqueYears)
            
            % The days for the particular year are isolated
            [r c] = find(year == uniqueYears(i));
            tmpDays = ones(numel(r), 1);
            
            % The first and last year are special, since it may not contain all litterfall days
            if i == 1 || i == numel(uniqueYears)
                
                % An array with the months of the current year is created
                [r c] = find(year == uniqueYears(i));
                tmpMonths = month(c);
                
                % The relative amount of C for these month is constructed
                uniqueMonths = unique(tmpMonths);
                tmpRel = NaN(size(tmpMonths));
                for j = 1:numel(uniqueMonths)
                    [r c] = find(tmpMonths == uniqueMonths(j));
                    % The number of days in this months is obtained
                    daysInThisMonth = numel(tmpMonths(tmpMonths == uniqueMonths(j)));
                    tmpRel(c) = relMonth(uniqueMonths(j))./daysInThisMonth;
                end
                
                % The inputs per day are formatted for the particular year
                inputsForThisYear = tmpRel.*C_input;
            
            % In all years but the first one and the last one
            elseif i > 1 && i < numel(uniqueYears)
                
                % An array with the dates of the current year (the account for leap years)
                datesInCurrentYear = datetime(['01/01/' num2str(uniqueYears(i))], 'InputFormat', 'dd/MM/yyyy'): ...
                    caldays(1):datetime(['31/12/' num2str(uniqueYears(i))], 'InputFormat', 'dd/MM/yyyy');
                % The months are obtained
                currentMonths = datesInCurrentYear.Month;
                
                % The relative amount of C for these month is constructed 
                uniqueMonths = unique(currentMonths);
                tmpRel = NaN(size(currentMonths));
                for j = 1:numel(uniqueMonths)
                    [r c] = find(currentMonths == uniqueMonths(j));
                    % The number of days in this months is obtained
                    daysInThisMonth = numel(currentMonths(currentMonths == uniqueMonths(j)));
                    tmpRel(c) = relMonth(uniqueMonths(j))./daysInThisMonth;
                end
                
                % The inputs per day are formatted for the particular year
                inputsForThisYear = tmpRel.*C_input;
                
            end
            
            % The calculated inputs per day are put in the right place in
            % the array with all simulateddays
            [r c] = find(year == uniqueYears(i));
            inputsPerDay(c) = inputsForThisYear;
            
        end

        % The output variable should contain the average litter C input per
        % model time step, so this has to be changed
        
        % First, the moving average is calculated, using the model timestep
        % as the window size
        if variableTimeStep == 0
            avgPerTimeStep = movmean(inputsPerDay, dt);
        elseif variableTimeStep == 1
            % This is done separately for the two time steps
            % The index of the first treatment timestep
            [r ind] = find(allDates == dates_treatmentRun(1));
            avgPerTimeStep1 = movmean(inputsPerDay(1:ind-1), d1);
            avgPerTimeStep2 = movmean(inputsPerDay(ind:end), d2);
            % Both are merged
            avgPerTimeStep = [avgPerTimeStep1; avgPerTimeStep2];
        end        
        
%         avgPerTimeStep = movmean(inputsPerDay, dt);
        
        % The indices of the modelled days are obtained
        [modelledDays c] = ismember(allDates, [dates_spinup dates_treatmentRun]);
        
        % Only the movmean values of the modelled days are retained
        avg_daily_C_input = avgPerTimeStep(modelledDays);       

    end
    
% If the soil C inputs are scaled relative to the litter C inputs
elseif n > 1
    
    % An array with the relative variation in litterfall C input
    % compared to the spinup amount is calculated
    relLitterInput = annualLitterCInputs./annualLitterCInputs(1);

    % The complete simulated time series is constructed
    allDates = date_start_spinup:caldays(1):endDate;
    % The years are isolated
    year = allDates.Year;
    uniqueYears = unique(year);

    % The annual soil C inputs are constructed
    annualCInputs = relLitterInput.*C_input;

    % This is converted to daily C inputs
    inputsPerDay = NaN(size(allDates,2),1);
        
    % If the inputs are constant throughout the year
    if seasonalCarbonInputs == 0
        
        for i = 1:numel(uniqueYears)
            % The indices of the current year are looked for
            [r c] = find(year == uniqueYears(i));
            % The number of total days in the current year
            potDays = numel(datetime(['01/01/' num2str(uniqueYears(i))], 'InputFormat', 'dd/MM/yyyy'): ...
                    caldays(1):datetime(['31/12/' num2str(uniqueYears(i))], 'InputFormat', 'dd/MM/yyyy'));
            % The daily C inputs
            inputsPerDay(c) = annualCInputs(i)./potDays;
            
        end
        
    % If the C inputs are distributed within the year
    elseif seasonalCarbonInputs == 1
        
        % The relative fraction of monthly C inputs
        relMonth = [0 0.001673936 0.018588437 0.05842778 0.155583436 ...
            0.205079066 0.176468904 0.174845489 0.10661491 0.082597127 ...
            0.020120914 0];
        % An id per month of the year
        relMonth_id = 1:1:12;
        
        % The complete simulated time series is constructed
        allDates = date_start_spinup:caldays(1):endDate;
        % The months are isolated
        month = allDates.Month;
        year = allDates.Year;
        % An array for the final C inputs in initialized
        inputsPerDay = NaN(size(allDates,2), 1);
        
        % A loop through all years to calculate the litter inputs
        uniqueYears = unique(allDates.Year);
        for i = 1:numel(uniqueYears)
            
            % The days for the particular year are isolated
            [r c] = find(year == uniqueYears(i));
            tmpDays = ones(numel(r), 1);
            
            % The first and last year are special, since it may not contain all litterfall days
            if i == 1 || i == numel(uniqueYears)
                
                % An array with the months of the current year is created
                [r c] = find(year == uniqueYears(i));
                tmpMonths = month(c);
                
                % The relative amount of C for these month is constructed
                uniqueMonths = unique(tmpMonths);
                tmpRel = NaN(size(tmpMonths));
                for j = 1:numel(uniqueMonths)
                    [r c] = find(tmpMonths == uniqueMonths(j));
                    % The number of days in this months is obtained
                    daysInThisMonth = numel(tmpMonths(tmpMonths == uniqueMonths(j)));
                    tmpRel(c) = relMonth(uniqueMonths(j))./daysInThisMonth;
                end
                
                % The inputs per day are formatted for the particular year
                inputsForThisYear = tmpRel.*annualCInputs(i);
            
            % In all years but the first one and the last one
            elseif i > 1 && i < numel(uniqueYears)
                
                % An array with the dates of the current year (the account for leap years)
                datesInCurrentYear = datetime(['01/01/' num2str(uniqueYears(i))], 'InputFormat', 'dd/MM/yyyy'): ...
                    caldays(1):datetime(['31/12/' num2str(uniqueYears(i))], 'InputFormat', 'dd/MM/yyyy');
                % The months are obtained
                currentMonths = datesInCurrentYear.Month;
                
                % The relative amount of C for these month is constructed 
                uniqueMonths = unique(currentMonths);
                tmpRel = NaN(size(currentMonths));
                for j = 1:numel(uniqueMonths)
                    [r c] = find(currentMonths == uniqueMonths(j));
                    % The number of days in this months is obtained
                    daysInThisMonth = numel(currentMonths(currentMonths == uniqueMonths(j)));
                    tmpRel(c) = relMonth(uniqueMonths(j))./daysInThisMonth;
                end
                
                % The inputs per day are formatted for the particular year
                inputsForThisYear = tmpRel.*annualCInputs(i);
                
            end
            
            % The calculated inputs per day are put in the right place in
            % the array with all simulateddays
            [r c] = find(year == uniqueYears(i));
            inputsPerDay(c) = inputsForThisYear;

        end
    
    end

    % The output variable should contain the average C input per
        % model time step, so this has to be changed
        
        % First, the moving average is calculated, using the model timestep
        % as the window size. 
        
        if variableTimeStep == 0
            avgPerTimeStep = movmean(inputsPerDay, dt);
        elseif variableTimeStep == 1
            % This is done separately for the two time steps
            % The index of the first treatment timestep
            [r ind] = find(allDates == dates_treatmentRun(1));
            avgPerTimeStep1 = movmean(inputsPerDay(1:ind-1), d1);
            avgPerTimeStep2 = movmean(inputsPerDay(ind:end), d2);
            % Both are merged
            avgPerTimeStep = [avgPerTimeStep1; avgPerTimeStep2];
        end
        
        % The indices of the modelled days are obtained
        [modelledDays c] = ismember(allDates, [dates_spinup dates_treatmentRun]);
        
        % Only the movmean values of the modelled days are retained
        avg_daily_C_input = avgPerTimeStep(modelledDays); 


end

