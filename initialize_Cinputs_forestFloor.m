% Function to initialize the total C inputs for the forest floor
function varargout = initialize_Cinputs_forestFloor(C_input, C_input_year, ...
            dates_spinup, dates_treatmentRun, dt, seasonalCarbonInputs, date_start_spinup, ...
            endDate, variableTimeStep)

%% The output arguments are defined

nOutputs = nargout;
varargout = cell(1,nOutputs);        

%% The function is run

% If a variable timestep is used, these different timesteps are isolated
if variableTimeStep == 1
    d1 = max(dt); % The largest timestep (used for the spinup)
    d2 = min(dt); % The smallest timestep (used for the treatment run)
end

% First check if a constant input or a time series is provided
n = numel(C_input);

% If just one number is provided, this means a constant C input for the entire run
if n == 1
    
    % If the inputs are constant throughout the year
    if seasonalCarbonInputs == 0
        
        % The calculation of the inputs depends on whether or not a
        % constant timestep is used throughout the simulation
        if variableTimeStep == 0 % The same timestep for the entire simulation is used
            
            % The number of timesteps is retrieved
            nT = numel([dates_spinup dates_treatmentRun]);
            avg_daily_C_input = ones(nT,1).*(C_input./365);
        
        elseif variableTimeStep == 1
            
            % The number of timesteps is retrieved
            nT = numel([dates_spinup dates_treatmentRun]);
            % The average daily C input is calculated
            dailyInput = (C_input./365);
            % The C input per timestep, depending on the size of the timestep, is calculated
            avg_daily_C_input = NaN(nT,1);
            % Timestep d1
            [r c] = find(dt == d1);
            avg_daily_C_input(r,1) = dailyInput;
            % Timestep d2
            [r c] = find(dt == d2);
            avg_daily_C_input(r,1) = dailyInput;
            
        end
    
    % If the litter inputs are occuring only between August and December
    elseif seasonalCarbonInputs == 1
        
        % The complete simulated time series is constructed
        allDates = date_start_spinup:caldays(1):endDate;
        % The months are isolated
        month = allDates.Month;
        year = allDates.Year;
        % The months in which litterfall takes place
        litterMonths = [8 9 10 11 12];
        % An array with 1 for days with litterfall, 0 for the other days
        [litterDays c] = ismember(month,litterMonths);
        % An array for the final C inputs in initialized
        inputsPerDay = NaN(size(allDates,2), 1);
        
        % A loop through all years to calculate the litter inputs
        uniqueYears = unique(allDates.Year);
        for i = 1:numel(uniqueYears)
            
            % The days for the particular year are isolated
            [r c] = find(year == uniqueYears(i));
            tmpDays = litterDays(c);
            
            % The first and last year are special, since it may not contain all litterfall days
            if i == 1 || i == numel(uniqueYears)
                
                % The potential number of litterfall days in the first year is calculated
                potDays = numel(datetime(['01/08/' num2str(uniqueYears(i))], 'InputFormat', 'dd/MM/yyyy'): ...
                    caldays(1):datetime(['31/12/' num2str(uniqueYears(i))], 'InputFormat', 'dd/MM/yyyy'));
                
                % The C inputs per day are calculated
                dailyInput = C_input/potDays;
                
                % The inputs per day are formatted for the particular year
                inputsForThisYear = tmpDays.*dailyInput;
            
            % In all years but the first one and the last one
            elseif i > 1 && i < numel(uniqueYears)
                
                % The number of days with litterfall in this year is obtained
                numDays = sum(tmpDays);
                
                % The C inputs per day are calculated
                dailyInput = C_input/numDays;
                
                % The inputs per day are formatted for the particular year
                inputsForThisYear = tmpDays.*dailyInput;
                
            end
            
            % The calculated inputs per day are put in the right place in
            % the array with all simulated days
            [r c] = find(year == uniqueYears(i));
            inputsPerDay(c) = inputsForThisYear;
            
        end

        % The output variable should contain the average litter C input per
        % model time step, so this has to be changed
        % This will depend whether a constant time step is used or not
        
            if variableTimeStep == 0 % The same timestep for the entire simulation is used
                
                % First, the moving average is calculated, using the model timestep
                % as the window size
                avgPerTimeStep = movmean(inputsPerDay, dt);

                % The indices of the modelled days are obtained
                [modelledDays c] = ismember(allDates, [dates_spinup dates_treatmentRun]);

                % Only the movmean values of the modelled days are retained
                avg_daily_C_input = avgPerTimeStep(modelledDays);

            elseif variableTimeStep == 1
                
                % The movmean for the first and second timesteps are calculated
                nd1 = numel(find(dt == d1));                                % The number of timesteps with length d1
                avgPerTimeStep_d1 = movmean(inputsPerDay(1:nd1*d1), d1);    % The average inputs per timestep when timestep is d1
                avgPerTimeStep_d2 = movmean(inputsPerDay(dt == d2), d2);    % The average inputs per timestep when timestep is d2
                avgPerTimeStep = [avgPerTimeStep_d1; avgPerTimeStep_d2];
                
                % The indices of the modelled days are obtained
                [modelledDays c] = ismember(allDates, [dates_spinup dates_treatmentRun]);

                % Only the movmean values of the modelled days are retained
                avg_daily_C_input = avgPerTimeStep(modelledDays);
                
            end

    end
    
% If a time series is provided
elseif n > 1
    
    % The number of columns that the C_input file has is determined
    % => 1 column: only data for the control treatment is provided
    % => 2 columns: data for the control and warmed run is provided
    [dummy nCol] = size(C_input);
    
    % The treatment years are isolated
    years_treatment = unique(dates_treatmentRun.Year);
    % Only litterfall inputs for the treatment years are retained
    [r c] = find(ismember(C_input_year, years_treatment));
    C_input = C_input(r,:);
    C_input_year = C_input_year(r);

    % The average of measured litterfall for the control run is calculated to be used
    % during the spinup and treatment years for which data is missing
    avgMeasuredInput = mean(C_input(:,1));

    % The complete simulated time series is constructed
    allDates = date_start_spinup:caldays(1):endDate;
    % The years for which the model is run are obtained
    uniqueYears = unique(allDates.Year)';

    % An array to store the annual litterfall is created
    annualInput = NaN(size(uniqueYears,1), nCol);

    % Get the indices of the treatment years for which input data is available
    [r c] = find(ismember(uniqueYears, C_input_year));

    % Fill in the measured data
    annualInput(r,:) = C_input;

    % In the other years, the average of measurements is filled
    annualInput(isnan(annualInput)) = avgMeasuredInput(1); % The value for the control run is used
        
    % If the inputs are constant throughout the year
    if seasonalCarbonInputs == 0
        
        % The complete simulated time series is constructed
        allDates = date_start_spinup:caldays(1):endDate;
        % The months are isolated
        year = allDates.Year;
        % An array for the final C inputs in initialized
        inputsPerDay = NaN(size(allDates,2), nCol);
        
        % A loop through all years to calculate the litter inputs
        uniqueYears = unique(allDates.Year);
        for i = 1:numel(uniqueYears)
            
            % The days for the particular year are isolated
            [r c] = find(year == uniqueYears(i));
            tmpDays = ones(numel(r), 1);
            
            % The first and last year are special, since it may not contain
            % all possible days of a year
            if i == 1 || i == numel(uniqueYears)
                
                % The number of days in the first or last year is calculated
                numDays = numel(tmpDays);
                
                % The C inputs per day are calculated
                dailyInput = annualInput(i,:)/numDays;
                
                % The inputs per day are formatted for the particular year
                inputsForThisYear = repmat(tmpDays,1,nCol).*dailyInput;
            
            % In all years but the first one and the last one
            elseif i > 1 && i < numel(uniqueYears)
                
                % The number of days in this year is calculated
                numDays = numel(tmpDays);
                
                % The C inputs per day are calculated
                dailyInput = annualInput(i,:)/numDays;
                
                % The inputs per day are formatted for the particular year
                inputsForThisYear = repmat(tmpDays,1,nCol).*dailyInput;
                
            end
            
            % The calculated inputs per day are put in the right place in
            % the array with all simulateddays
            [r c] = find(year == uniqueYears(i));
            inputsPerDay(c,:) = inputsForThisYear;
            
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
            avgPerTimeStep1 = movmean(inputsPerDay(1:ind-1,:), d1, 1);
            avgPerTimeStep2 = movmean(inputsPerDay(ind:end,:), d2, 1);
            % Both are merged
            avgPerTimeStep = [avgPerTimeStep1; avgPerTimeStep2];
        end
        
%         avgPerTimeStep = movmean(inputsPerDay, dt,1);
        
        % The indices of the modelled days are obtained
        [modelledDays c] = ismember(allDates, [dates_spinup dates_treatmentRun]);
        
        % Only the movmean values of the modelled days are retained
        avg_daily_C_input = avgPerTimeStep(modelledDays,:);
        
    % If the litter inputs are occuring only between August and December
    elseif seasonalCarbonInputs == 1
        
        % The complete simulated time series is constructed
        allDates = date_start_spinup:caldays(1):endDate;
        % The months are isolated
        month = allDates.Month;
        year = allDates.Year;
        % The months in which litterfall takes place
        litterMonths = [8 9 10 11 12];
        % An array with 1 for days with litterfall, 0 for the other days
        [litterDays c] = ismember(month,litterMonths);
        % An array for the final C inputs in initialized
        inputsPerDay = NaN(size(allDates,2), nCol);
        
        % A loop through all years to calculate the litter inputs
        uniqueYears = unique(allDates.Year);
        for i = 1:numel(uniqueYears)
            
            % The days for the particular year are isolated
            [r c] = find(year == uniqueYears(i));
            tmpDays = litterDays(c);
            
            % The first and last year are special, since it may not contain all litterfall days
            if i == 1 || i == numel(uniqueYears)
                
                % The potential number of litterfall days in the first year is calculated
                potDays = numel(datetime(['01/08/' num2str(uniqueYears(i))], 'InputFormat', 'dd/MM/yyyy'): ...
                    caldays(1):datetime(['31/12/' num2str(uniqueYears(i))], 'InputFormat', 'dd/MM/yyyy'));
                
                % The C inputs per day are calculated
                dailyInput = annualInput(i,:)./potDays;
                
                % The inputs per day are formatted for the particular year
                inputsForThisYear = repmat(tmpDays,nCol,1).*repmat(dailyInput',1,size(tmpDays,2));
            
            % In all years but the first one and the last one
            elseif i > 1 && i < numel(uniqueYears)
                
                % The number of days with litterfall in this year is obtained
                numDays = sum(tmpDays);
                
                % The C inputs per day are calculated
                dailyInput = annualInput(i,:)/numDays;
                
                % The inputs per day are formatted for the particular year
                inputsForThisYear = repmat(tmpDays,nCol,1).*repmat(dailyInput',1,size(tmpDays,2));
                
            end
            
            % The calculated inputs per day are put in the right place in
            % the array with all simulateddays
            [r c] = find(year == uniqueYears(i));
            inputsPerDay(c,:) = inputsForThisYear';
            
        end

        % The output variable should contain the average litter C input per
        % model time step, so this has to be changed
        
        % First, the moving average is calculated, using the model timestep
        % as the window size. 
        
        if variableTimeStep == 0
            avgPerTimeStep = movmean(inputsPerDay, dt);
        elseif variableTimeStep == 1
            % This is done separately for the two time steps
            % The index of the first treatment timestep
            [r ind] = find(allDates == dates_treatmentRun(1));
            avgPerTimeStep1 = movmean(inputsPerDay(1:ind-1,:), d1, 1);
            avgPerTimeStep2 = movmean(inputsPerDay(ind:end,:), d2, 1);
            % Both are merged
            avgPerTimeStep = [avgPerTimeStep1; avgPerTimeStep2];
        end

        % The indices of the modelled days are obtained
        [modelledDays c] = ismember(allDates, [dates_spinup dates_treatmentRun]);
        
        % Only the movmean values of the modelled days are retained
        avg_daily_C_input = avgPerTimeStep(modelledDays,:);
  
    end
    
end

%% The outputs are formatted
if n == 1
    
    varargout{1} = avg_daily_C_input;

elseif n > 1
    
    varargout{1} = avg_daily_C_input;
    varargout{2} = annualInput;
    
end


end