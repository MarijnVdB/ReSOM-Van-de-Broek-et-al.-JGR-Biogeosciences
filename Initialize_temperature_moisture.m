%% In this script, the temperature and moisture inputs for every timestep are loaded

% --------------------------
% Temperature
% --------------------------

%% The data is loaded

% ---------------------------------------------------------
% The historic temperature data (1 Jan 1965 - 30 July 1991)
% ---------------------------------------------------------

% This data only has to be loaded if temperature measurements are used
if useRealTemp == 1

    raw = readtable('Soil temperature 1978 - 2003.csv');

    dates_historic_temperature = datetime(raw.Date,'InputFormat','dd/MM/yyyy');
    dates_historic_temperature.Format = 'dd/MM/yyyy';

    historic_temperature_data = raw.Soil_temperature;

end

% ---------------------------------------------------------
% The historic average annual temperature curve
% ---------------------------------------------------------

day_historic_annualTemperatureCurve = xlsread('Annual curve soil temperature.xlsx', 'Sheet1', 'A2:A367');
month_historic_annualTemperatureCurve = xlsread('Annual curve soil temperature.xlsx', 'Sheet1', 'B2:B367');
historic_annualTemperatureCurve_dayAndMonth = [day_historic_annualTemperatureCurve  month_historic_annualTemperatureCurve];
clear day_historic_annualTemperatureCurve month_historic_annualTemperatureCurve

historic_annualTemperatureCurve_data = xlsread('Annual curve soil temperature.xlsx', 'Sheet1', 'C2:C367');

% -----------------------------------------------------------
% The measured temperature at 4 cm depth at the control plots
% -----------------------------------------------------------

% This data only has to be loaded if temperature measurements are used
if useRealTemp == 1
    
    raw = readtable('Soil temperatures control treatment.csv');
    dates_temperature_controlPlots = datetime(raw.Date,'InputFormat','dd/MM/yyyy');
    dates_temperature_controlPlots.Format = 'dd/MM/yyyy';

    temperature_controlPlots_data = raw.Soil_temperature;

end

% -----------------------------------------------------------
% The measured temperature at 4 cm depth at the warmed plots
% -----------------------------------------------------------

% This data only has to be loaded if temperature measurements are used
if useRealTemp == 1

    raw = readtable('Soil temperatures warming treatment.csv');
    dates_temperature_warmedPlots = datetime(raw.Date,'InputFormat','dd/MM/yyyy');
    dates_temperature_warmedPlots.Format = 'dd/MM/yyyy';

    temperature_warmedPlots_data = raw.Soil_temperature;

end

%% The temperature data is formatted according to the modelled time step

% -----------------------------------------------------------
% The temperature for the spin-up run is defined, before 1965
% -----------------------------------------------------------

% First check if articial or real temperature measurements are used

% If the temperature measurement are used
if useRealTemp == 1

    % As we only have temperature data from 1 Jan 1965 onwards,
    % it is checked if the spin-up run goes further back in time

    % First, the average time step for every day of the year for the avg
    % temperature curve is calculated
    % The annual T curve is replicated 3 times
    Tcurve_T_long = repmat(historic_annualTemperatureCurve_data,3,1);

    % If the spin-up run starts before 1965
    if date_start_spinup < dates_historic_temperature(1)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The temperature up to 1 Jan 1965 is calculated, based on the timestep
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % This is done using the annual average temperature curve
        dates_for_avg_annual_Tcurve = dates_spinup(dates_spinup < dates_historic_temperature(1));

        % For every date, the average temperature for the length of the time
        % step is calculated
        % For every date, the day and month are extracted
        spinupDates_day = dates_for_avg_annual_Tcurve.Day;
        spinupDates_month = dates_for_avg_annual_Tcurve.Month;
        spinupDates_dayAndMonth = [spinupDates_day' spinupDates_month'];
        clear spinupDates_day spinupDates_month
        
        % First we need to check if a variable or constant time step is used
        if variableTimeStep == 0
            
            % The moving average of the middle dates (1 year) is calculated
            Tcurve_movmean = movmean(Tcurve_T_long(numel(historic_annualTemperatureCurve_data)+1:numel(historic_annualTemperatureCurve_data)*2),dt);
            
            % For every spinup date, the movmean of the annual T curve is stored
            temperature_spinup_Tcurve = NaN(size(spinupDates_dayAndMonth,1),1);

            for i = 1:size(spinupDates_dayAndMonth,1)

                % The day and month of every timestep is matched with the Tcurve data
                [r, c] = ismember(spinupDates_dayAndMonth(i,:), historic_annualTemperatureCurve_dayAndMonth, 'rows');
                temperature_spinup_Tcurve(i,1) = Tcurve_movmean(c, 1);

            end
            
        elseif variableTimeStep == 1
            
            % The moving average of the middle dates (1 year) is calculated
            Tcurve_movmean1 = movmean(Tcurve_T_long(numel(historic_annualTemperatureCurve_data)+1:numel(historic_annualTemperatureCurve_data)*2),d1);
            Tcurve_movmean2 = movmean(Tcurve_T_long(numel(historic_annualTemperatureCurve_data)+1:numel(historic_annualTemperatureCurve_data)*2),d2);
            
            % The timesteps before the first date with temperature measurements are obtained
            [r, c] = find(dates_spinup < dates_historic_temperature(1));

            % Check if the small time step (d2) is present
            [r_smallTimeStep c_smallTimeStep] = find(dt_array_full(c) == d2);

            if ~isempty(r_smallTimeStep)     % If the small timestep is present

                % Find the indices for the large timesteps (d1)
                [r_largeTimeStep c_largeTimeStep] = find(dt_array_full(c) == d1);

                % For every spinup date, the movmean of the annual T curve is stored
                temperature_spinup_Tcurve1 = NaN(size(r_largeTimeStep,1),1);

                for i = 1:size(r_largeTimeStep,1)

                    % The day and month of every timestep is matched with the Tcurve data
                    tmp = spinupDates_dayAndMonth(r_largeTimeStep,:);
                    [r1, c1] = ismember(tmp(i,:), historic_annualTemperatureCurve_dayAndMonth, 'rows');
                    temperature_spinup_Tcurve1(i,1) = Tcurve_movmean1(c1, 1);

                end

                % The T curve data for the small timestep is stored
                temperature_spinup_Tcurve2 = NaN(size(r_smallTimeStep,1),1);

                for i = 1:size(r_smallTimeStep,1)

                    % The day and month of every timestep is matched with the Tcurve data
                    tmp = spinupDates_dayAndMonth(r_smallTimeStep,:);
                    [r1, c1] = ismember(tmp(i,:), historic_annualTemperatureCurve_dayAndMonth, 'rows');
                    temperature_spinup_Tcurve2(i,1) = Tcurve_movmean2(c1, 1);

                end

                % Both arrays are merged
                temperature_spinup_Tcurve = [temperature_spinup_Tcurve1; temperature_spinup_Tcurve2];
        
        else                % If the small timestep is not present
            
            Tcurve_movmean = movmean(Tcurve_T_long(numel(historic_annualTemperatureCurve_data)+1:numel(historic_annualTemperatureCurve_data)*2),d1);
            
            % For every spinup date, the movmean of the annual T curve is stored
            temperature_spinup_Tcurve = NaN(size(spinupDates_dayAndMonth,1),1);

            for i = 1:size(spinupDates_dayAndMonth,1)

                % The day and month of every timestep is matched with the Tcurve data
                [r, c] = ismember(spinupDates_dayAndMonth(i,:), historic_annualTemperatureCurve_dayAndMonth, 'rows');
                temperature_spinup_Tcurve(i,1) = Tcurve_movmean(c, 1);
            end

        end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The daily temperature data after 1965 is added
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if variableTimeStep == 0
        
            dates_for_daily_temperature_spinup = dates_spinup(dates_spinup >= dates_historic_temperature(1));

            % The average temperature per timestep is calculated
            daily_temp_spinup_movmean = movmean(historic_temperature_data, dt);

            % The dates are matched and the temperature variable is updated
            temperature_Tcurve = NaN(size(dates_for_daily_temperature_spinup,2),1);

            for i = 1:size(dates_for_daily_temperature_spinup, 2)

                [r, c] = ismember(dates_for_daily_temperature_spinup(1,i), dates_historic_temperature, 'rows');
                Tcurve(i,1) = daily_temp_spinup_movmean(c,1);

            end

            % The soil temperature for the spinup are merged
            temperature_spinup = [temperature_spinup_Tcurve; Tcurve];
            
        elseif variableTimeStep == 1
            
            dates_for_daily_temperature_spinup = dates_spinup(dates_spinup >= dates_historic_temperature(1));

            % The average temperature per timestep is calculated
            daily_temp_spinup_movmean1 = movmean(historic_temperature_data, d1);        % For the large timestep (d1)
            daily_temp_spinup_movmean2 = movmean(historic_temperature_data, d2);        % For the small timestep (d2)
            
            % Find the indices for the large and small timesteps
            [r c] = find(dates_spinup >= dates_historic_temperature(1));
            [r_largeTimeStep c_largeTimeStep] = find(dt_array_full(c) == d1);
            [r_smallTimeStep c_smallTimeStep] = find(dt_array_full(c) == d2);
            % The indices are converted to the indices in the full dt_array
            r_largeTimeStep = r_largeTimeStep + c(1)-1;
            r_smallTimeStep = r_smallTimeStep + c(1)-1;
            
            % The temperature for the large timesteps is stored
            Tcurve1 = NaN(size(r_largeTimeStep,1),1);

            for i = 1:size(r_largeTimeStep, 1)

                tmp = dates_spinup(1,r_largeTimeStep);
                [r1 c1] = ismember(tmp(1,i), dates_historic_temperature, 'rows');
                Tcurve1(i,1) = daily_temp_spinup_movmean1(c1,1);

            end
            
            % The temperature for the small timesteps is stored
            Tcurve2 = NaN(size(r_smallTimeStep,1),1);

            for i = 1:size(r_smallTimeStep, 1)

                tmp = dates_spinup(1,r_smallTimeStep);
                [r1, c1] = ismember(tmp(1,i), dates_historic_temperature, 'rows');
                Tcurve2(i,1) = daily_temp_spinup_movmean2(c1,1);

            end
            
            Tcurve = [Tcurve1; Tcurve2];            

            % The soil temperature for the spinup are merged
            temperature_spinup = [temperature_spinup_Tcurve; Tcurve];
 
        end

    else    % If the spin-up run starts after 1965
        
        if variableTimeStep == 0
            
            % The daily temperature data after 1965 is added
            dates_for_daily_temperature_spinup = dates_spinup;

            % The average temperature per timestep is calculated
            daily_temp_spinup_movmean = movmean(historic_temperature_data, dt);

            % The dates are matched and the temperature variable is updated
            temperature_historicT = NaN(size(dates_for_daily_temperature_spinup,2),1);

            for i = 1:size(dates_for_daily_temperature_spinup, 2)

                [r c] = ismember(dates_for_daily_temperature_spinup(1,i), dates_historic_temperature, 'rows');
                historicT(i,1) = daily_temp_spinup_movmean(c,1);

            end

        % The soil temperature for the spinup are merged
        temperature_spinup = historicT;
            
        elseif variableTimeStep == 1
            
            dates_for_daily_temperature_spinup = dates_spinup;

            % The average temperature per timestep is calculated
            daily_temp_spinup_movmean1 = movmean(historic_temperature_data, d1);        % For the large timestep (d1)
            daily_temp_spinup_movmean2 = movmean(historic_temperature_data, d2);        % For the small timestep (d2)
            
            % Find the indices for the large and small timesteps
            [r_largeTimeStep ~] = find(dt_array_full(1:numel(dates_spinup)) == d1);
            [r_smallTimeStep ~] = find(dt_array_full(1:numel(dates_spinup)) == d2);
            
            % The dates are matched and the temperature variable is updated
            temperature_historicT1 = NaN(size(r_largeTimeStep,2),1);
            
            for i = 1:size(r_largeTimeStep, 1)

                tmp = dates_spinup(1,r_largeTimeStep);
                [r1 c1] = ismember(tmp(1,i), dates_historic_temperature, 'rows');
                temperature_historicT1(i,1) = daily_temp_spinup_movmean1(c1,1);

            end
            
            % The temperature for the small timesteps is stored
            temperature_historicT2 = NaN(size(r_smallTimeStep,1),1);

            for i = 1:size(r_smallTimeStep, 1)

                tmp = dates_spinup(1,r_smallTimeStep);
                [r1 c1] = ismember(tmp(1,i), dates_historic_temperature, 'rows');
                temperature_historicT2(i,1) = daily_temp_spinup_movmean2(c1,1);

            end

            temperature_spinup = [temperature_historicT1; temperature_historicT2]; 
                        
        end

    end

    % --------------------------------------------------------------------
    % The temperature for the control and heated treatment are constructed
    % --------------------------------------------------------------------
    
    % If a variable timestep is used, the timestep for the treatment run
    % should always be the small timestep (d2)
    if variableTimeStep == 1
        dt = d2;
    end

    % Temperature for the control run

    % The daily temperature data after 1965 is added
    dates_for_daily_temperature_treatmentRun = dates_treatmentRun;

    % The average temperature per timestep is calculated
    daily_temp_spinup_movmean = movmean(temperature_controlPlots_data, dt);

    % The dates are matched and the temperature variable is updated
    temperature_controlRun = NaN(size(dates_for_daily_temperature_treatmentRun,2),1);

    for i = 1:size(dates_for_daily_temperature_treatmentRun, 2)

        [r c] = ismember(dates_for_daily_temperature_treatmentRun(1,i), dates_temperature_controlPlots, 'rows');
        temperature_controlRun(i,1) = daily_temp_spinup_movmean(c,1);

    end

    % Temperature for the warming run  

    % The daily temperature data after 1965 is added
    dates_for_daily_temperature_treatmentRun = dates_treatmentRun;

    % The average temperature per timestep is calculated
    daily_temp_spinup_movmean = movmean(temperature_warmedPlots_data, dt);

    % The dates are matched and the temperature variable is updated
    temperature_warmingRun = NaN(size(dates_for_daily_temperature_treatmentRun,2),1);

    for i = 1:size(dates_for_daily_temperature_treatmentRun, 2)

        [r c] = ismember(dates_for_daily_temperature_treatmentRun(1,i), dates_temperature_controlPlots, 'rows');
        temperature_warmingRun(i,1) = daily_temp_spinup_movmean(c,1);

    end  
    
    % -----------------------------------------------------------
    % Merging temperatures spinup and treatment run
    % -----------------------------------------------------------

    soilTemperature_spinupAndControl = [temperature_spinup; temperature_controlRun];
    soilTemperature_spinupAndWarming = [temperature_spinup; temperature_warmingRun];
    soilTemperature_warmingOnly = temperature_warmingRun;
    
% If the artificial temperature scenario is run
elseif useArtificialTemp == 1
    
    % -------------------------------------------------------------------
    % The temperature for the spinup and control treatment are calculated
    % -------------------------------------------------------------------
    
    % Tcurve_T_long is the annual temperature curve that will be used 
    Tcurve_T_long = repmat(historic_annualTemperatureCurve_data,3,1);

    % If the timestep is the same throughout the simulation
    if variableTimeStep == 0
        
        % The moving average of the middle dates (1 year) is calculated
        Tcurve_movmean = movmean(Tcurve_T_long(numel(historic_annualTemperatureCurve_data)+1:numel(historic_annualTemperatureCurve_data)*2),dt);

        % The temperature up to 1 Jan 1965 are calculated, based on the timestep
        % This is done using the annual average temperature curve
        dates_all_artificialT = dates_all;

        % For every date, the average temperature for the length of the time step is calculated
        % For every date, the day and month are extracted
        dates_all_artificialT_day = dates_all_artificialT.Day;
        dates_all_artificialT_month = dates_all_artificialT.Month;
        dates_all_artificialT_dayAndMonth = [dates_all_artificialT_day' dates_all_artificialT_month'];
        clear dates_all_artificialT_day dates_all_artificialT_month

        % For every spinup and control date, the movmean of the annual T curve is stored
        soilTemperature_spinupAndControl = NaN(size(dates_all_artificialT_dayAndMonth,1),1);

        for i = 1:size(dates_all_artificialT_dayAndMonth,1)

            % The day and month of every timestep is matched with the Tcurve data
            [r, c] = ismember(dates_all_artificialT_dayAndMonth(i,:), historic_annualTemperatureCurve_dayAndMonth, 'rows');
            soilTemperature_spinupAndControl(i,1) = Tcurve_movmean(c, 1);

        end
        
    % If the timestep varies throughout the simulation
    elseif variableTimeStep == 1
        
        % The moving average of the middle dates (1 year) is calculated
        Tcurve_movmean_d1 = movmean(Tcurve_T_long(numel(historic_annualTemperatureCurve_data)+1:numel(historic_annualTemperatureCurve_data)*2),d1);
        Tcurve_movmean_d2 = movmean(Tcurve_T_long(numel(historic_annualTemperatureCurve_data)+1:numel(historic_annualTemperatureCurve_data)*2),d2);

        % The temperature up to 1 Jan 1965 are calculated, based on the timestep
        % This is done using the annual average temperature curve
        dates_all_artificialT = dates_all;

        % For every date, the average temperature for the length of the time step is calculated
        % For every date, the day and month are extracted
        dates_all_artificialT_day = dates_all_artificialT.Day;
        dates_all_artificialT_month = dates_all_artificialT.Month;
        dates_all_artificialT_dayAndMonth = [dates_all_artificialT_day' dates_all_artificialT_month'];
        clear dates_all_artificialT_day dates_all_artificialT_month

        % For every spinup date, the movmean of the annual T curve is stored
        soilTemperature_spinupAndControl = NaN(size(dates_all_artificialT_dayAndMonth,1),1);

        for i = 1:size(dates_all_artificialT_dayAndMonth,1)

            % The day and month of every timestep is matched with the Tcurve data
            [r, c] = ismember(dates_all_artificialT_dayAndMonth(i,:), historic_annualTemperatureCurve_dayAndMonth, 'rows');
            
            % Depending of the dt value of that timestep, a different average temperature has to be stored
            currentdt = dt_array_full(i,1);
            if currentdt == d1
                soilTemperature_spinupAndControl(i,1) = Tcurve_movmean_d1(c, 1);
            elseif currentdt == d2
                soilTemperature_spinupAndControl(i,1) = Tcurve_movmean_d2(c, 1);
            end

        end
        
    end
    
    % -------------------------------------------------------------------
    % The temperature for the spinup and warming treatment are calculated
    % -------------------------------------------------------------------
    
    % The indices of the treatment run are looked for
    ind_warming = find(dates_all >= min(dates_treatmentRun));
    
    % The temperature for the treatment year are increased with dTemp
    soilTemperature_spinupAndWarming = soilTemperature_spinupAndControl;
    % The temperature increase can be a one-time step or a linear increase
    if linearTempIncrease == 0  % If there is a one-time step increase in temperature
        soilTemperature_spinupAndWarming(ind_warming) = soilTemperature_spinupAndControl(ind_warming) + dTemp;
    else                        % If the temperature increases linearly over the treatment run
        relativeIncrease = linspace(0,dTemp,size(ind_warming,2));
        soilTemperature_spinupAndWarming(ind_warming) = soilTemperature_spinupAndControl(ind_warming) + relativeIncrease';
    end       
    
    soilTemperature_warmingOnly = soilTemperature_spinupAndWarming(ind_warming);
    
    clear ind_warming ind_spinup
    
end
   
% Some working variables are removed

clear r c i daily_temp_spinup_movmean dates_for_daily_temperature_spinup dates_for_daily_temperature_treatmentRun...
    historic_annualTemperatureCurve_data historic_annualTemperatureCurve_dayAndMonth ...
    historicT Tcurve_movmean Tcurve_T_long ...
    temperature_historicT temperature_treatmentRun relativeIncrease

%% --------------------------
% Moisture
% --------------------------

%% The data is loaded

% -----------------------------------------------------------
% The measured soil moisture
% -----------------------------------------------------------

% This only needs to be done if  soil moisture is simulated
if simulateSoilMoisture == 1
    
    raw = readtable('Soil moisture control treatment Barre Woods.csv');
    soilMoisture_controlPlots_data = raw.Soil_moisture;

    raw = readtable('Soil moisture heated treatment Barre Woods.csv');
    soilMoisture_warmedPlots_data = raw.Soil_moisture;

    day_historic_annualMoistureCurve = xlsread('Annual curve soil moisture.xlsx', 'Sheet1', 'A2:A367');
    month_historic_annualMoistureCurve = xlsread('Annual curve soil moisture.xlsx', 'Sheet1', 'B2:B367');
    historic_annualMoistureCurve_dayAndMonth = [day_historic_annualMoistureCurve  month_historic_annualMoistureCurve];
    clear day_historic_annualMoistureCurve month_historic_annualMoistureCurve

    historic_annualMoistureCurve_data = xlsread('Annual curve soil moisture.xlsx', 'Sheet1', 'C2:C367');        

end

%% The data is formatted

% This only needs to be done if  soil moisture is simulated
if simulateSoilMoisture == 1
    
    % First, the average annual soil moisture curve is used to construct
    % soil moisture during the spinup run
    
    % The annual moisture curve is replicated 3 times
    Mcurve_M_long = repmat(historic_annualMoistureCurve_data,3,1);

    dates_for_avg_annual_Mcurve = dates_spinup;

    % For every date, the average temperature for the length of the time
    % step is calculated
    % For every date, the day and month are extracted
    spinupDates_day = dates_for_avg_annual_Mcurve.Day;
    spinupDates_month = dates_for_avg_annual_Mcurve.Month;
    spinupDates_dayAndMonth = [spinupDates_day' spinupDates_month'];
    clear spinupDates_day spinupDates_month
    
    % A different procedure is used, depending on if the timestep is
    % constant or variable
    if variableTimeStep == 0
        
        % The moving average of the middle dates (1 year) is calculated
        Mcurve_movmean = movmean(Mcurve_M_long(numel(historic_annualMoistureCurve_data)+1:numel(historic_annualMoistureCurve_data)*2),dt);

        % For every spinup date, the movmean of the annual curve is stored
        temperature_spinup_Tcurve = NaN(size(spinupDates_dayAndMonth,1),1);

        for i = 1:size(spinupDates_dayAndMonth,1)

            % The day and month of every timestep is matched with the Tcurve data
            [r, c] = ismember(spinupDates_dayAndMonth(i,:), historic_annualMoistureCurve_dayAndMonth, 'rows');
            moisture_spinup_Mcurve(i,1) = Mcurve_movmean(c, 1);

        end

        moisture_spinup = moisture_spinup_Mcurve;
    
    elseif variableTimeStep == 1
        
        % The moving average of the middle dates (1 year) is calculated
        Mcurve_movmean1 = movmean(Mcurve_M_long(numel(historic_annualMoistureCurve_data)+1:numel(historic_annualMoistureCurve_data)*2),d1);
        Mcurve_movmean2 = movmean(Mcurve_M_long(numel(historic_annualMoistureCurve_data)+1:numel(historic_annualMoistureCurve_data)*2),d2);
        
        % The indices of dates with large and small timesteps are obtained
        % The number of spinup time steps
        ns = numel(dates_spinup);
        [r_largeTimeStep, ~] = find(dt_array_full(1:ns) == d1);
        [r_smallTimeStep, ~] = find(dt_array_full(1:ns) == d2);
        
        % For every spinup date with a large timestep, the movmean of the annual curve is stored
        moisture_spinup_Mcurve1 = NaN(size(r_largeTimeStep,1),1);

        tmp = spinupDates_dayAndMonth(r_largeTimeStep,:);
        for i = 1:size(r_largeTimeStep,1)

            % The day and month of every timestep is matched with the Tcurve data
            [r, c] = ismember(tmp(i,:), historic_annualMoistureCurve_dayAndMonth, 'rows');
            moisture_spinup_Mcurve1(i,1) = Mcurve_movmean1(c, 1);

        end
        
        % For every spinup date with a small timestep, the movmean of the annual curve is stored
        moisture_spinup_Mcurve2 = NaN(size(r_smallTimeStep,1),1);

        tmp = spinupDates_dayAndMonth(r_smallTimeStep,:);
        for i = 1:size(r_smallTimeStep,1)

            % The day and month of every timestep is matched with the Tcurve data
            [r, c] = ismember(tmp(i,:), historic_annualMoistureCurve_dayAndMonth, 'rows');
            moisture_spinup_Mcurve2(i,1) = Mcurve_movmean2(c, 1);

        end

        moisture_spinup = [moisture_spinup_Mcurve1; moisture_spinup_Mcurve2];
         
    end
    
    % -----------------------------------------------------------
    % The moisture for the control and treatment are constructed
    % -----------------------------------------------------------
    
    if variableTimeStep == 1
        dt = d2;
    end
    
    % If soil moisture measurements are used
    if useRealMoisture == 1
    
        % Moisture for the control run

        % The daily temperature data after 1965 is added
        dates_for_daily_moisture_treatmentRun = dates_treatmentRun;

        % The average temperature per timestep is calculated
        daily_moisture_spinup_movmean = movmean(soilMoisture_controlPlots_data, dt);

        % The dates are matched and the moisture variable is updated
        moisture_controlRun = NaN(size(dates_for_daily_moisture_treatmentRun,2),1);

        for i = 1:size(dates_for_daily_moisture_treatmentRun, 2)

            [r c] = ismember(dates_for_daily_moisture_treatmentRun(1,i), dates_temperature_controlPlots, 'rows');
            moisture_controlRun(i,1) = daily_moisture_spinup_movmean(c,1);

        end

        % Temperature for the warming run  

        % The daily temperature data after 1965 is added
        dates_for_daily_moisture_treatmentRun = dates_treatmentRun;

        % The average temperature per timestep is calculated
        daily_moisture_spinup_movmean = movmean(soilMoisture_warmedPlots_data, dt);

        % The dates are matched and the temperature variable is updated
        moisture_warmingRun = NaN(size(dates_for_daily_moisture_treatmentRun,2),1);

        for i = 1:size(dates_for_daily_moisture_treatmentRun, 2)

            [r c] = ismember(dates_for_daily_moisture_treatmentRun(1,i), dates_temperature_controlPlots, 'rows');
            moisture_warmingRun(i,1) = daily_moisture_spinup_movmean(c,1);

        end  

        % -----------------------------------------------------------
        % Merging temperatures spinup and treatment run
        % -----------------------------------------------------------

        soilMoisture_spinupAndControl = [moisture_spinup; moisture_controlRun];
        soilMoisture_spinupAndWarming = [moisture_spinup; moisture_warmingRun];
        soilMoisture_warmingOnly = moisture_warmingRun;
    
    % If the spin-up soil moisture curve is used for the treatment runs (to serve as an 'artificial' soil moisture curve)
    elseif useRealMoisture == 0
        
        % ----------------------------------------------------
        % The moisture for the control treatment is calculated
        % ----------------------------------------------------

        % The annual moisture curve is replicated 3 times
        Mcurve_M_long = repmat(historic_annualMoistureCurve_data,3,1);

        % The moving average of the middle dates (1 year) is calculated
        Mcurve_movmean = movmean(Mcurve_M_long(numel(historic_annualMoistureCurve_data)+1:numel(historic_annualMoistureCurve_data)*2),dt);

        % The dates for the treatment run are isolated
        dates_treatment_artificialM = dates_treatmentRun;

        % For every date, the average temperature for the length of the time step is calculated
        % For every date, the day and month are extracted
        dates_all_artificialM_day = dates_treatment_artificialM.Day;
        dates_all_artificialM_month = dates_treatment_artificialM.Month;
        dates_all_artificialM_dayAndMonth = [dates_all_artificialM_day' dates_all_artificialM_month'];
        clear dates_all_artificialT_day dates_all_artificialT_month

        % For every treatment date, the movmean of the annual T curve is stored
        moisture_controlRun = NaN(size(dates_all_artificialM_dayAndMonth,1),1);

        for i = 1:size(dates_all_artificialM_dayAndMonth,1)

            % The day and month of every timestep is matched with the Tcurve data
            [r, c] = ismember(dates_all_artificialM_dayAndMonth(i,:), historic_annualMoistureCurve_dayAndMonth, 'rows');
            moisture_controlRun(i,1) = Mcurve_movmean(c, 1);

        end
        
        % ----------------------------------------------------
        % The moisture for the control treatment is calculated
        % ----------------------------------------------------
        
        % If no change is soil moisture is assumed with warming, the
        % moisture of the heated and control runs are equal
        if linearMoistChange == 0
            
            moisture_warmingRun = moisture_controlRun.*dMoist;
            
        % Soil moisture can also change linearly with time
        elseif linearMoistChange == 1 
            
            relativeChange = linspace(1,dMoist,size(moisture_controlRun,1));
            moisture_warmingRun = moisture_controlRun .* relativeChange';
            
        end
        
        % -----------------------------------------------------------
        % Merging temperatures spinup and treatment run
        % -----------------------------------------------------------

        soilMoisture_spinupAndControl = [moisture_spinup; moisture_controlRun];
        soilMoisture_spinupAndWarming = [moisture_spinup; moisture_warmingRun];
        soilMoisture_warmingOnly = moisture_warmingRun;
        
    end
    
end