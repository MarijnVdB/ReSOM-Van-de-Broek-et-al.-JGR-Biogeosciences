% ===========================
% The measured data is loaded
% ===========================

% ----------------------------------------------
% Measured CO2 fluxes
% ----------------------------------------------

warning('off')
raw = readtable('CO2 measurements per day.csv');    % The file containing measured CO2 fluxes
if strcmp(raw.Properties.VariableNames{1}(1:4), 'x___') == 1
    raw.Properties.VariableNames{1} = ...
        strrep(raw.Properties.VariableNames{1}, 'x___', '');
elseif strcmp(raw.Properties.VariableNames{1}(1:2), 'x_') == 1
    raw.Properties.VariableNames{1} = ...
        strrep(raw.Properties.VariableNames{1}, 'x_', '');
end
warning('on')

dates_CO2_measurements_control = datetime(raw.Date_control,'InputFormat','dd/MM/yyyy');
dates_CO2_measurements_control.Format = 'dd/MM/yyyy';

dates_CO2_measurements_warming = datetime(raw.Date_warming,'InputFormat','dd/MM/yyyy');
dates_CO2_measurements_warming.Format = 'dd/MM/yyyy';

% The measured heterotrophic CO2 flux is loaded
Measured_CO2_flux_control = raw.CO2_control;
Measured_CO2_flux_warming = raw.CO2_warming;

% The standard errors
Measured_CO2_flux_control_stdError = raw.stdError_control;
Measured_CO2_flux_warming_stdError = raw.stedError_warming;

clear raw

% If there are measurement at dates earlier than the start of the modelled
% treatment run, these are removed
[r, c] = find(dates_CO2_measurements_control < date_startWarming);
if ~isempty(r)
   dates_CO2_measurements_control(r) = [];
   dates_CO2_measurements_warming(r) = [];
   Measured_CO2_flux_control(r) = [];
   Measured_CO2_flux_warming(r) = [];
   Measured_CO2_flux_control_stdError(r) = [];
   Measured_CO2_flux_warming_stdError(r) = [];
end

% ----------------------------------------------
% Measured SOC data
% ----------------------------------------------

date_SOC_measurement = datetime(2011,11,01);    % This has to be adjusted if SOC has been sampled on a different date

% Organic horizon
Cstock_forestFloor_control = xlsread('SOC data.xlsx', 'Feuil1', 'B5');  
CstockStdError_forestFloor_control = xlsread('SOC data.xlsx', 'Feuil1', 'C5');

Cstock_forestFloor_warming = xlsread('SOC data.xlsx', 'Feuil1', 'B9');  
CstockStdError_forestFloor_warmed = xlsread('SOC data.xlsx', 'Feuil1', 'C9');  

% Soil organic carbon
Cstock_soil_control = xlsread('SOC data.xlsx', 'Feuil1', 'B16');  
CstockStdError_soil_control = xlsread('SOC data.xlsx', 'Feuil1', 'C16');  

Cstock_soil_warming = xlsread('SOC data.xlsx', 'Feuil1', 'B20');  
CstockStdError_soil_warmed = xlsread('SOC data.xlsx', 'Feuil1', 'C20');

% ----------------------------------------------
% Measured annual CO2 flux (optional, not used in this version)
% ----------------------------------------------

years_annualMeasuredCO2 = xlsread('Annual CO2 fluxes.xlsx', 'Sheet1', 'A2:A29');
measuredAnnualCO2_control = xlsread('Annual CO2 fluxes.xlsx', 'Sheet1', 'B2:B29');
measuredAnnualCO2_stDev_control = xlsread('Annual CO2 fluxes.xlsx', 'Sheet1', 'D2:D29');
measuredAnnualCO2_warming = xlsread('Annual CO2 fluxes.xlsx', 'Sheet1', 'C2:C29');
measuredAnnualCO2_stDev_warming = xlsread('Annual CO2 fluxes.xlsx', 'Sheet1', 'E2:E29');

% ----------------------------------------------
% The measured data is formatted for plotting
% ----------------------------------------------

% First check if the last day of treatment run is included in the measured
% data timeseries

% If the end date is before the last years with available data
if max(dates_treatmentRun) <= max(dates_CO2_measurements_control)
    
    datesWithData = dates_treatmentRun; % The dates for which data is available
    
% If the treatment run is longer than the last date for which data is available
else 
    
    [r c] = find(dates_treatmentRun < max(dates_CO2_measurements_control));
    datesWithData = dates_treatmentRun(:,c);
    
end

% The annual CO2 fluxes
simulatedYears = unique(datesWithData.Year);

years_annualMeasuredCO2_tmp = years_annualMeasuredCO2;
years_annualMeasuredCO2 = NaN(numel(simulatedYears), 1);

measuredAnnualCO2_control_tmp = measuredAnnualCO2_control;
measuredAnnualCO2_control = NaN(numel(simulatedYears), 1);

measuredAnnualCO2_stDev_control_tmp = measuredAnnualCO2_stDev_control;
measuredAnnualCO2_stDev_control = NaN(numel(simulatedYears), 1);

measuredAnnualCO2_warming_tmp = measuredAnnualCO2_warming;
measuredAnnualCO2_warming = NaN(numel(simulatedYears), 1);

measuredAnnualCO2_stDev_warming_tmp = measuredAnnualCO2_stDev_warming;
measuredAnnualCO2_stDev_warming = NaN(numel(simulatedYears), 1);

for i = 1:size(simulatedYears, 2)

    [r c] = ismember(simulatedYears(1,i), years_annualMeasuredCO2_tmp, 'rows');
    years_annualMeasuredCO2(i,1) = years_annualMeasuredCO2_tmp(c,1);
    measuredAnnualCO2_control(i,1) = measuredAnnualCO2_control_tmp(c,1);
    measuredAnnualCO2_stDev_control(i,1) = measuredAnnualCO2_stDev_control_tmp(c,1);
    measuredAnnualCO2_warming(i,1) = measuredAnnualCO2_warming_tmp(c,1);
    measuredAnnualCO2_stDev_warming(i,1) = measuredAnnualCO2_stDev_warming_tmp(c,1);    

end

% ----------------------------------------------
% Measurements are combined for calibration
% ----------------------------------------------

allMeasurements = struct();

allMeasurements.dates_treatmentRun = datesWithData;

allMeasurements.dates_CO2_measurements_control = dates_CO2_measurements_control;
allMeasurements.dates_CO2_measurements_warming = dates_CO2_measurements_warming;
allMeasurements.Measured_CO2_flux_control = Measured_CO2_flux_control;
allMeasurements.Measured_CO2_flux_warming = Measured_CO2_flux_warming;

allMeasurements.Cstock_forestFloor_control = Cstock_forestFloor_control;
allMeasurements.Cstock_soil_control = Cstock_soil_control;
allMeasurements.Cstock_forestFloor_warming = Cstock_forestFloor_warming;
allMeasurements.Cstock_soil_warming = Cstock_soil_warming;

allMeasurements.date_SOC_measurement = date_SOC_measurement;

allMeasurements.measuredAnnualCO2_control = measuredAnnualCO2_control;
allMeasurements.measuredAnnualCO2_warming = measuredAnnualCO2_warming;
allMeasurements.measuredAnnualCO2_stDev_control = measuredAnnualCO2_stDev_control;
allMeasurements.measuredAnnualCO2_stDev_warming = measuredAnnualCO2_stDev_warming;