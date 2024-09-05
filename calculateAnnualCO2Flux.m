function [annualCO2_control, annualCO2_warming] = calculateAnnualCO2Flux(dates, CO2_measurements, dt)

    dates_all = dates.dates_all;
    dates_spinup = dates.dates_spinup;
    
    numberOfSpinupTimesteps = numel(dates_spinup);
%     Cpools_litter_warming = [Cpools_litter_control(2:numberOfSpinupTimesteps+1,:); Cpools_litter_warmingOnly];
%     Cpools_soil_warming = [Cpools_soil_control(2:numberOfSpinupTimesteps+1,:); Cpools_soil_warmingOnly];
    
    % The model results after 1991 are isolated
    firstTreatmentDate = dates.dates_treatmentRun(1);
    [r c] = find(dates_all > firstTreatmentDate);
%     [r c] = find(dates_all > datetime(1990,07,01));
    treatmentDates = dates_all(1,c(1:end-1));
    all_CO2_control = CO2_measurements.CO2_rStrat_litter_control + CO2_measurements.CO2_kStrat_litter_control + CO2_measurements.CO2_soil_control;
    CO2_treatment_control = all_CO2_control(c(1:end-1),1);
    all_CO2_warming = CO2_measurements.CO2_rStrat_litter_warmed + CO2_measurements.CO2_kStrat_litter_warmed + CO2_measurements.CO2_soil_warmed;
    
    all_CO2_warming = [all_CO2_control(2:numberOfSpinupTimesteps+1,:); all_CO2_warming];
    
    
    CO2_treatment_warming = all_CO2_warming(c(1:end-1),1);
    
    % The daily fluxes are calculated
    % The flux per timestep (dt) is converted to the daily flux
    dailyCO2Flux_modDates_control = CO2_treatment_control./dt;
    dailyCO2Flux_modDates_warming = CO2_treatment_warming./dt;
    
    % These daily fluxes are interpolated to a daily timestep
    days_CO2_forPlotting = treatmentDates(1):caldays(1):treatmentDates(end);
    dailyCO2Flux_control = interp1(treatmentDates,dailyCO2Flux_modDates_control,days_CO2_forPlotting);
    dailyCO2Flux_warming = interp1(treatmentDates,dailyCO2Flux_modDates_warming,days_CO2_forPlotting);
    
    % The total annual flux is calculated
    allYears = days_CO2_forPlotting.Year;
    uniqueYears = unique(allYears);
    uniqueYears(uniqueYears<1991) = [];
    
    annualCO2_control = NaN(size(uniqueYears,2),1);
    annualCO2_warming = NaN(size(uniqueYears,2),1);
    
    % A loop through all years of warming
    for i = 1:size(uniqueYears,2)
        
        % Find the rownumbers for the specific year
        yy = uniqueYears(1,i);
        [r c] = find(allYears == yy);
        annualCO2_control(i,1) = sum(dailyCO2Flux_control(1,c));
        annualCO2_warming(i,1) = sum(dailyCO2Flux_warming(1,c));
        
    end



end

