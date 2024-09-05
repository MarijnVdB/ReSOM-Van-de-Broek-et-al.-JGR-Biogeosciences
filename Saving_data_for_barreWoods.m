%% This script is used to save the data for plotting purposes
% To save the model output, run this script after the model has run
%% The folder name is constructed

% The scenario is saved, to be included in the filename
if thermalAdapt == 0
    scenario = 'noAdapt';
elseif thermalAdapt == 1
    if enzymeRigidity == 1
        scenario = 'enzymeRigidity';
    elseif optimumDriven == 1
        scenario = 'optimumDriven';
    end
end

% The folder where the data has to be stored
mainFolder = 'Insert folder';

folder = [mainFolder '\' scenario];

%% The output variable is created

out = struct;

out.dates_all = dates_all;

out.Cpools_litter_control = Cpools_litter_control;
out.Cpools_litter_warming = Cpools_litter_warming;
out.Cpools_soil_control = Cpools_soil_control;
out.Cpools_soil_warming = Cpools_soil_warming;

out.totalCO2_control = totalCO2_control;
out.totalCO2_warmed = totalCO2_warmed;

out.histData = histData;
out.uniqueYears_CO2diff = uniqueYears_CO2diff;
out.diffAnnualCO2 = diffAnnualCO2;
out.diffCumul = diffCumul;
out.dates_after1991 = dates_after1991;

out.soilTemperature_spinupAndControl = soilTemperature_spinupAndControl;
out.soilTemperature_spinupAndWarming = soilTemperature_spinupAndWarming;

out.CO2flux_negRes_meas_control = CO2flux_negRes_meas_control;
out.CO2flux_negres_mod_control = CO2flux_negres_mod_control;
out.CO2flux_posRes_meas_control = CO2flux_posRes_meas_control;
out.CO2flux_posres_mod_control = CO2flux_posres_mod_control;

out.CO2flux_negRes_meas_warming = CO2flux_negRes_meas_warming;
out.CO2flux_negres_mod_warming = CO2flux_negres_mod_warming;
out.CO2flux_posRes_meas_warming = CO2flux_posRes_meas_warming;
out.CO2flux_posres_mod_warming = CO2flux_posres_mod_warming;

out.dt_array_full = dt_array_full;

out.CO2_rStrat_litter_control = CO2_rStrat_litter_control;
out.CO2_kStrat_litter_control = CO2_kStrat_litter_control;
out.CO2_soil_control = CO2_soil_control;
out.totalCO2_control = totalCO2_control;

out.CO2_rStrat_litter_warmed = CO2_rStrat_litter_warmed;
out.CO2_kStrat_litter_warmed = CO2_kStrat_litter_warmed;
out.CO2_soil_warmed = CO2_soil_warmed;
out.totalCO2_warmed = totalCO2_warmed;

out.CUE = CUE;

if thermalAdapt == 0
    out.Tref_MMRT = Tref_MMRT;
    out.par_MMRT = par_MMRT;
elseif thermalAdapt == 1
    if enzymeRigidity == 1
        out.alphaT = alphaT;
        out.betaT = betaT;
        out.alphaC = alphaC;
        out.betaC = betaC;
        out.nYears_thermalAdapt = nYears_thermalAdapt;
        out.par_MMRT = par_MMRT;
    elseif optimumDriven == 1
        out.alphaT = alphaT;
        out.betaT = betaT;
        out.nYears_thermalAdapt = nYears_thermalAdapt;
        out.par_MMRT = par_MMRT;
    elseif thermalBreadth == 1
        out.alphaC = alphaC;
        out.betaC = betaC;
        out.Topt_MMRT = Topt_MMRT;
        out.nYears_thermalAdapt = nYears_thermalAdapt;
        out.par_MMRT = par_MMRT;
    end
end

%% The file is saved

save([scenario '.mat'],'out')
