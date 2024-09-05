function varargout = run_INTERFACE(x_soil, x_litter, paramForFunction, soilTemperature, ...
    calibPar, allMeasurements, dates, soilMoisture, parameters)

% x_soil = x_in_soil;
% x_litter = x_in_litter;
% paramForFunction = paramForFunction;
% soilTemperature = inTemperature;
% calibPar = calibPar;
% allMeasurements = allMeasurements;
% dates = inDates;
% soilMoisture = inMoisture;

% Relict from a previous version
continueToWarming = 1;

%% The output arguments are defined
if calibPar.calibrationMode == 0
    nOutputs = nargout;
    varargout = cell(1,nOutputs);
elseif calibPar.calibrationMode == 1 
    nOutputs = 1;
    varargout = cell(1,nOutputs);
end

% Variable necessary to stop the run if things go wrong
stop = 0;

%% The temperature series are formatted

temperature_spinupAndControl = soilTemperature.spinupAndControl;
soilTemperature_spinupAndWarming = soilTemperature.soilTemperature_spinupAndWarming;
temperature_warming = soilTemperature.warming;

%% The moisture series are formatted

simulateSoilMoisture = soilMoisture.simulateSoilMoisture;
if simulateSoilMoisture == 1
    soilMoisture_spinupAndControl = soilMoisture.spinupAndControl;
    soilMoisture_spinupAndWarming = soilMoisture.soilMoisture_spinupAndWarming;
    soilMoisture_warmingOnly = soilMoisture.warming;
end

%% The dates are reformatted

dates_all = dates.dates_all;
dates_spinup = dates.dates_spinup;
dates_treatmentRun = dates.dates_treatmentRun;
date_start_spinup = dates.date_start_spinup;
date_startWarming = dates.date_startWarming;
endDate = dates.endDate;

%% The input parameters are retrieved from the structure

variableTimeStep = paramForFunction.variableTimeStep;
if variableTimeStep == 0
    dt = paramForFunction.dt;
end
if variableTimeStep == 1
    dt_array_full = paramForFunction.dt_array_full;
    dt_array_treatment = paramForFunction.dt_array_treatment;
end
thermalAdapt = paramForFunction.thermalAdapt;
if thermalAdapt == 0
    Topt = paramForFunction.Topt;
elseif thermalAdapt == 1
    nYears_thermalAdapt = paramForFunction.nYears_thermalAdapt;
    enzymeRigidity = paramForFunction.enzymeRigidity;
    optimumDriven = paramForFunction.optimumDriven;
    if optimumDriven == 1 || enzymeRigidity == 1
        alphaT = paramForFunction.alphaT;
        betaT = paramForFunction.betaT;
        Topt_MMRT_spinup = paramForFunction.Topt_MMRT_spinup;
    end
    if enzymeRigidity == 1
        alphaC = paramForFunction.alphaC;
        betaC = paramForFunction.betaC;
        Cref_MMRT_spinup = paramForFunction.Cref_MMRT_spinup;
    end
end
varyInputsInterAnnually = paramForFunction.varyInputsInterAnnually;
seasonalCarbonInputs = paramForFunction.seasonalCarbonInputs;
scaleSoilCInputs = paramForFunction.scaleSoilCInputs;
site = paramForFunction.site;
differentWarmingCarbonInputs = paramForFunction.differentWarmingCarbonInputs;
useRealTemp = paramForFunction.useRealTemp;
useRealMoisture = paramForFunction.useRealMoisture;
changeHeatedCInputs = paramForFunction.changeHeatedCInputs;
calib_moisture = paramForFunction.calib_moisture;
if changeHeatedCInputs == 1
    dCInput = paramForFunction.dCInput;
end

% -------------------------------------------------------------------------
%% The spin-up and control run
% -------------------------------------------------------------------------

% The parameters are defined

Parameter_initialisation;

% -------------------------------------------------------------------------
% Set up and run the model
% -------------------------------------------------------------------------

kend = size(temperature_spinupAndControl,1);            % Number of timesteps

TOUT_ctrl = zeros(kend+1,1);                            % Output array for timesteps
Cpools_litter_control = zeros(kend+1,length(x_litter)); % TEMP solution; size of the pools at every timestep
Cpools_litter_control(1,:) = x_litter;                  % The initial pool sizes are stored in x
Cpools_soil_control = zeros(kend+1,length(x_soil));     % TEMP solution; size of the pools at every timestep
Cpools_soil_control(1,:) = x_soil;                      % The initial pool sizes are stored in x

% Some globals to be able to export CUE values out of this function
global CUE_soil_control;
global CUE_litter_rStrat_control;
global CUE_litter_kStrat_control;

CUE_soil_control = NaN(kend,1);
CUE_litter_rStrat_control = NaN(kend,1);
CUE_litter_kStrat_control = NaN(kend,1);

for kk = 1 : kend       % A step for every timestep
    
    if calibPar.calibrationMode == 0
        
        % Display the number of time steps to go
        timestepsToGo = kend-kk;
        if rem(timestepsToGo,100) == 0
            disp(['Still ' num2str(timestepsToGo) ' timesteps to go'])
        end
    
    end
    
    % If a variable time step is used, the timestep is defined
    if variableTimeStep == 1
        dt = dt_array_full(kk);
    end
    
    t = (kk-0.5)*dt;
    
    % -----------------------------------------------------
    % Incorporate temperature effects on parameters
    % -----------------------------------------------------
    
    % =========
    % Litter
    % =========
    
    % Define physiological temperature response curve
    Tref_litter = Tref;          %reference temperature
    
    % The input temperature is converted to kelvin
    temp = 273.15 + temperature_spinupAndControl(kk,1);
        
    fref0_litter = temp/Tref_litter;        % To be used for the calculations of temperature-dependent parameters (for non-equilibrium reactions)
    Tinv_litter = 1/temp - 1/Tref_litter;   % To be used for the calculations of temperature-dependent parameters (for equilibrium reactions)
      
    n_biggest = 1;
    fref_litter = zeros(n_biggest,1);

    if thermalAdapt == 1    % If thermal adaptation takes place
        
        if variableTimeStep == 0
            nTimeSteps = ceil((365*nYears_thermalAdapt)/dt);
        elseif variableTimeStep == 1
            % The timesteps up to this date are isolated
            ts = dt_array_full(1:kk, 1);
            % This array is flipped and the cumulative sum is taken
            ts_rev = flipud(ts);
            ts_cum = cumsum(ts_rev);
            % The index of the value closest to the desired number of days is looked for
            [~, nTimeSteps] = min(abs(ts_cum - 365*nYears_thermalAdapt));
        end
                
        if optimumDriven == 1

            if kk < nTimeSteps+1       % If nYears have not yet passed
                Topt = Topt_MMRT_spinup;
%                 Tref_MMRT = 299.65;
                R = 8.314;
                par_MMRT(2) = (Topt*(-par_MMRT(1)*1e3-R) + par_MMRT(1)*1e3*Tref_MMRT)/1e3; % dH is calculated based on Topt (which is optimized)
                par_MMRT(3) = (0.0033 * (par_MMRT(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9)) (dataset Alster et al. 2018, only soil data used)
            else
                % The weighted mean temperature is calculated (because previous time steps may have a different length)
                temps = temperature_spinupAndControl(kk-nTimeSteps:kk,1);
                dts = ts(kk-nTimeSteps:kk,1);
                meanTemp = sum(temps.*dts)./sum(dts);
                Topt = alphaT*meanTemp + betaT + 273.15;
                
%                 Tref_MMRT = 299.65;
                R = 8.314;
                par_MMRT(2) = (Topt*(-par_MMRT(1)*1e3-R) + par_MMRT(1)*1e3*Tref_MMRT)/1e3; % dH is calculated based on Topt (which is optimized)
                par_MMRT(3) = (0.0033 * (par_MMRT(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9)) (dataset Alster et al. 2018, only soil data used)
            end
            
        elseif enzymeRigidity == 1
            
            if kk < nTimeSteps+1       % If nYears have not yet passed
                Topt = Topt_MMRT_spinup;
                par_MMRT(1) = Cref_MMRT_spinup;
%                 Tref_MMRT = Topt - 6;
                R = 8.314;
                par_MMRT(2) = (Topt*(-par_MMRT(1)*1e3-R) + par_MMRT(1)*1e3*Tref_MMRT)/1e3; % dH is calculated based on Topt (which is optimized)
                par_MMRT(3) = (0.0033 * (par_MMRT(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9)) (dataset Alster et al. 2018, only soil data used)
            else
                temps = temperature_spinupAndControl(kk-nTimeSteps:kk,1);
                dts = ts(kk-nTimeSteps:kk,1);
                meanTemp = sum(temps.*dts)./sum(dts);
                Topt = alphaT*meanTemp + betaT + 273.15;
                par_MMRT(1) = alphaC*meanTemp + betaC;
%                 Tref_MMRT = Topt - 6;
                R = 8.314;
                par_MMRT(2) = (Topt*(-par_MMRT(1)*1e3-R) + par_MMRT(1)*1e3*Tref_MMRT)/1e3; % dH is calculated based on Topt (which is optimized)
                par_MMRT(3) = (0.0033 * (par_MMRT(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9)) (dataset Alster et al. 2018, only soil data used)
            end 
        end
    end
    

    % The rate modifier based on MMRT is calculated
    MMRTfactor = MMRT(par_MMRT, Tref_MMRT, temp, Topt);
    
    % -----------------------------------------------------
    % Temperature-dependent parameter values are calculated
    % -----------------------------------------------------
    
    % ------------------------------
    % Litter microbes: r-strategists
    % ------------------------------

    % equilibrium processes - Arrhenius equation
    par_litter.par_mic_litter_rStrat.Kaff_monomer = par_mic_ref_litter_rStrat.Kaff_monomer.*exp(-TEa_litter.Ea_Kaff_monomer_micb(1,:)*Tinv_litter); % microbial doc affinity (K_BC)
    par_litter.par_mic_litter_rStrat.mr_micb = par_mic_ref_litter_rStrat.mr_micb.*exp(-TEa_litter.Ea_mr_micb(1)*Tinv_litter); % Microbial maintenance respiration

    % non-equilibrium and enzyme processes - Eyring's transition state theory and Murphy's equations
    par_litter.par_mic_litter_rStrat.Vmax_micb = par_mic_ref_litter_rStrat.Vmax_micb*MMRTfactor;    % maxium DOC uptake rate by microbe
    par_litter.par_mic_litter_rStrat.kappa_micb = par_mic_ref_litter_rStrat.kappa_micb*MMRTfactor;  %  Metabolic turnover rate for plastic microbe

    % ------------------------------
    % Litter microbes: K-strategists
    % ------------------------------
    
    % equilibrium processes - Arrhenius equation
    par_litter.par_mic_litter_kStrat.Kaff_monomer = par_mic_ref_litter_kStrat.Kaff_monomer.*exp(-TEa_litter.Ea_Kaff_monomer_micb(1,:)*Tinv_litter); % microbial doc affinity (K_BC)
    par_litter.par_mic_litter_kStrat.mr_micb = par_mic_ref_litter_kStrat.mr_micb.*exp(-TEa_litter.Ea_mr_micb(1)*Tinv_litter); % Microbial maintenance respiration

    % non-equilibrium and enzyme processes - Eyring's transition state theory and Murphy's equations
    par_litter.par_mic_litter_kStrat.Vmax_micb = par_mic_ref_litter_kStrat.Vmax_micb*MMRTfactor; % maxium DOC uptake rate by microbe
    par_litter.par_mic_litter_kStrat.kappa_micb = par_mic_ref_litter_kStrat.kappa_micb*MMRTfactor; %  Metabolic turnover rate for plastic microbe
    
    % ------------------------------
    % Enzymes from r-strategists
    % ------------------------------
    
    % equilibrium processes
    par_litter.par_enz_litter_rStrat.Kaff_ee_polymer = par_enz_ref_litter_rStrat.Kaff_ee_polymer.*exp(-TEa_litter.Ea_Kaff_ee_polymer_micb(1,1)*Tinv_litter); % enzyme affinity to polymer (K_ES)
    par_litter.par_enz_litter_rStrat.Kaff_ee_msurf = par_enz_ref_litter_rStrat.Kaff_ee_msurf.*exp(-TEa_litter.Ea_Kaff_ee_msurf(1,:)*Tinv_litter); % enzyme affinity for adsorptive surface (K_ME)

    % non-equilibrium and enzyme processes
    par_litter.par_enz_litter_rStrat.Vmax_ee = par_enz_ref_litter_rStrat.Vmax_ee*MMRTfactor; % maximum som degradation rate (V_Emax)
    
    % non-equilibrium and non-enzyme processes (mineral adsorption processes)
    par_litter.par_enz_litter_rStrat.Vmax_ads_enzyme = par_enz_ref_litter_rStrat.Vmax_ads_enzyme*fref0_litter.*exp(-TEa_litter.Ea_Vmax_ads_enzyme(1,:)*Tinv_litter); % maximum enzyme adsorption rate (1/day) (??)

    % ------------------------------
    % Enzymes from K-strategists
    % ------------------------------
    
    % equilibrium processes
    par_litter.par_enz_litter_kStrat.Kaff_ee_polymer = par_enz_ref_litter_kStrat.Kaff_ee_polymer.*exp(-TEa_litter.Ea_Kaff_ee_polymer_micb(1,1)*Tinv_litter); % enzyme affinity to polymer (K_ES)
    par_litter.par_enz_litter_kStrat.Kaff_ee_msurf = par_enz_ref_litter_kStrat.Kaff_ee_msurf.*exp(-TEa_litter.Ea_Kaff_ee_msurf(1,:)*Tinv_litter); % enzyme affinity for adsorptive surface (K_ME)

    % non-equilibrium and enzyme processes
    par_litter.par_enz_litter_kStrat.Vmax_ee = par_enz_ref_litter_kStrat.Vmax_ee*MMRTfactor; % maximum som degradation rate (V_Emax)
    
    % non-equilibrium and non-enzyme processes (mineral adsorption processes)
    par_litter.par_enz_litter_kStrat.Vmax_ads_enzyme = par_enz_ref_litter_kStrat.Vmax_ads_enzyme*fref0_litter.*exp(-TEa_litter.Ea_Vmax_ads_enzyme(1,:)*Tinv_litter); % maximum enzyme adsorption rate (1/day) (??)
    
    % ------------------------------
    % Parameters for surfaces
    % ------------------------------
    
    % equilibrium processes
    par_litter.par_surface_litter.Kaff_monomer = par_surface_ref_litter.Kaff_monomer.*exp(-TEa_litter.Ea_Kaff_monomer_msurf(1,:)*Tinv_litter); % adsorption surface doc affinity (K_MC) 

    % non-equilibrium and non-enzyme processes
    par_litter.par_surface_litter.Vmax_ads_monomer = par_surface_ref_litter.Vmax_ads_monomer*fref0_litter.*exp(-TEa_litter.Ea_Vmax_ads_monomer(1,:)*Tinv_litter); % maximum monomer adsorption rate    
    
    par_litter.par_surface_litter.adsmon_decay = par_surface_ref_litter.adsmon_decay*MMRTfactor; % Desorption of monomers
    par_litter.par_enz_litter_rStrat.adsenz_decay = par_enz_ref_litter_rStrat.adsenz_decay * MMRTfactor; % Desorption of enzymes rStrat
    
    par_litter.par_enz_litter_kStrat.adsenz_decay = par_enz_ref_litter_kStrat.adsenz_decay * MMRTfactor; % Desorption of enzymes kStrat
    
    % ------------------------------
    % C inputs for the current timestep are saved in
    % par_litter.input_litter
    % ------------------------------
    
    par_litter.input_litter.polymers = cpar_litter{1}(kk,:);
    
    % =========
    % Soil
    % =========
   
    % Define physiological temperature response curve
    Tref_soil = Tref;	%reference temperature

    temp = 273.15 + temperature_spinupAndControl(kk,1);
        
    fref0_soil = temp/Tref_soil;        % To be used for the calculations of temperature-dependent parameters (for non-equilibrium reactions)
    Tinv_soil = 1/temp - 1/Tref_soil;   % To be used for the calculations of temperature-dependent parameters (for equilibrium reactions)

    if thermalAdapt == 1    % If thermal adaptation takes place
        
        if variableTimeStep == 0
            nTimeSteps = ceil((365*nYears_thermalAdapt)/dt);
        elseif variableTimeStep == 1
            % The timesteps up to this date are isolated
            ts = dt_array_full(1:kk, 1);
            % This array is flipped and the cumulative sum is taken
            ts_rev = flipud(ts);
            ts_cum = cumsum(ts_rev);
            % The index of the value closest to the desired number of days is looked for
            [~, nTimeSteps] = min(abs(ts_cum - 365*nYears_thermalAdapt));
        end
            
        if optimumDriven == 1
            
            if kk < nTimeSteps+1       % If nYears have not yet passed
                Topt = Topt_MMRT_spinup;
                Tref_MMRT = Topt - 6;
                R = 8.314;
                par_MMRT(2) = (Topt*(-par_MMRT(1)*1e3-R) + par_MMRT(1)*1e3*Tref_MMRT)/1e3;
            else
                temps = temperature_spinupAndControl(kk-nTimeSteps:kk,1);
                dts = ts(kk-nTimeSteps:kk,1);
                meanTemp = sum(temps.*dts)./sum(dts);
                Topt = alphaT*meanTemp + betaT + 273.15;
                Tref_MMRT = Topt - 6;
                R = 8.314;
                par_MMRT(2) = (Topt*(-par_MMRT(1)*1e3-R) + par_MMRT(1)*1e3*Tref_MMRT)/1e3;
            end
            
        elseif enzymeRigidity == 1
            
            if kk < nTimeSteps+1       % If nYears have not yet passed
                Topt = Topt_MMRT_spinup;
                Tref_MMRT = Topt - 6;
                R = 8.314;
                par_MMRT(2) = (Topt*(-par_MMRT(1)*1e3-R) + par_MMRT(1)*1e3*Tref_MMRT)/1e3;
                par_MMRT(1) = Cref_MMRT_spinup;
            else
                temps = temperature_spinupAndControl(kk-nTimeSteps:kk,1);
                dts = ts(kk-nTimeSteps:kk,1);
                meanTemp = sum(temps.*dts)./sum(dts);
                Topt = alphaT*meanTemp + betaT + 273.15;
                par_MMRT(1) = alphaC*meanTemp + betaC;
                Tref_MMRT = Topt - 6;
                R = 8.314;
                par_MMRT(2) = (Topt*(-par_MMRT(1)*1e3-R) + par_MMRT(1)*1e3*Tref_MMRT)/1e3;
            end
        end
    end
    
    % The rate modifier based on MMRT is calculated
    MMRTfactor = MMRT(par_MMRT, Tref_MMRT, temp, Topt);
    
    % -----------------------------------------------------
    % Temperature-dependent parameter values are calculated
    % -----------------------------------------------------

    for i=1:n_microbep_soil
        % equilibrium processes - Arrhenius equation
        par_soil.par_mic_soil(i).Kaff_monomer = par_mic_ref_soil(i).Kaff_monomer.*exp(-TEa_soil.Ea_Kaff_monomer_micb(i,:)*Tinv_soil); % microbial doc affinity (K_BC)
        par_soil.par_mic_soil(i).mr_micb = par_mic_ref_soil(i).mr_micb.*exp(-TEa_soil.Ea_mr_micb(i)*Tinv_soil); % Microbial maintenance respiration
        
        % non-equilibrium and enzyme processes - Eyring's transition state theory and Murphy's equations
        par_soil.par_mic_soil(i).Vmax_micb = par_mic_ref_soil(i).Vmax_micb*MMRTfactor; % maxium DOC uptake rate by microbe
        par_soil.par_mic_soil(i).kappa_micb = par_mic_ref_soil(i).kappa_micb*MMRTfactor; %  Metabolic turnover rate for plastic microbe
    end
    
    for i=1:n_enzymes_soil
        % equilibrium processes
        par_soil.par_enz_soil(i).Kaff_ee_polymer = par_enz_ref_soil(i).Kaff_ee_polymer.*exp(-TEa_soil.Ea_Kaff_ee_polymer_micb(i,:)*Tinv_soil); % enzyme affinity to polymer (K_ES)
        par_soil.par_enz_soil(i).Kaff_ee_msurf = par_enz_ref_soil(i).Kaff_ee_msurf.*exp(-TEa_soil.Ea_Kaff_ee_msurf(i,:)*Tinv_soil); % enzyme affinity for adsorptive surface (K_ME)
        
        % non-equilibrium and enzyme processes
        par_soil.par_enz_soil(i).Vmax_ee = par_enz_ref_soil(i).Vmax_ee*MMRTfactor; % maximum som degradation rate (V_Emax)

        % non-equilibrium and non-enzyme processes (mineral adsorption processes)
        par_soil.par_enz_soil(i).Vmax_ads_enzyme = par_enz_ref_soil(i).Vmax_ads_enzyme*fref0_soil.*exp(-TEa_soil.Ea_Vmax_ads_enzyme(i,:)*Tinv_soil); % maximum enzyme adsorption rate (1/day) (??)
    end
    
    for i=1:n_surfaces_soil
        % equilibrium processes
        par_soil.par_surface_soil(i).Kaff_monomer = par_surface_ref_soil(i).Kaff_monomer.*exp(-TEa_soil.Ea_Kaff_monomer_msurf(i,:)*Tinv_soil); % adsorption surface doc affinity (K_MC) 

        % non-equilibrium and non-enzyme processes
        par_soil.par_surface_soil(i).Vmax_ads_monomer = par_surface_ref_soil(i).Vmax_ads_monomer*fref0_soil.*exp(-TEa_soil.Ea_Vmax_ads_monomer(i,:)*Tinv_soil); % maximum monomer adsorption rate
        par_soil.par_surface_soil(i).adsmon_decay = par_surface_ref_soil(i).adsmon_decay*MMRTfactor;
        par_soil.par_enz_soil(i).adsenz_decay = par_enz_ref_soil.adsenz_decay * MMRTfactor;
    end
    
    % ----------------------------------------
    % The C inputs are assigned
    % ----------------------------------------
    
    par_soil.input_soil.polymers = cpar_soil{1}(kk,:);
    par_soil.input_soil.monomers = cpar_soil{2}(kk,:);
    
    % ----------------------------------------
    % Soil moisture constraints are applied
    % ----------------------------------------
    
    if simulateSoilMoisture == 1

        % The water-filled pore space is calculated
        % The volumetric moisture content of this time step
        VM = soilMoisture_spinupAndControl(kk);
        wfps = VM/porosity;
        
        % The rate modifier is calculated
        rm_moisture = interp1(s12, fth, wfps);
        
        % The tortuosity is calculated
%         tort = (VM^(10/3))/(porosity^2);
        
        % ---------------------------------------------------
        % The parameters for the organic horizon are adjusted
        % ---------------------------------------------------
        
        % The maximum DOC uptake by microbes is modified
        % r-strategists
        par_litter.par_mic_litter_rStrat.Vmax_micb = par_litter.par_mic_litter_rStrat.Vmax_micb*rm_moisture; 
        % K-strategists
        par_litter.par_mic_litter_kStrat.Vmax_micb = par_litter.par_mic_litter_kStrat.Vmax_micb*rm_moisture;
        
        % The affinity parameter (K) for depolimerisation is adjusted
        % r-strategists
%         par_litter.par_enz_litter_rStrat.Kaff_ee_polymer = par_litter.par_enz_litter_rStrat(1).Kaff_ee_polymer/tort;
        par_litter.par_enz_litter_rStrat.Kaff_ee_polymer = par_litter.par_enz_litter_rStrat(1).Kaff_ee_polymer;
%         par_litter.par_enz_litter_rStrat.Kaff_ee_msurf = par_litter.par_enz_litter_rStrat.Kaff_ee_msurf/tort;
        par_litter.par_enz_litter_rStrat.Kaff_ee_msurf = par_litter.par_enz_litter_rStrat.Kaff_ee_msurf;
        % K-strategists
%         par_litter.par_enz_litter_kStrat.Kaff_ee_polymer = par_litter.par_enz_litter_kStrat(1).Kaff_ee_polymer/tort;
        par_litter.par_enz_litter_kStrat.Kaff_ee_polymer = par_litter.par_enz_litter_kStrat(1).Kaff_ee_polymer;
%         par_litter.par_enz_litter_kStrat.Kaff_ee_msurf = par_litter.par_enz_litter_kStrat.Kaff_ee_msurf/tort;
        par_litter.par_enz_litter_kStrat.Kaff_ee_msurf = par_litter.par_enz_litter_kStrat.Kaff_ee_msurf;
        
        % ----------------------------------------------
        % The parameters for the soil layer are adjusted
        % ----------------------------------------------
        
        % The maximum DOC uptake by microbes is modified
        par_soil.par_mic_soil.Vmax_micb = par_soil.par_mic_soil.Vmax_micb*rm_moisture;
        
        % The affinity parameter (K) for depolimerisation is adjusted
%         par_soil.par_enz_soil.Kaff_ee_polymer = par_soil.par_enz_soil.Kaff_ee_polymer/tort;
        par_soil.par_enz_soil.Kaff_ee_polymer = par_soil.par_enz_soil.Kaff_ee_polymer;
%         par_soil.par_enz_soil.Kaff_ee_msurf = par_soil.par_enz_soil.Kaff_ee_msurf/tort;
        par_soil.par_enz_soil.Kaff_ee_msurf = par_soil.par_enz_soil.Kaff_ee_msurf;
        
    end

    % -----------------------------------------------------
    % Run core program
    % -----------------------------------------------------

    % The litter C pools are calculated
    Cpools_litter_control(kk+1,:) = adptmbbks1(@ReSOM_litter,Cpools_litter_control(kk,:),TOUT_ctrl(kk)+dt/2,dt,par_litter);
    
    % During calibration, some parameter combinations can lead to no
    % changes in the C pools. In that case, the current run is ended
    if Cpools_litter_control(kk+1,:) == Cpools_litter_control(kk,:)
        stop = 1;
        break
    end
    
    % The soil C pools are calculated
    Cpools_soil_control(kk+1,:) = adptmbbks1(@many_bug_many_substrate,Cpools_soil_control(kk,:),TOUT_ctrl(kk)+dt/2,dt, par_soil);
    
    % During calibration, some parameter combinations can lead to no
    % changes in the C pools. In that case, the current run is ended
    if Cpools_soil_control(kk+1,:) == Cpools_soil_control(kk,:)
        stop = 1;
        break
    end
    
    % The time array is updated
    TOUT_ctrl(kk+1)=TOUT_ctrl(kk)+dt;
    
    % The CUE values are stored
    global value_CUE_soil
    CUE_soil_control(kk) = value_CUE_soil;
    
    global value_CUE_litter_rStrat
    CUE_litter_rStrat_control(kk) = value_CUE_litter_rStrat;
    
    global value_CUE_litter_kStrat
    CUE_litter_kStrat_control(kk) = value_CUE_litter_kStrat;

end

% -------------------------------------------------------------------------
%% The warming run
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Set up and run the model
% -------------------------------------------------------------------------

kend = size(temperature_warming,1);         % Number of timesteps

numberOfSpinupTimesteps = numel(dates_spinup);

TOUT_warming = zeros(kend+1,1);                                                     % Output array for timesteps
Cpools_litter_warming = zeros(kend+1,length(x_litter));                             % TEMP solution; size of the pools at every timestep
Cpools_litter_warming(1,:) = Cpools_litter_control(numberOfSpinupTimesteps+1,:);    % The initial pool sizes are stored in x
Cpools_soil_warming = zeros(kend+1,length(x_soil));                                 % TEMP solution; size of the pools at every timestep
Cpools_soil_warming(1,:) = Cpools_soil_control(numberOfSpinupTimesteps+1,:);        % The initial pool sizes are stored in x

% Some globals to be able to export CUE values out of this function
global CUE_soil_warming;
global CUE_litter_rStrat_warming;
global CUE_litter_kStrat_warming;

CUE_soil_warming = NaN(kend,1);
CUE_litter_rStrat_warming = NaN(kend,1);
CUE_litter_kStrat_warming = NaN(kend,1);

if continueToWarming == 1

for kk = 1 : kend       % A step for every timestep (not day)
    
    % If something went wrong during the spin-up run, the run is stopped
    if stop == 1
        break
    end
    
    % The number of time steps left are displayed
    if calibPar.calibrationMode == 0
        timestepsToGo = kend-kk;
        if rem(timestepsToGo,100) == 0
            disp(['Still ' num2str(timestepsToGo) ' timesteps to go'])
        end
    end
    
    % If a variable timestep is used, the time step is defined
    if variableTimeStep == 1
        dt = dt_array_treatment(kk);
    end
    
    t=(kk-0.5)*dt;
    
    % -----------------------------------------------------
    % Incorporate temperature effects on parameters
    % -----------------------------------------------------
    
    % =========
    % Litter
    % =========
    
    % Define physiological temperature response curve
    Tref_litter = Tref;          %reference temperature

    % The temperature is converted to Kelvin
    temp = 273.15 + temperature_warming(kk,1);
        
    fref0_litter = temp/Tref_litter;        % To be used for the calculations of temperature-dependent parameters (for non-equilibrium reactions)
    Tinv_litter = 1/temp - 1/Tref_litter;   % To be used for the calculations of temperature-dependent parameters (for equilibrium reactions)

    % The first date of the warming treatment is looked for 
    if kk == 1
        firstDate = dates_treatmentRun(1);
        [~, dateCounter] = ismember(firstDate, dates_all);
    end

    if thermalAdapt == 1    % If thermal adaptation takes place
        
        if variableTimeStep == 0
            nTimeSteps = ceil((365*nYears_thermalAdapt)/dt);
        elseif variableTimeStep == 1
            % The timesteps up to this date are isolated
            ts = dt_array_full(1:dateCounter, 1);
            % This array is flipped and the cumulative sum is taken
            ts_rev = flipud(ts);
            ts_cum = cumsum(ts_rev);
            % The index of the value closest to the desired number of days is looked for
            [~, nTimeSteps] = min(abs(ts_cum - 365*nYears_thermalAdapt));
        end

        if optimumDriven == 1
    %         Tref_MMRT = 273.15 +  mean(soilTemperature_spinupAndWarming(dateCounter-nTimeSteps:dateCounter,1)) + deltaTemp;
            temps = soilTemperature_spinupAndWarming(dateCounter-nTimeSteps:dateCounter,1);
            dts = ts(dateCounter-nTimeSteps:dateCounter,1);
            meanTemp = sum(temps.*dts)./sum(dts);
            Topt = alphaT*meanTemp + betaT + 273.15;
            Tref_MMRT = Topt - 6;
            R = 8.314;
            par_MMRT(2) = (Topt*(-par_MMRT(1)*1e3-R) + par_MMRT(1)*1e3*Tref_MMRT)/1e3;
        elseif enzymeRigidity == 1
            temps = soilTemperature_spinupAndWarming(dateCounter-nTimeSteps:dateCounter,1);
            dts = ts(dateCounter-nTimeSteps:dateCounter,1);
            meanTemp = sum(temps.*dts)./sum(dts);
            Topt = alphaT*meanTemp + betaT + 273.15;
            par_MMRT(1) = alphaC*meanTemp + betaC;
            Tref_MMRT = Topt - 6;
            R = 8.314;
            par_MMRT(2) = (Topt*(-par_MMRT(1)*1e3-R) + par_MMRT(1)*1e3*Tref_MMRT)/1e3;
        end
                    
    end

    % The MMRT factor is calculated
    MMRTfactor = MMRT(par_MMRT, Tref_MMRT, temp, Topt);
    
    % -----------------------------------------------------
    % Temperature-dependent parameter values are calculated
    % -----------------------------------------------------
    
    % ------------------------------
    % Litter microbes: r-strategists
    % ------------------------------
    
    % Some variables have to be declared again in the parfor loop
    par_mic_ref_litter_rStrat = par_mic_ref_litter_rStrat;
    TEa_litter = TEa_litter;
    par_mic_ref_litter_kStrat = par_mic_ref_litter_kStrat;
    par_enz_ref_litter_rStrat = par_enz_ref_litter_rStrat;
    par_enz_ref_litter_kStrat = par_enz_ref_litter_kStrat;
    par_surface_litter = par_surface_litter;
    
    par_litter = par_litter;

    % equilibrium processes - Arrhenius equation
    par_litter.par_mic_litter_rStrat.Kaff_monomer = par_mic_ref_litter_rStrat.Kaff_monomer.*exp(-TEa_litter.Ea_Kaff_monomer_micb(1,:)*Tinv_litter); % microbial doc affinity (K_BC)
    par_litter.par_mic_litter_rStrat.mr_micb = par_mic_ref_litter_rStrat.mr_micb.*exp(-TEa_litter.Ea_mr_micb(1)*Tinv_litter); % Microbial maintenance respiration

    % non-equilibrium and enzyme processes - Eyring's transition state theory and Murphy's equations
    par_litter.par_mic_litter_rStrat.Vmax_micb = par_mic_ref_litter_rStrat.Vmax_micb*MMRTfactor; % maxium DOC uptake rate by microbe
    par_litter.par_mic_litter_rStrat.kappa_micb = par_mic_ref_litter_rStrat.kappa_micb*MMRTfactor; %  Metabolic turnover rate for plastic microbe

    % ------------------------------
    % Litter microbes: K-strategists
    % ------------------------------
    
    % equilibrium processes - Arrhenius equation
    par_litter.par_mic_litter_kStrat.Kaff_monomer = par_mic_ref_litter_kStrat.Kaff_monomer.*exp(-TEa_litter.Ea_Kaff_monomer_micb(1,:)*Tinv_litter); % microbial doc affinity (K_BC)
    par_litter.par_mic_litter_kStrat.mr_micb = par_mic_ref_litter_kStrat(1).mr_micb.*exp(-TEa_litter.Ea_mr_micb(1)*Tinv_litter); % Microbial maintenance respiration

    % non-equilibrium and enzyme processes - Eyring's transition state theory and Murphy's equations
    par_litter.par_mic_litter_kStrat.Vmax_micb = par_mic_ref_litter_kStrat(1).Vmax_micb*MMRTfactor; % maxium DOC uptake rate by microbe
    par_litter.par_mic_litter_kStrat.kappa_micb = par_mic_ref_litter_kStrat(1).kappa_micb*MMRTfactor; %  Metabolic turnover rate for plastic microbe
    
    % ------------------------------
    % Enzymes from r-strategists
    % ------------------------------
    
    % equilibrium processes
    par_litter.par_enz_litter_rStrat.Kaff_ee_polymer = par_enz_ref_litter_rStrat.Kaff_ee_polymer.*exp(-TEa_litter.Ea_Kaff_ee_polymer_micb(1,1)*Tinv_litter); % enzyme affinity to polymer (K_ES)
    par_litter.par_enz_litter_rStrat.Kaff_ee_msurf = par_enz_ref_litter_rStrat.Kaff_ee_msurf.*exp(-TEa_litter.Ea_Kaff_ee_msurf(1,:)*Tinv_litter); % enzyme affinity for adsorptive surface (K_ME)

    % non-equilibrium and enzyme processes
    par_litter.par_enz_litter_rStrat.Vmax_ee = par_enz_ref_litter_rStrat.Vmax_ee*MMRTfactor; % maximum som degradation rate (V_Emax)
    
    % non-equilibrium and non-enzyme processes (mineral adsorption processes)
    par_litter.par_enz_litter_rStrat.Vmax_ads_enzyme = par_enz_ref_litter_rStrat.Vmax_ads_enzyme*fref0_litter.*exp(-TEa_litter.Ea_Vmax_ads_enzyme(1,:)*Tinv_litter); % maximum enzyme adsorption rate (1/day) (??)

    % ------------------------------
    % Enzymes from K-strategists
    % ------------------------------
    
    % equilibrium processes
    par_litter.par_enz_litter_kStrat.Kaff_ee_polymer = par_enz_ref_litter_kStrat.Kaff_ee_polymer.*exp(-TEa_litter.Ea_Kaff_ee_polymer_micb(1,1)*Tinv_litter); % enzyme affinity to polymer (K_ES)
    par_litter.par_enz_litter_kStrat.Kaff_ee_msurf = par_enz_ref_litter_kStrat.Kaff_ee_msurf.*exp(-TEa_litter.Ea_Kaff_ee_msurf(1,:)*Tinv_litter); % enzyme affinity for adsorptive surface (K_ME)

    % non-equilibrium and enzyme processes
    par_litter.par_enz_litter_kStrat.Vmax_ee = par_enz_ref_litter_kStrat.Vmax_ee*MMRTfactor; % maximum som degradation rate (V_Emax)
    
    % non-equilibrium and non-enzyme processes (mineral adsorption processes)
    par_litter.par_enz_litter_kStrat.Vmax_ads_enzyme = par_enz_ref_litter_kStrat.Vmax_ads_enzyme*fref0_litter.*exp(-TEa_litter.Ea_Vmax_ads_enzyme(1,:)*Tinv_litter); % maximum enzyme adsorption rate (1/day) (??)
    
    % ------------------------------
    % Parameters for surfaces
    % ------------------------------
    
    % equilibrium processes
    par_litter.par_surface_litter.Kaff_monomer = par_surface_ref_litter.Kaff_monomer.*exp(-TEa_litter.Ea_Kaff_monomer_msurf(1,:)*Tinv_litter); % adsorption surface doc affinity (K_MC) 

    % non-equilibrium and non-enzyme processes
    par_litter.par_surface_litter.Vmax_ads_monomer = par_surface_ref_litter.Vmax_ads_monomer*fref0_litter.*exp(-TEa_litter.Ea_Vmax_ads_monomer(1,:)*Tinv_litter); % maximum monomer adsorption rate    
    
    par_litter.par_surface_litter.adsmon_decay = par_surface_ref_litter.adsmon_decay * MMRTfactor; % Desorption of monomers
    par_litter.par_enz_litter_rStrat.adsenz_decay = par_enz_ref_litter_rStrat.adsenz_decay * MMRTfactor; % Desorption of enzymes rStrat
    
    par_litter.par_enz_litter_kStrat.adsenz_decay = par_enz_ref_litter_kStrat.adsenz_decay * MMRTfactor; % Desorption of enzymes kStrat
    
    % ------------------------------
    % C inputs for the current timestep are saved in
    % par_litter.input_litter
    % ------------------------------
    
    % Depending on whether or not separate litter C inputs for the heated
    % treatment are provided, different litter C inputs are used
    if heatedLitterInputs == 0 % No separate litter inputs for the heated treatment are provided
        par_litter.input_litter.polymers = cpar_litter{1}(dateCounter,:);
    elseif heatedLitterInputs == 1 % Separate litter inputs for the heated treatment are provided
        par_litter.input_litter.polymers = cpar_litter{17}(dateCounter,:);
    end
    
    % =========
    % Soil
    % =========
   
    % Define physiological temperature response curve
    Tref_soil = Tref;          %reference temperature
    
    % The temperatures are converted to kelvin
    temp = 273.15 + temperature_warming(kk,1);
        
    fref0_soil = temp/Tref_soil;        % To be used for the calculations of temperature-dependent parameters (for non-equilibrium reactions)
    Tinv_soil = 1/temp - 1/Tref_soil;   % To be used for the calculations of temperature-dependent parameters (for equilibrium reactions)

    if thermalAdapt == 1    % If thermal adaptation takes place
        
        if variableTimeStep == 0
            nTimeSteps = ceil((365*nYears_thermalAdapt)/dt);
        elseif variableTimeStep == 1
            % The timesteps up to this date are isolated
            ts = dt_array_full(1:dateCounter, 1);
            % This array is flipped and the cumulative sum is taken
            ts_rev = flipud(ts);
            ts_cum = cumsum(ts_rev);
            % The index of the value closest to the desired number of days is looked for
            [~, nTimeSteps] = min(abs(ts_cum - 365*nYears_thermalAdapt));
        end

        if optimumDriven == 1
            temps = soilTemperature_spinupAndWarming(dateCounter-nTimeSteps:dateCounter,1);
            dts = ts(dateCounter-nTimeSteps:dateCounter,1);
            meanTemp = sum(temps.*dts)./sum(dts);
            Topt = alphaT*meanTemp + betaT + 273.15;
            Tref_MMRT = Topt - 6;
            R = 8.314;
            par_MMRT(2) = (Topt*(-par_MMRT(1)*1e3-R) + par_MMRT(1)*1e3*Tref_MMRT)/1e3;
        elseif enzymeRigidity == 1
            temps = soilTemperature_spinupAndWarming(dateCounter-nTimeSteps:dateCounter,1);
            dts = ts(dateCounter-nTimeSteps:dateCounter,1);
            meanTemp = sum(temps.*dts)./sum(dts);
            Topt = alphaT*meanTemp + betaT + 273.15;
            par_MMRT(1) = alphaC*meanTemp + betaC;
            Tref_MMRT = Topt - 6;
            R = 8.314;
            par_MMRT(2) = (Topt*(-par_MMRT(1)*1e3-R) + par_MMRT(1)*1e3*Tref_MMRT)/1e3;
        end
    end

    % The MMRT factor is calculated
    MMRTfactor = MMRT(par_MMRT, Tref_MMRT, temp, Topt);
    
    % -----------------------------------------------------
    % Temperature-dependent parameter values are calculated
    % -----------------------------------------------------
    
    % Some variables have to be declared again in the parfor loop
    par_mic_ref_soil = par_mic_ref_soil;
    TEa_soil = TEa_soil;
    par_enz_ref_soil = par_enz_ref_soil;
    par_surface_ref_soil = par_surface_ref_soil;
    par_soil = par_soil;
    
    for i=1:n_microbep_soil
        % equilibrium processes - Arrhenius equation
        par_soil.par_mic_soil(i).Kaff_monomer = par_mic_ref_soil(i).Kaff_monomer.*exp(-TEa_soil.Ea_Kaff_monomer_micb(i,:)*Tinv_soil); % microbial doc affinity (K_BC)
        par_soil.par_mic_soil(i).mr_micb = par_mic_ref_soil(i).mr_micb.*exp(-TEa_soil.Ea_mr_micb(i)*Tinv_soil); % Microbial maintenance respiration
        
        % non-equilibrium and enzyme processes - Eyring's transition state theory and Murphy's equations
        par_soil.par_mic_soil(i).Vmax_micb = par_mic_ref_soil(i).Vmax_micb*MMRTfactor; % maxium DOC uptake rate by microbe
        par_soil.par_mic_soil(i).kappa_micb = par_mic_ref_soil(i).kappa_micb*MMRTfactor; %  Metabolic turnover rate for plastic microbe
    end
    
    for i=1:n_enzymes_soil
        % equilibrium processes
        par_soil.par_enz_soil(i).Kaff_ee_polymer = par_enz_ref_soil(i).Kaff_ee_polymer.*exp(-TEa_soil.Ea_Kaff_ee_polymer_micb(i,:)*Tinv_soil); % enzyme affinity to polymer (K_ES)
        par_soil.par_enz_soil(i).Kaff_ee_msurf = par_enz_ref_soil(i).Kaff_ee_msurf.*exp(-TEa_soil.Ea_Kaff_ee_msurf(i,:)*Tinv_soil); % enzyme affinity for adsorptive surface (K_ME)
        
        % non-equilibrium and enzyme processes
        par_soil.par_enz_soil(i).Vmax_ee = par_enz_ref_soil(i).Vmax_ee*MMRTfactor; % maximum som degradation rate (V_Emax)

        % non-equilibrium and non-enzyme processes (mineral adsorption processes)
        par_soil.par_enz_soil(i).Vmax_ads_enzyme = par_enz_ref_soil(i).Vmax_ads_enzyme*fref0_soil.*exp(-TEa_soil.Ea_Vmax_ads_enzyme(i,:)*Tinv_soil); % maximum enzyme adsorption rate (1/day) (??)
    end
    
    for i=1:n_surfaces_soil
        % equilibrium processes
        par_soil.par_surface_soil(i).Kaff_monomer = par_surface_ref_soil(i).Kaff_monomer.*exp(-TEa_soil.Ea_Kaff_monomer_msurf(i,:)*Tinv_soil); % adsorption surface doc affinity (K_MC) 

        % non-equilibrium and non-enzyme processes
        par_soil.par_surface_soil(i).Vmax_ads_monomer = par_surface_ref_soil(i).Vmax_ads_monomer*fref0_soil.*exp(-TEa_soil.Ea_Vmax_ads_monomer(i,:)*Tinv_soil); % maximum monomer adsorption rate
%         par_soil.par_surface_soil(i).Vmax_ads_monomer = par_surface_ref_soil(i).Vmax_ads_monomer*MMRTfactor; % maximum monomer adsorption rate
        par_soil.par_surface_soil(i).adsmon_decay = par_surface_ref_soil(i).adsmon_decay*MMRTfactor;
        par_soil.par_enz_soil(i).adsenz_decay = par_enz_ref_soil.adsenz_decay * MMRTfactor;
    end
    
    % ----------------------------------------
    % The C inputs are assigned
    % ----------------------------------------
    
    par_soil.input_soil.polymers = cpar_soil{16}(dateCounter,:);
    par_soil.input_soil.monomers = cpar_soil{17}(dateCounter,:);
    
    % ----------------------------------------
    % Soil moisture constraints are applied
    % ----------------------------------------
    
    if simulateSoilMoisture == 1

        % The water-filled pore space is calculated
        % The volumetric moisture content of this time step
        VM = soilMoisture_warmingOnly(kk);
        wfps = VM/porosity;
        
        % The rate modifier is calculated
        rm_moisture = interp1(s12, fth,wfps);
        
        % The tortuosity is calculated
%         tort = (VM^(10/3))/(porosity^2);
        
        % ------------------------------------------------------
        % The parameters for the forest floor layer are adjusted
        % ------------------------------------------------------
        
        % The maximum DOC uptake by microbes is modified
        % r-strategists
        par_litter.par_mic_litter_rStrat.Vmax_micb = par_litter.par_mic_litter_rStrat.Vmax_micb*rm_moisture; 
        % K-strategists
        par_litter.par_mic_litter_kStrat.Vmax_micb = par_litter.par_mic_litter_kStrat.Vmax_micb*rm_moisture;
        
        % The affinity parameter (K) for depolimerisation is adjusted
        % r-strategists
%         par_litter.par_enz_litter_rStrat.Kaff_ee_polymer = par_litter.par_enz_litter_rStrat(1).Kaff_ee_polymer/tort;
        par_litter.par_enz_litter_rStrat.Kaff_ee_polymer = par_litter.par_enz_litter_rStrat(1).Kaff_ee_polymer;
%         par_litter.par_enz_litter_rStrat.Kaff_ee_msurf = par_litter.par_enz_litter_rStrat.Kaff_ee_msurf/tort;
        par_litter.par_enz_litter_rStrat.Kaff_ee_msurf = par_litter.par_enz_litter_rStrat.Kaff_ee_msurf;
        % K-strategists
%         par_litter.par_enz_litter_kStrat.Kaff_ee_polymer = par_litter.par_enz_litter_kStrat(1).Kaff_ee_polymer/tort;
        par_litter.par_enz_litter_kStrat.Kaff_ee_polymer = par_litter.par_enz_litter_kStrat(1).Kaff_ee_polymer;
%         par_litter.par_enz_litter_kStrat.Kaff_ee_msurf = par_litter.par_enz_litter_kStrat.Kaff_ee_msurf/tort;
        par_litter.par_enz_litter_kStrat.Kaff_ee_msurf = par_litter.par_enz_litter_kStrat.Kaff_ee_msurf;
        
        % ----------------------------------------------
        % The parameters for the soil layer are adjusted
        % ----------------------------------------------
        
        % The maximum DOC uptake by microbes is modified
        par_soil.par_mic_soil.Vmax_micb = par_soil.par_mic_soil(1).Vmax_micb*rm_moisture;   

        % The affinity parameter (K) for depolimerisation is adjusted
%         par_soil.par_enz_soil.Kaff_ee_polymer = par_soil.par_enz_soil(1).Kaff_ee_polymer/tort;
        par_soil.par_enz_soil.Kaff_ee_polymer = par_soil.par_enz_soil(1).Kaff_ee_polymer;
%         par_soil.par_enz_soil(i).Kaff_ee_msurf = par_soil.par_enz_soil(i).Kaff_ee_msurf/tort;
        par_soil.par_enz_soil(i).Kaff_ee_msurf = par_soil.par_enz_soil(i).Kaff_ee_msurf;
        
    end
    
    % The date counter is updated
    dateCounter = dateCounter + 1;

    % -----------------------------------------------------
    % Run core program
    % -----------------------------------------------------
    
    % The litter C pools are calculated
    Cpools_litter_warming(kk+1,:) = adptmbbks1(@ReSOM_litter,Cpools_litter_warming(kk,:),TOUT_warming(kk)+dt/2,dt,par_litter);
       
    % The soil C pools are calculated
    Cpools_soil_warming(kk+1,:) = adptmbbks1(@many_bug_many_substrate,Cpools_soil_warming(kk,:),TOUT_warming(kk)+dt/2,dt, par_soil);
    
    % The time array is updated
    TOUT_warming(kk+1) = TOUT_warming(kk)+dt;

    % The CUE values are stored
    global value_CUE_soil
    CUE_soil_warming(kk) = value_CUE_soil;
    
    global value_CUE_litter_rStrat
    CUE_litter_rStrat_warming(kk) = value_CUE_litter_rStrat;
    
    global value_CUE_litter_kStrat
    CUE_litter_kStrat_warming(kk) = value_CUE_litter_kStrat;

end

%% If a calibration is run, the errors are calculated

% ===================================================
% The annual CO2 fluxes are calculated
% ===================================================

if calibPar.calibrationMode == 1
        
    % ===================================================
    % The annual CO2 fluxes are calculated
    % ===================================================
    
    CO2_measurements = struct();

    CO2_measurements.CO2_rStrat_litter_control = diff(Cpools_litter_control(2:end,14));
    CO2_measurements.CO2_kStrat_litter_control = diff(Cpools_litter_control(2:end,15));
    CO2_measurements.CO2_soil_control = diff(Cpools_soil_control(2:end,9));

    CO2_measurements.CO2_rStrat_litter_warmed = diff(Cpools_litter_warming(2:end,14));
    CO2_measurements.CO2_kStrat_litter_warmed = diff(Cpools_litter_warming(2:end,15));
    CO2_measurements.CO2_soil_warmed = diff(Cpools_soil_warming(2:end,9));

    [annualCO2_control, annualCO2_warming] = calculateAnnualCO2Flux(dates, CO2_measurements, dt);    

    % The error is calculated
    error_annualCO2_control = -nash_sutcliffe(allMeasurements.measuredAnnualCO2_control, annualCO2_control);
    error_annualCO2_control = error_annualCO2_control+1;

    error_annualCO2_warming = -nash_sutcliffe(allMeasurements.measuredAnnualCO2_warming, annualCO2_warming);
    error_annualCO2_warming = error_annualCO2_warming+1;
    
    % ===================================================
    % Errors for the control run
    % ===================================================

    % ---------------------------------------------------
    % The error of CO2 fluxes is calculated
    % ---------------------------------------------------
    
    if variableTimeStep == 0
        
        % For now, CO2 is not being evaluated when the timestep is not variable
        
    elseif variableTimeStep == 1
        
        % The model time steps for which CO2 measurements are available are extracted
        firstTreatmentDate = date_startWarming;
        [r c] = find(dates_all >= firstTreatmentDate);

        % The modelled CO2 fluxes at these dates are extracted
        CO2_cumul_control = Cpools_litter_control(c,14) + Cpools_litter_control(c,15) + Cpools_soil_control(c,9);
        CO2_mod_control = diff(CO2_cumul_control);

        % The indices of the modelled timeseries for which measured data is available are isolated
        % The measurement dates that are modelled are isolated
        tmpMeasDates = allMeasurements.dates_CO2_measurements_control(allMeasurements.dates_CO2_measurements_control < endDate);
        tmpMeas = allMeasurements.Measured_CO2_flux_control(allMeasurements.dates_CO2_measurements_control < endDate);
        
        % To calculate the error for the first 4 years only
%         tmpMeasDates = allMeasurements.dates_CO2_measurements_control(allMeasurements.dates_CO2_measurements_control < datetime(2007,12,31));
%         tmpMeas = allMeasurements.Measured_CO2_flux_control(allMeasurements.dates_CO2_measurements_control < datetime(2007,12,31));
        
        [r c] = find(ismember(dates_treatmentRun, tmpMeasDates));
        % r is adjusted because the derivative is taken, so the indices have shifted
        c = c-1;

        % The error is calculated
        CO2_mod_control = CO2_mod_control(c,1);
        error_timeStepCO2_control = -nash_sutcliffe(tmpMeas, CO2_mod_control);
        error_timeStepCO2_control = error_timeStepCO2_control+1;
        
    end
   
    % ---------------------------------------------------
    % The error of C stocks is calculated
    % ---------------------------------------------------

    % Forest floor carbon
    forestFloor_MAOC_meas_control  = allMeasurements.Cstock_forestFloor_control*0.058;
    forestFloor_metabolicC_meas_control = allMeasurements.Cstock_forestFloor_control*0.942*0.66;
    forestFloor_structC_meas_control  = allMeasurements.Cstock_forestFloor_control*0.942*0.34;

    % Find the timestep closest to the measurement date
    measDate = allMeasurements.date_SOC_measurement;
    c = find(dates_all > measDate, 1, 'first');
    diff1 = between(measDate, dates_all(c),{'days'});
    diff2 = between(dates_all(c-1),measDate,{'days'});

    if caldays(diff1) < caldays(diff2)
        matchedCol = c;
    elseif caldays(diff2) < caldays(diff1)
        matchedCol = c-1;
    else
        matchedCol = c;
    end

    Crows_forestFloor_MAOC = [7 12 13];
    Crows_forestFloor_metabolicC = [8];
    Crows_forestFloor_structC = [9];

    forestFloor_MAOC_mod_control = sum(Cpools_litter_control(matchedCol+1,Crows_forestFloor_MAOC));
    forestFloor_metabolicC_mod_control = sum(Cpools_litter_control(matchedCol+1,Crows_forestFloor_metabolicC));
    forestFloor_structC_mod_control = sum(Cpools_litter_control(matchedCol+1,Crows_forestFloor_structC));

    error_forestFloor_MAOC_control = sqrt(((forestFloor_MAOC_meas_control-forestFloor_MAOC_mod_control)/forestFloor_MAOC_meas_control)^2);
    error_forestFloor_metabolicC_control = sqrt(((forestFloor_metabolicC_meas_control-forestFloor_metabolicC_mod_control)/forestFloor_metabolicC_meas_control)^2);
    error_forestFloor_structC_control = sqrt(((forestFloor_structC_meas_control-forestFloor_structC_mod_control)/forestFloor_structC_meas_control)^2);

    % Soil organic carbon
    soil_C_meas_control = allMeasurements.Cstock_soil_control;

    soil_C_meas_MAOC_control = soil_C_meas_control.*0.58;
    Soil_C_meas_POC_control = soil_C_meas_control.*0.42;
    soil_C_mod_MAOC_control = sum(Cpools_soil_control(c+1,[5 8]));
    Soil_C_mod_POC_control = sum(Cpools_soil_control(c+1,[1 2 4 6 7]));

    % The error is calculated
    error_soil_MAOC_control = sqrt(((soil_C_meas_MAOC_control - soil_C_mod_MAOC_control)/soil_C_meas_MAOC_control)^2);
    error_soil_POC_control = sqrt(((Soil_C_meas_POC_control - Soil_C_mod_POC_control)/Soil_C_meas_POC_control)^2);
            
    % ===================================================
    % Errors for the warming run
    % ===================================================

    % ---------------------------------------------------
    % The error of CO2 fluxes is calculated
    % ---------------------------------------------------
    
    if variableTimeStep == 0
        
        % For now, CO2 is not being evaluated when the timestep is not variable

    elseif variableTimeStep == 1
        
       % The model time steps for which CO2 measurements are available are extracted
        firstTreatmentDate = date_startWarming;
        [r c] = find(dates_treatmentRun >= firstTreatmentDate);

        % The modelled CO2 fluxes at these dates are extracted
        CO2_cumul_warming = Cpools_litter_warming(c,14) + Cpools_litter_warming(c,15) + Cpools_soil_warming(c,9);
        CO2_mod_warming = diff(CO2_cumul_warming);

        % The indices of the modelled timeseries for which measured data is available are isolated
        % The measurement dates that are modelled are isolated
        tmpMeasDates = allMeasurements.dates_CO2_measurements_warming(allMeasurements.dates_CO2_measurements_warming < endDate);
        tmpMeas = allMeasurements.Measured_CO2_flux_warming(allMeasurements.dates_CO2_measurements_warming < endDate);
        
        % To calculate the error for the first 4 years only
%         tmpMeasDates = allMeasurements.dates_CO2_measurements_warming(allMeasurements.dates_CO2_measurements_warming < datetime(2007,12,31));
%         tmpMeas = allMeasurements.Measured_CO2_flux_warming(allMeasurements.dates_CO2_measurements_warming < datetime(2007,12,31));
        
        [r c] = find(ismember(dates_treatmentRun, tmpMeasDates));
        % r is adjusted because the derivative is taken, so the indices have shifted
        c = c-1;
        CO2_mod_warming = CO2_mod_warming(c,1);

        % The error is calculated
        error_timeStepCO2_warming = -nash_sutcliffe(tmpMeas, CO2_mod_warming);
        error_timeStepCO2_warming = error_timeStepCO2_warming+1; 
        
    end

    % ---------------------------------------------------
    % The error of C stocks is calculated
    % ---------------------------------------------------

    % Forest floor carbon
    forestFloor_MAOC_meas_warming  = allMeasurements.Cstock_forestFloor_warming*0.058;
    forestFloor_metabolicC_meas_warming = allMeasurements.Cstock_forestFloor_warming*0.942*0.66;
    forestFloor_structC_meas_warming  = allMeasurements.Cstock_forestFloor_warming*0.942*0.34;
    forestFloor_totalC_meas_warming  = allMeasurements.Cstock_forestFloor_warming;

    % Find the timestep closest to the measurement date
    measDate = allMeasurements.date_SOC_measurement;
    c = find(dates_treatmentRun > measDate, 1, 'first');
    % If c == 1, the measurements are taken on the first simulation day
    if c == 1
        matchedCol = c;
    else    % Otherwise, the closest date is looked for
        diff1 = between(measDate, dates_treatmentRun(c),{'days'});
        diff2 = between(dates_treatmentRun(c-1),measDate,{'days'});

        if caldays(diff1) < caldays(diff2)
            matchedCol = c;
        elseif caldays(diff2) < caldays(diff1)
            matchedCol = c-1;
        else
            matchedCol = c;
        end
    end

    Crows_forestFloor_MAOC = [7 12 13];
    Crows_forestFloor_metabolicC = [8];
    Crows_forestFloor_structC = [9];
    Crows_forestFloor_totalC = [1 2 3 4 6 7 8 9 10 11 12 13];

    forestFloor_MAOC_mod_warming = sum(Cpools_litter_warming(matchedCol+1,Crows_forestFloor_MAOC));
    forestFloor_metabolicC_mod_warming = sum(Cpools_litter_warming(matchedCol+1,Crows_forestFloor_metabolicC));
    forestFloor_structC_mod_warming = sum(Cpools_litter_warming(matchedCol+1,Crows_forestFloor_structC));
    forestFloor_totalC_mod_warming = sum(Cpools_litter_warming(matchedCol+1,Crows_forestFloor_totalC));

    error_forestFloor_MAOC_warming = sqrt(((forestFloor_MAOC_meas_warming-forestFloor_MAOC_mod_warming)/forestFloor_MAOC_meas_warming)^2);
    error_forestFloor_metabolicC_warming = sqrt(((forestFloor_metabolicC_meas_warming-forestFloor_metabolicC_mod_warming)/forestFloor_metabolicC_meas_warming)^2);
    error_forestFloor_structC_warming = sqrt(((forestFloor_structC_meas_warming-forestFloor_structC_mod_warming)/forestFloor_structC_meas_warming)^2);
    error_forestFloor_totalC_warming = sqrt(((forestFloor_totalC_meas_warming-forestFloor_totalC_mod_warming)/forestFloor_totalC_meas_warming)^2);

    % Soil organic carbon
    soil_C_meas_warming = allMeasurements.Cstock_soil_warming;
    Crows_soil = 1:1:8;

    soil_C_meas_MAOC_warming = soil_C_meas_warming.*0.58;
    Soil_C_meas_POC_warming = soil_C_meas_warming.*0.42;

    soil_C_mod_MAOC_warming = sum(Cpools_soil_warming(c+1,[5 8]));
    Soil_C_mod_POC_warming = sum(Cpools_soil_warming(c+1,[1 2 4 6 7]));
    Soil_C_mod_totalC_warming = sum(Cpools_soil_warming(c+1,[1 2 4 5 6 7 8]));

    % The error is calculated
    error_soil_MAOC_warming = sqrt(((soil_C_meas_MAOC_warming - soil_C_mod_MAOC_warming)/soil_C_meas_MAOC_warming)^2);
    error_soil_POC_warming = sqrt(((Soil_C_meas_POC_warming - Soil_C_mod_POC_warming)/Soil_C_meas_POC_warming)^2);
    error_soil_totalC_warming = sqrt(((soil_C_meas_warming - Soil_C_mod_totalC_warming)/soil_C_meas_warming)^2);
        
    % ---------------------------------------------------
    % The error in the difference in CO2 between control and warming
    % ---------------------------------------------------
    
    diffAnnualCO2 = annualCO2_warming - annualCO2_control;
    measuredDiff = allMeasurements.measuredAnnualCO2_warming - allMeasurements.measuredAnnualCO2_control;
    
    % NaNs are removed
    [r c] = find(isnan(measuredDiff));
    measuredDiff(r) = [];
    diffAnnualCO2(r) = [];
    
    error_annualCO2Diff = -nash_sutcliffe(diffAnnualCO2, measuredDiff);
    error_annualCO2Diff = error_annualCO2Diff+1;
    
    
    % ===================================================
    % The total error is calculated
    % ===================================================

    totalError = ... % Errors for control run
            error_timeStepCO2_control*2 + error_forestFloor_MAOC_control + error_forestFloor_metabolicC_control + ...
            error_forestFloor_structC_control + error_soil_MAOC_control + error_soil_POC_control + ...
            ... % Errors for the warming run
            error_timeStepCO2_warming*2;


end

end % End for: if continueToWarming == 1

%% The outputs are formatted, depending if a calibration mode is run or not

if calibPar.calibrationMode == 1            % A calibration mode is run
    
    varargout{1} = totalError;
    
else        % A 'normal run' is run
    
    % The CUE arrays are combined
    CUE = struct();
    CUE.CUE_soil_control = CUE_soil_control;
    CUE.CUE_litter_rStrat_control = CUE_litter_rStrat_control;
    CUE.CUE_litter_kStrat_control = CUE_litter_kStrat_control;
    CUE.CUE_soil_warming = CUE_soil_warming;
    CUE.CUE_litter_rStrat_warming = CUE_litter_rStrat_warming;
    CUE.CUE_litter_kStrat_warming = CUE_litter_kStrat_warming;
 
    
    varargout{1} = Cpools_soil_control;
    varargout{2} = Cpools_litter_control;
    varargout{3} = TOUT_ctrl;
    varargout{4} = Cpools_soil_warming;
    varargout{5} = Cpools_litter_warming;
    varargout{6} = TOUT_warming;
    varargout{7} = CUE;

end
