%% Initialisation of the parameters

%% General parameters

Tref = 273.15+8; %MAT             % The reference temperature

% Parameters for MMRT, units are kJ mol-1 K-1
par_MMRT = [-2.17;...   % Cp
            NaN; ...  % dH
            NaN];     % dS        

Tref_MMRT = 299.9; % Median Topt from dataset for soil from Alster et al. (2018) is 305.9. TREF (= T0) is assumed to be 6 K lower than this value
        
if thermalAdapt == 0
    R = 8.314;
    par_MMRT(2) = (Topt*(-par_MMRT(1)*1e3-R) + par_MMRT(1)*1e3*Tref_MMRT)/1e3; % dH is calculated based on Topt (which is optimized)
    par_MMRT(3) = (0.0033 * (par_MMRT(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9))
end

%% Parameters for the organic horizon
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Set number of pools
% -------------------------------------------------------------------------

n_polymers_litter = 2;                  % cellulose-hemicellulose and ligin
n_monomers_litter = 1;                  % dissolved organic carbon (compounds of varying complexity)
n_enzymes_litter = 2;                   % use a generic enzyme to minimize model parameters (BG/Phenox activities reflect corresponding degradation rates) 
n_microbep_litter = 2;                  % microbes can directly uptake DOC 
n_micc_litter = n_microbep_litter;      % Microbial reserve pool
n_surfaces_litter = 1;                  % Mineral surfaces
n_co2_litter = 2;                       % CO2 pools
n_enzymes_ads_litter = n_enzymes_litter*n_surfaces_litter;      % Adsorbed enzymes
n_monomers_ads_litter = n_monomers_litter*n_surfaces_litter;    % Adsorbed monomers
n_all = n_microbep_litter + n_micc_litter + n_enzymes_litter + n_polymers_litter + ...
    n_monomers_litter + n_surfaces_litter + n_co2_litter + n_enzymes_ads_litter + n_monomers_ads_litter;

% -------------------------------------------------------------------------
% Set key parameters
% -------------------------------------------------------------------------

% Inputs
if varyInputsInterAnnually == 0     % If a constant litter input is used for every year
    
    C_input = 178.9;          % Total C inputs [g C/yr]
    C_input_year = [];
    
    % In case of a variable timestep, the full array of timesteps is passed to the function
    if variableTimeStep == 1
        dt = dt_array_full;
    end
    
    % Output is the average C input per time step
    input_polymers = initialize_Cinputs_forestFloor(C_input, C_input_year, dates_spinup, ...
        dates_treatmentRun, dt, seasonalCarbonInputs, date_start_spinup, ...
            endDate, variableTimeStep);
    
elseif varyInputsInterAnnually == 1 % Every year has a different litter input
    
    % In case of a variable timestep, the full array of timesteps is passed to the function
    if variableTimeStep == 1
        dt = dt_array_full;
    end
    
    % The time series of litterfall is loaded
    warning('off')
    raw = readtable('Litterfall inputs.csv');
    warning('on')
    if strcmp(raw.Properties.VariableNames{1}(1), 'x') == 1
        raw.Properties.VariableNames{1} = ...
            strrep(raw.Properties.VariableNames{1}, 'x___', '');
    end
    C_input_year = table2array(raw(:,1));
    C_input_control = table2array(raw(:,2));
    if size(raw,2) == 3
        C_input_warming = table2array(raw(:,3));
    end

    if size(raw,2) == 2
        C_input = C_input_control;
    elseif size(raw,2) == 3
        C_input = [C_input_control C_input_warming];
    end
    
    % Depending on if separate litterfall data is present for the control and
    % warmed run or not, the inputs are processed differently
    % Output is the average C input per time step (output: 1 column if only control inputs
    % are provided, 2 colums if litter inputs for the control and warming
    % treatment are provided)
    [input_polymers, annualLitterCInputs] = initialize_Cinputs_forestFloor(C_input, C_input_year, dates_spinup, ...
        dates_treatmentRun, dt, seasonalCarbonInputs, date_start_spinup, ...
            endDate, variableTimeStep);
    
end

% A variable is created that contains info about if litter inputs are
% provided only for the control treatment, or for both the control and
% heated treatment

[dummy, inputCols] = size(input_polymers);
if inputCols == 1
    heatedLitterInputs = 0;
elseif inputCols == 2
    heatedLitterInputs = 1;
end

fmet = 0.71;                % The fraction of inputs going to metabolic litter (Magill et al. 2004)
if heatedLitterInputs == 0  % Litter inputs are only provided for the control treatment
    input_metabolic_polymers_control = input_polymers.*fmet;
    input_structural_polymers_control = input_polymers.*(1-fmet);
elseif heatedLitterInputs == 1 % Litter inputs are provided for both the control and heated treatment
    input_metabolic_polymers_control = input_polymers(:,1).*fmet;
    input_structural_polymers_control = input_polymers(:,1).*(1-fmet);
    input_metabolic_polymers_heated = input_polymers(:,2).*fmet;
    input_structural_polymers_heated = input_polymers(:,2).*(1-fmet);
end
input_monomers = 0;     % Monomer input [d-1]

% In the case of a run with artificial temperature and moisture, the user
% can also choose to change C inputs upon heating with a fixed portion
% compared to the control run
if useRealTemp == 0 && useRealMoisture == 0
    if changeHeatedCInputs == 1
        heatedLitterInputs = 1;
        input_metabolic_polymers_heated = input_metabolic_polymers_control.*dCInput;
        input_structural_polymers_heated = input_structural_polymers_control.*dCInput;
    end
end
    
% Microbial parameters: r-strategists
micr_maint_rate_Rstrat = 0.0231;                % Microbial maintenance rate [d-1; fraction of microbial biomass]
res_turnover_rate_Rstrat = 0.655;               % Reserve turnover rate (metabolic turnover rate of reserve) [d-1]
max_doc_uptake_Rstrat = 167.56;                 % Maximum doc uptake rate [d-1]
micr_death_rate_Rstrat = 0.059;                 % Microbial death rate [d-1]
max_enz_prod_Rstrat = 0.0013;                   % Maximum enzyme production rate [d-1]
max_micr_growth_Rstrat = 0.8760;                % Maximum microbial growth rate [d-1]
enz_decay_rate_Rstrat = 0.0061;                 % Enzyme decay rate [d-1]
max_poly_degradation_Rstrat = 36.4;             % Maximum som degradation rate [d-1]
Kaff_ee_polymer_Rstrat = 200;                   % Affinity for polymers

% Other parameters
max_enz_ads_Rstrat = 0.0093;                    % Maximum enzyme adsorption rate for r-strategists [d-1]
ads_enz_decay_rate_Rstrat = 0.05;%0.0552;             % Decay rate of adsorbed enzyme for r-strategists [d-1]
max_monomer_ads = 0.0041;                       % Maximum monomer adsorption rate [d-1]
decay_rate_ads_monomer = 1e-4;%5.8439e-4;             % Decay rate of adsorbed monomer [d-1]
enz_prod_analyt_solution = 5d-6*386.294367;     % Enzyme production rate used for analytic solution only [d-1]

% In a ga calibration run, the calibration values are used
if calibPar.calibrationMode == 1
    
    % Parameters for the forest floor
    res_turnover_rate_Rstrat = parameters(1);
    micr_death_rate_Rstrat = parameters(2);
    max_poly_degradation_Rstrat = parameters(3);
    max_enz_prod_Rstrat = parameters(4);
    max_doc_uptake_Rstrat = parameters(5);
    max_monomer_ads = parameters(6);

    % Parameters related to the 'standard' MMRT curve
    if thermalAdapt == 0
            
        par_MMRT(1) = parameters(15); % dC
        Topt = parameters(16); % Topt

        R = 8.314;
        par_MMRT(2) = (Topt*(-par_MMRT(1)*1e3-R) + par_MMRT(1)*1e3*Tref_MMRT)/1e3; % dH is calculated based on Topt (which is optimized)
        par_MMRT(3) = (0.0033 * (par_MMRT(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9)) (dataset Alster et al. 2018, only soil data used)
        
    elseif thermalAdapt == 1
        if enzymeRigidity == 1
            % Nothing has to happen here
        elseif optimumDriven == 1
            par_MMRT(1) = parameters(15); % dC
        end
    end

    % Parameters related to thermal adaptation of the MMRT curve
    if thermalAdapt == 1
        if enzymeRigidity == 1
            alphaT = parameters(15);
            betaT = parameters(16);
            alphaC = parameters(17);
            betaC = parameters(18);
            nYears_thermalAdapt = parameters(19);
        elseif optimumDriven == 1
            alphaT = parameters(16);
            betaT = parameters(17);
            nYears_thermalAdapt = parameters(18);
        end
    end
    
% In case a provided parameter set has to be used
elseif calibPar.useCalibrationData == 1
    
    xga = calibPar.calibValues;
    
    % Parameters for the forest floor
    res_turnover_rate_Rstrat = xga(1);
    micr_death_rate_Rstrat = xga(2);
    max_poly_degradation_Rstrat = xga(3);
    max_enz_prod_Rstrat = xga(4);
    max_doc_uptake_Rstrat = xga(5);
    max_monomer_ads = xga(6);

    if thermalAdapt == 0
        par_MMRT(1) = xga(15);
        Topt = xga(16);
        
        R = 8.314;
        par_MMRT(2) = (Topt*(-par_MMRT(1)*1e3-R) + par_MMRT(1)*1e3*Tref_MMRT)/1e3; % dH is calculated based on Topt (which is optimized)
        par_MMRT(3) = (0.0033 * (par_MMRT(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9)) (dataset Alster et al. 2018, only soil data used)
    elseif thermalAdapt == 1
        if enzymeRigidity == 1
            % Nothing has to happen here
        elseif optimumDriven == 1
            par_MMRT(1) = xga(15);
        end
    end

    if thermalAdapt == 1
        if enzymeRigidity == 1
            alphaT = xga(15);
            betaT = xga(16);
            alphaC = xga(17);
            betaC = xga(18);
            nYears_thermalAdapt = xga(19);
        elseif optimumDriven == 1
            alphaT = xga(16);
            betaT = xga(17);
            nYears_thermalAdapt = xga(18);
        end
    end 
end

% Parameters for the K-strategists in the organic horizon

% Microbial parameters: K-strategists
micr_maint_rate_Kstrat = 0.021;                             % Microbial maintenance rate [d-1; fraction of microbial biomass]
res_turnover_rate_Kstrat = res_turnover_rate_Rstrat;        % Reserve turnover rate (metabolic turnover rate of reserve) [d-1]
max_doc_uptake_Kstrat = max_doc_uptake_Rstrat;              % Maximum doc uptake rate [d-1]
micr_death_rate_Kstrat = 0.045;                             % Microbial death rate [d-1]
max_enz_prod_Kstrat = max_enz_prod_Rstrat;                  % Maximum enzyme production rate [d-1]
max_micr_growth_Kstrat = 0.535;                             % Maximum microbial growth rate [d-1]
enz_decay_rate_Kstrat = enz_decay_rate_Rstrat;              % Enzyme decay rate [d-1]
max_poly_degradation_Kstrat = max_poly_degradation_Rstrat;  % Maximum som degradation rate [d-1]
Kaff_ee_polymer_Kstrat = Kaff_ee_polymer_Rstrat;            % Affinity for polymers

if calibPar.calibrationMode == 1            % In case of a calibration mode
    
    micr_maint_rate_Kstrat = parameters(7);
    micr_death_rate_Kstrat = parameters(8);
    max_micr_growth_Kstrat = parameters(9);
   
elseif calibPar.useCalibrationData == 1    % If a provided parameter set has to be used as inputs
    
    micr_maint_rate_Kstrat = xga(7);
    micr_death_rate_Kstrat = xga(8);
    max_micr_growth_Kstrat = xga(9);
    
end

max_enz_ads_Kstrat = max_enz_ads_Rstrat;                % Maximum enzyme adsorption rate for K-strategists [d-1]
ads_enz_decay_rate_Kstrat = ads_enz_decay_rate_Rstrat;  % Decay rate of adsorbed enzyme for K-strategists [d-1]

% All parameters for the organic horizon combined
cpar_litter=[
    {[input_metabolic_polymers_control input_structural_polymers_control]};    % polymer inputs (metabolic ¦ structural), 1/day
    {input_monomers};                                          % 1 x nmonomers, monomer input, 1/day
    {[micr_maint_rate_Rstrat micr_maint_rate_Kstrat]};         % Microbial maintenance rate [d-1] (r-strat ¦ K-strat)
    {[res_turnover_rate_Rstrat res_turnover_rate_Kstrat]};     % Reserve turnover rate [d-1] (r-strat ¦ K-strat)
    {[max_doc_uptake_Rstrat max_doc_uptake_Kstrat]};           % Maximum doc uptake rate [d-1] (r-strat ¦ K-strat)
    
    {[micr_death_rate_Rstrat micr_death_rate_Kstrat]};         % Microbial death rate [d-1] (r-strat ¦ K-strat)
    {[max_enz_prod_Rstrat max_enz_prod_Kstrat]};               % Maximum enzyme production rate [d-1] (r-strat ¦ K-strat)
    {[max_micr_growth_Rstrat max_micr_growth_Kstrat]};         % Maximum microbial growth rate [d-1] (r-strat ¦ K-strat)
    {[enz_decay_rate_Rstrat enz_decay_rate_Kstrat]};           % Enzyme decay rate [d-1] (r-strat ¦ K-strat)
    {[max_poly_degradation_Rstrat max_poly_degradation_Kstrat]};  % Maximum som degradation rate [d-1]
    
    {[max_enz_ads_Rstrat max_enz_ads_Kstrat]};                 % Maximum enzyme adsorption rate [d-1] (r-strat ¦ K-strat)
    {[ads_enz_decay_rate_Rstrat ads_enz_decay_rate_Kstrat]};   % Decay rate of adsorbed enzyme [d-1] (r-strat ¦ K-strat)
    {[max_monomer_ads]};                                       % Maximum monomer adsorption rate [d-1]
    {[decay_rate_ads_monomer]};                                % Decay rate of adsorbed monomer [d-1]
    {enz_prod_analyt_solution};                                % Enzyme production rate used for analytic solution only [d-1]
    
    {[Kaff_ee_polymer_Rstrat Kaff_ee_polymer_Kstrat]};
    ];

% If separate litter C inputs for the heated treatment are present, these
% inputs are appended to the cpar_litter structure
if heatedLitterInputs == 1
    cpar_litter{end+1} = [input_metabolic_polymers_heated input_structural_polymers_heated];
end

% -------------------------------------------------------------------------
% Set id index - organic horizon
% -------------------------------------------------------------------------

vid_litter.microbep = 1:n_microbep_litter;                                      % Microbial biomass
vid_litter.micc = vid_litter.microbep(end)+(1:n_microbep_litter);               % Reserve biomass
vid_litter.surfaces = vid_litter.micc(end)+(1:n_surfaces_litter);               % Mineral surface sites
vid_litter.monomers = vid_litter.surfaces(end)+(1:n_monomers_litter);           % Monomeric C
vid_litter.monomers_ads = vid_litter.monomers(end)+(1:n_monomers_ads_litter);   % Adsorbed monomers
vid_litter.polymers = vid_litter.monomers_ads(end)+(1:n_polymers_litter);       % Polymeric C
vid_litter.enzymes = vid_litter.polymers(end)+(1:n_enzymes_litter);             % Enzymes
vid_litter.enzymes_ads = vid_litter.enzymes(end)+(1:n_enzymes_ads_litter);      % Adsorbed enzymes
vid_litter.co2 = vid_litter.enzymes_ads(end)+(1:n_co2_litter);                  % Produced CO2
vid_litter.cue = vid_litter.co2(end)+(1:n_microbep_litter);                     % Carbon use efficiency
vid_litter.defactoTurnover = vid_litter.cue(end)+(1:n_polymers_litter);         % Turnover
vid_litter.maintCO2 = vid_litter.defactoTurnover(end) + [1 2];
vid_litter.growthCO2 = vid_litter.maintCO2(end) + [1 2];
vid_litter.enzymeCO2 = vid_litter.growthCO2(end) + [1 2];
vid_litter.uptakeCO2 = vid_litter.enzymeCO2(end) + [1 2];

% -------------------------------------------------------------------------
% Set external input (should change to time series)
% -------------------------------------------------------------------------

% Actual inputs per time step for polymers are added in the loop through all time steps
input_litter.polymers = [];                 % Is defined in run_INTERFACE
input_litter.monomers = cpar_litter{2};

% -------------------------------------------------------------------------
% What used to be globals are now stored in a structure
% -------------------------------------------------------------------------

par_litter = struct();
par_litter.cpar_litter = cpar_litter;                       % a vector of considered calibrating parameters
par_litter.vid_litter = vid_litter;                         % id index

% ---------------
% Microbes
% ---------------

% Parameter for r-strategists in the litter layer
par_mic_litter_rStrat = set_microbe_par_default_litter_rStrat(1, 1, par_litter);

% Parameter for K-strategists in the litter layer
par_mic_litter_kStrat = set_microbe_par_default_litter_kStrat(1, 1, par_litter);

% ---------------
% Enzymes
% ---------------

% Parameters for enzymes from r-strategists in the litter layer
par_enz_litter_rStrat = set_enzyme_par_default_litter_rStrat(1, 1, par_litter);           % enzyme-related parameters

% Some variable values will depend on wether or not soil moisture is simulated
if simulateSoilMoisture == 1
    par_enz_litter_rStrat.Kaff_ee_polymer = 2;
    par_enz_litter_rStrat.Kaff_ee_msurf = 0.5;
elseif simulateSoilMoisture == 0
    par_enz_litter_rStrat.Kaff_ee_polymer = 200;
    par_enz_litter_rStrat.Kaff_ee_msurf = 50;
end

% Parameters for enzymes from K-strategists in the litter layer
par_enz_litter_kStrat = set_enzyme_par_default_litter_kStrat(1, 1, par_litter);           % enzyme-related parameters

% Some variable values will depend on wether or not soil moisture is simulated
if simulateSoilMoisture == 1
    par_enz_litter_kStrat.Kaff_ee_polymer = 2;
    par_enz_litter_kStrat.Kaff_ee_msurf = 0.5;
elseif simulateSoilMoisture == 0
    par_enz_litter_kStrat.Kaff_ee_polymer = 200;
    par_enz_litter_kStrat.Kaff_ee_msurf = 50;
end

par_ss_litter = set_substrate_par_default_litter(par_litter);         % substrate-related parameters

% ---------------
% Surfaces
% ---------------

par_surface_litter = set_msurface_par_default_litter(par_litter);     % mineral surface-related parameters
for i = 2:n_surfaces_litter
    par_surface_litter(i) = set_msurface_par_default_litter(par_litter);
end


% -------------------------------------------------------------------------
% Make a copy of those temperature dependent parameters
% -------------------------------------------------------------------------

par_mic_ref_litter_rStrat = par_mic_litter_rStrat;
par_mic_ref_litter_kStrat = par_mic_litter_kStrat;
par_enz_ref_litter_rStrat = par_enz_litter_rStrat;
par_enz_ref_litter_kStrat = par_enz_litter_kStrat;
par_surface_ref_litter = par_surface_litter;

% -------------------------------------------------------------------------
% Define activation energies
% -------------------------------------------------------------------------

TEa_litter = set_Ea_default_litter(par_litter);
TEa_litter.Ea_vmax_ee = 45000/8.31446; %temp sensitivity of decomp (currently default)

% -------------------------------------------------------------------------
% What used to be globals are now stored in a structure
% -------------------------------------------------------------------------

par_litter.par_mic_litter_rStrat = par_mic_litter_rStrat;   % each par represents a microbe
par_litter.par_mic_litter_kStrat = par_mic_litter_kStrat;   % each par represents a microbe
par_litter.par_enz_litter_rStrat = par_enz_litter_rStrat;   % each par represents a enzyme
par_litter.par_enz_litter_kStrat = par_enz_litter_kStrat;   % each par represents a enzyme
par_litter.par_ss_litter = par_ss_litter;                   % parameter structure for substrate
par_litter.input_litter = input_litter;                     % substrate input structure
par_litter.par_surface_litter = par_surface_litter;         % each par represents a mineral surface

% ------------------------------------------------------------------------
%% Parameters for the soil layer
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Set number of pools
% -------------------------------------------------------------------------

n_polymers_soil = 1;            % cellulose-hemicellulose and ligin
n_monomers_soil = 1;            % dissolved organic carbon (compounds of varying complexity)
n_enzymes_soil = 1;             % use a generic enzyme to minimize model parameters (BG/Phenox activities reflect corresponding degradation rates) 
n_microbep_soil = 1;            % microbes can directly uptake DOC 
n_micc_soil = n_microbep_soil;  % Microbial reserve pool
n_surfaces_soil = 1;            % Mineral surfaces
n_co2_soil = 1;                 % CO2 pool
n_enzymes_ads_soil = n_enzymes_soil*n_surfaces_soil; % consider every surface has all enzymes/monomers (e.g., 3enzymes x 2surfaces network: enz1-sur1,enz1-sur2,enz2-sur1,...,enz3-sur2)
n_monomers_ads_soil = n_monomers_soil*n_surfaces_soil;
n_all = n_microbep_soil + n_micc_soil + n_enzymes_soil + n_polymers_soil + ...
    n_monomers_soil + n_surfaces_soil + n_co2_soil + n_enzymes_ads_soil + n_monomers_ads_soil;

% -------------------------------------------------------------------------
% Set key parameters
% -------------------------------------------------------------------------

% Inputs
if scaleSoilCInputs == 0     % A constant C input is used for every year
    
    C_input = 332;          % Total C inputs [g C/yr]
    annualLitterCInputs = [];
    
    % Output is the average C input per time step
    inputSoil_control = initialize_Cinputs_soil(C_input, annualLitterCInputs,...
            dates_spinup, dates_treatmentRun, dt, date_start_spinup, ...
            date_startWarming, endDate, seasonalCarbonInputs, variableTimeStep);
    
    % If soil C inputs are the same for control and heated plots
    if differentWarmingCarbonInputs == 0
        inputSoil_warming = inputSoil_control;
    % If the soil C inputs for the heated plots are scaled relative to the
    % values for the control plots
    elseif differentWarmingCarbonInputs == 1
        
        % The annual relative change in soil C inputs for the heated plots is defined in an Excel file
        yrsDummy = xlsread('Change in soil C inputs for heated plots.xlsx', 'Sheet1', 'A2:A15');
        annualInputDiff = xlsread('Change in soil C inputs for heated plots.xlsx', 'Sheet1', ['B2:B' int2str(numel(yrsDummy)+1)]);
        clear yrsDummy
        
        % Find the row number for the C input array that are part of the treatment run
        [r c] = find(dates_all >= date_startWarming);
        % The unique years are isolated
        uniqueYearsWarming = unique((dates_all(1,c).Year));
        
        % First, the soil C inputs for the heated treatment are equal to the C inputs for the control treatment
        inputSoil_warming = inputSoil_control;
        
        % The soil C inputs for the heated treatment are calculated, by scaling the soil C inputs for the control treatment
        for i = 1:numel(uniqueYearsWarming)
            % Find the column numbers for the current year are isolated
            [r c] = find(dates_all.Year == uniqueYearsWarming(i));
            inputSoil_warming(c,1) = inputSoil_control(c,1).*annualInputDiff(i); 
        end
        clear uniqueYearsWarming
    end
    
    % In the case of a run with artificial temperature and moisture, the user
    % can also choose to change C inputs upon heating with a fixed portion
    % compared to the control run
    if useRealTemp == 0 && useRealMoisture == 0
        if changeHeatedCInputs == 1
            inputSoil_warming = inputSoil_control.*dCInput;
        end
    end
    
elseif scaleSoilCInputs == 1 % Every year has a different litter input
    
    C_input = 332;          % Total C inputs [g C/yr] (Finzi et al. 2020)
    
    % Output is the average C input per time step
    inputSoil_control = initialize_Cinputs_soil(C_input, annualLitterCInputs,...
            dates_spinup, dates_treatmentRun, dt, date_start_spinup, ...
            date_startWarming, endDate, seasonalCarbonInputs, variableTimeStep);
    
end

fractInputsToPoly = 0.79;                                               % Fraction of inputs to polymers (Finzi et al. 2020, Table 3)
% Soil C inputs for the control treatment
input_polymers_control = inputSoil_control.*fractInputsToPoly;          % Polymer input [d-1]
input_monomers_control = inputSoil_control.*(1-fractInputsToPoly);      % Monomer input [d-1]
% Soil C inputs for the warming treatment
input_polymers_warming = inputSoil_warming.*fractInputsToPoly;          % Polymer input [d-1]
input_monomers_warming = inputSoil_warming.*(1-fractInputsToPoly);      % Monomer input [d-1]

% Microbial parameters
micr_maint_rate = micr_maint_rate_Rstrat;           % Microbial maintenance rate [d-1; fraction of microbial biomass]
res_turnover_rate = 13.502;                         % Reserve turnover rate (metabolic turnover rate of reserve) [d-1]
max_doc_uptake = 304.24;                            % Maximum doc uptake rate [d-1]
micr_death_rate = 0.01;                             % Microbial death rate [d-1]
max_enz_prod = max_enz_prod_Rstrat;                 % Maximum enzyme production rate [d-1]
max_micr_growth = max_micr_growth_Rstrat;           % Maximum microbial growth rate [d-1]
enz_decay_rate = enz_decay_rate_Rstrat;             % Enzyme decay rate [d-1]
max_poly_degradation = 69.31;                       % Maximum som degradation rate [d-1]

max_enz_ads = max_enz_ads_Rstrat;                   % Maximum enzyme adsorption rate [d-1]
ads_enz_decay_rate = ads_enz_decay_rate_Rstrat;     % Decay rate of adsorbed enzyme [d-1]
max_monomer_ads = 0.01;                             % Maximum monomer adsorption rate [d-1]
decay_rate_ads_monomer = 1e-4;%1.6941e-4;                 % Decay rate of adsorbed monomer [d-1]
enz_prod_analyt_solution = 0.0019;                  % Enzyme production rate used for analytic solution only [d-1]

% In a calibration run, the calibration values are used
if calibPar.calibrationMode == 1
    
    res_turnover_rate = parameters(10);
    micr_death_rate = parameters(11);
    max_poly_degradation = parameters(12);
    max_doc_uptake = parameters(13);
    max_monomer_ads = parameters(14);
    
% If a provided set of input parameters is used
elseif calibPar.useCalibrationData == 1
    
    res_turnover_rate = xga(10);
    micr_death_rate = xga(11);
    max_poly_degradation = xga(12);
    max_doc_uptake = xga(13);
    max_monomer_ads = xga(14);

end

% The parameters are combined
cpar0_soil = [0 ...
        0 ...
        micr_maint_rate ...
        res_turnover_rate ...
        max_doc_uptake ...
        micr_death_rate ...
        max_enz_prod ...
        max_micr_growth ...
        enz_decay_rate ...
        max_poly_degradation ...
        max_enz_ads ...
        ads_enz_decay_rate ...
        max_monomer_ads ...
        decay_rate_ads_monomer ...
        enz_prod_analyt_solution]';

% A vector of considered calibrating parameters
cpar_soil=[
    {input_polymers_control};                       % 1 x npolymers, polymer input, 1/day
    {input_monomers_control};                       % 1 x nmonomers, monomer input, 1/day
    {cpar0_soil(3)};                                % 1 x nmicrobes, microbial maintenance rate (1/day)
    {cpar0_soil(4)};                                % 1 x nmicrobes, reserve turnover rate, 1/day
    {repmat(cpar0_soil(5),1,n_monomers_soil)};      % 1 x nmonomers, maximum doc uptake rate                   (1/day)
    
    {cpar0_soil(6)};                                % 1 x nmicrobes, microbial death rate       (1/day)
    {cpar0_soil(7)};                                % 1 x nenzymes, maximum enzyme production rate     (1/day) [inverted from analytic solution]
    {cpar0_soil(8)};                                % 1 x nmicrobes, maximum microbial growth rate (1/day) [inverted from analytic solution]
    {cpar0_soil(9)};                                % 1 x nenzymes, enzyme decay rate          (1/day)
    {repmat(cpar0_soil(10),1,n_polymers_soil)};     % 1 x npolymers, maximum som degradation rate              (1/day)
    
    {repmat(cpar0_soil(11),1,n_surfaces_soil)};     % 1 x nsurfaces, maximum enzyme adsorption rate (1/day)
    {repmat(cpar0_soil(12),1,n_surfaces_soil)};     % 1 x nsurfaces, decay rate of adsorbed enzyme (1/day)
    {repmat(cpar0_soil(13),1,n_monomers_soil)};     % 1 x nmonomers, maximum monomer adsorption rate (1/day)
    {repmat(cpar0_soil(14),1,n_monomers_soil)};     % 1 x nmonomers, decay rate of adsorbed monomer (1/day)
    {cpar0_soil(15)};                               % 1 x nmicrobes, enzyme production rate used for analytic solution only (1/day)
    {input_polymers_warming};
    {input_monomers_warming};
    ];

% -------------------------------------------------------------------------
% Set id index - C cycling
% -------------------------------------------------------------------------

vid_soil.microbep = 1:n_microbep_soil;                                      % Microbial biomass
vid_soil.micc = vid_soil.microbep(end)+(1:n_microbep_soil);                 % Reserve biomass
vid_soil.surfaces = vid_soil.micc(end)+(1:n_surfaces_soil);                 % Mineral surface sites
vid_soil.monomers = vid_soil.surfaces(end)+(1:n_monomers_soil);             % Monomeric C
vid_soil.monomers_ads = vid_soil.monomers(end)+(1:n_monomers_ads_soil);     % Adsorbed monomers
vid_soil.polymers = vid_soil.monomers_ads(end)+(1:n_polymers_soil);         % Polymeric C
vid_soil.enzymes = vid_soil.polymers(end)+(1:n_enzymes_soil);               % Enzymes
vid_soil.enzymes_ads = vid_soil.enzymes(end)+(1:n_enzymes_ads_soil);        % Adsorbed enzymes
vid_soil.co2 = vid_soil.enzymes_ads(end)+(1:n_co2_soil);                    % Produced CO2
vid_soil.cue = vid_soil.co2(end)+(1:n_microbep_soil);                       % Carbon use efficiency
vid_soil.defactoTurnover = vid_soil.cue(end)+(1:n_polymers_soil);           % Turnover
vid_soil.maintCO2 = vid_soil.defactoTurnover(end) + 1;
vid_soil.growthCO2 = vid_soil.maintCO2(end) + 1;
vid_soil.enzymeCO2 = vid_soil.growthCO2(end) + 1;
vid_soil.uptakeCO2 = vid_soil.enzymeCO2(end) + 1;

% -------------------------------------------------------------------------
% Set external input (should change to time series)
% -------------------------------------------------------------------------

% Actual inputs per time step are added in the loop through all time steps
input_soil.polymers = []; 
input_soil.monomers = [];

% -------------------------------------------------------------------------
% What used to be globals are now stored in a structure
% -------------------------------------------------------------------------

par_soil = struct();
par_soil.cpar_soil = cpar_soil;     % a vector of considered calibrating parameters                 
par_soil.vid_soil =  vid_soil;      % id index

% -------------------------------------------------------------------------
% Set key parameters (same for each function group of microbe/surface, 
% might change to consider different groups)
% -------------------------------------------------------------------------

par_mic_soil = set_microbe_par_default_soil(1, 1, par_soil);          % microbe-related parameters
for i = 2:n_microbep_soil
    par_mic_soil(i) = set_microbe_par_default_soil(1, 1, par_soil); 
end

par_enz_soil = set_enzyme_par_default_soil(1, 1, par_soil);           % enzyme-related parameters
for i = 2:n_enzymes_soil
    par_enz_soil(i) = set_enzyme_par_default_soil(1, 1, par_soil);
end

% Some variable values will depend on wether or not soil moisture is simulated
if simulateSoilMoisture == 1
    par_enz_soil.Kaff_ee_polymer = 2;
    par_enz_soil.Kaff_ee_msurf = 0.5;
elseif simulateSoilMoisture == 0
    par_enz_soil.Kaff_ee_polymer = 20d1;
    par_enz_soil.Kaff_ee_msurf = 5d1;
end

par_surface_soil = set_msurface_par_default_soil(par_soil);     % mineral surface-related parameters
for i = 2:n_surfaces_soil
    par_surface_soil(i) = set_msurface_par_default_soil(par_soil);
end

par_ss_soil = set_substrate_par_default_soil(par_soil);         % substrate-related parameters

% -------------------------------------------------------------------------
% Make a copy of those temperature dependent parameters
% -------------------------------------------------------------------------

par_surface_ref_soil = par_surface_soil;
par_mic_ref_soil = par_mic_soil;
par_enz_ref_soil = par_enz_soil;
par_enz_ref_soil.Vmax_ee = cpar0_soil(10);      % max decomp rate 

% -------------------------------------------------------------------------
% Define activation energies
% -------------------------------------------------------------------------

TEa_soil = set_Ea_default_soil(par_soil);
TEa_soil.Ea_vmax_ee = 45000/8.31446; %temp sensitivity of decomp (currently default)

% -------------------------------------------------------------------------
% What used to be globals are now stored in a structure
% -------------------------------------------------------------------------

% par_soil = struct();
par_soil.par_mic_soil = par_mic_soil;               % each par represents a microbe
par_soil.par_enz_soil = par_enz_soil;               % each par represents a enzyme
par_soil.par_surface_soil = par_surface_soil;       % each par represents a mineral surface
par_soil.par_ss_soil = par_ss_soil;                 % parameter structure for substrate
par_soil.input_soil = input_soil;                   % substrate input structure

%% The soil moisture response curve is determined, if this is simulated

if simulateSoilMoisture == 1

    porosity = 0.74; % Soil layer
    
    % The soil moisture value for optimal respiration is calculated
    theta_op = calculate_theta_op();
    theta_op = porosity * 0.65;
        
    % If theta_op is calibrated, some changes are made
    if calibPar.calibrationMode == 1 && calib_moisture == 1
        
        if thermalAdapt == 0
            theta_op = porosity * parameters(17);
        elseif thermalAdapt == 1
            if enzymeRigidity == 1
                theta_op = porosity * parameters(20);
            elseif optimumDriven == 1
                theta_op = porosity * parameters(19);
            end
        end
    
    elseif calibPar.useCalibrationData == 1 && calib_moisture == 1
        
        if thermalAdapt == 0
            theta_op = porosity * xga(17);
        elseif thermalAdapt == 1
            if enzymeRigidity == 1
                theta_op = porosity * xga(20);
            elseif optimumDriven == 1
                theta_op = porosity * xga(19);
            end
        end
        
    end
    
    % The moisture function is calculated
    [s12, fth] = function_fp(theta_op);
    
end





























