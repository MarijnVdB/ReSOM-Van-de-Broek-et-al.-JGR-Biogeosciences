% ReSOM script used in the article:

% Van de Broek M., Riley W. J., Tang J., Frey S. D., Schmidt M. W. I.
% Thermal adaptation of enzyme-mediated processes reduces simulated soil CO2 fluxes upon soil warming

% This script is based on earlier Matlab versions of the ReSOM model:
% - Core model by Jinyun Tang, in Tang and Riley, 2015 (https://doi.org/10.1038/nclimate2438)
% - Adaptation by Rose Abramoff, in Sulman et al., 2018 (https://doi.org/10.1007/s10533-018-0509-z)

%% Variables are initialised

% Clear command window, variables and close graphs
clc; clearvars; close all

% Start timer
tic

% Define which site is run
site = 'barreWoods';

% Previous path connections are removed
rmpath([pwd '/Measured Data Barre Woods'])

% The folders containing necessary data are added
addpath('Measured data Barre Woods');

% =================================================
% Define if a normal run or calibration takes place
% =================================================

calibPar = struct();                % A variable to store multiple other variables in
calibPar.calibrationMode = 0;       % 1 if a calibration takes place, 0 if a 'normal run' is ran
calibPar.useCalibrationData = 1;    % 1 if the output from a previous calibration is used, 0 if not

% If the output from a previous calibration is to be used
if calibPar.useCalibrationData == 1
    % The path to the output from the genetic algorithm
%     tmp = load('GA output/Barre Woods - noTau/xga_noAdapt_26-Apr-2024.mat');
    tmp = load('GA output/Barre Woods - noTau - 5 June 2024/xga_noAdapt_16-May-2024.mat');
%     tmp = load('GA output/Barre Woods - noTau/xga_optDriv_29-Apr-2024.mat');
    calibPar.calibValues = tmp.xga;
    clear tmp
end

% The output from a previous calibration cannot be used during a calibration run
if calibPar.calibrationMode == 1 & calibPar.useCalibrationData == 1
    error('Please verify the mode the model is run in!')
end

% ======================================================
% Define which temperature and moisture scenario is used
% ======================================================

% Temperature
useRealTemp = 1;                % If 1, use the measured temperature, if 0, an artificial annual temperature curve is used
useArtificialTemp = 0;          % If 1, the same temperature cycle is used every year, if 0, the measured temperature is used
if useArtificialTemp == 1
    dTemp = 3;                  % The temperature increase in the heated plots, this temperature is 
                                % reached at the end of the simulation period (linear increase)
    linearTempIncrease = 1;     % If 1, the temperature increases linear during the treatment run
end

% Moisture
useRealMoisture = 0;            % If 1, use the measured soil moisture, if 0, an artificial annual moisture curve is used
useArtificialMoisture = 1;      % If 1, the soil moisture in the heated treatment is changed relative to the control treatment, if 0, the measured soil moisture is used
if useArtificialMoisture == 1
    dMoist = 1;                 % The change in soil moisture, relative to control
    linearMoistChange = 0;      % If 1, the moisture changes linear during the treatment run
end

% Carbon inputs for the heated treatment
changeHeatedCInputs = 0;        % If 0, user-defined C inputs for the heated treatment are used,
                                % if 1, the C inputs for the heated treatment are scaled relative to the control
if changeHeatedCInputs == 1
    dCInput = 1.1;              % E.g. if C inputs are increased with 10 % relative to the control, the value should be 1.1
end
    
if useRealTemp == 1 && useArtificialTemp == 1
    error('Choose which temperature scenario is used!')
end

if useRealMoisture == 1 && useArtificialMoisture == 1
    error('Choose which temperature scenario is used!')
end

% ===========================================
% Define which soil moisture scenario is used
% ===========================================

% If 1, soil moisture constraints on respiration are simulated, if 0, not
simulateSoilMoisture = 1;

% If 1, optimal soil moisture is calibrated, if 0, not (only for the no
% thermal adaptation scenario)
calib_moisture = 0;

% ===========================================
% Define if thermal adaptation takes place for MMRT
% ===========================================

thermalAdapt = 0;               % If 1: thermal adaptation is simulated, if 0, not

% Choose which thermal adaptation hypothesis is applied
if thermalAdapt == 0            % If no thermal adaptation is simulated
    
    Topt = 309.65;
        
elseif thermalAdapt == 1        % If thermal adaptation is simulated
    
    enzymeRigidity = 1;         % If 1, the enzyme rigidity scenario is simulated
    optimumDriven = 0;          % If 1, the optimum driven scenario is simulated
    
    if enzymeRigidity == 1      % Parameters for the enzyme rigidity scenario
        alphaT = 0.7;
        betaT = 13 + 273.15;
        alphaC = 0.3;
        betaC = -10;
        Topt_MMRT_spinup = 273.15 + 50;     % The Tref_MMRT used during the first initial years of the spinup
        Cref_MMRT_spinup = -2.98;           % The dCp used during the first initial years of the spinup
        nYears_thermalAdapt = 1;            % The number of years over which the average T is calculated
    elseif optimumDriven == 1   % Parameters for the optimum driven scenario
        alphaT = 0.5;
        betaT = 28 + 173.15;
        Topt_MMRT_spinup = 273.15 + 50;     % The Tref_MMRT used during the first initial years of the spinup
        nYears_thermalAdapt = 1;            % The number of years over which the average T is calculated
    end
    
end

% ===========================================
% Define how C inputs are treated
% ===========================================

% 0 if a constant litter C input for every year is used
% 1 if a series of annual litter inputs is provided and is being used
varyInputsInterAnnually = 1;

% In the case of an 'artificial run' (no real temperature or moisture used)
if useRealTemp == 0 && useRealMoisture == 0
    varyInputsInterAnnually = 0;
end

% 0 if C inputs are linearly distributed over all time steps
% 1 if litter C inputs are concentrated in autumn, and soil C inputs are varied over the growing season
seasonalCarbonInputs = 1;

% 0 if soil C inputs are constant for every year
% 1 if soil C inputs are scaled based on the magnitude of litterfall
scaleSoilCInputs = 0;

% 0 if soil C inputs are not changed upon warming (only for soil, not litter)
% 1 if soil C inputs are changed upon warming (only for soil, not litter)
differentWarmingCarbonInputs = 0;
% These relative changes have to be defined per year in an Excel file (see script Parameter_initialisation.m)
% If 1, the relative change in C inputs is defined
% < 1 for a decrease, > 1 for an increase
if differentWarmingCarbonInputs == 0
    disp('Note: the same soil C inputs are used for the control and heated treatments.')
elseif differentWarmingCarbonInputs == 1
    disp(['Note: the soil C inputs for the warming run are changed relative to the control.'])
end

% A check if separate litter C inputs for the warming treatment are provided in the .csv file
warning('off')
raw = readtable('Litterfall inputs.csv');   % This is the file containing annual litter C inputs
warning('on')
C_input_year = table2array(raw(:,1));       % The years are retrieved
C_input_control = table2array(raw(:,2));    % The annual litter C inputs for the control are retrieved
if size(raw,2) == 3                         % If 3 columns are provided in the .csv file, this means separate litter C inputs for the heated treatment are provided
    C_input_warming = table2array(raw(:,3));
end

% A message is displayed based on the information provided in 'Litterfall inputs.csv'
if size(raw,2) == 2
    disp('Note: the same litter C inputs are used for the control and heated treatments.')
elseif size(raw,2) == 3
    disp('Note: separate litter C inputs are provided for the control and heated treatments.')
else
    disp('Something is wrong with the litter input file')
end
clear dummy

if changeHeatedCInputs == 1 && differentWarmingCarbonInputs == 1
    error('You cannot both change heated soil C inputs with a fixed fraction (changeHeatedCInputs == 1) and change soil C inputs based on the size of litter C inputs (differentWarmingCarbonInputs == 1) at the same time!')
end

% ===========================================
% The calendar dates are defined
% ===========================================

% The number of treatment years is defined
maxYearsWithData = 13;      % The number of years for which data is available
treatment_years = 13;       % The number of years the treatment is run

% Spinup for a defined number of years to steady state
date_startWarming = datetime(2003,05,21);                               % The date the warming started (yyyy/MM/dd)
spinup_years = 50;                                                      % Number of spin-up years
date_start_spinup = date_startWarming - calyears(spinup_years);         % The data the spinup is started
endDate_startPlusTen = date_startWarming + calyears(treatment_years);   % The last date of the treatment run
endDate_fullYear = ['31/12/' int2str(endDate_startPlusTen.Year)];       % This date is adjusted so it's the last day of that year
endDate = datetime(endDate_fullYear, 'InputFormat', 'dd/MM/yyyy', 'Format', 'dd/MM/yyyy');
 
% The number of simulated spinup and treatment days is calculated
numberOfSpinupDays = between(date_start_spinup, date_startWarming, 'days');
numberOfSpinupDays = split(numberOfSpinupDays, 'days');

numberOfTreatmentDays = between(date_startWarming, endDate, 'days');
numberOfTreatmentDays = split(numberOfTreatmentDays, 'days');

% The time step can be constant or variable
variableTimeStep = 1;
if variableTimeStep == 0        % In case the same timestep is used throughout the simulation
    dt = 40;
elseif variableTimeStep == 1    % In case the spinup has a different timestep untill x years before the treatment run starts (to increase runtime)
    d1 = 6;                     % The large timestep during most of the spinup (days)
    d2 = 1;                     % The small timestep, used during the x last spinup years and the treatment run (days)
    timeStepBuffer = 30;        % The number of years before the treatment run that the small time step has to be applied [years] (should be > 0)
end

% An array containing the dates of simulated timesteps is created
if variableTimeStep == 0
    if dt > 1
        dates_spinup = date_start_spinup+caldays(dt/2):caldays(dt):date_startWarming;
        dates_spinup.Format = 'dd/MM/yyyy';
        dates_treatmentRun = date_startWarming+caldays(dt/2):caldays(dt):endDate;
        dates_treatmentRun.Format = 'dd/MM/yyyy';
    elseif dt == 1
        dates_spinup = date_start_spinup:caldays(dt):date_startWarming-caldays(dt);
        dates_spinup.Format = 'dd/MM/yyyy';
        dates_treatmentRun = date_startWarming:caldays(dt):endDate;
        dates_treatmentRun.Format = 'dd/MM/yyyy';
    end
elseif variableTimeStep == 1
    datesLargeTimestep = date_start_spinup+caldays(d1/2):caldays(d1):date_startWarming-calyears(timeStepBuffer);
    datesLargeTimestep.Format = 'dd/MM/yyyy';
    datesSmallTimestep = datesLargeTimestep(end) + d1/2:d2:date_startWarming-caldays(d2);
    dates_spinup = [datesLargeTimestep datesSmallTimestep];
    
    dates_treatmentRun = date_startWarming:caldays(d2):endDate;
    dates_treatmentRun.Format = 'dd/MM/yyyy'; 
    
    % An array with the time step lengt per time step is created
    dt_array_full = ones(size(datesLargeTimestep,2),1).*d1;
    dt_array_full = [dt_array_full; ones(size(datesSmallTimestep,2),1).*d2];
    dt_array_full = [dt_array_full; ones(size(dates_treatmentRun,2),1).*d2];
    
    dt_array_treatment = ones(size(dates_treatmentRun,2),1).*d2;
end
    
dates_all = [dates_spinup dates_treatmentRun];  % All simulated dates are combined

% ==============================
% The measured data is loaded
% ==============================

Load_measured_data;

% ==============================================
% The temperature for every time step is defined
% ==============================================

Initialize_temperature_moisture;

%% Cluster properties for calibration
% Only in case the calibration needs to run parallel (if not, comment this section out)
% Should be adjusted depending on the cluster to be used

if calibPar.calibrationMode == 1
   
    % The names of the available cluster profiles are loaded
    allProfiles = parallel.clusterProfiles;
    
    % Decide if the Euler cluster profile is used (if the job is submitted
    % there), or the local cluster on the virtual machine
    if ismember('EulerLSF8h',allProfiles) % The Euler cluster is detected
        
        batch_job = parcluster('EulerLSF8h');   % To use the ETH Euler cluster
        batch_job.SubmitArguments = '-W 72:00 -R "rusage[mem=2700]"';
        pool = parpool(batch_job, 42);
        fprintf('The parallel loop started at... %s\n', datestr(now,'HH:MM:SS.FFF'))
        
    else % The local cluster on a virtual or local machine is used
        
        batch_job = parpool('local',90);       % To run a local job using the cores of the current machine
        fprintf('The parallel loop started at... %s\n', datestr(now,'HH:MM:SS.FFF'))
        
    end
end

%% The model is ran for the spinup and control scenario

% During the spin-up run, x_in (initial pool sizes) are defined in the function

% ===========================================
% The initial pool sizes are defined
% ===========================================
% See parameters vid_litter and vid_soil (in Parameter_initialisation.m)
% for the names of the different pools

% Initial pool sizes for the organic horizon
x_in_litter = [0.36, 0.46, 1.39, 1.79, 638.35, 18.77, 94.55, 1002.33, ...
                575.79, 0.501, 0.36, 0.043, 0.035, 0, 0, ...
                1000, 1000, 0.5, 0.5, 0, 0, 0, ...
                67.028, 10.27, 9.19, 2688.95, 2721.00];
            
% Initial pool sizes for the mineral soil
x_in_soil = [1.19, 0.28, 5895.75, 54.73, 2646.23, 1974.79, 1.81, 0.205, ...
                0, 1000, 1, 0, 0, 0, 0];

% ===========================================
% The temperature data are combined
% ===========================================

inTemperature = struct();
inTemperature.spinupAndControl = soilTemperature_spinupAndControl;
inTemperature.soilTemperature_spinupAndWarming = soilTemperature_spinupAndWarming;
inTemperature.warming = soilTemperature_warmingOnly;

% ===========================================
% The soil moisture data are combined
% ===========================================

inMoisture = struct();
inMoisture.simulateSoilMoisture = simulateSoilMoisture;
if simulateSoilMoisture == 1
    inMoisture.spinupAndControl = soilMoisture_spinupAndControl;
    inMoisture.soilMoisture_spinupAndWarming = soilMoisture_spinupAndWarming;
    inMoisture.warming = soilMoisture_warmingOnly;
    inMoisture.simulateSoilMoisture = simulateSoilMoisture;
end

% ===========================================
% The dates are combined
% ===========================================

inDates = struct();
inDates.dates_all = dates_all;
inDates.dates_spinup = dates_spinup;
inDates.dates_treatmentRun = dates_treatmentRun;
inDates.date_start_spinup = date_start_spinup;
inDates.date_startWarming = date_startWarming;
inDates.endDate = endDate;

% ===========================================
% Parameters for the function (run_interface.m) are grouped
% ===========================================

paramForFunction = struct();
if variableTimeStep == 0
    paramForFunction.dt = dt;
end
paramForFunction.variableTimeStep = variableTimeStep;
if variableTimeStep == 1
    paramForFunction.dt_array_full = dt_array_full;
    paramForFunction.dt_array_treatment = dt_array_treatment;
end
paramForFunction.thermalAdapt = thermalAdapt;
if thermalAdapt == 0
    paramForFunction.Topt = Topt;
elseif thermalAdapt == 1
    paramForFunction.nYears_thermalAdapt = nYears_thermalAdapt;
    paramForFunction.enzymeRigidity = enzymeRigidity;
    paramForFunction.optimumDriven = optimumDriven;
    if optimumDriven == 1 || enzymeRigidity == 1
        paramForFunction.alphaT = alphaT;
        paramForFunction.betaT = betaT;
        paramForFunction.Topt_MMRT_spinup = Topt_MMRT_spinup;
    end
    if enzymeRigidity == 1
        paramForFunction.alphaC = alphaC;
        paramForFunction.betaC = betaC;
        paramForFunction.Cref_MMRT_spinup = Cref_MMRT_spinup;
    end
end
paramForFunction.varyInputsInterAnnually = varyInputsInterAnnually;
paramForFunction.seasonalCarbonInputs = seasonalCarbonInputs;
paramForFunction.scaleSoilCInputs = scaleSoilCInputs;
paramForFunction.site = site;
paramForFunction.differentWarmingCarbonInputs = differentWarmingCarbonInputs;
paramForFunction.useRealTemp = useRealTemp;
paramForFunction.useRealMoisture = useRealMoisture;
paramForFunction.changeHeatedCInputs = changeHeatedCInputs;
paramForFunction.calib_moisture = calib_moisture;
if changeHeatedCInputs == 1
    paramForFunction.dCInput = dCInput;
end

% ===========================================
%% The model is run
% ===========================================

% Run ReSOM
if calibPar.calibrationMode == 0        % A 'normal run' is run

    % The function run_interface.m is run
    [Cpools_soil_control, Cpools_litter_control, TOUT_control, ...
        Cpools_soil_warmingOnly, Cpools_litter_warmingOnly, TOUT_warmingOnly, CUE] = run_INTERFACE(...
        x_in_soil, x_in_litter, paramForFunction, inTemperature, calibPar, ...
        allMeasurements, inDates, inMoisture);

    % The spinup outputs are added to the warming run outputs
    numberOfSpinupTimesteps = numel(dates_spinup);
    Cpools_litter_warming = [Cpools_litter_control(1:numberOfSpinupTimesteps,:); Cpools_litter_warmingOnly];
    Cpools_soil_warming = [Cpools_soil_control(1:numberOfSpinupTimesteps,:); Cpools_soil_warmingOnly];
    Tout_warming = TOUT_control;
    
elseif calibPar.calibrationMode == 1        % A calibration using the genetic algorithm is run
    
    % The calibration parameters are set
    Set_calibration_parameters;
    
    % The objective function is defined
    objective = @(parameters) run_INTERFACE(x_in_soil, x_in_litter, paramForFunction, ...
        inTemperature, calibPar, allMeasurements, inDates, inMoisture, parameters);
    
    rng default % to reproduce the optimization (same random numbers)
    
    % Different options are defined depending on if an initial population
    % is defined
    if exist('initialPop', 'var')
        gaoptions = gaoptimset('UseParallel',true,'Generations', 30, 'StallGenLimit', 15,'Display','iter', 'PlotFcns',{@gaplotbestf,@gaplotstopping}, ...
                'PopulationSize', 180, 'MutationFcn', {'mutationadaptfeasible'}, 'InitialPopulation', initialPop);
        disp('This is the GA with initial parameter values')
    else
        gaoptions = gaoptimset('UseParallel',true,'Generations', 30, 'StallGenLimit', 15,'Display','iter', 'PlotFcns',{@gaplotbestf,@gaplotstopping}, ...
                'PopulationSize', 180, 'MutationFcn', {'mutationadaptfeasible'});
        disp('This is the GA without initial parameter values')
    end

    % The genetic algorithm is run and the optimal parameter
    % combinations are stored in the variable 'xga'
    if exist('A', 'var')
        xga = ga(objective,numberOfParametersToCalibrate,A,b,[],[],lb,ub,[],gaoptions);
    else
        xga = ga(objective,numberOfParametersToCalibrate,[],[],[],[],lb,ub,[],gaoptions);
    end
end

%% Plotting
if calibPar.calibrationMode == 0
    plotting;
end

%% If a calibration is run, the output is saved

 % The thermal adaptation scenario is stored as a string to provide
% this information to the output variables
if thermalAdapt == 0
    scen = 'noAdapt';
elseif thermalAdapt == 1
    if enzymeRigidity == 1
        scen = 'enzRig';                
    elseif optimumDriven == 1
        scen = 'optDriv';
    end
end

if calibPar.calibrationMode == 1
    
    c = date;
    save(['xga_' c], 'xga')
    save(['xga_' scen '_' c], 'xga')
    fprintf('... And ended at... %s\n', datestr(now,'HH:MM:SS.FFF'))
    
end
%% End timer
toc
