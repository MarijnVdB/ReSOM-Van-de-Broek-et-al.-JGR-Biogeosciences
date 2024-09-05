%% This script is used to plot the results of the 100-year run with linearly increasing soil temperature

clc; clearvars; close all

%%

% --------------------------------------
% Different thermal adaptation scenarios
% --------------------------------------

%% The data is loaded


mainFolder = '/Users/vamarijn/Documents/ETH/Papers/19. ReSOM manuscript/7. Second round of revisions JGR/ReSOM - March 2023/Data artificial runs/100 years - 3 degrees - June2024';

% -------------------
% No adaptation
% -------------------

load([mainFolder '/noAdapt.mat']); noAdaptData = out; clear out

% -------------------
% Optimum driven runs
% -------------------

load([mainFolder '/optimumDriven.mat']); optimumDrivenData = out; clear out

% --------------------
% Enzyme rigidity runs
% --------------------

load([mainFolder '/enzymeRigidity.mat']); enzRigData = out; clear out

clear mainFolder

%% The measured data is loaded

site = 'barreWoods';

rmpath([pwd '/Measured Data Barre Woods'])
addpath('Measured data Barre Woods')

% The number of treatment years is defined
maxYearsWithData = 13;      % The number of years for which data is available
treatment_years = 13;       % The number of years the treatment is run

% spinup for a defined number of years to steady state
continueOn = 0; % 
date_startWarming = datetime(2003,05,21);       % The date the warming started
spinup_years = 30;%100;     % Number of spin-up years
date_start_spinup = date_startWarming - calyears(spinup_years);
endDate_startPlusTen = date_startWarming + calyears(treatment_years);
endDate_fullYear = ['31/12/' int2str(endDate_startPlusTen.Year)];
endDate = datetime(endDate_fullYear, 'InputFormat', 'dd/MM/yyyy', 'Format', 'dd/MM/yyyy');     

d1 = 40;        % The large timestep (days)
d2 = 1;         % The small timestep (days)
timeStepBuffer = 10;  

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

dates_all = [dates_spinup dates_treatmentRun];

Load_measured_data;

%% The OC pools are plotted for the FOREST FLOOR

% The colors are defined
color_monomers_litter = [140,107,177]./255;
color_adsorbedC_litter = [140,107,177]./255;
color_metabolic_litter = [254,196,79]./255;
color_structural_litter = [217,95,14]./255;
color_totalCarbon_litter = 'k';
color_measuredPOC = [153,52,4]./255;

lineWidth = 2;
startWarming = datetime('01/01/2000', 'format', 'dd/MM/yyyy');

% =====================
% The data is formatted
% =====================

% Forest floor - No adaptation - control
adsorbedC_ff_noAdapt_control = noAdaptData.Cpools_litter_control(2:end,7) + noAdaptData.Cpools_litter_control(2:end,12) + noAdaptData.Cpools_litter_control(2:end,13);
metabolic_ff_noAdapt_control = noAdaptData.Cpools_litter_control(2:end,8);
structural_ff_noAdapt_control = noAdaptData.Cpools_litter_control(2:end,9);
enz_rStrat_ff_noAdapt_control = noAdaptData.Cpools_litter_control(2:end,10);
enz_kStrat_ff_noAdapt_control = noAdaptData.Cpools_litter_control(2:end,11);
totalCarbon_ff_noAdapt_control = sum(noAdaptData.Cpools_litter_control(2:end,[1:4 6:13]),2);

% Forest floor - No adaptation - warming
adsorbedC_ff_noAdapt_warming = noAdaptData.Cpools_litter_warming(2:end,7) + noAdaptData.Cpools_litter_warming(2:end,12) + noAdaptData.Cpools_litter_warming(2:end,13);
metabolic_ff_noAdapt_warming = noAdaptData.Cpools_litter_warming(2:end,8);
structural_ff_noAdapt_warming = noAdaptData.Cpools_litter_warming(2:end,9);
enz_rStrat_ff_noAdapt_warming = noAdaptData.Cpools_litter_warming(2:end,10);
enz_kStrat_ff_noAdapt_warming = noAdaptData.Cpools_litter_warming(2:end,11);
totalCarbon_ff_noAdapt_warming = sum(noAdaptData.Cpools_litter_warming(2:end,[1:4 6:13]),2);

% Forest floor - Optimum driven - control
adsorbedC_ff_optimumDriv_control = optimumDrivenData.Cpools_litter_control(2:end,7) + optimumDrivenData.Cpools_litter_control(2:end,12) + optimumDrivenData.Cpools_litter_control(2:end,13);
metabolic_ff_optimumDriv_control = optimumDrivenData.Cpools_litter_control(2:end,8);
structural_ff_optimumDriv_control = optimumDrivenData.Cpools_litter_control(2:end,9);
enz_rStrat_ff_optimumDriv_control = optimumDrivenData.Cpools_litter_control(2:end,10);
enz_kStrat_ff_optimumDriv_control = optimumDrivenData.Cpools_litter_control(2:end,11);
totalCarbon_ff_optimumDriv_control = sum(optimumDrivenData.Cpools_litter_control(2:end,[1:4 6:13]),2);

% Forest floor - Optimum driven - warming
adsorbedC_ff_optimumDriv_warming = optimumDrivenData.Cpools_litter_warming(2:end,7) + optimumDrivenData.Cpools_litter_warming(2:end,12) + optimumDrivenData.Cpools_litter_warming(2:end,13);
metabolic_ff_optimumDriv_warming = optimumDrivenData.Cpools_litter_warming(2:end,8);
structural_ff_optimumDriv_warming = optimumDrivenData.Cpools_litter_warming(2:end,9);
enz_rStrat_ff_optimumDriv_warming = optimumDrivenData.Cpools_litter_warming(2:end,10);
enz_kStrat_ff_optimumDriv_warming = optimumDrivenData.Cpools_litter_warming(2:end,11);
totalCarbon_ff_optimumDriv_warming = sum(optimumDrivenData.Cpools_litter_warming(2:end,[1:4 6:13]),2);

% % Forest floor - Enzyme rigidity - control
adsorbedC_ff_enzRig_control = enzRigData.Cpools_litter_control(2:end,7) + enzRigData.Cpools_litter_control(2:end,12) + enzRigData.Cpools_litter_control(2:end,13);
metabolic_ff_enzRig_control = enzRigData.Cpools_litter_control(2:end,8);
structural_ff_enzRig_control = enzRigData.Cpools_litter_control(2:end,9);
enz_rStrat_ff_enzRig_control = enzRigData.Cpools_litter_control(2:end,10);
enz_kStrat_ff_enzRig_control = enzRigData.Cpools_litter_control(2:end,11);
totalCarbon_ff_enzRig_control = sum(enzRigData.Cpools_litter_control(2:end,[1:4 6:13]),2);

% Forest floor - Enzyme rigidity - warming
adsorbedC_ff_enzRig_warming = enzRigData.Cpools_litter_warming(2:end,7) + enzRigData.Cpools_litter_warming(2:end,12) + enzRigData.Cpools_litter_warming(2:end,13);
metabolic_ff_enzRig_warming = enzRigData.Cpools_litter_warming(2:end,8);
structural_ff_enzRig_warming = enzRigData.Cpools_litter_warming(2:end,9);
enz_rStrat_ff_enzRig_warming = enzRigData.Cpools_litter_warming(2:end,10);
enz_kStrat_ff_enzRig_warming = enzRigData.Cpools_litter_warming(2:end,11);
totalCarbon_ff_enzRig_warming = sum(enzRigData.Cpools_litter_warming(2:end,[1:4 6:13]),2);

% =====================
% The plots are created
% =====================

dotSize = 8;
yMin = 0;
yMax = 2500;

scale = 2;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 12*scale], 'color', [1 1 1])

% --------------------------------------
% Forest floor - no adaptation - control
% --------------------------------------
subplot(3,2,1)
hold on

dates_all = optimumDrivenData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_noAdapt_control,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_noAdapt_control,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_noAdapt_control, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_noAdapt_control, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2105', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - no adaptation')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Forest floor - Optimum driven - control
% --------------------------------------
subplot(3,2,3)
hold on

dates_all = optimumDrivenData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_optimumDriv_control, 'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_optimumDriv_control,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_optimumDriv_control, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_optimumDriv_control, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - optimum driven')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Forest floor - enzyme rigidity - control
% --------------------------------------
subplot(3,2,5)
hold on

dates_all = optimumDrivenData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_enzRig_control,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_enzRig_control,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_enzRig_control, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_enzRig_control, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - enzyme rigidity')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Forest floor - no adaptation - warming
% --------------------------------------
subplot(3,2,2)
hold on

dates_all = optimumDrivenData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_noAdapt_warming,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_noAdapt_warming,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_noAdapt_warming, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_noAdapt_warming, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - no adaptation')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Forest floor - Optimum driven - warming
% --------------------------------------
subplot(3,2,4)
hold on

dates_all = optimumDrivenData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_optimumDriv_warming,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_optimumDriv_warming,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_optimumDriv_warming, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_optimumDriv_warming, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - optimum driven')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Forest floor - enzyme rigidity - warming
% --------------------------------------
subplot(3,2,6)
hold on

dates_all = optimumDrivenData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_enzRig_warming,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_enzRig_warming,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_enzRig_warming, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_enzRig_warming, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - enzyme rigidity')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
text(0.38,0.98,'Forest floor organic carbon stocks','fontweight','bold','fontsize',14)

% Letters
text(0.04, 0.95, '(a)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.65, '(c)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.35, '(e)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.95, '(b)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.65, '(d)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.35, '(f)', 'FontSize', 14, 'FontWeight', 'bold')

% The legend
l = legend([p4, p3, p2, p1], 'Total C', 'Structural C', 'Metabolic C', 'Adsorbed C', 'Orientation', 'Horizontal', 'FontSIze', 14);
pos = l.Position;
pos(1) = 0.3; % x
pos(2) = 0.02; % y
l.Position = pos;
legend('boxoff')

% =========================================================================

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Forest floor carbon stocks.tiff' -r300 -nocrop -transparent
% =========================================================================

%% The OC pools are plotted for the SOIL

% The colors are defined
color_monomers_soil = [140,107,177]./255;
color_adsorbedC_soil = [140,107,177]./255;
color_polymers_soil = [153,52,4]./255;
color_totalCarbon_soil = 'k';

lineWidth = 2;
startWarming = datetime('01/01/2000', 'format', 'dd/MM/yyyy');

% =====================
% The data is formatted
% =====================

% Soil - No adaptation - control
adsorbedC_soil_noAdapt_control = noAdaptData.Cpools_soil_control(2:end,5) + noAdaptData.Cpools_soil_control(2:end,8);
monomers_soil_noAdapt_control = noAdaptData.Cpools_soil_control(2:end,4);
polymers_soil_noAdapt_control = noAdaptData.Cpools_soil_control(2:end,6);
totalCarbon_soil_noAdapt_control = sum(noAdaptData.Cpools_soil_control(2:end,[1 2 4 5 6 7 8]),2);

% Soil - No adaptation - warming
adsorbedC_soil_noAdapt_warming = noAdaptData.Cpools_soil_warming(2:end,5) + noAdaptData.Cpools_soil_warming(2:end,8);
monomers_soil_noAdapt_warming = noAdaptData.Cpools_soil_warming(2:end,4);
polymers_soil_noAdapt_warming = noAdaptData.Cpools_soil_warming(2:end,6);
totalCarbon_soil_noAdapt_warming = sum(noAdaptData.Cpools_soil_warming(2:end,[1 2 4 5 6 7 8]),2);

% Soil - Optimum driven - control
adsorbedC_soil_optimumDriv_control = optimumDrivenData.Cpools_soil_control(2:end,5) + optimumDrivenData.Cpools_soil_control(2:end,8);
monomers_soil_optimumDriv_control = optimumDrivenData.Cpools_soil_control(2:end,4);
polymers_soil_optimumDriv_control = optimumDrivenData.Cpools_soil_control(2:end,6);
totalCarbon_soil_optimumDriv_control = sum(optimumDrivenData.Cpools_soil_control(2:end,[1 2 4 5 6 7 8]),2);

% Soil - Optimum driven - warming
adsorbedC_soil_optimumDriv_warming = optimumDrivenData.Cpools_soil_warming(2:end,5) + optimumDrivenData.Cpools_soil_warming(2:end,8);
monomers_soil_optimumDriv_warming = optimumDrivenData.Cpools_soil_warming(2:end,4);
polymers_soil_optimumDriv_warming = optimumDrivenData.Cpools_soil_warming(2:end,6);
totalCarbon_soil_optimumDriv_warming = sum(optimumDrivenData.Cpools_soil_warming(2:end,[1 2 4 5 6 7 8]),2);

% Soil - Enzyme rigidity - control
adsorbedC_soil_enzRig_control = enzRigData.Cpools_soil_control(2:end,5) + enzRigData.Cpools_soil_control(2:end,8);
monomers_soil_enzRig_control = enzRigData.Cpools_soil_control(2:end,4);
polymers_soil_enzRig_control = enzRigData.Cpools_soil_control(2:end,6);
totalCarbon_soil_enzRig_control = sum(enzRigData.Cpools_soil_control(2:end,[1 2 4 5 6 7 8]),2);

% Soil - Enzyme rigidity - warming
adsorbedC_soil_enzRig_warming = enzRigData.Cpools_soil_warming(2:end,5) + enzRigData.Cpools_soil_warming(2:end,8);
monomers_soil_enzRig_warming = enzRigData.Cpools_soil_warming(2:end,4);
polymers_soil_enzRig_warming = enzRigData.Cpools_soil_warming(2:end,6);
totalCarbon_soil_enzRig_warming = sum(enzRigData.Cpools_soil_warming(2:end,[1 2 4 5 6 7 8]),2);

% =====================
% The plots are created
% =====================

dotSize = 8;
yMin = 0;
yMax = 6000;

scale = 2;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 12*scale], 'color', [1 1 1])

% --------------------------------------
% Soil - no adaptation - control
% --------------------------------------
subplot(3,2,1)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_noAdapt_control,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_noAdapt_control, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_noAdapt_control, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2105', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - no adaptation')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - Optimum driven - control
% --------------------------------------
subplot(3,2,3)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_optimumDriv_control,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_optimumDriv_control, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_optimumDriv_control, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - optimum driven')

set(gca, 'FontSize', 12, 'Ytick', [yMin:1000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - Enzyme rigidity - control
% --------------------------------------
subplot(3,2,5)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_enzRig_control,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_enzRig_control, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_enzRig_control, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - enzyme rigidity')

set(gca, 'FontSize', 12, 'Ytick', [yMin:1000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - no adaptation - warming
% --------------------------------------
subplot(3,2,2)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_noAdapt_warming,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_noAdapt_warming, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_noAdapt_warming, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - no adaptation')

set(gca, 'FontSize', 12, 'Ytick', [yMin:1000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - Optimum driven - warming
% --------------------------------------
subplot(3,2,4)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_optimumDriv_warming,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_optimumDriv_warming, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_optimumDriv_warming, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - optimum driven')

set(gca, 'FontSize', 12, 'Ytick', [yMin:1000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil -Enzyme rigidity - warming
% --------------------------------------
subplot(3,2,6)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_enzRig_warming,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_enzRig_warming, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_enzRig_warming, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - enzyme rigidity')

set(gca, 'FontSize', 12, 'Ytick', [yMin:1000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
text(0.41,0.98,'Soil organic carbon stocks','fontweight','bold','fontsize',14)

% Letters
text(0.04, 0.95, '(a)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.65, '(c)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.35, '(e)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.95, '(b)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.65, '(d)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.35, '(f)', 'FontSize', 14, 'FontWeight', 'bold')

% The legend
l = legend([p3, p2, p1], 'Total C', 'Polymeric C', 'Adsorbed C', 'Orientation', 'Horizontal', 'FontSIze', 14);
pos = l.Position;
pos(1) = 0.3; % x
pos(2) = 0.02; % y
l.Position = pos;
legend('boxoff')
% =========================================================================

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Soil carbon stocks.tiff' -r300 -nocrop -transparent
% =========================================================================

%% Plotting total, forest floor and soil C stocks

% ---------------------
% The data is formatted
% ---------------------

% The modelled total C stocks
totalC_noAdapt_control = totalCarbon_ff_noAdapt_control + totalCarbon_soil_noAdapt_control;
totalC_noAdapt_warming = totalCarbon_ff_noAdapt_warming + totalCarbon_soil_noAdapt_warming;

totalC_optimumDriv_control = totalCarbon_ff_optimumDriv_control + totalCarbon_soil_optimumDriv_control;
totalC_optimumDriv_warming = totalCarbon_ff_optimumDriv_warming + totalCarbon_soil_optimumDriv_warming;

totalC_enzRig_control = totalCarbon_ff_enzRig_control + totalCarbon_soil_enzRig_control;
totalC_enzRig_warming = totalCarbon_ff_enzRig_warming + totalCarbon_soil_enzRig_warming;

% Measurements
totalC_noAdapt_control_meas = Cstock_forestFloor_control + Cstock_soil_control;
totalC_noAdapt_heated_meas = Cstock_forestFloor_warming + Cstock_soil_warming;

totalC_noAdapt_control_meas_stdError = sqrt((CstockStdError_forestFloor_control^2) + (CstockStdError_soil_control^2));
totalC_noAdapt_heated_meas_stdError = sqrt((CstockStdError_forestFloor_warmed^2) + (CstockStdError_soil_warmed^2));

% =====================
% The plots are created
% =====================

color_total = 'k'; 
color_ff = [35,139,69]./255;
color_soil = [204,76,2]./255;

dotSize = 8;
yMin = 0;
yMax = 8000;

scale = 2;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 12*scale], 'color', [1 1 1])

% --------------------------------------
% Soil - no adaptation - control
% --------------------------------------
subplot(3,2,1)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_noAdapt_control,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_noAdapt_control, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_noAdapt_control, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2105', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - no adaptation')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - Optimum driven - control
% --------------------------------------
subplot(3,2,3)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_optimumDriv_control,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_optimumDriv_control, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_optimumDriv_control, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - optimum driven')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - Enzyme rigidity - control
% --------------------------------------
subplot(3,2,5)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_enzRig_control,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_enzRig_control, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_enzRig_control, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - enzyme rigidity')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - no adaptation - warming
% --------------------------------------
subplot(3,2,2)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_noAdapt_warming,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_noAdapt_warming, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_noAdapt_warming, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - no adaptation')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - Optimum driven - warming
% --------------------------------------
subplot(3,2,4)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_optimumDriv_warming,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_optimumDriv_warming, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_optimumDriv_warming, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - optimum driven')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil -Enzyme rigidity - warming
% --------------------------------------
subplot(3,2,6)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_enzRig_warming,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_enzRig_warming, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_enzRig_warming, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - enzyme rigidity')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
text(0.41,0.98,'Organic carbon stocks','fontweight','bold','fontsize',14)

% Letters
text(0.04, 0.95, '(a)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.65, '(c)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.35, '(e)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.95, '(b)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.65, '(d)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.35, '(f)', 'FontSize', 14, 'FontWeight', 'bold')

% The legend
l = legend([p3, p2, p1], 'Total C', 'Soil organic C', 'Forest floor C', 'Orientation', 'Horizontal', 'FontSIze', 14);
pos = l.Position;
pos(1) = 0.3; % x
pos(2) = 0.02; % y
l.Position = pos;
legend('boxoff')
% =========================================================================

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Total carbon stocks.tiff' -r300 -nocrop -transparent
% =========================================================================

%% Cumulative CO2 emissions are plotted

% The colors are defined
c_noAdapt = [31,120,180]./255;
c_optDriv = [51,160,44]./255;
c_enzRig = [255,127,0]./255;

yMin = 0;
yMax = 1600;

lineWidth = 3;
scale = 2;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [6 6 30*scale 12*scale], 'color', [1 1 1])
hold on

date1 = datetime('21/05/2003', 'format', 'dd/MM/yyyy');
date2 = datetime('31/12/2105', 'format', 'dd/MM/yyyy');

p1 = plot(noAdaptData.dates_after1991, noAdaptData.diffCumul, 'color', c_noAdapt, 'linewidth', lineWidth);
p2 = plot(optimumDrivenData.dates_after1991, optimumDrivenData.diffCumul, 'color', c_optDriv, 'linewidth', lineWidth);
p3 = plot(enzRigData.dates_after1991, enzRigData.diffCumul, 'color', c_enzRig, 'linewidth', lineWidth);

lgd = legend([p1 p3 p2], 'No adaptation', 'Enzyme rigidity', 'Optimum driven', 'location', 'northwest');%, 'FontSize', 15)
% lgd = legend([p1 p2 p4], 'No adaptation', 'Optimum driven', 'Measured', 'location', 'northwest');%, 'FontSize', 15)
set(lgd,'FontSize',20);
legend('boxoff')

xlim([date1 date2])

xlabel('Year')
ylabel('g CO_{2}-C')
title('Difference in cumulative CO_{2} emisions between control and heated')

set(gca, 'FontSize', 22, 'Ytick', [yMin:200:yMax])

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Cumulative differences in CO2 fluxes_100year.tiff' -r300 -nocrop -transparent
% =========================================================================

%% Cumulative relative CO2 emissions compared to initial SOC are plotted

% The amount of litter + soil OC at the start of warming is calculated
rowNum = find(dates_all == datetime(2003,05,21));
Cinit_noAdapt = totalC_noAdapt_control(rowNum);
Cinit_optimumDriv = totalC_optimumDriv_control(rowNum);
Cinit_enzRig = totalC_enzRig_control(rowNum);

% The colors are defined
c_noAdapt = [31,120,180]./255;
c_optDriv = [51,160,44]./255;
c_enzRig = [255,127,0]./255;

yMin = 0;
yMax = 25;

lineWidth = 2;
scale = 2;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 8*scale 8*scale], 'color', [1 1 1])
hold on

date1 = datetime('21/05/2003', 'format', 'dd/MM/yyyy');
date2 = datetime('31/12/2105', 'format', 'dd/MM/yyyy');

p1 = plot(noAdaptData.dates_after1991, noAdaptData.diffCumul/Cinit_noAdapt*100, 'color', c_noAdapt, 'linewidth', lineWidth);
p2 = plot(optimumDrivenData.dates_after1991, optimumDrivenData.diffCumul/Cinit_optimumDriv*100, 'color', c_optDriv, 'linewidth', lineWidth);
p3 = plot(enzRigData.dates_after1991, enzRigData.diffCumul/Cinit_enzRig*100, 'color', c_enzRig, 'linewidth', lineWidth);

lgd = legend([p1 p3 p2], 'No adaptation', 'Enzyme rigidity', 'Optimum driven', 'location', 'northwest');%, 'FontSize', 14)
% lgd = legend([p1 p2 p4], 'No adaptation', 'Optimum driven', 'Measured', 'location', 'northwest');%, 'FontSize', 15)
set(lgd,'FontSize',14);
legend('boxoff')

xlim([date1 date2])
ylim([0 20.5])

xlabel('Years after the initiation of soil warming')
ylabel('Relative SOC loss (%)')

set(gca, 'FontSize', 14, 'Ytick', [yMin:5:yMax], 'Xtick', [date1:calyears(20):date2], 'xticklabels', {'0','20','40','60','80', '100'}, 'box', 'on');

title([' Relative SOC loss as a consequence of '; 'linear 3 °C soil warming over 100 years'], 'FontSize', 14)

pos = lgd.Position;
pos(1) = 0.15; % x
pos(2) = 0.78; % y
lgd.Position = pos;
legend('boxoff')

% =========================================================================
% Export the figure
addpath('export_fig')
export_fig 'Relative SOC loss_100year.tiff' -r300 -nocrop -transparent
% =========================================================================

%% Annual CO2 emissions

% The colors are defined
c_noAdapt = [31,120,180]./255;
c_optDriv = [51,160,44]./255;
c_enzRig = [255,127,0]./255;
color_diff = [.5 .5 .5];%[1,102,94]./255;

lineWidth = .5;
scale = 2;

yMin = -50;
yMax = 20;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 80*scale 12*scale], 'color', [1 1 1])
hold on

date1 = datetime('01/01/1989', 'format', 'dd/MM/yyyy');
date2 = datetime('31/12/2105', 'format', 'dd/MM/yyyy');
p0 = plot([2003 2105],[0 0], '--', 'color', [.2 .2 .2]);

% Let's add a datapoint for 2002, when control and heated fluxes are equal
yrs = [noAdaptData.uniqueYears_CO2diff(1)-1 noAdaptData.uniqueYears_CO2diff];

% Modelled
p1 = plot(yrs, [0; noAdaptData.diffAnnualCO2], '--o', 'color', c_noAdapt, 'linewidth', lineWidth, 'MarkerSize', 6, 'MarkerEdgeColor', c_noAdapt,'MarkerFaceColor', c_noAdapt);
p2 = plot(yrs, [0; optimumDrivenData.diffAnnualCO2], '--o', 'color', c_optDriv, 'linewidth', lineWidth, 'MarkerSize', 6, 'MarkerEdgeColor', c_optDriv,'MarkerFaceColor', c_optDriv);
p3 = plot(yrs, [0; enzRigData.diffAnnualCO2], '--o', 'color', c_enzRig, 'linewidth', lineWidth, 'MarkerSize', 6, 'MarkerEdgeColor', c_enzRig,'MarkerFaceColor', c_enzRig);

% p1 = scatter(noAdaptData.uniqueYears_CO2diff,noAdaptData.diffAnnualCO2, 'MarkerFaceColor', c_noAdapt);
% p2 = scatter(optimumDrivenData.uniqueYears_CO2diff,optimumDrivenData.diffAnnualCO2, 'MarkerFaceColor', c_optDriv);
     
lgd = legend([p1 p3 p2], 'No adaptation', 'Enzyme rigidity', 'Optimum driven', 'location', 'northwest');%, 'FontSize', 15)
% lgd = legend([b1 p1 p2], 'Measured (Melillo et al., 2011)', 'No adaptation', 'Optimum driven', 'location', 'southwest');%, 'FontSize', 15)
set(lgd,'FontSize',20);
legend('boxoff')

xlim([2001 2105])
ylim([0 yMax])
set(gca, 'FontSize', 22, 'Ytick', [yMin:5:yMax])

xlabel('Year')
ylabel('g CO_{2}-C yr^{-1}')
title('Difference in annual CO_{2} emissions between control and heated')

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Annual differences in CO2 fluxes.tiff' -r300 -nocrop -transparent
% =========================================================================

%% The daily CO2 fluxes are plotted for total CO2

% ------------------
% The data is loaded
% ------------------

totalCO2_noAdapt_control = noAdaptData.totalCO2_control;
totalCO2_noAdapt_warmed = noAdaptData.totalCO2_warmed;

totalCO2_optimumDriv_control = optimumDrivenData.totalCO2_control;
totalCO2_optimumDriv_warmed = optimumDrivenData.totalCO2_warmed;

totalCO2_enzRig_control = enzRigData.totalCO2_control;
totalCO2_enzRig_warmed = enzRigData.totalCO2_warmed;

dates_all = optimumDrivenData.dates_all;

% --------------------
% The data are plotted
% --------------------

% The colors are defined
c_noAdapt = [31,120,180]./255;
c_optDriv = [51,160,44]./255;
c_enzRig = [255,127,0]./255;
c_meas = [215,48,31]./255;

dotSize = 20;
lineWidth = 1;
scale = 2;

yMin = 0;
yMax = 5.5;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 10*scale], 'color', [1 1 1])
hold on

% -----------------
% The control plots
% -----------------
subplot(2,1,1)
hold on

% The modelled CO2 is plotted
p1 = plot(dates_all,totalCO2_noAdapt_control, 'color', c_noAdapt, 'LineWidth', lineWidth);
p2 = plot(dates_all,totalCO2_optimumDriv_control, 'color', c_optDriv, 'LineWidth', lineWidth);
p3 = plot(dates_all,totalCO2_enzRig_control, 'color', c_enzRig, 'LineWidth', lineWidth);

% The start of warming is indicated
plot([date_startWarming date_startWarming],[0 yMax], ':', 'color', 'k', 'lineWidth', 2);

% Formatting
xlabel('Years')
ylabel('CO_{2}-C (g m^{-2} day{-1})')
title('Control treatment')

date1 = datetime('01/01/2000', 'InputFormat', 'dd/MM/yyyy');
date2 = datetime('31/12/2023', 'InputFormat', 'dd/MM/yyyy');

xlim([date1 date2])
ylim([yMin yMax])

set(gca, 'FontSize', 12, 'Ytick', [yMin:1:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.61; % y
pos(3) = 0.9; % width
pos(4) = 0.34; % height
set(gca,'position',pos)

% ----------------
% The heated plots
% ----------------
subplot(2,1,2)
hold on

% The modelled CO2 is plotted
p1 = plot(dates_all, totalCO2_noAdapt_warmed, 'color', c_noAdapt, 'LineWidth', lineWidth);
p2 = plot(dates_all, totalCO2_optimumDriv_warmed, 'color', c_optDriv, 'LineWidth', lineWidth);
p3 = plot(dates_all, totalCO2_enzRig_warmed, 'color', c_enzRig, 'LineWidth', lineWidth);

% The start of warming is indicated
plot([date_startWarming date_startWarming],[0 yMax], ':', 'color', 'k', 'lineWidth', 2);

% Formatting
xlabel('Years')
ylabel('CO_{2}-C (g m^{-2} day{-1})')
title('Heated treatment')

date1 = datetime('01/01/2000', 'InputFormat', 'dd/MM/yyyy');
date2 = datetime('31/12/2023', 'InputFormat', 'dd/MM/yyyy');

xlim([date1 date2])
ylim([yMin yMax])

set(gca, 'FontSize', 12, 'Ytick', [yMin:1:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.14; % y
pos(3) = 0.9; % width
pos(4) = 0.34; % height
set(gca,'position',pos)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% Letters
text(0.03, 0.97, '(a)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.03, 0.50, '(b)', 'FontSize', 14, 'FontWeight', 'bold')

% % The legend
l = legend([p1, p2, p3], 'No adaptation', 'Optimum driven', 'Enzyme rigidity', 'Orientation', 'Horizontal', 'FontSIze', 14);
pos = l.Position;
pos(1) = 0.3; % x
pos(2) = 0.02; % y
l.Position = pos;
legend('boxoff')
% =========================================================================

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Daily CO2 flux timeseries.tiff' -r300 -nocrop -transparent
% =========================================================================

%%

clc; clearvars; close all

% -------------------------------------------
% Different C inputs for the heated treatment
% -------------------------------------------

%% The data is loaded


mainFolder = '/Users/vamarijn/Documents/ETH/DEEP - C/Codes/ReSOM Matlab/16. ReSOM - Matlab - new MMRT formlation -max1 - newCO2Data/model/Data Barre Woods/Data artificial run/100 years - C inputs';

% -------------------
% No adaptation
% -------------------

load([mainFolder '/noAdapt_CinputConstant.mat']); CinputConstant = out; clear out

% -------------------
% Optimum driven runs
% -------------------

load([mainFolder '/noAdapt_CinputMinus10perc.mat']); CinputMinus10perc = out; clear out

% --------------------
% Enzyme rigidity runs
% --------------------

load([mainFolder '/noAdapt_CinputPlus10perc.mat']); CinputPlus10perc = out; clear out

clear mainFolder

%% The measured data is loaded

site = 'barreWoods';

rmpath([pwd '/Measured Data Barre Woods'])
addpath('Measured data Barre Woods')

% The number of treatment years is defined
maxYearsWithData = 13;      % The number of years for which data is available
treatment_years = 13;       % The number of years the treatment is run

% spinup for a defined number of years to steady state
continueOn = 0; % 
date_startWarming = datetime(2003,05,21);       % The date the warming started
spinup_years = 30;%100;     % Number of spin-up years
date_start_spinup = date_startWarming - calyears(spinup_years);
endDate_startPlusTen = date_startWarming + calyears(treatment_years);
endDate_fullYear = ['31/12/' int2str(endDate_startPlusTen.Year)];
endDate = datetime(endDate_fullYear, 'InputFormat', 'dd/MM/yyyy', 'Format', 'dd/MM/yyyy');     

d1 = 40;        % The large timestep (days)
d2 = 1;         % The small timestep (days)
timeStepBuffer = 10;  

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

dates_all = [dates_spinup dates_treatmentRun];

Load_measured_data;

%% The OC pools are plotted for the FOREST FLOOR

% The colors are defined
color_monomers_litter = [140,107,177]./255;
color_adsorbedC_litter = [140,107,177]./255;
color_metabolic_litter = [254,196,79]./255;
color_structural_litter = [217,95,14]./255;
color_totalCarbon_litter = 'k';
color_measuredPOC = [153,52,4]./255;

lineWidth = 2;
startWarming = datetime('01/01/2000', 'format', 'dd/MM/yyyy');

% =====================
% The data is formatted
% =====================

% Forest floor - No adaptation - control
adsorbedC_ff_CinputConstant_control = CinputConstant.Cpools_litter_control(2:end,7) + CinputConstant.Cpools_litter_control(2:end,12) + CinputConstant.Cpools_litter_control(2:end,13);
metabolic_ff_CinputConstant_control = CinputConstant.Cpools_litter_control(2:end,8);
structural_ff_CinputConstant_control = CinputConstant.Cpools_litter_control(2:end,9);
enz_rStrat_ff_CinputConstant_control = CinputConstant.Cpools_litter_control(2:end,10);
enz_kStrat_ff_CinputConstant_control = CinputConstant.Cpools_litter_control(2:end,11);
totalCarbon_ff_CinputConstant_control = sum(CinputConstant.Cpools_litter_control(2:end,[1:4 6:13]),2);

% Forest floor - No adaptation - warming
adsorbedC_ff_CinputConstant_warming = CinputConstant.Cpools_litter_warming(2:end,7) + CinputConstant.Cpools_litter_warming(2:end,12) + CinputConstant.Cpools_litter_warming(2:end,13);
metabolic_ff_CinputConstant_warming = CinputConstant.Cpools_litter_warming(2:end,8);
structural_ff_CinputConstant_warming = CinputConstant.Cpools_litter_warming(2:end,9);
enz_rStrat_ff_CinputConstant_warming = CinputConstant.Cpools_litter_warming(2:end,10);
enz_kStrat_ff_CinputConstant_warming = CinputConstant.Cpools_litter_warming(2:end,11);
totalCarbon_ff_CinputConstant_warming = sum(CinputConstant.Cpools_litter_warming(2:end,[1:4 6:13]),2);

% Forest floor - Optimum driven - control
adsorbedC_ff_CinputMinus10perc_control = CinputMinus10perc.Cpools_litter_control(2:end,7) + CinputMinus10perc.Cpools_litter_control(2:end,12) + CinputMinus10perc.Cpools_litter_control(2:end,13);
metabolic_ff_CinputMinus10perc_control = CinputMinus10perc.Cpools_litter_control(2:end,8);
structural_ff_CinputMinus10perc_control = CinputMinus10perc.Cpools_litter_control(2:end,9);
enz_rStrat_ff_CinputMinus10perc_control = CinputMinus10perc.Cpools_litter_control(2:end,10);
enz_kStrat_ff_CinputMinus10perc_control = CinputMinus10perc.Cpools_litter_control(2:end,11);
totalCarbon_ff_CinputMinus10perc_control = sum(CinputMinus10perc.Cpools_litter_control(2:end,[1:4 6:13]),2);

% Forest floor - Optimum driven - warming
adsorbedC_ff_CinputMinus10perc_warming = CinputMinus10perc.Cpools_litter_warming(2:end,7) + CinputMinus10perc.Cpools_litter_warming(2:end,12) + CinputMinus10perc.Cpools_litter_warming(2:end,13);
metabolic_ff_CinputMinus10perc_warming = CinputMinus10perc.Cpools_litter_warming(2:end,8);
structural_ff_CinputMinus10perc_warming = CinputMinus10perc.Cpools_litter_warming(2:end,9);
enz_rStrat_ff_CinputMinus10perc_warming = CinputMinus10perc.Cpools_litter_warming(2:end,10);
enz_kStrat_ff_CinputMinus10perc_warming = CinputMinus10perc.Cpools_litter_warming(2:end,11);
totalCarbon_ff_CinputMinus10perc_warming = sum(CinputMinus10perc.Cpools_litter_warming(2:end,[1:4 6:13]),2);

% % Forest floor - Enzyme rigidity - control
adsorbedC_ff_CinputPlus10perc_control = CinputPlus10perc.Cpools_litter_control(2:end,7) + CinputPlus10perc.Cpools_litter_control(2:end,12) + CinputPlus10perc.Cpools_litter_control(2:end,13);
metabolic_ff_CinputPlus10perc_control = CinputPlus10perc.Cpools_litter_control(2:end,8);
structural_ff_CinputPlus10perc_control = CinputPlus10perc.Cpools_litter_control(2:end,9);
enz_rStrat_ff_CinputPlus10perc_control = CinputPlus10perc.Cpools_litter_control(2:end,10);
enz_kStrat_ff_CinputPlus10perc_control = CinputPlus10perc.Cpools_litter_control(2:end,11);
totalCarbon_ff_CinputPlus10perc_control = sum(CinputPlus10perc.Cpools_litter_control(2:end,[1:4 6:13]),2);

% Forest floor - Enzyme rigidity - warming
adsorbedC_ff_CinputPlus10perc_warming = CinputPlus10perc.Cpools_litter_warming(2:end,7) + CinputPlus10perc.Cpools_litter_warming(2:end,12) + CinputPlus10perc.Cpools_litter_warming(2:end,13);
metabolic_ff_CinputPlus10perc_warming = CinputPlus10perc.Cpools_litter_warming(2:end,8);
structural_ff_CinputPlus10perc_warming = CinputPlus10perc.Cpools_litter_warming(2:end,9);
enz_rStrat_ff_CinputPlus10perc_warming = CinputPlus10perc.Cpools_litter_warming(2:end,10);
enz_kStrat_ff_CinputPlus10perc_warming = CinputPlus10perc.Cpools_litter_warming(2:end,11);
totalCarbon_ff_CinputPlus10perc_warming = sum(CinputPlus10perc.Cpools_litter_warming(2:end,[1:4 6:13]),2);

% =====================
% The plots are created
% =====================

dotSize = 8;
yMin = 0;
yMax = 2500;

scale = 2;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 12*scale], 'color', [1 1 1])

% --------------------------------------
% Forest floor - CinputConstant - control
% --------------------------------------
subplot(3,2,1)
hold on

dates_all = CinputPlus10perc.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_CinputConstant_control,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_CinputConstant_control,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_CinputConstant_control, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_CinputConstant_control, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2105', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Constant inputs')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Forest floor - CinputMinus10perc - control
% --------------------------------------
subplot(3,2,3)
hold on

% dates_all = CinputMinus10perc.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_CinputMinus10perc_control, 'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_CinputMinus10perc_control,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_CinputMinus10perc_control, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_CinputMinus10perc_control, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Inputs -10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Forest floor - CinputPlus10perc - control
% --------------------------------------
subplot(3,2,5)
hold on

% dates_all = noAdaptData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_CinputPlus10perc_control,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_CinputPlus10perc_control,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_CinputPlus10perc_control, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_CinputPlus10perc_control, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Inputs +10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Forest floor - CinputConstant - warming
% --------------------------------------
subplot(3,2,2)
hold on

% dates_all = noAdaptData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_CinputConstant_warming,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_CinputConstant_warming,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_CinputConstant_warming, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_CinputConstant_warming, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Constant inputs')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Forest floor - CinputMinus10perc - warming
% --------------------------------------
subplot(3,2,4)
hold on

% dates_all = optimumDrivenData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_CinputMinus10perc_warming,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_CinputMinus10perc_warming,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_CinputMinus10perc_warming, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_CinputMinus10perc_warming, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Inputs -10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Forest floor - CinputPlus10perc - warming
% --------------------------------------
subplot(3,2,6)
hold on

% dates_all = noAdaptData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_CinputPlus10perc_warming,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_CinputPlus10perc_warming,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_CinputPlus10perc_warming, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_CinputPlus10perc_warming, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Inputs +10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
text(0.38,0.98,'Forest floor organic carbon stocks','fontweight','bold','fontsize',14)

% Letters
text(0.04, 0.95, '(a)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.65, '(c)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.35, '(e)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.95, '(b)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.65, '(d)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.35, '(f)', 'FontSize', 14, 'FontWeight', 'bold')

% The legend
l = legend([p4, p3, p2, p1], 'Total C', 'Structural C', 'Metabolic C', 'Adsorbed C', 'Orientation', 'Horizontal', 'FontSIze', 14);
pos = l.Position;
pos(1) = 0.3; % x
pos(2) = 0.02; % y
l.Position = pos;
legend('boxoff')

% =========================================================================

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Forest floor carbon stocks.tiff' -r300 -nocrop -transparent
% =========================================================================

%% The OC pools are plotted for the SOIL

% The colors are defined
color_monomers_soil = [140,107,177]./255;
color_adsorbedC_soil = [140,107,177]./255;
color_polymers_soil = [153,52,4]./255;
color_totalCarbon_soil = 'k';

lineWidth = 2;
startWarming = datetime('01/01/2000', 'format', 'dd/MM/yyyy');

% =====================
% The data is formatted
% =====================

% Soil - No adaptation - control
adsorbedC_soil_CinputConstant_control = CinputConstant.Cpools_soil_control(2:end,5) + CinputConstant.Cpools_soil_control(2:end,8);
monomers_soil_CinputConstant_control = CinputConstant.Cpools_soil_control(2:end,4);
polymers_soil_CinputConstant_control = CinputConstant.Cpools_soil_control(2:end,6);
totalCarbon_soil_CinputConstant_control = sum(CinputConstant.Cpools_soil_control(2:end,[1 2 4 5 6 7 8]),2);

% Soil - No adaptation - warming
adsorbedC_soil_CinputConstant_warming = CinputConstant.Cpools_soil_warming(2:end,5) + CinputConstant.Cpools_soil_warming(2:end,8);
monomers_soil_CinputConstant_warming = CinputConstant.Cpools_soil_warming(2:end,4);
polymers_soil_CinputConstant_warming = CinputConstant.Cpools_soil_warming(2:end,6);
totalCarbon_soil_CinputConstant_warming = sum(CinputConstant.Cpools_soil_warming(2:end,[1 2 4 5 6 7 8]),2);

% Soil - Optimum driven - control
adsorbedC_soil_CinputMinus10perc_control = CinputMinus10perc.Cpools_soil_control(2:end,5) + CinputMinus10perc.Cpools_soil_control(2:end,8);
monomers_soil_CinputMinus10perc_control = CinputMinus10perc.Cpools_soil_control(2:end,4);
polymers_soil_CinputMinus10perc_control = CinputMinus10perc.Cpools_soil_control(2:end,6);
totalCarbon_soil_CinputMinus10perc_control = sum(CinputMinus10perc.Cpools_soil_control(2:end,[1 2 4 5 6 7 8]),2);

% Soil - Optimum driven - warming
adsorbedC_soil_CinputMinus10perc_warming = CinputMinus10perc.Cpools_soil_warming(2:end,5) + CinputMinus10perc.Cpools_soil_warming(2:end,8);
monomers_soil_CinputMinus10perc_warming = CinputMinus10perc.Cpools_soil_warming(2:end,4);
polymers_soil_CinputMinus10perc_warming = CinputMinus10perc.Cpools_soil_warming(2:end,6);
totalCarbon_soil_CinputMinus10perc_warming = sum(CinputMinus10perc.Cpools_soil_warming(2:end,[1 2 4 5 6 7 8]),2);

% Soil - Enzyme rigidity - control
adsorbedC_soil_CinputPlus10perc_control = CinputPlus10perc.Cpools_soil_control(2:end,5) + CinputPlus10perc.Cpools_soil_control(2:end,8);
monomers_soil_CinputPlus10perc_control = CinputPlus10perc.Cpools_soil_control(2:end,4);
polymers_soil_CinputPlus10perc_control = CinputPlus10perc.Cpools_soil_control(2:end,6);
totalCarbon_soil_CinputPlus10perc_control = sum(CinputPlus10perc.Cpools_soil_control(2:end,[1 2 4 5 6 7 8]),2);

% Soil - Enzyme rigidity - warming
adsorbedC_soil_CinputPlus10perc_warming = CinputPlus10perc.Cpools_soil_warming(2:end,5) + CinputPlus10perc.Cpools_soil_warming(2:end,8);
monomers_soil_CinputPlus10perc_warming = CinputPlus10perc.Cpools_soil_warming(2:end,4);
polymers_soil_CinputPlus10perc_warming = CinputPlus10perc.Cpools_soil_warming(2:end,6);
totalCarbon_soil_CinputPlus10perc_warming = sum(CinputPlus10perc.Cpools_soil_warming(2:end,[1 2 4 5 6 7 8]),2);

% =====================
% The plots are created
% =====================

dotSize = 8;
yMin = 0;
yMax = 6000;

scale = 2;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 12*scale], 'color', [1 1 1])

% --------------------------------------
% Soil - CinputConstant - control
% --------------------------------------
subplot(3,2,1)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_CinputConstant_control,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_CinputConstant_control, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_CinputConstant_control, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2105', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Constant inputs')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - CinputMinus10perc - control
% --------------------------------------
subplot(3,2,3)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_CinputMinus10perc_control,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_CinputMinus10perc_control, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_CinputMinus10perc_control, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Inputs -10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:1000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - CinputPlus10perc - control
% --------------------------------------
subplot(3,2,5)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_CinputPlus10perc_control,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_CinputPlus10perc_control, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_CinputPlus10perc_control, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Inputs +10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:1000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - CinputConstant - warming
% --------------------------------------
subplot(3,2,2)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_CinputConstant_warming,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_CinputConstant_warming, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_CinputConstant_warming, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Constant inputs')

set(gca, 'FontSize', 12, 'Ytick', [yMin:1000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - CinputMinus10perc - warming
% --------------------------------------
subplot(3,2,4)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_CinputMinus10perc_warming,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_CinputMinus10perc_warming, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_CinputMinus10perc_warming, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Inputs -10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:1000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - CinputPlus10perc - warming
% --------------------------------------
subplot(3,2,6)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_CinputPlus10perc_warming,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_CinputPlus10perc_warming, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_CinputPlus10perc_warming, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Inputs +10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:1000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
text(0.41,0.98,'Soil organic carbon stocks','fontweight','bold','fontsize',14)

% Letters
text(0.04, 0.95, '(a)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.65, '(c)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.35, '(e)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.95, '(b)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.65, '(d)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.35, '(f)', 'FontSize', 14, 'FontWeight', 'bold')

% The legend
l = legend([p3, p2, p1], 'Total C', 'Polymeric C', 'Adsorbed C', 'Orientation', 'Horizontal', 'FontSIze', 14);
pos = l.Position;
pos(1) = 0.3; % x
pos(2) = 0.02; % y
l.Position = pos;
legend('boxoff')
% =========================================================================

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Soil carbon stocks.tiff' -r300 -nocrop -transparent
% =========================================================================

%% Plotting total, forest floor and soil C stocks

% ---------------------
% The data is formatted
% ---------------------

% The modelled total C stocks
totalC_CinputConstant_control = totalCarbon_ff_CinputConstant_control + totalCarbon_soil_CinputConstant_control;
totalC_CinputConstant_warming = totalCarbon_ff_CinputConstant_warming + totalCarbon_soil_CinputConstant_warming;

totalC_CinputMinus10perc_control = totalCarbon_ff_CinputMinus10perc_control + totalCarbon_soil_CinputMinus10perc_control;
totalC_CinputMinus10perc_warming = totalCarbon_ff_CinputMinus10perc_warming + totalCarbon_soil_CinputMinus10perc_warming;

totalC_CinputPlus10perc_control = totalCarbon_ff_CinputPlus10perc_control + totalCarbon_soil_CinputPlus10perc_control;
totalC_CinputPlus10perc_warming = totalCarbon_ff_CinputPlus10perc_warming + totalCarbon_soil_CinputPlus10perc_warming;

% Measurements
totalC_noAdapt_control_meas = Cstock_forestFloor_control + Cstock_soil_control;
totalC_noAdapt_heated_meas = Cstock_forestFloor_warming + Cstock_soil_warming;

totalC_noAdapt_control_meas_stdError = sqrt((CstockStdError_forestFloor_control^2) + (CstockStdError_soil_control^2));
totalC_noAdapt_heated_meas_stdError = sqrt((CstockStdError_forestFloor_warmed^2) + (CstockStdError_soil_warmed^2));

% =====================
% The plots are created
% =====================

color_total = 'k'; 
color_ff = [35,139,69]./255;
color_soil = [204,76,2]./255;

dotSize = 8;
yMin = 0;
yMax = 8000;

scale = 2;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 12*scale], 'color', [1 1 1])

% --------------------------------------
% Soil - CinputConstant - control
% --------------------------------------
subplot(3,2,1)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_CinputConstant_control,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_CinputConstant_control, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_CinputConstant_control, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2105', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Constant inputs')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - CinputMinus10perc - control
% --------------------------------------
subplot(3,2,3)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_CinputMinus10perc_control,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_CinputMinus10perc_control, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_CinputMinus10perc_control, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Inputs -10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - CinputPlus10perc - control
% --------------------------------------
subplot(3,2,5)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_CinputPlus10perc_control,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_CinputPlus10perc_control, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_CinputPlus10perc_control, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Inputs +10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - CinputConstant - warming
% --------------------------------------
subplot(3,2,2)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_CinputConstant_warming,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_CinputConstant_warming, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_CinputConstant_warming, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Constant inputs')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - CinputMinus10perc - warming
% --------------------------------------
subplot(3,2,4)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_CinputMinus10perc_warming,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_CinputMinus10perc_warming, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_CinputMinus10perc_warming, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Inputs -10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - CinputPlus10perc - warming
% --------------------------------------
subplot(3,2,6)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_CinputPlus10perc_warming,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_CinputPlus10perc_warming, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_CinputPlus10perc_warming, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Inputs +10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
text(0.41,0.98,'Organic carbon stocks','fontweight','bold','fontsize',14)

% Letters
text(0.04, 0.95, '(a)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.65, '(c)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.35, '(e)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.95, '(b)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.65, '(d)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.35, '(f)', 'FontSize', 14, 'FontWeight', 'bold')

% The legend
l = legend([p3, p2, p1], 'Total C', 'Soil organic C', 'Forest floor C', 'Orientation', 'Horizontal', 'FontSIze', 14);
pos = l.Position;
pos(1) = 0.3; % x
pos(2) = 0.02; % y
l.Position = pos;
legend('boxoff')
% =========================================================================

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Total carbon stocks.tiff' -r300 -nocrop -transparent
% =========================================================================

%% Cumulative CO2 emissions are plotted

% The colors are defined
c_noAdapt = [31,120,180]./255;
c_optDriv = [51,160,44]./255;
c_enzRig = [255,127,0]./255;

yMin = -4000;
yMax = 7000;

lineWidth = 2;
scale = 2;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 8*scale 8*scale], 'color', [1 1 1])
hold on

p1 = plot(CinputConstant.dates_after1991, CinputConstant.diffCumul, 'color', c_noAdapt, 'linewidth', lineWidth);
p2 = plot(CinputMinus10perc.dates_after1991, CinputMinus10perc.diffCumul, 'color', c_optDriv, 'linewidth', lineWidth);
p3 = plot(CinputPlus10perc.dates_after1991, CinputPlus10perc.diffCumul, 'color', c_enzRig, 'linewidth', lineWidth);

lgd = legend([p1 p3 p2], 'Constant inputs', 'Inputs +10 %', 'Inputs -10 %', 'location', 'northwest');%, 'FontSize', 15)
% lgd = legend([p1 p2 p4], 'No adaptation', 'Optimum driven', 'Measured', 'location', 'northwest');%, 'FontSize', 15)
set(lgd,'FontSize',12);
legend('boxoff')

ylim([yMin yMax])

xlabel('Year')
ylabel('g CO_{2}-C')
title('Difference in cumulative CO_{2} emisions between control and heated')

set(gca, 'FontSize', 12, 'Ytick', [yMin:1000:yMax])

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Cumulative differences in CO2 fluxes.tiff' -r300 -nocrop -transparent
% =========================================================================

%% Annual CO2 emissions

% The colors are defined
c_noAdapt = [31,120,180]./255;
c_optDriv = [51,160,44]./255;
c_enzRig = [255,127,0]./255;
color_diff = [.5 .5 .5];%[1,102,94]./255;

lineWidth = .5;
scale = 2;

yMin = -50;
yMax = 70;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 8*scale], 'color', [1 1 1])
hold on

date1 = datetime('01/01/1989', 'format', 'dd/MM/yyyy');
date2 = datetime('31/12/2018', 'format', 'dd/MM/yyyy');
p0 = plot([2003 2105],[0 0], '--', 'color', [.2 .2 .2]);

% Let's add a datapoint for 2002, when control and heated fluxes are equal
yrs = [CinputConstant.uniqueYears_CO2diff(1)-1 CinputConstant.uniqueYears_CO2diff];

% Modelled
p1 = plot(yrs, [0; CinputConstant.diffAnnualCO2], '--o', 'color', c_noAdapt, 'linewidth', lineWidth, 'MarkerSize', 6, 'MarkerEdgeColor', c_noAdapt,'MarkerFaceColor', c_noAdapt);
p2 = plot(yrs, [0; CinputMinus10perc.diffAnnualCO2], '--o', 'color', c_optDriv, 'linewidth', lineWidth, 'MarkerSize', 6, 'MarkerEdgeColor', c_optDriv,'MarkerFaceColor', c_optDriv);
p3 = plot(yrs, [0; CinputPlus10perc.diffAnnualCO2], '--o', 'color', c_enzRig, 'linewidth', lineWidth, 'MarkerSize', 6, 'MarkerEdgeColor', c_enzRig,'MarkerFaceColor', c_enzRig);

% p1 = scatter(noAdaptData.uniqueYears_CO2diff,noAdaptData.diffAnnualCO2, 'MarkerFaceColor', c_noAdapt);
% p2 = scatter(optimumDrivenData.uniqueYears_CO2diff,optimumDrivenData.diffAnnualCO2, 'MarkerFaceColor', c_optDriv);
     
lgd = legend([p1 p3 p2], 'Constant inputs', 'Inputs +10 %', 'Inputs -10 %', 'location', 'northwest');%, 'FontSize', 15)
% lgd = legend([b1 p1 p2], 'Measured (Melillo et al., 2011)', 'No adaptation', 'Optimum driven', 'location', 'southwest');%, 'FontSize', 15)
set(lgd,'FontSize',12);
legend('boxoff')

ylim([yMin yMax])

xlim([2001 2105])
% ylim([0 yMax])
set(gca, 'FontSize', 12, 'Ytick', [yMin:10:yMax])

xlabel('Year')
ylabel('g CO_{2}-C yr^{-1}')
title('Difference in annual CO_{2} emissions between control and heated')

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Annual differences in CO2 fluxes.tiff' -r300 -nocrop -transparent
% =========================================================================

%%

clc; clearvars; close all

% -------------------------------------------
% Different moisture for the heated treatment
% -------------------------------------------

%% The data is loaded


mainFolder = '/Users/vamarijn/Documents/ETH/DEEP - C/Codes/ReSOM Matlab/16. ReSOM - Matlab - new MMRT formlation -max1 - newCO2Data/model/Data Barre Woods/Data artificial run/100 years - Moisture';

% -------------------
% No adaptation
% -------------------

load([mainFolder '/noAdapt_moistureConstant.mat']); moistureConstant = out; clear out

% -------------------
% Optimum driven runs
% -------------------

load([mainFolder '/noAdapt_moistureMinus10perc.mat']); moistureMinus10perc = out; clear out

% --------------------
% Enzyme rigidity runs
% --------------------

load([mainFolder '/noAdapt_moisturePlus10perc.mat']); moisturePlus10perc = out; clear out

clear mainFolder

%% The measured data is loaded

site = 'barreWoods';

rmpath([pwd '/Measured Data Barre Woods'])
addpath('Measured data Barre Woods')

% The number of treatment years is defined
maxYearsWithData = 13;      % The number of years for which data is available
treatment_years = 13;       % The number of years the treatment is run

% spinup for a defined number of years to steady state
continueOn = 0; % 
date_startWarming = datetime(2003,05,21);       % The date the warming started
spinup_years = 30;%100;     % Number of spin-up years
date_start_spinup = date_startWarming - calyears(spinup_years);
endDate_startPlusTen = date_startWarming + calyears(treatment_years);
endDate_fullYear = ['31/12/' int2str(endDate_startPlusTen.Year)];
endDate = datetime(endDate_fullYear, 'InputFormat', 'dd/MM/yyyy', 'Format', 'dd/MM/yyyy');     

d1 = 40;        % The large timestep (days)
d2 = 1;         % The small timestep (days)
timeStepBuffer = 10;  

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

dates_all = [dates_spinup dates_treatmentRun];

Load_measured_data;

%% The OC pools are plotted for the FOREST FLOOR

% The colors are defined
color_monomers_litter = [140,107,177]./255;
color_adsorbedC_litter = [140,107,177]./255;
color_metabolic_litter = [254,196,79]./255;
color_structural_litter = [217,95,14]./255;
color_totalCarbon_litter = 'k';
color_measuredPOC = [153,52,4]./255;

lineWidth = 2;
startWarming = datetime('01/01/2000', 'format', 'dd/MM/yyyy');

% =====================
% The data is formatted
% =====================

% Forest floor - No adaptation - control
adsorbedC_ff_moistureConstant_control = moistureConstant.Cpools_litter_control(2:end,7) + moistureConstant.Cpools_litter_control(2:end,12) + moistureConstant.Cpools_litter_control(2:end,13);
metabolic_ff_moistureConstant_control = moistureConstant.Cpools_litter_control(2:end,8);
structural_ff_moistureConstant_control = moistureConstant.Cpools_litter_control(2:end,9);
enz_rStrat_ff_moistureConstant_control = moistureConstant.Cpools_litter_control(2:end,10);
enz_kStrat_ff_moistureConstant_control = moistureConstant.Cpools_litter_control(2:end,11);
totalCarbon_ff_moistureConstant_control = sum(moistureConstant.Cpools_litter_control(2:end,[1:4 6:13]),2);

% Forest floor - No adaptation - warming
adsorbedC_ff_moistureConstant_warming = moistureConstant.Cpools_litter_warming(2:end,7) + moistureConstant.Cpools_litter_warming(2:end,12) + moistureConstant.Cpools_litter_warming(2:end,13);
metabolic_ff_moistureConstant_warming = moistureConstant.Cpools_litter_warming(2:end,8);
structural_ff_moistureConstant_warming = moistureConstant.Cpools_litter_warming(2:end,9);
enz_rStrat_ff_moistureConstant_warming = moistureConstant.Cpools_litter_warming(2:end,10);
enz_kStrat_ff_moistureConstant_warming = moistureConstant.Cpools_litter_warming(2:end,11);
totalCarbon_ff_moistureConstant_warming = sum(moistureConstant.Cpools_litter_warming(2:end,[1:4 6:13]),2);

% Forest floor - Optimum driven - control
adsorbedC_ff_moistureMinus10perc_control = moistureMinus10perc.Cpools_litter_control(2:end,7) + moistureMinus10perc.Cpools_litter_control(2:end,12) + moistureMinus10perc.Cpools_litter_control(2:end,13);
metabolic_ff_moistureMinus10perc_control = moistureMinus10perc.Cpools_litter_control(2:end,8);
structural_ff_moistureMinus10perc_control = moistureMinus10perc.Cpools_litter_control(2:end,9);
enz_rStrat_ff_moistureMinus10perc_control = moistureMinus10perc.Cpools_litter_control(2:end,10);
enz_kStrat_ff_moistureMinus10perc_control = moistureMinus10perc.Cpools_litter_control(2:end,11);
totalCarbon_ff_moistureMinus10perc_control = sum(moistureMinus10perc.Cpools_litter_control(2:end,[1:4 6:13]),2);

% Forest floor - Optimum driven - warming
adsorbedC_ff_moistureMinus10perc_warming = moistureMinus10perc.Cpools_litter_warming(2:end,7) + moistureMinus10perc.Cpools_litter_warming(2:end,12) + moistureMinus10perc.Cpools_litter_warming(2:end,13);
metabolic_ff_moistureMinus10perc_warming = moistureMinus10perc.Cpools_litter_warming(2:end,8);
structural_ff_moistureMinus10perc_warming = moistureMinus10perc.Cpools_litter_warming(2:end,9);
enz_rStrat_ff_moistureMinus10perc_warming = moistureMinus10perc.Cpools_litter_warming(2:end,10);
enz_kStrat_ff_moistureMinus10perc_warming = moistureMinus10perc.Cpools_litter_warming(2:end,11);
totalCarbon_ff_moistureMinus10perc_warming = sum(moistureMinus10perc.Cpools_litter_warming(2:end,[1:4 6:13]),2);

% % Forest floor - Enzyme rigidity - control
adsorbedC_ff_moisturePlus10perc_control = moisturePlus10perc.Cpools_litter_control(2:end,7) + moisturePlus10perc.Cpools_litter_control(2:end,12) + moisturePlus10perc.Cpools_litter_control(2:end,13);
metabolic_ff_moisturePlus10perc_control = moisturePlus10perc.Cpools_litter_control(2:end,8);
structural_ff_moisturePlus10perc_control = moisturePlus10perc.Cpools_litter_control(2:end,9);
enz_rStrat_ff_moisturePlus10perc_control = moisturePlus10perc.Cpools_litter_control(2:end,10);
enz_kStrat_ff_moisturePlus10perc_control = moisturePlus10perc.Cpools_litter_control(2:end,11);
totalCarbon_ff_moisturePlus10perc_control = sum(moisturePlus10perc.Cpools_litter_control(2:end,[1:4 6:13]),2);

% Forest floor - Enzyme rigidity - warming
adsorbedC_ff_moisturePlus10perc_warming = moisturePlus10perc.Cpools_litter_warming(2:end,7) + moisturePlus10perc.Cpools_litter_warming(2:end,12) + moisturePlus10perc.Cpools_litter_warming(2:end,13);
metabolic_ff_moisturePlus10perc_warming = moisturePlus10perc.Cpools_litter_warming(2:end,8);
structural_ff_moisturePlus10perc_warming = moisturePlus10perc.Cpools_litter_warming(2:end,9);
enz_rStrat_ff_moisturePlus10perc_warming = moisturePlus10perc.Cpools_litter_warming(2:end,10);
enz_kStrat_ff_moisturePlus10perc_warming = moisturePlus10perc.Cpools_litter_warming(2:end,11);
totalCarbon_ff_moisturePlus10perc_warming = sum(moisturePlus10perc.Cpools_litter_warming(2:end,[1:4 6:13]),2);

% =====================
% The plots are created
% =====================

dotSize = 8;
yMin = 0;
yMax = 2500;

scale = 2;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 12*scale], 'color', [1 1 1])

% --------------------------------------
% Forest floor - moistureConstant - control
% --------------------------------------
subplot(3,2,1)
hold on

dates_all = moisturePlus10perc.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_moistureConstant_control,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_moistureConstant_control,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_moistureConstant_control, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_moistureConstant_control, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2105', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Constant moisture')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Forest floor - moistureMinus10perc - control
% --------------------------------------
subplot(3,2,3)
hold on

% dates_all = moistureMinus10perc.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_moistureMinus10perc_control, 'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_moistureMinus10perc_control,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_moistureMinus10perc_control, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_moistureMinus10perc_control, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Moisture -10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Forest floor - moisturePlus10perc - control
% --------------------------------------
subplot(3,2,5)
hold on

% dates_all = noAdaptData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_moisturePlus10perc_control,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_moisturePlus10perc_control,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_moisturePlus10perc_control, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_moisturePlus10perc_control, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Moisture +10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Forest floor - moistureConstant - warming
% --------------------------------------
subplot(3,2,2)
hold on

% dates_all = noAdaptData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_moistureConstant_warming,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_moistureConstant_warming,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_moistureConstant_warming, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_moistureConstant_warming, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Constant moisture')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Forest floor - moistureMinus10perc - warming
% --------------------------------------
subplot(3,2,4)
hold on

% dates_all = optimumDrivenData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_moistureMinus10perc_warming,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_moistureMinus10perc_warming,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_moistureMinus10perc_warming, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_moistureMinus10perc_warming, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Moisture -10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Forest floor - moisturePlus10perc - warming
% --------------------------------------
subplot(3,2,6)
hold on

% dates_all = noAdaptData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_moisturePlus10perc_warming,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_moisturePlus10perc_warming,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_moisturePlus10perc_warming, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_moisturePlus10perc_warming, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Moisture +10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
text(0.38,0.98,'Forest floor organic carbon stocks','fontweight','bold','fontsize',14)

% Letters
text(0.04, 0.95, '(a)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.65, '(c)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.35, '(e)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.95, '(b)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.65, '(d)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.35, '(f)', 'FontSize', 14, 'FontWeight', 'bold')

% The legend
l = legend([p4, p3, p2, p1], 'Total C', 'Structural C', 'Metabolic C', 'Adsorbed C', 'Orientation', 'Horizontal', 'FontSIze', 14);
pos = l.Position;
pos(1) = 0.3; % x
pos(2) = 0.02; % y
l.Position = pos;
legend('boxoff')

% =========================================================================

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Forest floor carbon stocks.tiff' -r300 -nocrop -transparent
% =========================================================================

%% The OC pools are plotted for the SOIL

% The colors are defined
color_monomers_soil = [140,107,177]./255;
color_adsorbedC_soil = [140,107,177]./255;
color_polymers_soil = [153,52,4]./255;
color_totalCarbon_soil = 'k';

lineWidth = 2;
startWarming = datetime('01/01/2000', 'format', 'dd/MM/yyyy');

% =====================
% The data is formatted
% =====================

% Soil - No adaptation - control
adsorbedC_soil_moistureConstant_control = moistureConstant.Cpools_soil_control(2:end,5) + moistureConstant.Cpools_soil_control(2:end,8);
monomers_soil_moistureConstant_control = moistureConstant.Cpools_soil_control(2:end,4);
polymers_soil_moistureConstant_control = moistureConstant.Cpools_soil_control(2:end,6);
totalCarbon_soil_moistureConstant_control = sum(moistureConstant.Cpools_soil_control(2:end,[1 2 4 5 6 7 8]),2);

% Soil - No adaptation - warming
adsorbedC_soil_moistureConstant_warming = moistureConstant.Cpools_soil_warming(2:end,5) + moistureConstant.Cpools_soil_warming(2:end,8);
monomers_soil_moistureConstant_warming = moistureConstant.Cpools_soil_warming(2:end,4);
polymers_soil_moistureConstant_warming = moistureConstant.Cpools_soil_warming(2:end,6);
totalCarbon_soil_moistureConstant_warming = sum(moistureConstant.Cpools_soil_warming(2:end,[1 2 4 5 6 7 8]),2);

% Soil - Optimum driven - control
adsorbedC_soil_moistureMinus10perc_control = moistureMinus10perc.Cpools_soil_control(2:end,5) + moistureMinus10perc.Cpools_soil_control(2:end,8);
monomers_soil_moistureMinus10perc_control = moistureMinus10perc.Cpools_soil_control(2:end,4);
polymers_soil_moistureMinus10perc_control = moistureMinus10perc.Cpools_soil_control(2:end,6);
totalCarbon_soil_moistureMinus10perc_control = sum(moistureMinus10perc.Cpools_soil_control(2:end,[1 2 4 5 6 7 8]),2);

% Soil - Optimum driven - warming
adsorbedC_soil_moistureMinus10perc_warming = moistureMinus10perc.Cpools_soil_warming(2:end,5) + moistureMinus10perc.Cpools_soil_warming(2:end,8);
monomers_soil_moistureMinus10perc_warming = moistureMinus10perc.Cpools_soil_warming(2:end,4);
polymers_soil_moistureMinus10perc_warming = moistureMinus10perc.Cpools_soil_warming(2:end,6);
totalCarbon_soil_moistureMinus10perc_warming = sum(moistureMinus10perc.Cpools_soil_warming(2:end,[1 2 4 5 6 7 8]),2);

% Soil - Enzyme rigidity - control
adsorbedC_soil_moisturePlus10perc_control = moisturePlus10perc.Cpools_soil_control(2:end,5) + moisturePlus10perc.Cpools_soil_control(2:end,8);
monomers_soil_moisturePlus10perc_control = moisturePlus10perc.Cpools_soil_control(2:end,4);
polymers_soil_moisturePlus10perc_control = moisturePlus10perc.Cpools_soil_control(2:end,6);
totalCarbon_soil_moisturePlus10perc_control = sum(moisturePlus10perc.Cpools_soil_control(2:end,[1 2 4 5 6 7 8]),2);

% Soil - Enzyme rigidity - warming
adsorbedC_soil_moisturePlus10perc_warming = moisturePlus10perc.Cpools_soil_warming(2:end,5) + moisturePlus10perc.Cpools_soil_warming(2:end,8);
monomers_soil_moisturePlus10perc_warming = moisturePlus10perc.Cpools_soil_warming(2:end,4);
polymers_soil_moisturePlus10perc_warming = moisturePlus10perc.Cpools_soil_warming(2:end,6);
totalCarbon_soil_moisturePlus10perc_warming = sum(moisturePlus10perc.Cpools_soil_warming(2:end,[1 2 4 5 6 7 8]),2);

% =====================
% The plots are created
% =====================

dotSize = 8;
yMin = 0;
yMax = 6000;

scale = 2;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 12*scale], 'color', [1 1 1])

% --------------------------------------
% Soil - moistureConstant - control
% --------------------------------------
subplot(3,2,1)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_moistureConstant_control,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_moistureConstant_control, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_moistureConstant_control, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2105', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Constant moisture')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - moistureMinus10perc - control
% --------------------------------------
subplot(3,2,3)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_moistureMinus10perc_control,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_moistureMinus10perc_control, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_moistureMinus10perc_control, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Moisture -10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:1000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - moisturePlus10perc - control
% --------------------------------------
subplot(3,2,5)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_moisturePlus10perc_control,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_moisturePlus10perc_control, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_moisturePlus10perc_control, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Moisture +10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:1000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - moistureConstant - warming
% --------------------------------------
subplot(3,2,2)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_moistureConstant_warming,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_moistureConstant_warming, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_moistureConstant_warming, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Constant moisture')

set(gca, 'FontSize', 12, 'Ytick', [yMin:1000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - moistureMinus10perc - warming
% --------------------------------------
subplot(3,2,4)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_moistureMinus10perc_warming,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_moistureMinus10perc_warming, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_moistureMinus10perc_warming, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Moisture -10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:1000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - moisturePlus10perc - warming
% --------------------------------------
subplot(3,2,6)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_moisturePlus10perc_warming,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_moisturePlus10perc_warming, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_moisturePlus10perc_warming, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Moisture +10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:1000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
text(0.41,0.98,'Soil organic carbon stocks','fontweight','bold','fontsize',14)

% Letters
text(0.04, 0.95, '(a)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.65, '(c)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.35, '(e)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.95, '(b)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.65, '(d)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.35, '(f)', 'FontSize', 14, 'FontWeight', 'bold')

% The legend
l = legend([p3, p2, p1], 'Total C', 'Polymeric C', 'Adsorbed C', 'Orientation', 'Horizontal', 'FontSIze', 14);
pos = l.Position;
pos(1) = 0.3; % x
pos(2) = 0.02; % y
l.Position = pos;
legend('boxoff')
% =========================================================================

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Soil carbon stocks.tiff' -r300 -nocrop -transparent
% =========================================================================

%% Plotting total, forest floor and soil C stocks

% ---------------------
% The data is formatted
% ---------------------

% The modelled total C stocks
totalC_moistureConstant_control = totalCarbon_ff_moistureConstant_control + totalCarbon_soil_moistureConstant_control;
totalC_moistureConstant_warming = totalCarbon_ff_moistureConstant_warming + totalCarbon_soil_moistureConstant_warming;

totalC_moistureMinus10perc_control = totalCarbon_ff_moistureMinus10perc_control + totalCarbon_soil_moistureMinus10perc_control;
totalC_moistureMinus10perc_warming = totalCarbon_ff_moistureMinus10perc_warming + totalCarbon_soil_moistureMinus10perc_warming;

totalC_moisturePlus10perc_control = totalCarbon_ff_moisturePlus10perc_control + totalCarbon_soil_moisturePlus10perc_control;
totalC_moisturePlus10perc_warming = totalCarbon_ff_moisturePlus10perc_warming + totalCarbon_soil_moisturePlus10perc_warming;

% Measurements
totalC_noAdapt_control_meas = Cstock_forestFloor_control + Cstock_soil_control;
totalC_noAdapt_heated_meas = Cstock_forestFloor_warming + Cstock_soil_warming;

totalC_noAdapt_control_meas_stdError = sqrt((CstockStdError_forestFloor_control^2) + (CstockStdError_soil_control^2));
totalC_noAdapt_heated_meas_stdError = sqrt((CstockStdError_forestFloor_warmed^2) + (CstockStdError_soil_warmed^2));

% =====================
% The plots are created
% =====================

color_total = 'k'; 
color_ff = [35,139,69]./255;
color_soil = [204,76,2]./255;

dotSize = 8;
yMin = 0;
yMax = 8000;

scale = 2;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 12*scale], 'color', [1 1 1])

% --------------------------------------
% Soil - moistureConstant - control
% --------------------------------------
subplot(3,2,1)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_moistureConstant_control,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_moistureConstant_control, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_moistureConstant_control, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2105', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Constant moisture')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - moistureMinus10perc - control
% --------------------------------------
subplot(3,2,3)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_moistureMinus10perc_control,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_moistureMinus10perc_control, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_moistureMinus10perc_control, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Moisture -10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - moisturePlus10perc - control
% --------------------------------------
subplot(3,2,5)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_moisturePlus10perc_control,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_moisturePlus10perc_control, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_moisturePlus10perc_control, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Control - Moisture +10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - moistureConstant - warming
% --------------------------------------
subplot(3,2,2)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_moistureConstant_warming,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_moistureConstant_warming, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_moistureConstant_warming, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Constant moisture')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.72; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - moistureMinus10perc - warming
% --------------------------------------
subplot(3,2,4)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_moistureMinus10perc_warming,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_moistureMinus10perc_warming, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_moistureMinus10perc_warming, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Moisture -10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.42; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - moisturePlus10perc - warming
% --------------------------------------
subplot(3,2,6)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_moisturePlus10perc_warming,'color', color_ff, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),totalCarbon_soil_moisturePlus10perc_warming, 'color', color_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalC_moisturePlus10perc_warming, 'color', color_total, 'LineWidth', lineWidth); % Adsorbed enzymes

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Heated - Moisture +10 %')

set(gca, 'FontSize', 12, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.55; % x
pos(2) = 0.12; % y
pos(3) = 0.42; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
text(0.41,0.98,'Organic carbon stocks','fontweight','bold','fontsize',14)

% Letters
text(0.04, 0.95, '(a)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.65, '(c)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.04, 0.35, '(e)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.95, '(b)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.65, '(d)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.53, 0.35, '(f)', 'FontSize', 14, 'FontWeight', 'bold')

% The legend
l = legend([p3, p2, p1], 'Total C', 'Soil organic C', 'Forest floor C', 'Orientation', 'Horizontal', 'FontSIze', 14);
pos = l.Position;
pos(1) = 0.3; % x
pos(2) = 0.02; % y
l.Position = pos;
legend('boxoff')
% =========================================================================

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Total carbon stocks.tiff' -r300 -nocrop -transparent
% =========================================================================

%% Cumulative CO2 emissions are plotted

% The colors are defined
c_noAdapt = [31,120,180]./255;
c_optDriv = [51,160,44]./255;
c_enzRig = [255,127,0]./255;

yMin = -100;
yMax = 1600;

lineWidth = 2;
scale = 2;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 8*scale 8*scale], 'color', [1 1 1])
hold on

p1 = plot(moistureConstant.dates_after1991, moistureConstant.diffCumul, 'color', c_noAdapt, 'linewidth', lineWidth);
p2 = plot(moistureMinus10perc.dates_after1991, moistureMinus10perc.diffCumul, 'color', c_optDriv, 'linewidth', lineWidth);
p3 = plot(moisturePlus10perc.dates_after1991, moisturePlus10perc.diffCumul, 'color', c_enzRig, 'linewidth', lineWidth);

lgd = legend([p1 p3 p2], 'Constant moisture', 'Moisture +10 %', 'Moisture -10 %', 'location', 'northwest');%, 'FontSize', 15)
% lgd = legend([p1 p2 p4], 'No adaptation', 'Optimum driven', 'Measured', 'location', 'northwest');%, 'FontSize', 15)
set(lgd,'FontSize',12);
legend('boxoff')

ylim([yMin yMax])

xlabel('Year')
ylabel('g CO_{2}-C')
title('Difference in cumulative CO_{2} emisions between control and heated')

set(gca, 'FontSize', 12, 'Ytick', [yMin:100:yMax])

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Cumulative differences in CO2 fluxes.tiff' -r300 -nocrop -transparent
% =========================================================================

%% Annual CO2 emissions

% The colors are defined
c_noAdapt = [31,120,180]./255;
c_optDriv = [51,160,44]./255;
c_enzRig = [255,127,0]./255;
color_diff = [.5 .5 .5];%[1,102,94]./255;

lineWidth = .5;
scale = 2;

yMin = -20;
yMax = 30;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 8*scale], 'color', [1 1 1])
hold on

p0 = plot([2003 2105],[0 0], '--', 'color', [.2 .2 .2]);

% Let's add a datapoint for 2002, when control and heated fluxes are equal
yrs = [moistureConstant.uniqueYears_CO2diff(1)-1 moistureConstant.uniqueYears_CO2diff];

% Modelled
p1 = plot(yrs, [0; moistureConstant.diffAnnualCO2], '--o', 'color', c_noAdapt, 'linewidth', lineWidth, 'MarkerSize', 6, 'MarkerEdgeColor', c_noAdapt,'MarkerFaceColor', c_noAdapt);
p2 = plot(yrs, [0; moistureMinus10perc.diffAnnualCO2], '--o', 'color', c_optDriv, 'linewidth', lineWidth, 'MarkerSize', 6, 'MarkerEdgeColor', c_optDriv,'MarkerFaceColor', c_optDriv);
p3 = plot(yrs, [0; moisturePlus10perc.diffAnnualCO2], '--o', 'color', c_enzRig, 'linewidth', lineWidth, 'MarkerSize', 6, 'MarkerEdgeColor', c_enzRig,'MarkerFaceColor', c_enzRig);

% p1 = scatter(noAdaptData.uniqueYears_CO2diff,noAdaptData.diffAnnualCO2, 'MarkerFaceColor', c_noAdapt);
% p2 = scatter(optimumDrivenData.uniqueYears_CO2diff,optimumDrivenData.diffAnnualCO2, 'MarkerFaceColor', c_optDriv);
     
lgd = legend([p1 p3 p2], 'Constant moisture', 'Moisture +10 %', 'Moisture -10 %', 'location', 'northwest');%, 'FontSize', 15)
% lgd = legend([b1 p1 p2], 'Measured (Melillo et al., 2011)', 'No adaptation', 'Optimum driven', 'location', 'southwest');%, 'FontSize', 15)
set(lgd,'FontSize',12);
legend('boxoff')

ylim([yMin yMax])

xlim([2001 2105])
% ylim([0 yMax])
set(gca, 'FontSize', 12, 'Ytick', [yMin:10:yMax])

xlabel('Year')
ylabel('g CO_{2}-C yr^{-1}')
title('Difference in annual CO_{2} emissions between control and heated')

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Annual differences in CO2 fluxes.tiff' -r300 -nocrop -transparent
% =========================================================================

















