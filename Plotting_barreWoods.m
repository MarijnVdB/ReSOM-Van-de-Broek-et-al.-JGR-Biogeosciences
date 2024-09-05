%% This script is used to create the plots used in the manuscript

clc; clearvars; close all

%% The data is loaded

mainFolder = 'Data Barre Woods/Output 5 June 2024/';

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
adsorbedC_ff_optimumDriv_warming = optimumDrivenData.Cpools_litter_warming(2:end,7) + noAdaptData.Cpools_litter_warming(2:end,12) + noAdaptData.Cpools_litter_warming(2:end,13);
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

dates_all = noAdaptData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_noAdapt_control,'color', color_monomers_litter, 'LineWidth', lineWidth-1);
p2 = plot(datenum(dates_all),metabolic_ff_noAdapt_control,'color', color_metabolic_litter, 'LineWidth', lineWidth-1);
p3 = plot(datenum(dates_all),structural_ff_noAdapt_control, 'color', color_structural_litter, 'LineWidth', lineWidth-1);
p4 = plot(datenum(dates_all),totalCarbon_ff_noAdapt_control, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth);

% Measurements
d = datenum(date_SOC_measurement);
s1 = errorbar(d, Cstock_forestFloor_control, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');
s2 = errorbar(d, Cstock_forestFloor_control*0.07, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_monomers_litter);
s3 = errorbar(d, Cstock_forestFloor_control*0.93*0.66, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_metabolic_litter);
s4 = errorbar(d, Cstock_forestFloor_control*0.93*0.34, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_structural_litter);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

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

% Measurements
s1 = errorbar(d, Cstock_forestFloor_control, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');
s2 = errorbar(d, Cstock_forestFloor_control*0.07, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_monomers_litter);
s3 = errorbar(d, Cstock_forestFloor_control*0.93*0.66, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_metabolic_litter);
s4 = errorbar(d, Cstock_forestFloor_control*0.93*0.34, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_structural_litter);

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

dates_all = noAdaptData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_enzRig_control,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_enzRig_control,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_enzRig_control, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_enzRig_control, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

% Measurements
s1 = errorbar(d, Cstock_forestFloor_control, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');
s2 = errorbar(d, Cstock_forestFloor_control*0.07, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_monomers_litter);
s3 = errorbar(d, Cstock_forestFloor_control*0.93*0.66, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_metabolic_litter);
s4 = errorbar(d, Cstock_forestFloor_control*0.93*0.34, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_structural_litter);

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

dates_all = noAdaptData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_noAdapt_warming,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_noAdapt_warming,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_noAdapt_warming, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_noAdapt_warming, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

% Measurements
s1 = errorbar(d, Cstock_forestFloor_warming, CstockStdError_forestFloor_warmed, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');

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

% Measurements
s1 = errorbar(d, Cstock_forestFloor_warming, CstockStdError_forestFloor_warmed, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');

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

dates_all = noAdaptData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_enzRig_warming,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_enzRig_warming,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_enzRig_warming, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_enzRig_warming, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

% Measurements
s1 = errorbar(d, Cstock_forestFloor_warming, CstockStdError_forestFloor_warmed, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');

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

%% The OC pools are plotted for the FOREST FLOOR - COMBINED

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
adsorbedC_ff_optimumDriv_warming = optimumDrivenData.Cpools_litter_warming(2:end,7) + noAdaptData.Cpools_litter_warming(2:end,12) + noAdaptData.Cpools_litter_warming(2:end,13);
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
% Forest floor - no adaptation
% --------------------------------------
subplot(3,1,1)
hold on

dates_all = noAdaptData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_noAdapt_control,'color', color_monomers_litter, 'LineWidth', lineWidth-1);
p2 = plot(datenum(dates_all),metabolic_ff_noAdapt_control,'color', color_metabolic_litter, 'LineWidth', lineWidth-1);
p3 = plot(datenum(dates_all),structural_ff_noAdapt_control, 'color', color_structural_litter, 'LineWidth', lineWidth-1);
p4 = plot(datenum(dates_all),totalCarbon_ff_noAdapt_control, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth);

p5 = plot(datenum(dates_all),adsorbedC_ff_noAdapt_warming,'color', color_monomers_litter, 'LineWidth', lineWidth-1, 'LineStyle', '--');
p6 = plot(datenum(dates_all),metabolic_ff_noAdapt_warming,'color', color_metabolic_litter, 'LineWidth', lineWidth-1, 'LineStyle', '--');
p7 = plot(datenum(dates_all),structural_ff_noAdapt_warming, 'color', color_structural_litter, 'LineWidth', lineWidth-1, 'LineStyle', '--');
p8 = plot(datenum(dates_all),totalCarbon_ff_noAdapt_warming, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth, 'LineStyle', '--');

% Measurements
d = datenum(date_SOC_measurement);
s1 = errorbar(d, Cstock_forestFloor_control, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');
s2 = errorbar(d, Cstock_forestFloor_control*0.07, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_monomers_litter);
s3 = errorbar(d, Cstock_forestFloor_control*0.93*0.66, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_metabolic_litter);
s4 = errorbar(d, Cstock_forestFloor_control*0.93*0.34, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_structural_litter);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2002', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('No thermal adaptation')

set(gca, 'FontSize', 14, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.72; % y
pos(3) = 0.9; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Forest floor - Optimum driven
% --------------------------------------
subplot(3,1,2)
hold on

dates_all = optimumDrivenData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_optimumDriv_control, 'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_optimumDriv_control,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_optimumDriv_control, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_optimumDriv_control, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

p5 = plot(datenum(dates_all),adsorbedC_ff_optimumDriv_warming, 'color', color_monomers_litter, 'LineWidth', lineWidth-1, 'LineStyle', '--'); % Polymers
p6 = plot(datenum(dates_all),metabolic_ff_optimumDriv_warming,'color', color_metabolic_litter, 'LineWidth', lineWidth-1, 'LineStyle', '--'); % Enzymes
p7 = plot(datenum(dates_all),structural_ff_optimumDriv_warming, 'color', color_structural_litter, 'LineWidth', lineWidth-1, 'LineStyle', '--'); % Adsorbed enzymes
p8 = plot(datenum(dates_all),totalCarbon_ff_optimumDriv_warming, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth, 'LineStyle', '--'); % Adsorbed enzymes

% Measurements
s1 = errorbar(d, Cstock_forestFloor_control, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');
s2 = errorbar(d, Cstock_forestFloor_control*0.07, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_monomers_litter);
s3 = errorbar(d, Cstock_forestFloor_control*0.93*0.66, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_metabolic_litter);
s4 = errorbar(d, Cstock_forestFloor_control*0.93*0.34, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_structural_litter);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Optimum driven adaptation')

set(gca, 'FontSize', 14, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.43; % y
pos(3) = 0.9; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Forest floor - enzyme rigidity
% --------------------------------------
subplot(3,1,3)
hold on

dates_all = noAdaptData.dates_all;

% Model results
p1 = plot(datenum(dates_all),adsorbedC_ff_enzRig_control,'color', color_monomers_litter, 'LineWidth', lineWidth-1); % Polymers
p2 = plot(datenum(dates_all),metabolic_ff_enzRig_control,'color', color_metabolic_litter, 'LineWidth', lineWidth-1); % Enzymes
p3 = plot(datenum(dates_all),structural_ff_enzRig_control, 'color', color_structural_litter, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p4 = plot(datenum(dates_all),totalCarbon_ff_enzRig_control, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

p5 = plot(datenum(dates_all),adsorbedC_ff_enzRig_warming,'color', color_monomers_litter, 'LineWidth', lineWidth-1, 'LineStyle', '--'); % Polymers
p6 = plot(datenum(dates_all),metabolic_ff_enzRig_warming,'color', color_metabolic_litter, 'LineWidth', lineWidth-1, 'LineStyle', '--'); % Enzymes
p7 = plot(datenum(dates_all),structural_ff_enzRig_warming, 'color', color_structural_litter, 'LineWidth', lineWidth-1, 'LineStyle', '--'); % Adsorbed enzymes
p8 = plot(datenum(dates_all),totalCarbon_ff_enzRig_warming, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth, 'LineStyle', '--'); % Adsorbed enzymes

% Measurements
s1 = errorbar(d, Cstock_forestFloor_control, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');
s2 = errorbar(d, Cstock_forestFloor_control*0.07, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_monomers_litter);
s3 = errorbar(d, Cstock_forestFloor_control*0.93*0.66, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_metabolic_litter);
s4 = errorbar(d, Cstock_forestFloor_control*0.93*0.34, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_structural_litter);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 3000], ':', 'color', 'k', 'lineWidth', 1);

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Enzyme rigidity adaptation')

set(gca, 'FontSize', 14, 'Ytick', [yMin:500:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.13; % y
pos(3) = 0.9; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
text(0.38,0.98,'Organic horizon organic carbon stocks','fontweight','bold','fontsize',18)

% Letters
text(0.01, 0.94, '(A)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.01, 0.65, '(B)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.01, 0.35, '(C)', 'FontSize', 14, 'FontWeight', 'bold')

% The legend
l = legend([p4, p3, p2, p1, p8, p7, p6, p5], 'Total C - control', 'Structural C - control', 'Metabolic C - control', 'Adsorbed C - control', ...
            'Total C - heated', 'Structural C - heated', 'Metabolic C - heated', 'Adsorbed C - heated', 'Orientation', 'Horizontal', 'FontSize', 14, 'NumColumns', 4);
pos = l.Position;
pos(1) = 0.21; % x
pos(2) = 0.01; % y
l.Position = pos;
% legend('boxoff')
legend('box','off')

% =========================================================================

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Forest floor carbon stocks - combined.tiff' -r300 -nocrop -transparent
% =========================================================================

% ----------------------------
% The relative decrease in SOC
% ----------------------------

date1 = datetime('01/01/2000', 'format', 'dd/MM/yyyy');
date2 = datetime('31/12/2002', 'format', 'dd/MM/yyyy');

date3 = datetime('01/01/2014', 'format', 'dd/MM/yyyy');
date4 = datetime('31/12/2016', 'format', 'dd/MM/yyyy');

r1 = find(dates_all == date1);
r2 = find(dates_all == date2);

r3 = find(dates_all == date3);
r4 = find(dates_all == date4);

% No adaptation
C1 = mean(adsorbedC_ff_noAdapt_warming(r1:r2));
C2 = mean(adsorbedC_ff_noAdapt_warming(r3:r4));
dC_absorbed_noAdapt = (C1 - C2)/C1;

C3 = mean(metabolic_ff_noAdapt_warming(r1:r2));
C4 = mean(metabolic_ff_noAdapt_warming(r3:r4));
dC_metabolic_noAdapt = (C3 - C4)/C3;

C5 = mean(structural_ff_noAdapt_warming(r1:r2));
C6 = mean(structural_ff_noAdapt_warming(r3:r4));
dC_structural_noAdapt = (C5 - C6)/C5;

% Optimum driven
C1 = mean(adsorbedC_ff_optimumDriv_warming(r1:r2));
C2 = mean(adsorbedC_ff_optimumDriv_warming(r3:r4));
dC_absorbed_optDriv = (C1 - C2)/C1;

C3 = mean(metabolic_ff_optimumDriv_warming(r1:r2));
C4 = mean(metabolic_ff_optimumDriv_warming(r3:r4));
dC_metabolic_optDriv = (C3 - C4)/C3;

C5 = mean(structural_ff_optimumDriv_warming(r1:r2));
C6 = mean(structural_ff_optimumDriv_warming(r3:r4));
dC_structural_optDriv = (C5 - C6)/C5;

% Enzyme rigidity
C1 = mean(adsorbedC_ff_enzRig_warming(r1:r2));
C2 = mean(adsorbedC_ff_enzRig_warming(r3:r4));
dC_absorbed_enzRig = (C1 - C2)/C1;

C3 = mean(metabolic_ff_enzRig_warming(r1:r2));
C4 = mean(metabolic_ff_enzRig_warming(r3:r4));
dC_metabolic_enzRig = (C3 - C4)/C3;

C5 = mean(structural_ff_enzRig_warming(r1:r2));
C6 = mean(structural_ff_enzRig_warming(r3:r4));
dC_structural_enzRig = (C5 - C6)/C5;

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
adsorbedC_soil_optimumDriv_control = optimumDrivenData.Cpools_soil_control(2:end,5) + noAdaptData.Cpools_soil_control(2:end,8);
monomers_soil_optimumDriv_control = optimumDrivenData.Cpools_soil_control(2:end,4);
polymers_soil_optimumDriv_control = optimumDrivenData.Cpools_soil_control(2:end,6);
totalCarbon_soil_optimumDriv_control = sum(optimumDrivenData.Cpools_soil_control(2:end,[1 2 4 5 6 7 8]),2);

% Soil - Optimum driven - warming
adsorbedC_soil_optimumDriv_warming = optimumDrivenData.Cpools_soil_warming(2:end,5) + noAdaptData.Cpools_soil_warming(2:end,8);
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
% The standard error for POC and MAOC is calculated
% =====================

stdError_data = 0.163;

stdError_MAOC = sqrt((CstockStdError_soil_control/Cstock_soil_control)^2+(stdError_data/0.58)^2) * Cstock_soil_control.*0.58;
stdError_POC = sqrt((CstockStdError_soil_control/Cstock_soil_control)^2+(stdError_data/0.42)^2) * Cstock_soil_control.*0.42;

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

% Measurements
d = datenum(date_SOC_measurement);
s1 = errorbar(d, Cstock_soil_control, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');
s2 = errorbar(d, Cstock_soil_control.*0.58, stdError_MAOC, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_adsorbedC_soil);
s3 = errorbar(d, Cstock_soil_control.*0.42, stdError_POC, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_polymers_soil);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

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

% Measurements
s1 = errorbar(d, Cstock_soil_control, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');
s2 = errorbar(d, Cstock_soil_control.*0.57, stdError_MAOC, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_adsorbedC_soil);
s3 = errorbar(d, Cstock_soil_control.*0.43, stdError_POC, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_polymers_soil);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

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

% Measurements
s1 = errorbar(d, Cstock_soil_control, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');
s2 = errorbar(d, Cstock_soil_control.*0.57, stdError_MAOC, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_adsorbedC_soil);
s3 = errorbar(d, Cstock_soil_control.*0.43, stdError_POC, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_polymers_soil);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

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

% Measurements
s1 = errorbar(d, Cstock_soil_warming, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

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

% Measurements
s1 = errorbar(d, Cstock_soil_warming, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

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

% Measurements
s1 = errorbar(d, Cstock_soil_warming, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

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
pos(1) = 0.35; % x
pos(2) = 0.02; % y
l.Position = pos;
legend('boxoff')
% =========================================================================

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Soil carbon stocks.tiff' -r300 -nocrop -transparent
% =========================================================================

%% The OC pools are plotted for the SOIL - COMBINED

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
adsorbedC_soil_optimumDriv_control = optimumDrivenData.Cpools_soil_control(2:end,5) + noAdaptData.Cpools_soil_control(2:end,8);
monomers_soil_optimumDriv_control = optimumDrivenData.Cpools_soil_control(2:end,4);
polymers_soil_optimumDriv_control = optimumDrivenData.Cpools_soil_control(2:end,6);
totalCarbon_soil_optimumDriv_control = sum(optimumDrivenData.Cpools_soil_control(2:end,[1 2 4 5 6 7 8]),2);

% Soil - Optimum driven - warming
adsorbedC_soil_optimumDriv_warming = optimumDrivenData.Cpools_soil_warming(2:end,5) + noAdaptData.Cpools_soil_warming(2:end,8);
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
% The standard error for POC and MAOC is calculated
% =====================

stdError_data = 0.163;

stdError_MAOC = sqrt((CstockStdError_soil_control/Cstock_soil_control)^2+(stdError_data/0.58)^2) * Cstock_soil_control.*0.58;
stdError_POC = sqrt((CstockStdError_soil_control/Cstock_soil_control)^2+(stdError_data/0.42)^2) * Cstock_soil_control.*0.42;

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
% Soil - no adaptation
% --------------------------------------
subplot(3,1,1)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_noAdapt_control,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_noAdapt_control, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_noAdapt_control, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

p4 = plot(datenum(dates_all),adsorbedC_soil_noAdapt_warming,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1, 'LineStyle', '--'); % Enzymes
p5 = plot(datenum(dates_all),polymers_soil_noAdapt_warming, 'color', color_polymers_soil, 'LineWidth', lineWidth-1, 'LineStyle', '--'); % Adsorbed enzymes
p6 = plot(datenum(dates_all),totalCarbon_soil_noAdapt_warming, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth, 'LineStyle', '--'); % Adsorbed enzymes

% Measurements
d = datenum(date_SOC_measurement);
s1 = errorbar(d, Cstock_soil_control, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');
s2 = errorbar(d, Cstock_soil_control.*0.58, stdError_MAOC, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_adsorbedC_soil);
s3 = errorbar(d, Cstock_soil_control.*0.42, stdError_POC, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_polymers_soil);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2002', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('No thermal adaptation')

set(gca, 'FontSize', 14, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.72; % y
pos(3) = 0.9; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - Optimum driven
% --------------------------------------
subplot(3,1,2)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_optimumDriv_control,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_optimumDriv_control, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_optimumDriv_control, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

p4 = plot(datenum(dates_all),adsorbedC_soil_optimumDriv_warming,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1, 'LineStyle', '--'); % Enzymes
p5 = plot(datenum(dates_all),polymers_soil_optimumDriv_warming, 'color', color_polymers_soil, 'LineWidth', lineWidth-1, 'LineStyle', '--'); % Adsorbed enzymes
p6 = plot(datenum(dates_all),totalCarbon_soil_optimumDriv_warming, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth, 'LineStyle', '--'); % Adsorbed enzymes

% Measurements
s1 = errorbar(d, Cstock_soil_control, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');
s2 = errorbar(d, Cstock_soil_control.*0.57, stdError_MAOC, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_adsorbedC_soil);
s3 = errorbar(d, Cstock_soil_control.*0.43, stdError_POC, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_polymers_soil);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2002', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Optimum driven adaptation')

set(gca, 'FontSize', 14, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.43; % y
pos(3) = 0.9; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - Enzyme rigidity - control
% --------------------------------------
subplot(3,1,3)
hold on

% Model results
p1 = plot(datenum(dates_all),adsorbedC_soil_enzRig_control,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1); % Enzymes
p2 = plot(datenum(dates_all),polymers_soil_enzRig_control, 'color', color_polymers_soil, 'LineWidth', lineWidth-1); % Adsorbed enzymes
p3 = plot(datenum(dates_all),totalCarbon_soil_enzRig_control, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

p4 = plot(datenum(dates_all),adsorbedC_soil_enzRig_warming,'color', color_adsorbedC_soil, 'LineWidth', lineWidth-1, 'LineStyle', '--'); % Enzymes
p5 = plot(datenum(dates_all),polymers_soil_enzRig_warming, 'color', color_polymers_soil, 'LineWidth', lineWidth-1, 'LineStyle', '--'); % Adsorbed enzymes
p6 = plot(datenum(dates_all),totalCarbon_soil_enzRig_warming, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth, 'LineStyle', '--'); % Adsorbed enzymes

% Measurements
s1 = errorbar(d, Cstock_soil_control, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', 'k');
s2 = errorbar(d, Cstock_soil_control.*0.57, stdError_MAOC, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_adsorbedC_soil);
s3 = errorbar(d, Cstock_soil_control.*0.43, stdError_POC, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_polymers_soil);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 6000], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2002', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

xlabel('Year')
ylabel('g C m^{-2}')
title('Enzyme rigidity adaptation')

set(gca, 'FontSize', 14, 'Ytick', [yMin:2000:yMax])

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.13; % y
pos(3) = 0.9; % width
pos(4) = 0.19; % height
set(gca,'position',pos)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
text(0.38,0.98,'Soil organic carbon stocks','fontweight','bold','fontsize',18)

% Letters
text(0.01, 0.94, '(A)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.01, 0.65, '(B)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.01, 0.35, '(C)', 'FontSize', 14, 'FontWeight', 'bold')

% The legend
l = legend([p3, p2, p1, p6, p5, p4], 'Total C - control', 'Polymeric C - control', 'Adsorbed C - control', ...
            'Total C - heated', 'Polymeric C - heated', 'Adsorbed C - heated', 'Orientation', 'Horizontal', 'FontSize', 14, 'NumColumns', 3);
pos = l.Position;
pos(1) = 0.3; % x
pos(2) = 0.01; % y
l.Position = pos;
% legend('boxoff')
legend('box','off')

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Soil carbon stocks - combined.tiff' -r300 -nocrop -transparent
% =========================================================================

% ----------------------------
% The relative decrease in SOC
% ----------------------------

date1 = datetime('01/01/2000', 'format', 'dd/MM/yyyy');
date2 = datetime('31/12/2002', 'format', 'dd/MM/yyyy');

date3 = datetime('01/01/2014', 'format', 'dd/MM/yyyy');
date4 = datetime('31/12/2016', 'format', 'dd/MM/yyyy');

r1 = find(dates_all == date1);
r2 = find(dates_all == date2);

r3 = find(dates_all == date3);
r4 = find(dates_all == date4);

% % No adaptation
% C1 = mean(adsorbedC_soil_noAdapt_warming(r1:r2));
% C2 = mean(adsorbedC_soil_noAdapt_warming(r3:r4));
% dC_absorbed_noAdapt = (C2 - C1)/C1;
% 
% C3 = mean(polymers_soil_noAdapt_warming(r1:r2));
% C4 = mean(polymers_soil_noAdapt_warming(r3:r4));
% dC_polymeric_noAdapt = (C4 - C3)/C3;
% 
% % Optimum driven
% C1 = mean(adsorbedC_soil_optimumDriv_warming(r1:r2));
% C2 = mean(adsorbedC_soil_optimumDriv_warming(r3:r4));
% dC_absorbed_optimumDriv = (C2 - C1)/C1;
% 
% C3 = mean(polymers_soil_optimumDriv_warming(r1:r2));
% C4 = mean(polymers_soil_optimumDriv_warming(r3:r4));
% dC_polymeric_optimumDriv = (C4 - C3)/C3;
% 
% % Enzyme rigidity
% C1 = mean(adsorbedC_soil_enzRig_warming(r1:r2));
% C2 = mean(adsorbedC_soil_enzRig_warming(r3:r4));
% dC_absorbed_enzRig = (C2 - C1)/C1;
% 
% C3 = mean(polymers_soil_enzRig_warming(r1:r2));
% C4 = mean(polymers_soil_enzRig_warming(r3:r4));
% dC_polymeric_enzRig = (C4 - C3)/C3;

% No adaptation
C1 = mean(adsorbedC_soil_noAdapt_control(r3:r4));
C2 = mean(adsorbedC_soil_noAdapt_warming(r3:r4));
dC_absorbed_noAdapt = (C2 - C1)/C1;

C3 = mean(polymers_soil_noAdapt_control(r3:r4));
C4 = mean(polymers_soil_noAdapt_warming(r3:r4));
dC_polymeric_noAdapt = (C4 - C3)/C3;

% Optimum driven
C1 = mean(adsorbedC_soil_optimumDriv_control(r3:r4));
C2 = mean(adsorbedC_soil_optimumDriv_warming(r3:r4));
dC_absorbed_optimumDriv = (C2 - C1)/C1;

C3 = mean(polymers_soil_optimumDriv_control(r3:r4));
C4 = mean(polymers_soil_optimumDriv_warming(r3:r4));
dC_polymeric_optimumDriv = (C4 - C3)/C3;

% Enzyme rigidity
C1 = mean(adsorbedC_soil_enzRig_control(r3:r4));
C2 = mean(adsorbedC_soil_enzRig_warming(r3:r4));
dC_absorbed_enzRig = (C2 - C1)/C1;

C3 = mean(polymers_soil_enzRig_control(r3:r4));
C4 = mean(polymers_soil_enzRig_warming(r3:r4));
dC_polymeric_enzRig = (C4 - C3)/C3;

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

% Measurements
d = datenum(date_SOC_measurement);
s1 = errorbar(d, Cstock_forestFloor_control, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_ff);
s2 = errorbar(d, Cstock_soil_control, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_soil);
s3 = errorbar(d, totalC_noAdapt_control_meas, totalC_noAdapt_control_meas_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_total);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

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

% Measurements
d = datenum(date_SOC_measurement);
s1 = errorbar(d, Cstock_forestFloor_control, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_ff);
s2 = errorbar(d, Cstock_soil_control, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_soil);
s3 = errorbar(d, totalC_noAdapt_control_meas, totalC_noAdapt_control_meas_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_total);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

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

% Measurements
d = datenum(date_SOC_measurement);
s1 = errorbar(d, Cstock_forestFloor_control, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_ff);
s2 = errorbar(d, Cstock_soil_control, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_soil);
s3 = errorbar(d, totalC_noAdapt_control_meas, totalC_noAdapt_control_meas_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_total);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

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

% Measurements
d = datenum(date_SOC_measurement);
s1 = errorbar(d, Cstock_forestFloor_warming, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_ff);
s2 = errorbar(d, Cstock_soil_warming, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_soil);
s3 = errorbar(d, totalC_noAdapt_heated_meas, totalC_noAdapt_control_meas_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_total);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

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

% Measurements
d = datenum(date_SOC_measurement);
s1 = errorbar(d, Cstock_forestFloor_warming, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_ff);
s2 = errorbar(d, Cstock_soil_warming, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_soil);
s3 = errorbar(d, totalC_noAdapt_heated_meas, totalC_noAdapt_control_meas_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_total);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

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

% Measurements
d = datenum(date_SOC_measurement);
s1 = errorbar(d, Cstock_forestFloor_warming, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_ff);
s2 = errorbar(d, Cstock_soil_warming, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_soil);
s3 = errorbar(d, totalC_noAdapt_heated_meas, totalC_noAdapt_control_meas_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_total);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], ':', 'color', 'k', 'lineWidth', 1);

date1 = datenum('01/01/2000', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

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

%% Plotting total, forest floor and soil C stocks - COMBINED

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
set(hfig, 'Units', 'centimeter', 'Position', [15 15 9*scale 12*scale], 'color', [1 1 1])

% --------------------------------------
% Soil - no adaptation
% --------------------------------------
subplot(3,1,1)
hold on

% Model results - Control
p1 = plot(datenum(dates_all),totalCarbon_ff_noAdapt_control,'color', color_ff, 'LineWidth', lineWidth-1);
p2 = plot(datenum(dates_all),totalCarbon_soil_noAdapt_control, 'color', color_soil, 'LineWidth', lineWidth-1);
p3 = plot(datenum(dates_all),totalC_noAdapt_control, 'color', color_total, 'LineWidth', lineWidth);

% Model results - warming
p4 = plot(datenum(dates_all),totalCarbon_ff_noAdapt_warming,'color', color_ff, 'LineWidth', lineWidth-1, 'linestyle', '--');
p5 = plot(datenum(dates_all),totalCarbon_soil_noAdapt_warming, 'color', color_soil, 'LineWidth', lineWidth-1, 'linestyle', '--');
p6 = plot(datenum(dates_all),totalC_noAdapt_warming, 'color', color_total, 'LineWidth', lineWidth, 'linestyle', '--');

% Measurements
d = datenum(date_SOC_measurement);
s1 = errorbar(d, Cstock_forestFloor_control, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_ff);
s2 = errorbar(d, Cstock_soil_control, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_soil);
s3 = errorbar(d, totalC_noAdapt_control_meas, totalC_noAdapt_control_meas_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_total);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], '--', 'color', 'k', 'lineWidth', 1);

d1 = datenum('01/01/2002', 'dd/MM/yyyy');
d2 = datenum('01/01/2004', 'dd/MM/yyyy');
d3 = datenum('01/01/2006', 'dd/MM/yyyy');
d4 = datenum('01/01/2008', 'dd/MM/yyyy');
d5 = datenum('01/01/2010', 'dd/MM/yyyy');
d6 = datenum('01/01/2012', 'dd/MM/yyyy');
d7 = datenum('01/01/2014', 'dd/MM/yyyy');
d8 = datenum('01/01/2016', 'dd/MM/yyyy');

xticks = datenum([d1 d2 d3 d4 d5 d6 d7 d8]);
set(gca, 'xtick', xticks);

date1 = datenum('01/01/2002', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits', 'keepticks')
ylim([yMin yMax])

xlabel('Year')
ylabel('gC m^{-2}')
title('No thermal adaptation')

set(gca, 'FontSize', 14, 'Ytick', [yMin:2000:yMax], 'box', 'on')

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.13; % x
pos(2) = 0.75; % y
pos(3) = 0.84; % width
pos(4) = 0.21; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - Optimum driven
% --------------------------------------
subplot(3,1,2)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_optimumDriv_control,'color', color_ff, 'LineWidth', lineWidth-1);
p2 = plot(datenum(dates_all),totalCarbon_soil_optimumDriv_control, 'color', color_soil, 'LineWidth', lineWidth-1);
p3 = plot(datenum(dates_all),totalC_optimumDriv_control, 'color', color_total, 'LineWidth', lineWidth);

p4 = plot(datenum(dates_all),totalCarbon_ff_optimumDriv_warming,'color', color_ff, 'LineWidth', lineWidth-1, 'linestyle', '--');
p5 = plot(datenum(dates_all),totalCarbon_soil_optimumDriv_warming, 'color', color_soil, 'LineWidth', lineWidth-1, 'linestyle', '--');
p6 = plot(datenum(dates_all),totalC_optimumDriv_warming, 'color', color_total, 'LineWidth', lineWidth, 'linestyle', '--');

% Measurements
d = datenum(date_SOC_measurement);
s1 = errorbar(d, Cstock_forestFloor_control, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_ff);
s2 = errorbar(d, Cstock_soil_control, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_soil);
s3 = errorbar(d, totalC_noAdapt_control_meas, totalC_noAdapt_control_meas_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_total);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], '--', 'color', 'k', 'lineWidth', 1);

xticks = datenum([d1 d2 d3 d4 d5 d6 d7 d8]);
set(gca, 'xtick', xticks);

date1 = datenum('01/01/2002', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits', 'keepticks')
ylim([yMin yMax])

xlabel('Year')
ylabel('gC m^{-2}')
title('Optimum driven adaptation')

set(gca, 'FontSize', 14, 'Ytick', [yMin:2000:yMax], 'box', 'on')

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.13; % x
pos(2) = 0.44; % y
pos(3) = 0.84; % width
pos(4) = 0.21; % height
set(gca,'position',pos)

% --------------------------------------
% Soil - Enzyme rigidity
% --------------------------------------
subplot(3,1,3)
hold on

% Model results
p1 = plot(datenum(dates_all),totalCarbon_ff_enzRig_control,'color', color_ff, 'LineWidth', lineWidth-1);
p2 = plot(datenum(dates_all),totalCarbon_soil_enzRig_control, 'color', color_soil, 'LineWidth', lineWidth-1);
p3 = plot(datenum(dates_all),totalC_enzRig_control, 'color', color_total, 'LineWidth', lineWidth);

p4 = plot(datenum(dates_all),totalCarbon_ff_enzRig_warming,'color', color_ff, 'LineWidth', lineWidth-1, 'linestyle', '--');
p5 = plot(datenum(dates_all),totalCarbon_soil_enzRig_warming, 'color', color_soil, 'LineWidth', lineWidth-1, 'linestyle', '--');
p6 = plot(datenum(dates_all),totalC_enzRig_warming, 'color', color_total, 'LineWidth', lineWidth, 'linestyle', '--');

% Measurements
d = datenum(date_SOC_measurement);
s1 = errorbar(d, Cstock_forestFloor_control, CstockStdError_forestFloor_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_ff);
s2 = errorbar(d, Cstock_soil_control, CstockStdError_soil_control, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_soil);
s3 = errorbar(d, totalC_noAdapt_control_meas, totalC_noAdapt_control_meas_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', color_total);

plot([datenum(date_startWarming) datenum(date_startWarming)],[0 yMax], '--', 'color', 'k', 'lineWidth', 1);

xticks = datenum([d1 d2 d3 d4 d5 d6 d7 d8]);
set(gca, 'xtick', xticks);

date1 = datenum('01/01/2002', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits', 'keepticks')
ylim([yMin yMax])

xlabel('Year')
ylabel('gC m^{-2}')
title('Enzyme rigidity adaptation')

set(gca, 'FontSize', 14, 'Ytick', [yMin:2000:yMax], 'box', 'on')

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.13; % x
pos(2) = 0.14; % y
pos(3) = 0.84; % width
pos(4) = 0.21; % height
set(gca,'position',pos)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
% text(0.37,0.98,'Organic carbon stocks','fontweight','bold','fontsize',14)

% Letters
text(0.02, 0.98, '(A)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.02, 0.67, '(B)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.02, 0.37, '(C)', 'FontSize', 14, 'FontWeight', 'bold')

% The legend
l = legend([p3, p2, p1, p6, p5, p4], 'Total C - control', 'Mineral soil C - control', 'Organic horizon C - control', ...
            'Total C - heated', 'Mineral soil C - heated', 'Organic horizon C - heated', 'Orientation', 'Horizontal', 'FontSize', 12,'NumColumns',3);
pos = l.Position;
pos(1) = 0.04; % x
pos(2) = 0.02; % y
l.Position = pos;
% legend('boxoff')
legend('box','off')
% =========================================================================

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Total carbon stocks - combined.tiff' -r300 -nocrop -transparent
% =========================================================================

% ----------------------------
% The relative decrease in SOC
% ----------------------------

date1 = datetime('01/01/2000', 'format', 'dd/MM/yyyy');
date2 = datetime('31/12/2002', 'format', 'dd/MM/yyyy');

date3 = datetime('01/01/2014', 'format', 'dd/MM/yyyy');
date4 = datetime('31/12/2016', 'format', 'dd/MM/yyyy');

r1 = find(dates_all == date1);
r2 = find(dates_all == date2);

r3 = find(dates_all == date3);
r4 = find(dates_all == date4);

% % No adaptation
% C1 = mean(totalCarbon_ff_noAdapt_warming(r1:r2));
% C2 = mean(totalCarbon_ff_noAdapt_warming(r3:r4));
% dC_ff_noAdapt = (C1 - C2)/C1;
% 
% C3 = mean(totalCarbon_soil_noAdapt_warming(r1:r2));
% C4 = mean(totalCarbon_soil_noAdapt_warming(r3:r4));
% dC_soil_noAdapt = (C3 - C4)/C3;
% 
% C5 = mean(totalCarbon_ff_noAdapt_warming(r1:r2) + totalCarbon_soil_noAdapt_warming(r1:r2));
% C6 = mean(totalCarbon_ff_noAdapt_warming(r3:r4) + totalCarbon_soil_noAdapt_warming(r3:r4));
% dC_tot_noAdapt = (C5 - C6)/C5;
% 
% % Optimum driven
% C1 = mean(totalCarbon_ff_optimumDriv_warming(r1:r2));
% C2 = mean(totalCarbon_ff_optimumDriv_warming(r3:r4));
% dC_ff_optimumDriv = (C1 - C2)/C1;
% 
% C3 = mean(totalCarbon_soil_optimumDriv_warming(r1:r2));
% C4 = mean(totalCarbon_soil_optimumDriv_warming(r3:r4));
% dC_soil_optimumDriv = (C3 - C4)/C3;
% 
% C5 = mean(totalCarbon_ff_optimumDriv_warming(r1:r2) + totalCarbon_soil_optimumDriv_warming(r1:r2));
% C6 = mean(totalCarbon_ff_optimumDriv_warming(r3:r4) + totalCarbon_soil_optimumDriv_warming(r3:r4));
% dC_tot_optimumDriv = (C5 - C6)/C5;
% 
% % Enzyme rigidity
% C1 = mean(totalCarbon_ff_enzRig_warming(r1:r2));
% C2 = mean(totalCarbon_ff_enzRig_warming(r3:r4));
% dC_ff_enzRig = (C1 - C2)/C1;
% 
% C3 = mean(totalCarbon_soil_enzRig_warming(r1:r2));
% C4 = mean(totalCarbon_soil_enzRig_warming(r3:r4));
% dC_soil_enzRig = (C3 - C4)/C3;
% 
% C5 = mean(totalCarbon_ff_enzRig_warming(r1:r2) + totalCarbon_soil_enzRig_warming(r1:r2));
% C6 = mean(totalCarbon_ff_enzRig_warming(r3:r4) + totalCarbon_soil_enzRig_warming(r3:r4));
% dC_tot_enzRig = (C5 - C6)/C5;
% 
% dC_tot_thermAdapt = mean([dC_tot_optimumDriv dC_tot_enzRig]);

% No adaptation
C1 = mean(totalCarbon_ff_noAdapt_control(r3:r4));
C2 = mean(totalCarbon_ff_noAdapt_warming(r3:r4));
dC_ff_noAdapt = -(C1 - C2)/C1;

C3 = mean(totalCarbon_soil_noAdapt_control(r3:r4));
C4 = mean(totalCarbon_soil_noAdapt_warming(r3:r4));
dC_soil_noAdapt = -(C3 - C4)/C3;

C5 = mean(totalCarbon_ff_noAdapt_control(r3:r4) + totalCarbon_soil_noAdapt_control(r3:r4));
C6 = mean(totalCarbon_ff_noAdapt_warming(r3:r4) + totalCarbon_soil_noAdapt_warming(r3:r4));
dC_tot_noAdapt = -(C5 - C6)/C5;

% Optimum driven
C1 = mean(totalCarbon_ff_optimumDriv_control(r3:r4));
C2 = mean(totalCarbon_ff_optimumDriv_warming(r3:r4));
dC_ff_optimumDriv = -(C1 - C2)/C1;

C3 = mean(totalCarbon_soil_optimumDriv_control(r3:r4));
C4 = mean(totalCarbon_soil_optimumDriv_warming(r3:r4));
dC_soil_optimumDriv = -(C3 - C4)/C3;

C5 = mean(totalCarbon_ff_optimumDriv_control(r3:r4) + totalCarbon_soil_optimumDriv_control(r3:r4));
C6 = mean(totalCarbon_ff_optimumDriv_warming(r3:r4) + totalCarbon_soil_optimumDriv_warming(r3:r4));
dC_tot_optimumDriv = -(C5 - C6)/C5;

% Enzyme rigidity
C1 = mean(totalCarbon_ff_enzRig_control(r3:r4));
C2 = mean(totalCarbon_ff_enzRig_warming(r3:r4));
dC_ff_enzRig = -(C1 - C2)/C1;

C3 = mean(totalCarbon_soil_enzRig_control(r3:r4));
C4 = mean(totalCarbon_soil_enzRig_warming(r3:r4));
dC_soil_enzRig = -(C3 - C4)/C3;

C5 = mean(totalCarbon_ff_enzRig_control(r3:r4) + totalCarbon_soil_enzRig_control(r3:r4));
C6 = mean(totalCarbon_ff_enzRig_warming(r3:r4) + totalCarbon_soil_enzRig_warming(r3:r4));
dC_tot_enzRig = -(C5 - C6)/C5;

dC_tot_thermAdapt = mean([dC_tot_optimumDriv dC_tot_enzRig]);

%% Cumulative CO2 emissions are plotted

% The colors are defined
c_noAdapt = [31,120,180]./255;
c_optDriv = [51,160,44]./255;
c_enzRig = [255,127,0]./255;

yMin = 0;
yMax = 1600;

lineWidth = 2;
scale = 2;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 8*scale], 'color', [1 1 1])
hold on

p1 = plot(noAdaptData.dates_after1991, noAdaptData.diffCumul, 'color', c_noAdapt, 'linewidth', lineWidth);
p2 = plot(optimumDrivenData.dates_after1991, optimumDrivenData.diffCumul, 'color', c_optDriv, 'linewidth', lineWidth);
p3 = plot(enzRigData.dates_after1991, enzRigData.diffCumul, 'color', c_enzRig, 'linewidth', lineWidth);

lgd = legend([p1 p3 p2], 'No thermal adaptation', 'Enzyme rigidity', 'Optimum driven', 'location', 'northwest');%, 'FontSize', 15)
% lgd = legend([p1 p2 p4], 'No adaptation', 'Optimum driven', 'Measured', 'location', 'northwest');%, 'FontSize', 15)
set(lgd,'FontSize',14);
legend('boxoff')

set(gca, 'FontSize', 14, 'Ytick', [yMin:200:yMax], 'box', 'on')

xlabel('Year')
ylabel('Cumulative g CO_{2}-C m^{-2}')
title('Difference in simulated cumulative CO_{2} fluxes between the control and heated treatment', 'FontSize', 16)

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

lineWidth = 1;
scale = 2;

yMin = -50;
yMax = 350;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 8*scale 8*scale], 'color', [1 1 1])
hold on

date1 = datetime('01/01/1989', 'format', 'dd/MM/yyyy');
date2 = datetime('31/12/2018', 'format', 'dd/MM/yyyy');
p0 = plot([2001 2018],[0 0], '--', 'color', [.2 .2 .2]);

% measuredDiff = measuredAnnualCO2_warming - measuredAnnualCO2_control;
% b1 = bar(years_annualMeasuredCO2, measuredDiff, .4);
% set(b1(1), 'FaceColor', color_diff)
% set(b1(1), 'FaceAlpha', .2)

% Let's add a datapoint for 2002, when control and heated fluxes are equal
yrs = [noAdaptData.uniqueYears_CO2diff(1)-1 noAdaptData.uniqueYears_CO2diff];

% Modelled
p1 = plot(yrs, [0; noAdaptData.diffAnnualCO2], '-o', 'color', c_noAdapt, 'linewidth', lineWidth, 'MarkerSize', 6, 'MarkerEdgeColor', c_noAdapt,'MarkerFaceColor', c_noAdapt);
p2 = plot(yrs, [0; optimumDrivenData.diffAnnualCO2], '-o', 'color', c_optDriv, 'linewidth', lineWidth, 'MarkerSize', 6, 'MarkerEdgeColor', c_optDriv,'MarkerFaceColor', c_optDriv);
p3 = plot(yrs, [0; enzRigData.diffAnnualCO2], '-o', 'color', c_enzRig, 'linewidth', lineWidth, 'MarkerSize', 6, 'MarkerEdgeColor', c_enzRig,'MarkerFaceColor', c_enzRig);

% p1 = scatter(noAdaptData.uniqueYears_CO2diff,noAdaptData.diffAnnualCO2, 'MarkerFaceColor', c_noAdapt);
% p2 = scatter(optimumDrivenData.uniqueYears_CO2diff,optimumDrivenData.diffAnnualCO2, 'MarkerFaceColor', c_optDriv);

% Measured
% The standard deviations is calculated
% stDev_diff = sqrt((allMeasurements.measuredAnnualCO2_stDev_control.^2) + (allMeasurements.measuredAnnualCO2_stDev_warming.^2));
% % p3 = plot(years_annualMeasuredCO2, measuredDiff, 'color', color_diff, 'linewidth', 1.5);
% e1 = errorbar(years_annualMeasuredCO2, measuredDiff, stDev_diff./sqrt(6),...
%     'MarkerSize', 4, 'MarkerFaceColor', 'w', 'Color', [.3 .3 .3], 'MarkerEdgeColor', color_diff, 'LineStyle','none');
        
% lgd = legend([b1 p1 p3 p2], 'Calculated (Melillo et al., 2011)', 'No adaptation', 'Enzyme rigidity', 'Optimum driven', 'location', 'northeast');%, 'FontSize', 15)
lgd = legend([p1 p3 p2], 'No thermal adaptation', 'Enzyme rigidity', 'Optimum driven', 'location', 'northeast');%, 'FontSize', 15)
set(lgd,'FontSize',14);
legend('boxoff')

xlim([2001 2017])
ylim([-20 yMax])
set(gca, 'FontSize', 14, 'Ytick', [yMin:50:yMax], 'box', 'on');

xlabel('Year')
ylabel('\DeltaCO_{2} (gC m^{-2} yr^{-1})')
% title('Difference in the simulated annual CO_{2} flux between the control and heated treatment', 'FontSize', 16)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% Letters
text(0.05, 0.97, '(A)', 'FontSize', 14, 'FontWeight', 'bold')
% =========================================================================

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Annual differences in CO2 fluxes.tiff' -r300 -nocrop -transparent
% =========================================================================

dCO2_5to13_noAdapt = mean(noAdaptData.diffAnnualCO2(6:14));
dCO2_5to13_optDriv = mean(optimumDrivenData.diffAnnualCO2(6:14));
dCO2_5to13_enzRig = mean(enzRigData.diffAnnualCO2(6:14));

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
% c_meas = [215,48,31]./255;
c_meas = 'k';

dotSize = 5;
lineWidth = 1;
scale = 2;

yMin = 0;
yMax = 6;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 8*scale], 'color', [1 1 1])
hold on

% -----------------
% The control plots
% -----------------
subplot(2,1,1)
hold on

% The modelled CO2 is plotted
d = datenum(dates_all);
p1 = plot(d,totalCO2_noAdapt_control, 'color', c_noAdapt, 'LineWidth', lineWidth);
p2 = plot(d,totalCO2_optimumDriv_control, 'color', c_optDriv, 'LineWidth', lineWidth);
p3 = plot(d,totalCO2_enzRig_control, 'color', c_enzRig, 'LineWidth', lineWidth);

% The measurements are plotted
d = datenum(dates_CO2_measurements_control);
p4 = errorbar(d, Measured_CO2_flux_control, Measured_CO2_flux_control_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', c_meas, 'MarkerEdgeColor', c_meas, 'LineStyle', 'none');

% The start of warming is indicated
d1 = datenum(date_startWarming);
d2 = datenum(date_startWarming);
plot([d1 d2],[0 yMax], '--', 'color', 'k', 'lineWidth', 2);

date1 = datenum('01/01/2002', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

% Formatting
xlabel('Years')
ylabel('F_{CO_{2}} (gC m^{-2} d^{-1})')
title('Control treatment')

% date1 = datetime('01/01/2000', 'InputFormat', 'dd/MM/yyyy');
% date2 = datetime('31/12/2016', 'InputFormat', 'dd/MM/yyyy');
% 
% xlim([date1 date2])
% ylim([yMin yMax])

set(gca, 'FontSize', 14, 'Ytick', [yMin:1:yMax], 'box', 'on');

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.63; % y
pos(3) = 0.9; % width
pos(4) = 0.31; % height
set(gca,'position',pos)

% ----------------
% The heated plots
% ----------------
subplot(2,1,2)
hold on

% The modelled CO2 is plotted
d = datenum(dates_all);
p1 = plot(d, totalCO2_noAdapt_warmed, 'color', c_noAdapt, 'LineWidth', lineWidth);
p2 = plot(d, totalCO2_optimumDriv_warmed, 'color', c_optDriv, 'LineWidth', lineWidth);
p3 = plot(d, totalCO2_enzRig_warmed, 'color', c_enzRig, 'LineWidth', lineWidth);

% The measurements are plotted
d = datenum(dates_CO2_measurements_control);
p4 = errorbar(d, Measured_CO2_flux_warming, Measured_CO2_flux_warming_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', c_meas, 'MarkerEdgeColor', c_meas, 'LineStyle', 'none');

% The start of warming is indicated
d1 = datenum(date_startWarming);
d2 = datenum(date_startWarming);
plot([d1 d2],[0 yMax], '--', 'color', 'k', 'lineWidth', 2);

date1 = datenum('01/01/2002', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

% Formatting
xlabel('Years')
ylabel('F_{CO_{2}} (gC m^{-2} d^{-1})')
title('Heated treatment')

% date1 = datetime('01/01/2000', 'InputFormat', 'dd/MM/yyyy');
% date2 = datetime('31/12/2016', 'InputFormat', 'dd/MM/yyyy');
% 
% xlim([date1 date2])
% ylim([yMin yMax])

set(gca, 'FontSize', 14, 'Ytick', [yMin:1:yMax], 'box', 'on');

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.06; % x
pos(2) = 0.17; % y
pos(3) = 0.9; % width
pos(4) = 0.32; % height
set(gca,'position',pos)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% Letters
text(0.02, 0.98, '(A)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.02, 0.53, '(B)', 'FontSize', 14, 'FontWeight', 'bold')

% % The legend
l = legend([p1, p2, p3], 'No thermal adaptation', 'Optimum driven', 'Enzyme rigidity', 'Orientation', 'Horizontal', 'FontSIze', 14);
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

%% The daily fluxes with a moving average

% % ------------------
% % The data is loaded
% % ------------------
% 
% totalCO2_noAdapt_control = noAdaptData.totalCO2_control;
% totalCO2_noAdapt_warmed = noAdaptData.totalCO2_warmed;
% 
% totalCO2_optimumDriv_control = optimumDrivenData.totalCO2_control;
% totalCO2_optimumDriv_warmed = optimumDrivenData.totalCO2_warmed;
% 
% totalCO2_enzRig_control = enzRigData.totalCO2_control;
% totalCO2_enzRig_warmed = enzRigData.totalCO2_warmed;
% 
% dates_all = optimumDrivenData.dates_all;
% 
% % --------------------
% % The data are plotted
% % --------------------
% 
% % The colors are defined
% c_noAdapt = [31,120,180]./255;
% c_optDriv = [51,160,44]./255;
% c_enzRig = [255,127,0]./255;
% % c_meas = [215,48,31]./255;
% c_meas = 'k';
% 
% dotSize = 5;
% lineWidth = 1;
% scale = 2;
% 
% yMin = 0;
% yMax = 6;
% 
% hfig = figure;
% set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 8*scale], 'color', [1 1 1])
% hold on
% 
% % -----------------
% % The control plots
% % -----------------
% subplot(2,1,1)
% hold on
% 
% % The moving averages are calculated
% movMean_noAdapt_control = movmean(totalCO2_noAdapt_control,10);
% movMean_optimumDriv_control = movmean(totalCO2_optimumDriv_control,10);
% movMean_enzRig_control = movmean(totalCO2_enzRig_control,10);
% 
% % The modelled CO2 is plotted
% d = datenum(dates_all);
% % p1 = plot(d,totalCO2_noAdapt_control, 'color', c_noAdapt, 'LineWidth', lineWidth);
% % p2 = plot(d,totalCO2_optimumDriv_control, 'color', c_optDriv, 'LineWidth', lineWidth);
% % p3 = plot(d,totalCO2_enzRig_control, 'color', c_enzRig, 'LineWidth', lineWidth);
% 
% p1 = plot(d,movMean_noAdapt_control, 'color', c_noAdapt, 'LineWidth', lineWidth);
% p2 = plot(d,movMean_optimumDriv_control, 'color', c_optDriv, 'LineWidth', lineWidth);
% p3 = plot(d,movMean_enzRig_control, 'color', c_enzRig, 'LineWidth', lineWidth);
% 
% date1 = datenum('01/01/2003', 'dd/MM/yyyy');
% date2 = datenum('31/12/2008', 'dd/MM/yyyy');
% xlim([date1 date2])
% 
% % The measurements are plotted
% d = datenum(dates_CO2_measurements_control);
% p4 = errorbar(d, Measured_CO2_flux_control, Measured_CO2_flux_control_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', c_meas, 'MarkerEdgeColor', c_meas, 'LineStyle', 'none');
% 
% % The start of warming is indicated
% d1 = datenum(date_startWarming);
% d2 = datenum(date_startWarming);
% plot([d1 d2],[0 yMax], '--', 'color', 'k', 'lineWidth', 2);
% 
% % date1 = datenum('01/01/2002', 'dd/MM/yyyy');
% % date2 = datenum('31/12/2016', 'dd/MM/yyyy');
% % 
% % xlim([date1 date2])
% datetick('x', 'yyyy','keeplimits')
% ylim([yMin yMax])
% 
% % Formatting
% xlabel('Years')
% ylabel('F_{CO_{2}} (gC m^{-2} d^{-1})')
% title('Control treatment')
% 
% % date1 = datetime('01/01/2000', 'InputFormat', 'dd/MM/yyyy');
% % date2 = datetime('31/12/2016', 'InputFormat', 'dd/MM/yyyy');
% % 
% % xlim([date1 date2])
% % ylim([yMin yMax])
% 
% set(gca, 'FontSize', 14, 'Ytick', [yMin:1:yMax], 'box', 'on');
% 
% % The location of the plot is optimized
% ax = gca;
% pos = ax.Position;
% pos(1) = 0.06; % x
% pos(2) = 0.63; % y
% pos(3) = 0.9; % width
% pos(4) = 0.31; % height
% set(gca,'position',pos)
% 
% % ----------------
% % The heated plots
% % ----------------
% subplot(2,1,2)
% hold on
% 
% % The moving averages are calculated
% movMean_noAdapt_warmed = movmean(totalCO2_noAdapt_warmed,10);
% movMean_optimumDriv_warmed = movmean(totalCO2_optimumDriv_warmed,10);
% movMean_enzRig_warmed = movmean(totalCO2_enzRig_warmed,10);
% 
% % The modelled CO2 is plotted
% d = datenum(dates_all);
% % p1 = plot(d, totalCO2_noAdapt_warmed, 'color', c_noAdapt, 'LineWidth', lineWidth);
% % p2 = plot(d, totalCO2_optimumDriv_warmed, 'color', c_optDriv, 'LineWidth', lineWidth);
% % p3 = plot(d, totalCO2_enzRig_warmed, 'color', c_enzRig, 'LineWidth', lineWidth);
% 
% p1 = plot(d, movMean_noAdapt_warmed, 'color', c_noAdapt, 'LineWidth', lineWidth);
% p2 = plot(d, movMean_optimumDriv_warmed, 'color', c_optDriv, 'LineWidth', lineWidth);
% p3 = plot(d, movMean_enzRig_warmed, 'color', c_enzRig, 'LineWidth', lineWidth);
% 
% date1 = datenum('01/01/2003', 'dd/MM/yyyy');
% date2 = datenum('31/12/2008', 'dd/MM/yyyy');
% xlim([date1 date2])
% 
% % The measurements are plotted
% d = datenum(dates_CO2_measurements_control);
% p4 = errorbar(d, Measured_CO2_flux_warming, Measured_CO2_flux_warming_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', c_meas, 'MarkerEdgeColor', c_meas, 'LineStyle', 'none');
% 
% % The start of warming is indicated
% d1 = datenum(date_startWarming);
% d2 = datenum(date_startWarming);
% plot([d1 d2],[0 yMax], '--', 'color', 'k', 'lineWidth', 2);
% 
% % date1 = datenum('01/01/2002', 'dd/MM/yyyy');
% % date2 = datenum('31/12/2016', 'dd/MM/yyyy');
% % 
% % xlim([date1 date2])
% datetick('x', 'yyyy','keeplimits')
% ylim([yMin yMax])
% 
% % Formatting
% xlabel('Years')
% ylabel('F_{CO_{2}} (gC m^{-2} d^{-1})')
% title('Heated treatment')
% 
% % date1 = datetime('01/01/2000', 'InputFormat', 'dd/MM/yyyy');
% % date2 = datetime('31/12/2016', 'InputFormat', 'dd/MM/yyyy');
% % 
% % xlim([date1 date2])
% % ylim([yMin yMax])
% 
% set(gca, 'FontSize', 14, 'Ytick', [yMin:1:yMax], 'box', 'on');
% 
% % The location of the plot is optimized
% ax = gca;
% pos = ax.Position;
% pos(1) = 0.06; % x
% pos(2) = 0.17; % y
% pos(3) = 0.9; % width
% pos(4) = 0.32; % height
% set(gca,'position',pos)
% 
% % =========================================================================
% % A new invisible axis is created
% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% axes(ax1)
% xlim([0 1])
% ylim([0 1])
% 
% % Letters
% text(0.02, 0.98, '(A)', 'FontSize', 14, 'FontWeight', 'bold')
% text(0.02, 0.53, '(B)', 'FontSize', 14, 'FontWeight', 'bold')
% 
% % % The legend
% l = legend([p1, p2, p3], 'No thermal adaptation', 'Optimum driven', 'Enzyme rigidity', 'Orientation', 'Horizontal', 'FontSIze', 14);
% pos = l.Position;
% pos(1) = 0.3; % x
% pos(2) = 0.02; % y
% l.Position = pos;
% legend('boxoff')

%% The difference in simulated CO2 is plotted versus measurements

% The differences per treatment are calculated
diff_CO2_noAdapt = totalCO2_noAdapt_warmed - totalCO2_noAdapt_control;
diff_CO2_optimumDriven = totalCO2_optimumDriv_warmed - totalCO2_optimumDriv_control;
diff_CO2_enzRig = totalCO2_enzRig_warmed - totalCO2_enzRig_control;

% Some outliers are removed
val = 1.2;

[r c] = find(diff_CO2_noAdapt > val);
diff_CO2_noAdapt(r,c) = NaN;

[r c] = find(diff_CO2_optimumDriven > val);
diff_CO2_optimumDriven(r,c) = NaN;

[r c] = find(diff_CO2_enzRig > val);
diff_CO2_enzRig(r,c) = NaN;

dates_all = optimumDrivenData.dates_all;

% --------
% Plotting
% --------

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 8*scale], 'color', [1 1 1])
hold on

yMin = -1.1;
yMax = 1.25;
dotSize = 40;

% A line to indicate y = 0
yline(0, '--')

d = datenum(dates_all);
p1 = plot(d,diff_CO2_noAdapt, 'color', c_noAdapt, 'LineWidth', lineWidth);
p2 = plot(d,diff_CO2_optimumDriven, 'color', c_optDriv, 'LineWidth', lineWidth);
p3 = plot(d,diff_CO2_enzRig, 'color', c_enzRig, 'LineWidth', lineWidth);

% The start of warming is indicated
% d1 = datenum(date_startWarming);
% d2 = datenum(date_startWarming);
% plot([d1 d2],[yMin yMax], '--', 'color', 'k', 'lineWidth', 2);

date1 = datenum('01/01/2002', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

xlim([d1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

% Formatting
xlabel('Year')
ylabel('\DeltaCO_{2} (gC m^{-2} d^{-1})')
% title('Control treatment')

d = datenum(dates_CO2_measurements_control);
p5 = errorbar(d, Measured_CO2_flux_warming - Measured_CO2_flux_control, sqrt(Measured_CO2_flux_control_stdError.^2 + Measured_CO2_flux_warming_stdError.^2), 'color', [.5 .5 .5], 'LineStyle', 'none');
p4 = scatter(d, Measured_CO2_flux_warming - Measured_CO2_flux_control, dotSize, 'k', 'filled');
% p5 = errorbar(d, Measured_CO2_flux_warming - Measured_CO2_flux_control, sqrt(Measured_CO2_flux_control_stdError.^2 + Measured_CO2_flux_warming_stdError.^2), 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', c_meas, 'MarkerEdgeColor', c_meas, 'LineStyle', 'none');

set(gca, 'FontSize', 14, 'Ytick', [-1;0;1], 'box', 'on');

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% Letters
text(0.1, 0.97, '(B)', 'FontSize', 14, 'FontWeight', 'bold')

% % The legend
l = legend([p1, p2, p3], 'No thermal adaptation', 'Optimum driven', 'Enzyme rigidity', 'Orientation', 'Horizontal', 'FontSize', 14);
pos = l.Position;
pos(1) = 0.14; % x
pos(2) = 0.13; % y
l.Position = pos;
% legend('boxoff')
% =========================================================================

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'daily differences in CO2 fluxes.tiff' -r300 -nocrop -transparent
% =========================================================================

% ----------------------
% Same plot with movmean
% ----------------------

k = 10;

diff_CO2_noAdapt_mm = movmean(diff_CO2_noAdapt, k);
diff_CO2_optimumDriven_mm = movmean(diff_CO2_optimumDriven, k);
diff_CO2_enzRig_mm = movmean(diff_CO2_enzRig, k);

% dates_all = optimumDrivenData.dates_all;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 8*scale], 'color', [1 1 1])
hold on

yMin = -1.1;
yMax = 1.25;
dotSize = 40;

% A line to indicate y = 0
yline(0, '--')

lineWidth = 2;

d = datenum(dates_all);
p1 = plot(d,diff_CO2_noAdapt_mm, 'color', c_noAdapt, 'LineWidth', lineWidth);
p2 = plot(d,diff_CO2_optimumDriven_mm, 'color', c_optDriv, 'LineWidth', lineWidth);
p3 = plot(d,diff_CO2_enzRig_mm, 'color', c_enzRig, 'LineWidth', lineWidth);

% The start of warming is indicated
% d1 = datenum(date_startWarming);
% d2 = datenum(date_startWarming);
% plot([d1 d2],[yMin yMax], '--', 'color', 'k', 'lineWidth', 2);

date1 = datenum('01/01/2002', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

xlim([d1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

% Formatting
xlabel('Year')
ylabel('\DeltaCO_{2} (gC m^{-2} d^{-1})')
% title('Control treatment')

d = datenum(dates_CO2_measurements_control);
p5 = errorbar(d, Measured_CO2_flux_warming - Measured_CO2_flux_control, sqrt(Measured_CO2_flux_control_stdError.^2 + Measured_CO2_flux_warming_stdError.^2), 'color', [.5 .5 .5], 'LineStyle', 'none');
p4 = scatter(d, Measured_CO2_flux_warming - Measured_CO2_flux_control, dotSize, 'k', 'filled');
% p5 = errorbar(d, Measured_CO2_flux_warming - Measured_CO2_flux_control, sqrt(Measured_CO2_flux_control_stdError.^2 + Measured_CO2_flux_warming_stdError.^2), 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', c_meas, 'MarkerEdgeColor', c_meas, 'LineStyle', 'none');

set(gca, 'FontSize', 14, 'Ytick', [-1;0;1], 'box', 'on');

% ----------------------
% Scatterplots
% ----------------------

% The data for the first 4 years are isolated
% The modeled CO2 on measurement days is extracted

[r, c] = find(ismember(dates_all, dates_CO2_measurements_control));

diff_CO2_noAdapt_measDays = diff_CO2_noAdapt(c);
diff_CO2_optDriv_measDays = diff_CO2_optimumDriven(c);
diff_CO2_enzRig_measDays = diff_CO2_enzRig(c);

% The data for the first 4 years is extracted
initialYears = 2003:2006;
allYears = dates_CO2_measurements_control.Year;
[r_init, c] = find(ismember(allYears, initialYears));
[r_later, c] = find(~ismember(allYears, initialYears));

% The modelled data for the first 4 years
diff_CO2_noAdapt_measDays_initialYears = diff_CO2_noAdapt_measDays(r_init);
diff_CO2_optDriv_measDays_initialYears = diff_CO2_optDriv_measDays(r_init);
diff_CO2_enzRig_measDays_initialYears = diff_CO2_enzRig_measDays(r_init);

% The modelled data for the later years
diff_CO2_noAdapt_measDays_laterYears = diff_CO2_noAdapt_measDays;
diff_CO2_optDriv_measDays_laterYears = diff_CO2_optDriv_measDays;
diff_CO2_enzRig_measDays_laterYears = diff_CO2_enzRig_measDays;

diff_CO2_noAdapt_measDays_laterYears(r_init) = [];
diff_CO2_optDriv_measDays_laterYears(r_init) = [];
diff_CO2_enzRig_measDays_laterYears(r_init) = [];

% The measured data for the first 4 years
diff_CO2_measured = Measured_CO2_flux_warming - Measured_CO2_flux_control;
diff_CO2_measured_initialYears = diff_CO2_measured(r_init);

% The measured data for the later years
diff_CO2_measured_laterYears = Measured_CO2_flux_warming - Measured_CO2_flux_control;
diff_CO2_measured_laterYears(r_init) = [];

% Some outliers are removed
val = 1.2;

[r c] = find(diff_CO2_noAdapt > val);
diff_CO2_noAdapt(r,c) = NaN;

[r c] = find(diff_CO2_optimumDriven > val);
diff_CO2_optimumDriven(r,c) = NaN;

[r c] = find(diff_CO2_enzRig > val);
diff_CO2_enzRig(r,c) = NaN;

% === The figure ===

c_firstYear = [227,26,28]./255;
c_later = [.2 .2 .2];

nParam_noAdapt = 16;
nParam_optDriv = 18;
nParam_enzRig = 19;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [15 15 18*scale 6*scale], 'color', [1 1 1])

% --- No adaptation ---
subplot(1,3,1)
hold on

% The regression line
r = isnan(diff_CO2_noAdapt_measDays_initialYears);
diff_CO2_measured_initialYears_r = diff_CO2_measured_initialYears;
diff_CO2_noAdapt_measDays_initialYears_r = diff_CO2_noAdapt_measDays_initialYears;
diff_CO2_measured_initialYears_r(r) = [];
diff_CO2_noAdapt_measDays_initialYears_r(r) = [];

val = 1.2;
[r c] = find(diff_CO2_measured_initialYears_r > val);
diff_CO2_measured_initialYears_r(r) = NaN;
diff_CO2_measured_initialYears_r(r) = [];
diff_CO2_noAdapt_measDays_initialYears_r(r) = [];

P = polyfit(diff_CO2_measured_initialYears_r,diff_CO2_noAdapt_measDays_initialYears_r,1);
x0 = min(diff_CO2_measured_initialYears_r) ; 
x1 = max(diff_CO2_measured_initialYears_r) ;
xi = linspace(x0,x1) ;
yi = P(1)*xi+P(2);
hold on
% plot(xi, yi, '--', 'color', color_regr, 'linewidth', 2) ;

% ----
% Bias
% ----
bias_4_noAdapt = sum(diff_CO2_noAdapt_measDays_initialYears_r - diff_CO2_measured_initialYears_r) / numel(diff_CO2_measured_initialYears_r);
bias_4_noAdapt = round(bias_4_noAdapt*100)/100;

bias_laterYears_noAdapt = sum(diff_CO2_noAdapt_measDays_laterYears - diff_CO2_measured_laterYears) / numel(diff_CO2_measured_laterYears);
bias_laterYears_noAdapt = round(bias_laterYears_noAdapt*100)/100;

all_meas = [diff_CO2_measured_initialYears_r; diff_CO2_measured_laterYears];
all_mod = [diff_CO2_noAdapt_measDays_initialYears_r; diff_CO2_noAdapt_measDays_laterYears];
bias_all_noAdapt = sum(all_mod - all_meas) / numel(all_meas);
bias_all_noAdapt = round(bias_all_noAdapt*100)/100;

% Histogram of biases
% histogram(diff_CO2_noAdapt_measDays_initialYears_r - diff_CO2_measured_initialYears_r, 20)
% qqplot(diff_CO2_noAdapt_measDays_initialYears_r - diff_CO2_measured_initialYears_r)

% histogram(all_mod - all_meas, 10)
% qqplot(all_mod - all_meas)

% Null hypothesis that the data in x comes from a normal distribution with mean equal to zero
% The result h is 1 if the test rejects the null hypothesis at the 5 significance level
h1 = ttest(diff_CO2_noAdapt_measDays_initialYears_r - diff_CO2_measured_initialYears_r)
h2 = ttest(all_mod - all_meas)

% The value h = 1 indicates that the test rejects the null hypothesis of zero median.
p1 = signrank(diff_CO2_noAdapt_measDays_initialYears_r - diff_CO2_measured_initialYears_r)
p2 = signrank(all_mod - all_meas)

% ----
% RMSE
% ----
RMSE4_noAdapt = sqrt(mean((diff_CO2_noAdapt_measDays_initialYears_r - diff_CO2_measured_initialYears_r).^2));
RMSE4_noAdapt = round(RMSE4_noAdapt*100)/100;

RMSE_laterYears_noAdapt = sqrt(mean((diff_CO2_noAdapt_measDays_laterYears - diff_CO2_measured_laterYears).^2));
RMSE_laterYears_noAdapt = round(RMSE_laterYears_noAdapt*100)/100;

RMSE_all_noAdapt = sqrt(mean((all_mod - all_meas).^2));
RMSE_all_noAdapt = round(RMSE_all_noAdapt*100)/100;

% ---
% AIC
% ---

% First the SSR is calculated
SSR_4_noAdapt = sum((diff_CO2_noAdapt_measDays_initialYears_r - diff_CO2_measured_initialYears_r).^2);
SSR_laterYears_noAdapt = sum((diff_CO2_noAdapt_measDays_laterYears - diff_CO2_measured_laterYears).^2);
SSR_all_noAdapt = sum((all_mod - all_meas).^2);

% AIC
AIC_4_noAdapt = log(SSR_4_noAdapt/numel(diff_CO2_measured_initialYears_r)) + (2 * nParam_noAdapt);
AIC_laterYears_noAdapt = log(SSR_laterYears_noAdapt/numel(diff_CO2_measured_laterYears)) + (2 * nParam_noAdapt);
AIC_all_noAdapt = log(SSR_all_noAdapt/numel(all_meas)) + (2 * nParam_noAdapt);

AIC_4_noAdapt = round(AIC_4_noAdapt*100)/100;
AIC_laterYears_noAdapt = round(AIC_laterYears_noAdapt*100)/100;
AIC_all_noAdapt = round(AIC_all_noAdapt*100)/100;

xlim([-1 2])
ylim([-1 2])

scatter(diff_CO2_measured_laterYears, diff_CO2_noAdapt_measDays_laterYears, 50, 'MarkerFaceColor', c_later, 'MarkerEdgeColor', 'none');
scatter(diff_CO2_measured_initialYears_r, diff_CO2_noAdapt_measDays_initialYears_r, 50, 'MarkerFaceColor', c_firstYear, 'MarkerEdgeColor', 'none');

% The 1:1 line
hline1 = refline(1,0);
hline1.Color = [.5 .5 .5];
hline1.LineStyle = '--';

xlim([-1 1.5])
ylim([-1 1.5])

set(gca, 'FontSize', 12, 'Ytick', [yMin:1:yMax], 'box', 'on')

xlabel('Measured \DeltaCO_{2} (gC m^{-2} d^{-1})')
ylabel('Modelled \DeltaCO_{2} (gC m^{-2} d^{-1})')

title('No thermal adaptation')

% Adding the measures
text(0.5, -0.4, "Years")
text(1.1, -0.4, "1-4", 'HorizontalAlignment', 'right')
text(1.2, -0.4, "1-13", 'HorizontalAlignment', 'left')

plot([0.8 1.45], [-0.5 -0.5], 'color', [0 0 0]) % Horizontal line
plot([1.15 1.15], [-0.35 -0.48], 'color', [0 0 0]) % Vertical line top
plot([1.15 1.15], [-0.98 -0.52], 'color', [0 0 0]) % Vertical line bottom

text(0.5, -0.6, "MAE:")
text(1.1, -0.6, append(num2str(bias_4_noAdapt), "*"), 'HorizontalAlignment', 'right')
text(1.2, -0.6, append(num2str(bias_all_noAdapt), "*"), 'HorizontalAlignment', 'left')

text(0.5, -0.75, "RMSE:")
text(1.1, -0.75, num2str(RMSE4_noAdapt), 'HorizontalAlignment', 'right')
text(1.2, -0.75, num2str(RMSE_all_noAdapt), 'HorizontalAlignment', 'left')

text(0.5, -0.9, "AIC:")
text(1.1, -0.9, num2str(AIC_4_noAdapt), 'HorizontalAlignment', 'right')
text(1.2, -0.9, num2str(AIC_all_noAdapt), 'HorizontalAlignment', 'left')

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.05; % x
pos(2) = 0.23; % y
pos(3) = 0.25; % width
pos(4) = 0.69; % height
set(gca,'position',pos)

% --- Optimum driven ---
subplot(1,3,2)
hold on

% The regression line
r = isnan(diff_CO2_optDriv_measDays_initialYears);
diff_CO2_measured_initialYears_r = diff_CO2_measured_initialYears;
diff_CO2_optDriv_measDays_initialYears_r = diff_CO2_optDriv_measDays_initialYears;
diff_CO2_measured_initialYears_r(r) = [];
diff_CO2_optDriv_measDays_initialYears_r(r) = [];

val = 1.2;
[r c] = find(diff_CO2_measured_initialYears_r > val);
diff_CO2_measured_initialYears_r(r) = NaN;
diff_CO2_measured_initialYears_r(r) = [];
diff_CO2_optDriv_measDays_initialYears_r(r) = [];

P = polyfit(diff_CO2_measured_initialYears_r,diff_CO2_optDriv_measDays_initialYears_r,1);
x0 = min(diff_CO2_measured_initialYears_r) ; 
x1 = max(diff_CO2_measured_initialYears_r) ;
xi = linspace(x0,x1) ;
yi = P(1)*xi+P(2);
hold on
% plot(xi, yi, '--', 'color', color_regr, 'linewidth', 2) ;

% ----
% Bias
% ----
bias_4_optDriv = sum(diff_CO2_optDriv_measDays_initialYears_r - diff_CO2_measured_initialYears_r) / numel(diff_CO2_measured_initialYears_r);
bias_4_optDriv = round(bias_4_optDriv*100)/100;

bias_laterYears_optDriv = sum(diff_CO2_optDriv_measDays_laterYears - diff_CO2_measured_laterYears) / numel(diff_CO2_measured_laterYears);
bias_laterYears_optDriv = round(bias_laterYears_optDriv*100)/100;

all_meas = [diff_CO2_measured_initialYears_r; diff_CO2_measured_laterYears];
all_mod = [diff_CO2_optDriv_measDays_initialYears_r; diff_CO2_optDriv_measDays_laterYears];
bias_all_optDriv = sum(all_mod - all_meas) / numel(all_meas);
bias_all_optDriv = round(bias_all_optDriv*100)/100;

% Histogram of biases
% histogram(diff_CO2_optDriv_measDays_initialYears_r - diff_CO2_measured_initialYears_r, 8)
% qqplot(diff_CO2_optDriv_measDays_initialYears_r - diff_CO2_measured_initialYears_r)

% histogram(all_mod - all_meas, 10)
% qqplot(all_mod - all_meas)

% Null hypothesis that the data in x comes from a normal distribution with mean equal to zero
% The result h is 1 if the test rejects the null hypothesis at the 5 significance level
h1 = ttest(diff_CO2_optDriv_measDays_initialYears_r - diff_CO2_measured_initialYears_r)
h2 = ttest(all_mod - all_meas)

% The value h = 1 indicates that the test rejects the null hypothesis of zero median.
p1 = signrank(diff_CO2_optDriv_measDays_initialYears_r - diff_CO2_measured_initialYears_r)
p2 = signrank(all_mod - all_meas)

% ----
% RMSE
% ----
RMSE4_optDriv = sqrt(mean((diff_CO2_optDriv_measDays_initialYears_r - diff_CO2_measured_initialYears_r).^2));
RMSE4_optDriv = round(RMSE4_optDriv*100)/100;

RMSE_laterYears_optDriv = sqrt(mean((diff_CO2_optDriv_measDays_laterYears - diff_CO2_measured_laterYears).^2));
RMSE_laterYears_optDriv = round(RMSE_laterYears_optDriv*100)/100;

RMSE_all_optDriv = sqrt(mean((all_mod - all_meas).^2));
RMSE_all_optDriv = round(RMSE_all_optDriv*100)/100;

% ---
% AIC
% ---

% First the SSR is calculated
SSR_4_optDriv = sum((diff_CO2_optDriv_measDays_initialYears_r - diff_CO2_measured_initialYears_r).^2);
SSR_laterYears_optDriv = sum((diff_CO2_optDriv_measDays_laterYears - diff_CO2_measured_laterYears).^2);
SSR_all_optDriv = sum((all_mod - all_meas).^2);

% AIC
AIC_4_optDriv = log(SSR_4_optDriv/numel(diff_CO2_measured_initialYears_r)) + (2 * nParam_optDriv);
AIC_laterYears_optDriv = log(SSR_laterYears_optDriv/numel(diff_CO2_measured_laterYears)) + (2 * nParam_optDriv);
AIC_all_optDriv = log(SSR_all_optDriv/numel(all_meas)) + (2 * nParam_optDriv);

AIC_4_optDriv = round(AIC_4_optDriv*100)/100;
AIC_laterYears_optDriv = round(AIC_laterYears_optDriv*100)/100;
AIC_all_optDriv = round(AIC_all_optDriv*100)/100;

xlim([-1 2])
ylim([-1 2])

scatter(diff_CO2_measured_laterYears, diff_CO2_optDriv_measDays_laterYears, 50, 'MarkerFaceColor', c_later, 'MarkerEdgeColor', 'none');
scatter(diff_CO2_measured_initialYears_r, diff_CO2_optDriv_measDays_initialYears_r, 50, 'MarkerFaceColor', c_firstYear, 'MarkerEdgeColor', 'none');

% The 1:1 line
hline1 = refline(1,0);
hline1.Color = [.5 .5 .5];
hline1.LineStyle = '--';

xlim([-1 1.5])
ylim([-1 1.5])

set(gca, 'FontSize', 12, 'Ytick', [yMin:1:yMax], 'box', 'on')

xlabel('Measured \DeltaCO_{2} (gC m^{-2} d^{-1})')
ylabel('Modelled \DeltaCO_{2} (gC m^{-2} d^{-1})')

title('Optimum driven')

% Adding the measures
text(0.5, -0.4, "Years")
text(1.1, -0.4, "1-4", 'HorizontalAlignment', 'right')
text(1.2, -0.4, "1-13", 'HorizontalAlignment', 'left')

plot([0.8 1.45], [-0.5 -0.5], 'color', [0 0 0]) % Horizontal line
plot([1.15 1.15], [-0.35 -0.48], 'color', [0 0 0]) % Vertical line top
plot([1.15 1.15], [-0.98 -0.52], 'color', [0 0 0]) % Vertical line bottom

text(0.5, -0.6, "MAE:")
text(1.1, -0.6, num2str(bias_4_optDriv), 'HorizontalAlignment', 'right')
text(1.2, -0.6, append(num2str(bias_all_optDriv), "*"), 'HorizontalAlignment', 'left')

text(0.5, -0.75, "RMSE:")
text(1.1, -0.75, num2str(RMSE4_optDriv), 'HorizontalAlignment', 'right')
text(1.2, -0.75, num2str(RMSE_all_optDriv), 'HorizontalAlignment', 'left')

text(0.5, -0.9, "AIC:")
text(1.1, -0.9, num2str(AIC_4_optDriv), 'HorizontalAlignment', 'right')
text(1.2, -0.9, num2str(AIC_all_optDriv), 'HorizontalAlignment', 'left')

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.38; % x
pos(2) = 0.23; % y
pos(3) = 0.25; % width
pos(4) = 0.69; % height
set(gca,'position',pos)

% --- Enzyme rigidity ---
subplot(1,3,3)
hold on

% The regression line
r = isnan(diff_CO2_enzRig_measDays_initialYears);
diff_CO2_measured_initialYears_r = diff_CO2_measured_initialYears;
diff_CO2_enzRig_measDays_initialYears_r = diff_CO2_enzRig_measDays_initialYears;
diff_CO2_measured_initialYears_r(r) = [];
diff_CO2_enzRig_measDays_initialYears_r(r) = [];

val = 1.2;
[r c] = find(diff_CO2_measured_initialYears_r > val);
diff_CO2_measured_initialYears_r(r) = NaN;
diff_CO2_measured_initialYears_r(r) = [];
diff_CO2_enzRig_measDays_initialYears_r(r) = [];

P = polyfit(diff_CO2_measured_initialYears_r,diff_CO2_enzRig_measDays_initialYears_r,1);
x0 = min(diff_CO2_measured_initialYears_r) ; 
x1 = max(diff_CO2_measured_initialYears_r) ;
xi = linspace(x0,x1) ;
yi = P(1)*xi+P(2);
hold on
% plot(xi, yi, '--', 'color', color_regr, 'linewidth', 2) ;

% ----
% Bias
% ----
bias_4_enzRig = sum(diff_CO2_enzRig_measDays_initialYears_r - diff_CO2_measured_initialYears_r) / numel(diff_CO2_measured_initialYears_r);
bias_4_enzRig = round(bias_4_enzRig*100)/100;

bias_laterYears_enzRig = sum(diff_CO2_enzRig_measDays_laterYears - diff_CO2_measured_laterYears) / numel(diff_CO2_measured_laterYears);
bias_laterYears_enzRig = round(bias_laterYears_enzRig*100)/100;

all_meas = [diff_CO2_measured_initialYears_r; diff_CO2_measured_laterYears];
all_mod = [diff_CO2_enzRig_measDays_initialYears_r; diff_CO2_enzRig_measDays_laterYears];
bias_all_enzRig = sum(all_mod - all_meas) / numel(all_meas);
bias_all_enzRig = round(bias_all_enzRig*100)/100;

% Histogram of biases
% histogram(diff_CO2_enzRig_measDays_initialYears_r - diff_CO2_measured_initialYears_r, 8)
% qqplot(diff_CO2_enzRig_measDays_initialYears_r - diff_CO2_measured_initialYears_r)

% histogram(all_mod - all_meas, 10)
% qqplot(all_mod - all_meas)

% Null hypothesis that the data in x comes from a normal distribution with mean equal to zero
% The result h is 1 if the test rejects the null hypothesis at the 5 significance level
h1 = ttest(diff_CO2_enzRig_measDays_initialYears_r - diff_CO2_measured_initialYears_r)
h2 = ttest(all_mod - all_meas)

% The value h = 1 indicates that the test rejects the null hypothesis of zero median.
p1 = signrank(diff_CO2_enzRig_measDays_initialYears_r - diff_CO2_measured_initialYears_r)
p2 = signrank(all_mod - all_meas)

% ----
% RMSE
% ----
RMSE4_enzRig = sqrt(mean((diff_CO2_enzRig_measDays_initialYears_r - diff_CO2_measured_initialYears_r).^2));
RMSE4_enzRig = round(RMSE4_enzRig*100)/100;

RMSE_laterYears_enzRig = sqrt(mean((diff_CO2_enzRig_measDays_laterYears - diff_CO2_measured_laterYears).^2));
RMSE_laterYears_enzRig = round(RMSE_laterYears_enzRig*100)/100;

RMSE_all_enzRig = sqrt(mean((all_mod - all_meas).^2));
RMSE_all_enzRig = round(RMSE_all_enzRig*100)/100;

% ---
% AIC
% ---

% First the SSR is calculated
SSR_4_enzRig = sum((diff_CO2_enzRig_measDays_initialYears_r - diff_CO2_measured_initialYears_r).^2);
SSR_laterYears_enzRig = sum((diff_CO2_enzRig_measDays_laterYears - diff_CO2_measured_laterYears).^2);
SSR_all_enzRig = sum((all_mod - all_meas).^2);

% AIC
AIC_4_enzRig = log(SSR_4_enzRig/numel(diff_CO2_measured_initialYears_r)) + (2 * nParam_enzRig);
AIC_laterYears_enzRig = log(SSR_laterYears_enzRig/numel(diff_CO2_measured_laterYears)) + (2 * nParam_enzRig);
AIC_all_enzRig = log(SSR_all_enzRig/numel(all_meas)) + (2 * nParam_enzRig);

AIC_4_enzRig = round(AIC_4_enzRig*100)/100;
AIC_laterYears_enzRig = round(AIC_laterYears_enzRig*100)/100;
AIC_all_enzRig = round(AIC_all_enzRig*100)/100;

xlim([-1 2])
ylim([-1 2])

s1 = scatter(diff_CO2_measured_laterYears, diff_CO2_enzRig_measDays_laterYears, 50, 'MarkerFaceColor', c_later, 'MarkerEdgeColor', 'none');
s2 = scatter(diff_CO2_measured_initialYears_r, diff_CO2_enzRig_measDays_initialYears_r, 50, 'MarkerFaceColor', c_firstYear, 'MarkerEdgeColor', 'none');

% The 1:1 line
hline1 = refline(1,0);
hline1.Color = [.5 .5 .5];
hline1.LineStyle = '--';

xlim([-1 1.5])
ylim([-1 1.5])

set(gca, 'FontSize', 12, 'Ytick', [yMin:1:yMax], 'box', 'on')

xlabel('Measured \DeltaCO_{2} (gC m^{-2} d^{-1})')
ylabel('Modelled \DeltaCO_{2} (gC m^{-2} d^{-1})')

title('Enzyme rigidity')

% Adding the measures
text(0.5, -0.4, "Years")
text(1.1, -0.4, "1-4", 'HorizontalAlignment', 'right')
text(1.2, -0.4, "1-13", 'HorizontalAlignment', 'left')

plot([0.8 1.45], [-0.5 -0.5], 'color', [0 0 0]) % Horizontal line
plot([1.15 1.15], [-0.35 -0.48], 'color', [0 0 0]) % Vertical line top
plot([1.15 1.15], [-0.98 -0.52], 'color', [0 0 0]) % Vertical line bottom

text(0.5, -0.6, "MAE:")
text(1.1, -0.6, append(num2str(bias_4_enzRig), "*"), 'HorizontalAlignment', 'right')
text(1.2, -0.6, num2str(bias_all_enzRig), 'HorizontalAlignment', 'left')

text(0.5, -0.75, "RMSE:")
text(1.1, -0.75, num2str(RMSE4_enzRig), 'HorizontalAlignment', 'right')
text(1.2, -0.75, num2str(RMSE_all_enzRig), 'HorizontalAlignment', 'left')

text(0.5, -0.9, "AIC:")
text(1.1, -0.9, num2str(AIC_4_enzRig), 'HorizontalAlignment', 'right')
text(1.2, -0.9, num2str(AIC_all_enzRig), 'HorizontalAlignment', 'left')

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.71; % x
pos(2) = 0.23; % y
pos(3) = 0.25; % width
pos(4) = 0.69; % height
set(gca,'position',pos)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
% text(0.41,0.98,'Organic carbon stocks','fontweight','bold','fontsize',14)

% Letters
text(0.015, 0.96, '(C)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.345, 0.96, '(D)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.675, 0.96, '(E)', 'FontSize', 14, 'FontWeight', 'bold')

% The legend
l = legend([s2, s1], 'Data for first 4 years of warming', 'Data for years 5 - 13', 'Orientation', 'Horizontal', 'FontSIze', 14);
pos = l.Position;
pos(1) = 0.3; % x
pos(2) = 0.02; % y
l.Position = pos;
% =========================================================================

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Scatter dCO2 per treatment.tiff' -r300 -nocrop -transparent
% =========================================================================


%% Daily CO2 emissions are plotted per treatment

c_control_mod = [54,144,192]./255;
c_control_meas = [2,56,88]./255;
c_warming_mod = [239,101,72]./255;
c_warming_meas = [179,0,0]./255;

dotSize = 5;
lineWidth = 1;
scale = 2;

yMin = 0;
yMax = 6;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 18*scale 18*scale], 'color', [1 1 1])

% ---------------------
% No thermal adaptation
% ---------------------

subplot(3,1,1)
hold on

% Modelled
d = datenum(dates_all);
p1 = plot(d,totalCO2_noAdapt_control, 'color', c_control_mod, 'LineWidth', lineWidth);
p2 = plot(d, totalCO2_noAdapt_warmed, 'color', c_warming_mod, 'LineWidth', lineWidth);

% Measured
d = datenum(dates_CO2_measurements_control);
p3 = errorbar(d, Measured_CO2_flux_control, Measured_CO2_flux_control_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', c_control_meas, 'MarkerEdgeColor', c_control_meas, 'LineStyle', 'none');
p4 = errorbar(d, Measured_CO2_flux_warming, Measured_CO2_flux_warming_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', c_warming_meas, 'MarkerEdgeColor', c_warming_meas, 'LineStyle', 'none');

% The start of warming is indicated
d1 = datenum(date_startWarming);
d2 = datenum(date_startWarming);
plot([d1 d2],[0 yMax], ':', 'color', 'k', 'lineWidth', 2);

date1 = datenum('01/01/2002', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

% Formatting
xlabel('Years')
ylabel('CO_{2}-C (g m^{-2} day{-1})')
title('No thermal adaptation')

set(gca, 'FontSize', 14, 'Ytick', [yMin:1:yMax])

lgd = legend([p1, p2], 'Control', 'Heated', 'Location', 'northeast');
legend('box','off');

% ---------------------
% Optimum driven
% ---------------------

subplot(3,1,2)
hold on

% Modelled
d = datenum(dates_all);
p1 = plot(d,totalCO2_optimumDriv_control, 'color', c_control_mod, 'LineWidth', lineWidth);
p2 = plot(d, totalCO2_optimumDriv_warmed, 'color', c_warming_mod, 'LineWidth', lineWidth);

% Measured
d = datenum(dates_CO2_measurements_control);
p3 = errorbar(d, Measured_CO2_flux_control, Measured_CO2_flux_control_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', c_control_meas, 'MarkerEdgeColor', c_control_meas, 'LineStyle', 'none');
p4 = errorbar(d, Measured_CO2_flux_warming, Measured_CO2_flux_warming_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', c_warming_meas, 'MarkerEdgeColor', c_warming_meas, 'LineStyle', 'none');

% The start of warming is indicated
d1 = datenum(date_startWarming);
d2 = datenum(date_startWarming);
plot([d1 d2],[0 yMax], ':', 'color', 'k', 'lineWidth', 2);

date1 = datenum('01/01/2002', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

% Formatting
xlabel('Years')
ylabel('CO_{2}-C (g m^{-2} day{-1})')
title('Optimum driven')

set(gca, 'FontSize', 14, 'Ytick', [yMin:1:yMax])

lgd = legend([p1, p2], 'Control', 'Heated', 'Location', 'northeast');
legend('box','off');

% ---------------------
% Enzyme rigidity
% ---------------------

subplot(3,1,3)
hold on

% Modelled
d = datenum(dates_all);
p1 = plot(d,totalCO2_enzRig_control, 'color', c_control_mod, 'LineWidth', lineWidth);
p2 = plot(d, totalCO2_enzRig_warmed, 'color', c_warming_mod, 'LineWidth', lineWidth);

% Measured
d = datenum(dates_CO2_measurements_control);
p3 = errorbar(d, Measured_CO2_flux_control, Measured_CO2_flux_control_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', c_control_meas, 'MarkerEdgeColor', c_control_meas, 'LineStyle', 'none');
p4 = errorbar(d, Measured_CO2_flux_warming, Measured_CO2_flux_warming_stdError, 'k', 'marker', 'o', 'MarkerSize', dotSize, 'MarkerFaceColor', c_warming_meas, 'MarkerEdgeColor', c_warming_meas, 'LineStyle', 'none');

% The start of warming is indicated
d1 = datenum(date_startWarming);
d2 = datenum(date_startWarming);
plot([d1 d2],[0 yMax], ':', 'color', 'k', 'lineWidth', 2);

date1 = datenum('01/01/2002', 'dd/MM/yyyy');
date2 = datenum('31/12/2016', 'dd/MM/yyyy');

xlim([date1 date2])
datetick('x', 'yyyy','keeplimits')
ylim([yMin yMax])

% Formatting
xlabel('Years')
ylabel('CO_{2}-C (g m^{-2} day{-1})')
title('Enzyme rigidity')

set(gca, 'FontSize', 14, 'Ytick', [yMin:1:yMax])

lgd = legend([p1, p2], 'Control', 'Heated', 'Location', 'northeast');
legend('box','off');

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% Letters
text(0.09, 0.95, '(A)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.09, 0.64, '(B)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.09, 0.34, '(C)', 'FontSize', 14, 'FontWeight', 'bold')

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Daily CO2 flux timeseries_per treatment.tiff' -r300 -nocrop -transparent
% =========================================================================

% Relative differences in the CO2 flux between control and heated for the 
% 2nd and 3rd winter (to compare results to Contosta et al. 2011)

% 2nd winter: 2004 - 2005
date1 = datetime('01/12/2004', 'InputFormat', 'dd/MM/yyyy');
date2 = datetime('31/03/2005', 'InputFormat', 'dd/MM/yyyy');

[r c] = find(dates_all >= date1 & dates_all <= date2);

winterCO2_2005_noAdapt_control = sum(totalCO2_noAdapt_control(c));
winterCO2_2005_noAdapt_warmed = sum(totalCO2_noAdapt_warmed(c));

winterCO2_2005_optDriv_control = sum(totalCO2_optimumDriv_control(c));
winterCO2_2005_optDriv_warmed = sum(totalCO2_optimumDriv_warmed(c));

winterCO2_2005_enzRig_control = sum(totalCO2_enzRig_control(c));
winterCO2_2005_enzRig_warmed = sum(totalCO2_enzRig_warmed(c));

relDiff_2005_noAdapt = (winterCO2_2005_noAdapt_warmed - winterCO2_2005_noAdapt_control)/winterCO2_2005_noAdapt_control * 100;
relDiff_2005_optDriv = (winterCO2_2005_optDriv_warmed - winterCO2_2005_optDriv_control)/winterCO2_2005_optDriv_control * 100;
relDiff_2005_enzRig = (winterCO2_2005_enzRig_warmed - winterCO2_2005_enzRig_control)/winterCO2_2005_enzRig_control * 100;

% 3rd winter: 2004 - 2005
date1 = datetime('01/12/2005', 'InputFormat', 'dd/MM/yyyy');
date2 = datetime('31/03/2006', 'InputFormat', 'dd/MM/yyyy');

[r c] = find(dates_all >= date1 & dates_all <= date2);

winterCO2_2006_noAdapt_control = sum(totalCO2_noAdapt_control(c));
winterCO2_2006_noAdapt_warmed = sum(totalCO2_noAdapt_warmed(c));

winterCO2_2006_optDriv_control = sum(totalCO2_optimumDriv_control(c));
winterCO2_2006_optDriv_warmed = sum(totalCO2_optimumDriv_warmed(c));

winterCO2_2006_enzRig_control = sum(totalCO2_enzRig_control(c));
winterCO2_2006_enzRig_warmed = sum(totalCO2_enzRig_warmed(c));

relDiff_2006_noAdapt = (winterCO2_2006_noAdapt_warmed - winterCO2_2006_noAdapt_control)/winterCO2_2006_noAdapt_control * 100;
relDiff_2006_optDriv = (winterCO2_2006_optDriv_warmed - winterCO2_2006_optDriv_control)/winterCO2_2006_optDriv_control * 100;
relDiff_2006_enzRig = (winterCO2_2006_enzRig_warmed - winterCO2_2006_enzRig_control)/winterCO2_2006_enzRig_control * 100;

% The average for both years
relDiff_noAdapt = mean([relDiff_2005_noAdapt relDiff_2006_noAdapt]);
relDiff_optDriv = mean([relDiff_2005_optDriv relDiff_2006_optDriv]);
relDiff_enzRig = mean([relDiff_2005_enzRig relDiff_2006_enzRig]);

%% Scatterplot of modelled versus measured CO2 per timestep, for different years

% The data for the first 4 years are isolated
% Model results
% The modeled CO2 on measurement days is extracted
dates_all = optimumDrivenData.dates_all;

[r, c] = find(ismember(dates_all, dates_CO2_measurements_control));

mod_CO2_measDays_noAdapt_control = totalCO2_noAdapt_control(c);
mod_CO2_measDays_optimumDriv_control = totalCO2_optimumDriv_control(c);
mod_CO2_measDays_enzRig_control = totalCO2_enzRig_control(c);

[r, c] = find(ismember(dates_all, dates_CO2_measurements_warming));

mod_CO2_measDays_noAdapt_warmed = totalCO2_noAdapt_warmed(c);
mod_CO2_measDays_optimumDriv_warmed = totalCO2_optimumDriv_warmed(c);
mod_CO2_measDays_enzRig_warmed = totalCO2_enzRig_warmed(c);

% The data for the first 4 years is extracted
initialYears = 2003:2006;
allYears = dates_CO2_measurements_control.Year;
[r_init, c] = find(ismember(allYears, initialYears));
[r_later, c] = find(~ismember(allYears, initialYears));

mod_CO2_noAdapt_initialYears_control = mod_CO2_measDays_noAdapt_control(r_init);
mod_CO2_optimumDriv_initialYears_control = mod_CO2_measDays_optimumDriv_control(r_init);
mod_CO2_enzRig_initialYears_control = mod_CO2_measDays_enzRig_control(r_init);

mod_CO2_noAdapt_initialYears_warmed = mod_CO2_measDays_noAdapt_warmed(r_init);
mod_CO2_optimumDriv_initialYears_warmed = mod_CO2_measDays_optimumDriv_warmed(r_init);
mod_CO2_enzRig_initialYears_warmed = mod_CO2_measDays_enzRig_warmed(r_init);

% The data from the remaining years is extracted
mod_CO2_noAdapt_laterYears_control = mod_CO2_measDays_noAdapt_control(r_later);
mod_CO2_optimumDriv_laterYears_control = mod_CO2_measDays_optimumDriv_control(r_later);
mod_CO2_enzRig_laterYears_control = mod_CO2_measDays_enzRig_control(r_later);

mod_CO2_noAdapt_laterYears_warmed = mod_CO2_measDays_noAdapt_warmed(r_later);
mod_CO2_optimumDriv_laterYears_warmed = mod_CO2_measDays_optimumDriv_warmed(r_later);
mod_CO2_enzRig_laterYears_warmed = mod_CO2_measDays_enzRig_warmed(r_later);

% Measurements
meas_CO2_firstForYears_control = Measured_CO2_flux_control(r_init);
meas_CO2_firstForYears_warmed = Measured_CO2_flux_warming(r_init);

meas_CO2_laterYears_control = Measured_CO2_flux_control(r_later);
meas_CO2_laterYears_warmed = Measured_CO2_flux_warming(r_later);

% --------------------
% The data are plotted
% --------------------

dotSize_firstYears = 50;
dotSize_later = 50;
scale = 2;

yMin = 0;
yMax = 5;

c_firstYeart = [227,26,28]./255;
% c_firstYeart = [35,139,69]./255;
c_later = [.2 .2 .2];
% c_later = [.7 .7 .7];
% c_later = 'k';

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [15 15 18*scale 10*scale], 'color', [1 1 1])

% The number of free parameters for the different scenarios
nParam_noAdapt = 16;
nParam_optDriv = 18;
nParam_enzRig = 19;

% --------------------------------
% The control plots - no adapation
% --------------------------------
subplot(2,3,1)
hold on

scatter(meas_CO2_laterYears_control, mod_CO2_noAdapt_laterYears_control, dotSize_later, 'MarkerFaceColor', c_later, 'MarkerEdgeColor', 'none');
scatter(meas_CO2_firstForYears_control, mod_CO2_noAdapt_initialYears_control, dotSize_firstYears, 'MarkerFaceColor', c_firstYeart, 'MarkerEdgeColor', 'none');

xlim([yMin yMax])
ylim([yMin yMax])

% The 1:1 line
hline1 = refline(1,0);
hline1.Color = [.5 .5 .5];
hline1.LineStyle = '--';

% RMSE
RMSE4 = sqrt(mean((meas_CO2_firstForYears_control - mod_CO2_noAdapt_initialYears_control).^2));
RMSE4 = round(RMSE4*100)/100;

RMSE_lastYears_noAdapt_control = sqrt(mean((meas_CO2_laterYears_control - mod_CO2_noAdapt_laterYears_control).^2));

allMeas = [meas_CO2_firstForYears_control; meas_CO2_laterYears_control];
allMod = [mod_CO2_noAdapt_initialYears_control; mod_CO2_noAdapt_laterYears_control];
RMSEall = sqrt(mean((allMeas - allMod).^2));
RMSEall = round(RMSEall*100)/100;
% Tall = text(0.4, 4, ['RMSE_{all} = ' num2str(RMSEall)]);
% set(Tall, 'fontsize', 12);
% Tall = text(0.4, 4.5, ['RMSE_{4 | all} = ' num2str(RMSE4) ' | ' num2str(RMSEall)]);
% set(Tall, 'fontsize', 12);

% Bias
bias4 = sum(mod_CO2_noAdapt_initialYears_control - meas_CO2_firstForYears_control) / numel(meas_CO2_firstForYears_control);
bias4 = round(bias4*100)/100;
biasall = sum(allMod - allMeas) / numel(allMeas);
biasall = round(biasall*100)/100;
% Tall = text(0.4, 4, ['Bias_{4 | all} = ' num2str(bias4) ' | ' num2str(biasall)]);
% set(Tall, 'fontsize', 12);
bias_laterYears_noAdapt_control = sum(mod_CO2_noAdapt_laterYears_control - meas_CO2_laterYears_control) / numel(meas_CO2_laterYears_control);

% Histogram of biases
% histogram(mod_CO2_noAdapt_initialYears_control - meas_CO2_firstForYears_control, 8)
% qqplot(mod_CO2_noAdapt_initialYears_control - meas_CO2_firstForYears_control)

% histogram(allMod - allMeas, 10)
% qqplot(allMod - allMeas)

% Null hypothesis that the data in x comes from a normal distribution with mean equal to zero
% The result h is 1 if the test rejects the null hypothesis at the 5 significance level
h1 = ttest(mod_CO2_noAdapt_initialYears_control - meas_CO2_firstForYears_control)
h2 = ttest(allMod - allMeas)

% The value h = 1 indicates that the test rejects the null hypothesis of zero median.
p1 = signrank(mod_CO2_noAdapt_initialYears_control - meas_CO2_firstForYears_control)
p2 = signrank(allMod - allMeas)

% Slope
pFit = polyfit(meas_CO2_firstForYears_control, mod_CO2_noAdapt_initialYears_control,1);
slope_noAdapt_control_firstYears = pFit(1);

pFit = polyfit(meas_CO2_laterYears_control, mod_CO2_noAdapt_laterYears_control,1);
slope_noAdapt_control_laterYears = pFit(1);

pFit = polyfit(allMeas, allMod,1);
slope_noAdapt_control_allYears = pFit(1);

% ---
% AIC
% ---

% First the SSR is calculated
SSR_4_ctrl_noAdapt = sum((meas_CO2_firstForYears_control - mod_CO2_noAdapt_initialYears_control).^2);
SSR_end_ctrl_noAdapt = sum((meas_CO2_laterYears_control - mod_CO2_noAdapt_laterYears_control).^2);
SSR_all_ctrl_noAdapt = sum((allMeas - allMod).^2);

% AIC
AIC_4_ctrl_noAdapt = log(SSR_4_ctrl_noAdapt/numel(meas_CO2_firstForYears_control)) + (2 * nParam_noAdapt);
AIC_end_ctrl_noAdapt = log(SSR_end_ctrl_noAdapt/numel(meas_CO2_laterYears_control)) + (2 * nParam_noAdapt);
AIC_all_ctrl_noAdapt = log(SSR_all_ctrl_noAdapt/numel(allMeas)) + (2 * nParam_noAdapt);

AIC_4_ctrl_noAdapt = round(AIC_4_ctrl_noAdapt*100)/100;
AIC_end_ctrl_noAdapt = round(AIC_end_ctrl_noAdapt*100)/100;
AIC_all_ctrl_noAdapt = round(AIC_all_ctrl_noAdapt*100)/100;

set(gca, 'FontSize', 12, 'Ytick', [yMin:1:yMax], 'box', 'on')

xlabel('Measured F_{CO_{2}} (gC m^{-2} d^{-1})')
ylabel('Modelled F_{CO_{2}} (gC m^{-2} d^{-1})')
title('Control - no thermal adaptation')

% Adding the measures
text(3, 1.5, "Years")
text(4.2, 1.5, "1-4", 'HorizontalAlignment', 'right')
text(4.35, 1.5, "1-13", 'HorizontalAlignment', 'left')

plot([3.75 4.95], [1.3 1.3], 'color', [0 0 0]) % Horizontal line
plot([4.28 4.28], [1.35 1.6], 'color', [0 0 0]) % Vertical line top
plot([4.28 4.28], [0.2 1.25], 'color', [0 0 0]) % Vertical line bottom

text(3, 1.1, "MAE:")
text(4.2, 1.1, append(num2str(bias4), "*"), 'HorizontalAlignment', 'right')
text(4.35, 1.1, append(num2str(biasall), "*"), 'HorizontalAlignment', 'left')

text(3, 0.75, "RMSE:")
text(4.2, 0.75, num2str(RMSE4), 'HorizontalAlignment', 'right')
text(4.35, 0.75, num2str(RMSEall), 'HorizontalAlignment', 'left')

text(3, 0.4, "AIC:")
text(4.2, 0.4, num2str(AIC_4_ctrl_noAdapt), 'HorizontalAlignment', 'right')
text(4.35, 0.4, num2str(AIC_all_ctrl_noAdapt), 'HorizontalAlignment', 'left')

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.05; % x
pos(2) = 0.62; % y
pos(3) = 0.25; % width
pos(4) = 0.32; % height
set(gca,'position',pos)

% --------------------------------
% The control plots - Optimum driven
% --------------------------------
subplot(2,3,2)
hold on

scatter(meas_CO2_laterYears_control, mod_CO2_optimumDriv_laterYears_control, dotSize_later, 'MarkerFaceColor', c_later, 'MarkerEdgeColor', 'none');
scatter(meas_CO2_firstForYears_control, mod_CO2_optimumDriv_initialYears_control, dotSize_firstYears, 'MarkerFaceColor', c_firstYeart, 'MarkerEdgeColor', 'none');

xlim([yMin yMax])
ylim([yMin yMax])

% The 1:1 line
hline1 = refline(1,0);
hline1.Color = [.5 .5 .5];
hline1.LineStyle = '--';

% RMSE
RMSE4 = sqrt(mean((meas_CO2_firstForYears_control - mod_CO2_optimumDriv_initialYears_control).^2));
RMSE4 = round(RMSE4*100)/100;

RMSE_lastYears_optDriv_control = sqrt(mean((meas_CO2_laterYears_control - mod_CO2_optimumDriv_laterYears_control).^2));

allMeas = [meas_CO2_firstForYears_control; meas_CO2_laterYears_control];
allMod = [mod_CO2_optimumDriv_initialYears_control; mod_CO2_optimumDriv_laterYears_control];
RMSEall = sqrt(mean((allMeas - allMod).^2));
RMSEall = round(RMSEall*100)/100;
% Tall = text(0.4, 4.5, ['RMSE_{4 | all} = ' num2str(RMSE4) ' | ' num2str(RMSEall)]);
% set(Tall, 'fontsize', 12);

% Bias
bias4 = sum(mod_CO2_optimumDriv_initialYears_control - meas_CO2_firstForYears_control) / numel(meas_CO2_firstForYears_control);
bias4 = round(bias4*100)/100;
biasall = sum(allMod - allMeas) / numel(allMeas);
biasall = round(biasall*100)/100;
% Tall = text(0.4, 4, ['Bias_{4 | all} = ' num2str(bias4) ' | ' num2str(biasall)]);
% set(Tall, 'fontsize', 12);
bias_laterYears_optDriv_control = sum(mod_CO2_optimumDriv_laterYears_control - meas_CO2_laterYears_control) / numel(meas_CO2_laterYears_control);

% Histogram of biases
% histogram(mod_CO2_optimumDriv_initialYears_control - meas_CO2_firstForYears_control, 8)
% qqplot(mod_CO2_optimumDriv_initialYears_control - meas_CO2_firstForYears_control)

% histogram(allMod - allMeas, 10)
% qqplot(allMod - allMeas)

% Null hypothesis that the data in x comes from a normal distribution with mean equal to zero
% The result h is 1 if the test rejects the null hypothesis at the 5 significance level
h1 = ttest(mod_CO2_optimumDriv_initialYears_control - meas_CO2_firstForYears_control)
h2 = ttest(allMod - allMeas)

% The value h = 1 indicates that the test rejects the null hypothesis of zero median.
p1 = signrank(mod_CO2_optimumDriv_initialYears_control - meas_CO2_firstForYears_control)
p2 = signrank(allMod - allMeas)

% Slope
pFit = polyfit(meas_CO2_firstForYears_control, mod_CO2_optimumDriv_initialYears_control,1);
slope_optDriv_control_firstYears = pFit(1);

pFit = polyfit(meas_CO2_laterYears_control, mod_CO2_optimumDriv_laterYears_control,1);
slope_optDriv_control_laterYears = pFit(1);

pFit = polyfit(allMeas, allMod,1);
slope_optDriv_control_allYears = pFit(1);

% ---
% AIC
% ---

% First the SSR is calculated
SSR_4_ctrl_optDriv = sum((meas_CO2_firstForYears_control - mod_CO2_optimumDriv_initialYears_control).^2);
SSR_end_ctrl_optDriv = sum((meas_CO2_laterYears_control - mod_CO2_optimumDriv_laterYears_control).^2);
SSR_all_ctrl_optDriv = sum((allMeas - allMod).^2);

% AIC
AIC_4_ctrl_optDriv = log(SSR_4_ctrl_optDriv/numel(meas_CO2_firstForYears_control)) + (2 * nParam_optDriv);
AIC_end_ctrl_optDriv = log(SSR_end_ctrl_optDriv/numel(meas_CO2_laterYears_control)) + (2 * nParam_optDriv);
AIC_all_ctrl_optDriv = log(SSR_all_ctrl_optDriv/numel(allMeas)) + (2 * nParam_optDriv);

AIC_4_ctrl_optDriv = round(AIC_4_ctrl_optDriv*100)/100;
AIC_end_ctrl_optDriv = round(AIC_end_ctrl_optDriv*100)/100;
AIC_all_ctrl_optDriv = round(AIC_all_ctrl_optDriv*100)/100;

set(gca, 'FontSize', 12, 'Ytick', [yMin:1:yMax], 'box', 'on')

xlabel('Measured F_{CO_{2}} (gC m^{-2} d^{-1})')
ylabel('Modelled F_{CO_{2}} (gC m^{-2} d^{-1})')
title('Control - optimum driven')

% Adding the measures
text(3, 1.5, "Years")
text(4.2, 1.5, "1-4", 'HorizontalAlignment', 'right')
text(4.35, 1.5, "1-13", 'HorizontalAlignment', 'left')

plot([3.75 4.95], [1.3 1.3], 'color', [0 0 0]) % Horizontal line
plot([4.28 4.28], [1.35 1.6], 'color', [0 0 0]) % Vertical line top
plot([4.28 4.28], [0.2 1.25], 'color', [0 0 0]) % Vertical line bottom

text(3, 1.1, "MAE:")
text(4.2, 1.1, append(num2str(bias4), "*"), 'HorizontalAlignment', 'right')
text(4.35, 1.1, append(num2str(biasall), "*"), 'HorizontalAlignment', 'left')

text(3, 0.75, "RMSE:")
text(4.2, 0.75, num2str(RMSE4), 'HorizontalAlignment', 'right')
text(4.35, 0.75, num2str(RMSEall), 'HorizontalAlignment', 'left')

text(3, 0.4, "AIC:")
text(4.2, 0.4, num2str(AIC_4_ctrl_optDriv), 'HorizontalAlignment', 'right')
text(4.35, 0.4, num2str(AIC_all_ctrl_optDriv), 'HorizontalAlignment', 'left')

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.38; % x
pos(2) = 0.62; % y
pos(3) = 0.25; % width
pos(4) = 0.32; % height
set(gca,'position',pos)

% --------------------------------
% The control plots - Enzyme rigidity
% --------------------------------
subplot(2,3,3)
hold on

scatter(meas_CO2_laterYears_control, mod_CO2_enzRig_laterYears_control, dotSize_later, 'MarkerFaceColor', c_later, 'MarkerEdgeColor', 'none');
scatter(meas_CO2_firstForYears_control, mod_CO2_enzRig_initialYears_control, dotSize_firstYears, 'MarkerFaceColor', c_firstYeart, 'MarkerEdgeColor', 'none');

xlim([yMin yMax])
ylim([yMin yMax])

% The 1:1 line
hline1 = refline(1,0);
hline1.Color = [.5 .5 .5];
hline1.LineStyle = '--';

% RMSE
RMSE4 = sqrt(mean((meas_CO2_firstForYears_control - mod_CO2_enzRig_initialYears_control).^2));
RMSE4 = round(RMSE4*100)/100;

RMSE_lastYears_enzRig_control = sqrt(mean((meas_CO2_laterYears_control - mod_CO2_enzRig_laterYears_control).^2));

allMeas = [meas_CO2_firstForYears_control; meas_CO2_laterYears_control];
allMod = [mod_CO2_enzRig_initialYears_control; mod_CO2_enzRig_laterYears_control];
RMSEall = sqrt(mean((allMeas - allMod).^2));
RMSEall = round(RMSEall*100)/100;
% Tall = text(0.4, 4.5, ['RMSE_{4 | all} = ' num2str(RMSE4) ' | ' num2str(RMSEall)]);
% set(Tall, 'fontsize', 12);

% Bias
bias4 = sum(mod_CO2_enzRig_initialYears_control - meas_CO2_firstForYears_control) / numel(meas_CO2_firstForYears_control);
bias4 = round(bias4*100)/100;
biasall = sum(allMod - allMeas) / numel(allMeas);
biasall = round(biasall*100)/100;
% Tall = text(0.4, 4, ['Bias_{4 | all} = ' num2str(bias4) ' | ' num2str(biasall)]);
% set(Tall, 'fontsize', 12);
bias_laterYears_enzRig_control = sum(mod_CO2_enzRig_laterYears_control - meas_CO2_laterYears_control) / numel(meas_CO2_laterYears_control);

% Histogram of biases
% histogram(mod_CO2_enzRig_initialYears_control - meas_CO2_firstForYears_control, 8)
% qqplot(mod_CO2_enzRig_initialYears_control - meas_CO2_firstForYears_control)

% histogram(allMod - allMeas, 10)
% qqplot(allMod - allMeas)

% Null hypothesis that the data in x comes from a normal distribution with mean equal to zero
% The result h is 1 if the test rejects the null hypothesis at the 5 significance level
h1 = ttest(mod_CO2_enzRig_initialYears_control - meas_CO2_firstForYears_control)
h2 = ttest(allMod - allMeas)

% The value h = 1 indicates that the test rejects the null hypothesis of zero median.
p1 = signrank(mod_CO2_enzRig_initialYears_control - meas_CO2_firstForYears_control)
p2 = signrank(allMod - allMeas)

% Slope
pFit = polyfit(meas_CO2_firstForYears_control, mod_CO2_enzRig_initialYears_control,1);
slope_enzRig_control_firstYears = pFit(1);

pFit = polyfit(meas_CO2_laterYears_control, mod_CO2_enzRig_laterYears_control,1);
slope_enzRig_control_laterYears = pFit(1);

pFit = polyfit(allMeas, allMod,1);
slope_enzRig_control_allYears = pFit(1);

% ---
% AIC
% ---

% First the SSR is calculated
SSR_4_ctrl_enzRig = sum((meas_CO2_firstForYears_control - mod_CO2_enzRig_initialYears_control).^2);
SSR_end_ctrl_enzRig = sum((meas_CO2_laterYears_control - mod_CO2_enzRig_laterYears_control).^2);
SSR_all_ctrl_enzRig = sum((allMeas - allMod).^2);

% AIC
AIC_4_ctrl_enzRig = log(SSR_4_ctrl_enzRig/numel(meas_CO2_firstForYears_control)) + (2 * nParam_enzRig);
AIC_end_ctrl_enzRig = log(SSR_end_ctrl_enzRig/numel(meas_CO2_laterYears_control)) + (2 * nParam_enzRig);
AIC_all_ctrl_enzRig = log(SSR_all_ctrl_enzRig/numel(allMeas)) + (2 * nParam_enzRig);

AIC_4_ctrl_enzRig = round(AIC_4_ctrl_enzRig*100)/100;
AIC_end_ctrl_enzRig = round(AIC_end_ctrl_enzRig*100)/100;
AIC_all_ctrl_enzRig = round(AIC_all_ctrl_enzRig*100)/100;

set(gca, 'FontSize', 12, 'Ytick', [yMin:1:yMax], 'box', 'on')

xlabel('Measured F_{CO_{2}} (gC m^{-2} d^{-1})')
ylabel('Modelled F_{CO_{2}} (gC m^{-2} d^{-1})')
title('Control - enzyme rigidity')

% Adding the measures
text(3, 1.5, "Years")
text(4.2, 1.5, "1-4", 'HorizontalAlignment', 'right')
text(4.35, 1.5, "1-13", 'HorizontalAlignment', 'left')

plot([3.75 4.95], [1.3 1.3], 'color', [0 0 0]) % Horizontal line
plot([4.28 4.28], [1.35 1.6], 'color', [0 0 0]) % Vertical line top
plot([4.28 4.28], [0.2 1.25], 'color', [0 0 0]) % Vertical line bottom

text(3, 1.1, "MAE:")
text(4.2, 1.1, append(num2str(bias4), "*"), 'HorizontalAlignment', 'right')
text(4.35, 1.1, append(num2str(biasall), "*"), 'HorizontalAlignment', 'left')

text(3, 0.75, "RMSE:")
text(4.2, 0.75, num2str(RMSE4), 'HorizontalAlignment', 'right')
text(4.35, 0.75, num2str(RMSEall), 'HorizontalAlignment', 'left')

text(3, 0.4, "AIC:")
text(4.2, 0.4, num2str(AIC_4_ctrl_enzRig), 'HorizontalAlignment', 'right')
text(4.35, 0.4, num2str(AIC_all_ctrl_enzRig), 'HorizontalAlignment', 'left')

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.71; % x
pos(2) = 0.62; % y
pos(3) = 0.25; % width
pos(4) = 0.32; % height
set(gca,'position',pos)

% --------------------------------
% The heated plots - no adapation
% --------------------------------
subplot(2,3,4)
hold on

scatter(meas_CO2_laterYears_warmed, mod_CO2_noAdapt_laterYears_warmed, dotSize_later, 'MarkerFaceColor', c_later, 'MarkerEdgeColor', 'none');
scatter(meas_CO2_firstForYears_warmed, mod_CO2_noAdapt_initialYears_warmed, dotSize_firstYears, 'MarkerFaceColor', c_firstYeart, 'MarkerEdgeColor', 'none');

xlim([yMin yMax])
ylim([yMin yMax])

% The 1:1 line
hline1 = refline(1,0);
hline1.Color = [.5 .5 .5];
hline1.LineStyle = '--';

% RMSE
RMSE4 = sqrt(mean((meas_CO2_firstForYears_warmed - mod_CO2_noAdapt_initialYears_warmed).^2));
RMSE4 = round(RMSE4*100)/100;

RMSE_lastYears_noAdapt_warmed = sqrt(mean((meas_CO2_laterYears_warmed - mod_CO2_noAdapt_laterYears_warmed).^2));

allMeas = [meas_CO2_firstForYears_warmed; meas_CO2_laterYears_warmed];
allMod = [mod_CO2_noAdapt_initialYears_warmed; mod_CO2_noAdapt_laterYears_warmed];
RMSEall = sqrt(mean((allMeas - allMod).^2));
RMSEall = round(RMSEall*100)/100;
% Tall = text(0.4, 4.5, ['RMSE_{4 | all} = ' num2str(RMSE4) ' | ' num2str(RMSEall)]);
% set(Tall, 'fontsize', 12);

% Bias
bias4 = sum(mod_CO2_noAdapt_initialYears_warmed - meas_CO2_firstForYears_warmed) / numel(meas_CO2_firstForYears_warmed);
bias4 = round(bias4*100)/100;
biasall = sum(allMod - allMeas) / numel(allMeas);
biasall = round(biasall*100)/100;
% Tall = text(0.4, 4, ['Bias_{4 | all} = ' num2str(bias4) ' | ' num2str(biasall)]);
% set(Tall, 'fontsize', 12);
bias_laterYears_noAdapt_warmed = sum(mod_CO2_noAdapt_laterYears_warmed - meas_CO2_laterYears_warmed) / numel(meas_CO2_laterYears_warmed);

% Histogram of biases
% histogram(mod_CO2_noAdapt_initialYears_warmed - meas_CO2_firstForYears_warmed, 8)
% qqplot(mod_CO2_noAdapt_initialYears_warmed - meas_CO2_firstForYears_warmed)

% histogram(allMod - allMeas, 10)
% qqplot(allMod - allMeas)

% Null hypothesis that the data in x comes from a normal distribution with mean equal to zero
% The result h is 1 if the test rejects the null hypothesis at the 5 significance level
h1 = ttest(mod_CO2_noAdapt_initialYears_warmed - meas_CO2_firstForYears_warmed)
h2 = ttest(allMod - allMeas)

% The value h = 1 indicates that the test rejects the null hypothesis of zero median.
p1 = signrank(mod_CO2_noAdapt_initialYears_warmed - meas_CO2_firstForYears_warmed)
p2 = signrank(allMod - allMeas)

% Slope
pFit = polyfit(meas_CO2_firstForYears_warmed, mod_CO2_noAdapt_initialYears_warmed,1);
slope_noAdapt_warmed_firstYears = pFit(1);

pFit = polyfit(meas_CO2_laterYears_warmed, mod_CO2_noAdapt_laterYears_warmed,1);
slope_noAdapt_warmed_laterYears = pFit(1);

pFit = polyfit(allMeas, allMod,1);
slope_noAdapt_warmed_allYears = pFit(1);

% ---
% AIC
% ---

% First the SSR is calculated
SSR_4_heated_noAdapt = sum((meas_CO2_firstForYears_warmed - mod_CO2_noAdapt_initialYears_warmed).^2);
SSR_end_heated_noAdapt = sum((meas_CO2_laterYears_warmed - mod_CO2_noAdapt_laterYears_warmed).^2);
SSR_all_heated_noAdapt = sum((allMeas - allMod).^2);

% AIC
AIC_4_heated_noAdapt = log(SSR_4_heated_noAdapt/numel(meas_CO2_firstForYears_warmed)) + (2 * nParam_noAdapt);
AIC_end_heated_noAdapt = log(SSR_end_heated_noAdapt/numel(meas_CO2_laterYears_warmed)) + (2 * nParam_noAdapt);
AIC_all_heated_noAdapt = log(SSR_all_heated_noAdapt/numel(allMeas)) + (2 * nParam_noAdapt);

AIC_4_heated_noAdapt = round(AIC_4_heated_noAdapt*100)/100;
AIC_end_heated_noAdapt = round(AIC_end_heated_noAdapt*100)/100;
AIC_all_heated_noAdapt = round(AIC_all_heated_noAdapt*100)/100;

%BIC
BIC_4_heated_noAdapt = log(SSR_4_heated_noAdapt/numel(meas_CO2_firstForYears_warmed)) + ((2*numel(meas_CO2_firstForYears_warmed)*(nParam_noAdapt+1))/(numel(meas_CO2_firstForYears_warmed) - nParam_noAdapt - 2));
BIC_end_heated_noAdapt = log(SSR_end_heated_noAdapt/numel(meas_CO2_laterYears_warmed)) + ((2*numel(meas_CO2_laterYears_warmed)*(nParam_noAdapt+1))/(numel(meas_CO2_laterYears_warmed) - nParam_noAdapt - 2));
BIC_all_heated_noAdapt = log(SSR_all_heated_noAdapt/numel(allMeas)) + ((2*numel(allMeas)*(nParam_noAdapt+1))/(numel(allMeas) - nParam_noAdapt - 2));

set(gca, 'FontSize', 12, 'Ytick', [yMin:1:yMax], 'box', 'on')

xlabel('Measured F_{CO_{2}} (gC m^{-2} d^{-1})')
ylabel('Modelled F_{CO_{2}} (gC m^{-2} d^{-1})')
title('Heated - no thermal adaptation')

% Adding the measures
text(3, 1.5, "Years")
text(4.2, 1.5, "1-4", 'HorizontalAlignment', 'right')
text(4.35, 1.5, "1-13", 'HorizontalAlignment', 'left')

plot([3.75 4.95], [1.3 1.3], 'color', [0 0 0]) % Horizontal line
plot([4.28 4.28], [1.35 1.6], 'color', [0 0 0]) % Vertical line top
plot([4.28 4.28], [0.2 1.25], 'color', [0 0 0]) % Vertical line bottom

text(3, 1.1, "MAE:")
text(4.2, 1.1, append(num2str(bias4)), 'HorizontalAlignment', 'right')
text(4.35, 1.1, append(num2str(biasall), "*"), 'HorizontalAlignment', 'left')

text(3, 0.75, "RMSE:")
text(4.2, 0.75, num2str(RMSE4), 'HorizontalAlignment', 'right')
text(4.35, 0.75, num2str(RMSEall), 'HorizontalAlignment', 'left')

text(3, 0.4, "AIC:")
text(4.2, 0.4, num2str(AIC_4_heated_noAdapt), 'HorizontalAlignment', 'right')
text(4.35, 0.4, num2str(AIC_all_heated_noAdapt), 'HorizontalAlignment', 'left')

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.05; % x
pos(2) = 0.15; % y
pos(3) = 0.25; % width
pos(4) = 0.33; % height
set(gca,'position',pos)

% --------------------------------
% The heated plots - Optimum driven
% --------------------------------
subplot(2,3,5)
hold on

scatter(meas_CO2_laterYears_warmed, mod_CO2_optimumDriv_laterYears_warmed, dotSize_later, 'MarkerFaceColor', c_later, 'MarkerEdgeColor', 'none');
scatter(meas_CO2_firstForYears_warmed, mod_CO2_optimumDriv_initialYears_warmed, dotSize_firstYears, 'MarkerFaceColor', c_firstYeart, 'MarkerEdgeColor', 'none');

xlim([yMin yMax])
ylim([yMin yMax])

% The 1:1 line
hline1 = refline(1,0);
hline1.Color = [.5 .5 .5];
hline1.LineStyle = '--';

% RMSE
RMSE4 = sqrt(mean((meas_CO2_firstForYears_warmed - mod_CO2_optimumDriv_initialYears_warmed).^2));
RMSE4 = round(RMSE4*100)/100;

RMSE_lastYears_optDriv_warmed = sqrt(mean((meas_CO2_laterYears_warmed - mod_CO2_optimumDriv_laterYears_warmed).^2));

allMeas = [meas_CO2_firstForYears_warmed; meas_CO2_laterYears_warmed];
allMod = [mod_CO2_optimumDriv_initialYears_warmed; mod_CO2_optimumDriv_laterYears_warmed];
RMSEall = sqrt(mean((allMeas - allMod).^2));
RMSEall = round(RMSEall*100)/100;
% Tall = text(0.4, 4.5, ['RMSE_{4 | all} = ' num2str(RMSE4) ' | ' num2str(RMSEall)]);
% set(Tall, 'fontsize', 12);

% Bias
bias4 = sum(mod_CO2_optimumDriv_initialYears_warmed - meas_CO2_firstForYears_warmed) / numel(meas_CO2_firstForYears_warmed);
bias4 = round(bias4*100)/100;
biasall = sum(allMod - allMeas) / numel(allMeas);
biasall = round(biasall*100)/100;
% Tall = text(0.4, 4, ['Bias_{4 | all} = ' num2str(bias4) ' | ' num2str(biasall)]);
% set(Tall, 'fontsize', 12);
bias_laterYears_optDriv_warmed = sum(mod_CO2_optimumDriv_laterYears_warmed - meas_CO2_laterYears_warmed) / numel(meas_CO2_laterYears_warmed);

% Histogram of biases
% histogram(mod_CO2_optimumDriv_initialYears_warmed - meas_CO2_firstForYears_warmed, 8)
% qqplot(mod_CO2_optimumDriv_initialYears_warmed - meas_CO2_firstForYears_warmed)

% histogram(allMod - allMeas, 10)
% qqplot(allMod - allMeas)

% Null hypothesis that the data in x comes from a normal distribution with mean equal to zero
% The result h is 1 if the test rejects the null hypothesis at the 5 significance level
h1 = ttest(mod_CO2_optimumDriv_initialYears_warmed - meas_CO2_firstForYears_warmed)
h2 = ttest(allMod - allMeas)

% The value h = 1 indicates that the test rejects the null hypothesis of zero median.
p1 = signrank(mod_CO2_optimumDriv_initialYears_warmed - meas_CO2_firstForYears_warmed)
p2 = signrank(allMod - allMeas)

% Slope
pFit = polyfit(meas_CO2_firstForYears_warmed, mod_CO2_optimumDriv_initialYears_warmed,1);
slope_optDriv_warmed_firstYears = pFit(1);

pFit = polyfit(meas_CO2_laterYears_warmed, mod_CO2_optimumDriv_laterYears_warmed,1);
slope_optDriv_warmed_laterYears = pFit(1);

pFit = polyfit(allMeas, allMod,1);
slope_optDriv_warmed_allYears = pFit(1);

% ---
% AIC
% ---

% First the SSR is calculated
SSR_4_heated_optDriv = sum((meas_CO2_firstForYears_warmed - mod_CO2_optimumDriv_initialYears_warmed).^2);
SSR_end_heated_optDriv = sum((meas_CO2_laterYears_warmed - mod_CO2_optimumDriv_laterYears_warmed).^2);
SSR_all_heated_optDriv = sum((allMeas - allMod).^2);

% AIC
AIC_4_heated_optDriv = log(SSR_4_heated_optDriv/numel(meas_CO2_firstForYears_warmed)) + (2 * nParam_optDriv);
AIC_end_heated_optDriv = log(SSR_end_heated_optDriv/numel(meas_CO2_laterYears_warmed)) + (2 * nParam_optDriv);
AIC_all_heated_optDriv = log(SSR_all_heated_optDriv/numel(allMeas)) + (2 * nParam_optDriv);

AIC_4_heated_optDriv = round(AIC_4_heated_optDriv*100)/100;
AIC_end_heated_optDriv = round(AIC_end_heated_optDriv*100)/100;
AIC_all_heated_optDriv = round(AIC_all_heated_optDriv*100)/100;

%BIC
BIC_4_heated_optDriv = log(SSR_4_heated_optDriv/numel(meas_CO2_firstForYears_warmed)) + ((2*numel(meas_CO2_firstForYears_warmed)*(nParam_optDriv+1))/(numel(meas_CO2_firstForYears_warmed) - nParam_optDriv - 2));
BIC_end_heated_optDriv = log(SSR_end_heated_optDriv/numel(meas_CO2_laterYears_warmed)) + ((2*numel(meas_CO2_laterYears_warmed)*(nParam_optDriv+1))/(numel(meas_CO2_laterYears_warmed) - nParam_optDriv - 2));
BIC_all_optDriv = log(SSR_all_heated_optDriv/numel(allMeas)) + ((2*numel(allMeas)*(nParam_optDriv+1))/(numel(allMeas) - nParam_optDriv - 2));

set(gca, 'FontSize', 12, 'Ytick', [yMin:1:yMax], 'box', 'on')

xlabel('Measured F_{CO_{2}} (gC m^{-2} d^{-1})')
ylabel('Modelled F_{CO_{2}} (gC m^{-2} d^{-1})')
title('Heated - optimum driven')

% Adding the measures
text(3, 1.5, "Years")
text(4.2, 1.5, "1-4", 'HorizontalAlignment', 'right')
text(4.35, 1.5, "1-13", 'HorizontalAlignment', 'left')

plot([3.75 4.95], [1.3 1.3], 'color', [0 0 0]) % Horizontal line
plot([4.28 4.28], [1.35 1.6], 'color', [0 0 0]) % Vertical line top
plot([4.28 4.28], [0.2 1.25], 'color', [0 0 0]) % Vertical line bottom

text(3, 1.1, "MAE:")
text(4.2, 1.1, append(num2str(bias4), "*"), 'HorizontalAlignment', 'right')
text(4.35, 1.1, append(num2str(biasall), "*"), 'HorizontalAlignment', 'left')

text(3, 0.75, "RMSE:")
text(4.2, 0.75, num2str(RMSE4), 'HorizontalAlignment', 'right')
text(4.35, 0.75, num2str(RMSEall), 'HorizontalAlignment', 'left')

text(3, 0.4, "AIC:")
text(4.2, 0.4, num2str(AIC_4_heated_optDriv), 'HorizontalAlignment', 'right')
text(4.35, 0.4, num2str(AIC_all_heated_optDriv), 'HorizontalAlignment', 'left')

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.38; % x
pos(2) = 0.15; % y
pos(3) = 0.25; % width
pos(4) = 0.33; % height
set(gca,'position',pos)

% --------------------------------
% The heated plots - Enzyme rigidity
% --------------------------------
subplot(2,3,6)
hold on

s1 = scatter(meas_CO2_laterYears_warmed, mod_CO2_enzRig_laterYears_warmed, dotSize_later, 'MarkerFaceColor', c_later, 'MarkerEdgeColor', 'none');
s2 = scatter(meas_CO2_firstForYears_warmed, mod_CO2_enzRig_initialYears_warmed, dotSize_firstYears, 'MarkerFaceColor', c_firstYeart, 'MarkerEdgeColor', 'none');

xlim([yMin yMax])
ylim([yMin yMax])

% The 1:1 line
hline1 = refline(1,0);
hline1.Color = [.5 .5 .5];
hline1.LineStyle = '--';

% RMSE
RMSE4 = sqrt(mean((meas_CO2_firstForYears_warmed - mod_CO2_enzRig_initialYears_warmed).^2));
RMSE4 = round(RMSE4*100)/100;

RMSE_lastYears_enzRig_warmed = sqrt(mean((meas_CO2_laterYears_warmed - mod_CO2_enzRig_laterYears_warmed).^2));

allMeas = [meas_CO2_firstForYears_warmed; meas_CO2_laterYears_warmed];
allMod = [mod_CO2_enzRig_initialYears_warmed; mod_CO2_enzRig_laterYears_warmed];
RMSEall = sqrt(mean((allMeas - allMod).^2));
RMSEall = round(RMSEall*100)/100;
% Tall = text(0.4, 4.5, ['RMSE_{4 | all} = ' num2str(RMSE4) ' | ' num2str(RMSEall)]);
% set(Tall, 'fontsize', 12);

% Bias
bias4 = sum(mod_CO2_enzRig_initialYears_warmed - meas_CO2_firstForYears_warmed) / numel(meas_CO2_firstForYears_warmed);
bias4 = round(bias4*100)/100;
biasall = sum(allMod - allMeas) / numel(allMeas);
biasall = round(biasall*100)/100;
% Tall = text(0.4, 4, ['Bias_{4 | all} = ' num2str(bias4) ' | ' num2str(biasall)]);
% set(Tall, 'fontsize', 12);
bias_laterYears_enzRig_warmed = sum(mod_CO2_enzRig_laterYears_warmed - meas_CO2_laterYears_warmed) / numel(meas_CO2_laterYears_warmed);

% Histogram of biases
% histogram(mod_CO2_enzRig_initialYears_warmed - meas_CO2_firstForYears_warmed, 8)
% qqplot(mod_CO2_enzRig_initialYears_warmed - meas_CO2_firstForYears_warmed)

% histogram(allMod - allMeas, 10)
% qqplot(allMod - allMeas)

% Null hypothesis that the data in x comes from a normal distribution with mean equal to zero
% The result h is 1 if the test rejects the null hypothesis at the 5 significance level
h1 = ttest(mod_CO2_enzRig_initialYears_warmed - meas_CO2_firstForYears_warmed)
h2 = ttest(allMod - allMeas)

% The value h = 1 indicates that the test rejects the null hypothesis of zero median.
p1 = signrank(mod_CO2_enzRig_initialYears_warmed - meas_CO2_firstForYears_warmed)
p2 = signrank(allMod - allMeas)

% Slope
pFit = polyfit(meas_CO2_firstForYears_warmed, mod_CO2_enzRig_initialYears_warmed,1);
slope_enzRig_warmed_firstYears = pFit(1);

pFit = polyfit(meas_CO2_laterYears_warmed, mod_CO2_enzRig_laterYears_warmed,1);
slope_enzRig_warmed_laterYears = pFit(1);

pFit = polyfit(allMeas, allMod,1);
slope_enzRig_warmed_allYears = pFit(1);

% ---
% AIC
% ---

% First the SSR is calculated
SSR_4_heated_enzRig = sum((meas_CO2_firstForYears_warmed - mod_CO2_enzRig_initialYears_warmed).^2);
SSR_end_heated_enzRig = sum((meas_CO2_laterYears_warmed - mod_CO2_enzRig_laterYears_warmed).^2);
SSR_all_heated_enzRig = sum((allMeas - allMod).^2);

% AIC
AIC_4_heated_enzRig = log(SSR_4_heated_enzRig/numel(meas_CO2_firstForYears_warmed)) + (2 * nParam_enzRig);
AIC_end_heated_enzRig = log(SSR_end_heated_enzRig/numel(meas_CO2_laterYears_warmed)) + (2 * nParam_enzRig);
AIC_all_heated_enzRig = log(SSR_all_heated_enzRig/numel(allMeas)) + (2 * nParam_enzRig);

AIC_4_heated_enzRig = round(AIC_4_heated_enzRig*100)/100;
AIC_end_heated_enzRig = round(AIC_end_heated_enzRig*100)/100;
AIC_all_heated_enzRig = round(AIC_all_heated_enzRig*100)/100;

set(gca, 'FontSize', 12, 'Ytick', [yMin:1:yMax], 'box', 'on')

xlabel('Measured F_{CO_{2}} (gC m^{-2} d^{-1})')
ylabel('Modelled F_{CO_{2}} (gC m^{-2} d^{-1})')
title('Heated - enzyme rigidity')

% Adding the measures
text(3, 1.5, "Years")
text(4.2, 1.5, "1-4", 'HorizontalAlignment', 'right')
text(4.35, 1.5, "1-13", 'HorizontalAlignment', 'left')

plot([3.75 4.95], [1.3 1.3], 'color', [0 0 0]) % Horizontal line
plot([4.28 4.28], [1.35 1.6], 'color', [0 0 0]) % Vertical line top
plot([4.28 4.28], [0.2 1.25], 'color', [0 0 0]) % Vertical line bottom

text(3, 1.1, "MAE:")
text(4.2, 1.1, append(num2str(bias4), "*"), 'HorizontalAlignment', 'right')
text(4.35, 1.1, append(num2str(biasall), "*"), 'HorizontalAlignment', 'left')

text(3, 0.75, "RMSE:")
text(4.2, 0.75, num2str(RMSE4), 'HorizontalAlignment', 'right')
text(4.35, 0.75, num2str(RMSEall), 'HorizontalAlignment', 'left')

text(3, 0.4, "AIC:")
text(4.2, 0.4, num2str(AIC_4_heated_enzRig), 'HorizontalAlignment', 'right')
text(4.35, 0.4, num2str(AIC_all_heated_enzRig), 'HorizontalAlignment', 'left')

% The location of the plot is optimized
ax = gca;
pos = ax.Position;
pos(1) = 0.71; % x
pos(2) = 0.15; % y
pos(3) = 0.25; % width
pos(4) = 0.33; % height
set(gca,'position',pos)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
% text(0.41,0.98,'Organic carbon stocks','fontweight','bold','fontsize',14)

% Letters
text(0.015, 0.96, '(A)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.345, 0.96, '(B)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.675, 0.96, '(C)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.015, 0.5, '(D)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.345, 0.5, '(E)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.675, 0.5, '(F)', 'FontSize', 14, 'FontWeight', 'bold')

% The legend
l = legend([s2, s1], 'Data for first 4 years of warming', 'Data for years 5 - 13', 'Orientation', 'Horizontal', 'FontSIze', 14);
pos = l.Position;
pos(1) = 0.3; % x
pos(2) = 0.02; % y
l.Position = pos;
% =========================================================================

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Measured versus modelled CO2.tiff' -r300 -nocrop -transparent
% =========================================================================

%% Plotting the evolution of the errors over the years

% The unique years
years_perYear = unique(allYears);

% Vectors are created to store the error measures per year
RMSE_control_perYear = NaN(numel(years_perYear),3); % Columns are: noAdapt, enzRig and optDriv
RMSE_heated_perYear = NaN(numel(years_perYear),3);
NS_control_perYear = NaN(numel(years_perYear),3);
NS_heated_perYear = NaN(numel(years_perYear),3);
bias_control_perYear = NaN(numel(years_perYear),3);
bias_heated_perYear = NaN(numel(years_perYear),3);

ii = 1;

% A loop over all years to store the error measures
for ii = 1:numel(years_perYear)
    
    % The current year
    currentYear = years_perYear(ii);
    
    % The data for the control run for the current year
    [r, dummy] = find(ismember(allYears, currentYear));
    meas_control_currentYear = Measured_CO2_flux_control(r);
    mod_noAdapt_control = mod_CO2_measDays_noAdapt_control(r);
    mod_optDriv_control = mod_CO2_measDays_optimumDriv_control(r);
    mod_enzRig_control = mod_CO2_measDays_enzRig_control(r);
    
    % The data for the heated run for the current year
    [r, dummy] = find(ismember(allYears, currentYear));
    meas_heated_currentYear = Measured_CO2_flux_warming(r);
    mod_noAdapt_heated = mod_CO2_measDays_noAdapt_warmed(r);
    mod_optDriv_heated = mod_CO2_measDays_optimumDriv_warmed(r);
    mod_enzRig_heated = mod_CO2_measDays_enzRig_warmed(r);
    
    % The error are calculated and stored
    
    %RMSE
    RMSE_control_perYear(ii,1) = sqrt(mean((meas_control_currentYear - mod_noAdapt_control).^2));
    RMSE_control_perYear(ii,2) = sqrt(mean((meas_control_currentYear - mod_optDriv_control).^2));
    RMSE_control_perYear(ii,3) = sqrt(mean((meas_control_currentYear - mod_enzRig_control).^2));
    
    RMSE_heated_perYear(ii,1) = sqrt(mean((meas_heated_currentYear - mod_noAdapt_heated).^2));
    RMSE_heated_perYear(ii,2) = sqrt(mean((meas_heated_currentYear - mod_optDriv_heated).^2));
    RMSE_heated_perYear(ii,3) = sqrt(mean((meas_heated_currentYear - mod_enzRig_heated).^2));
    
    % Nash–Sutcliffe efficiency
    NS_control_perYear(ii,1) = nash_sutcliffe(meas_control_currentYear, mod_noAdapt_control);
    NS_control_perYear(ii,2) = nash_sutcliffe(meas_control_currentYear, mod_optDriv_control);
    NS_control_perYear(ii,3) = nash_sutcliffe(meas_control_currentYear, mod_enzRig_control);
    
    NS_heated_perYear(ii,1) = nash_sutcliffe(meas_heated_currentYear, mod_noAdapt_heated);
    NS_heated_perYear(ii,2) = nash_sutcliffe(meas_heated_currentYear, mod_optDriv_heated);
    NS_heated_perYear(ii,3) = nash_sutcliffe(meas_heated_currentYear, mod_enzRig_heated);
    
    % Bias
    bias_control_perYear(ii,1) = sum(mod_noAdapt_control - meas_control_currentYear) / numel(meas_control_currentYear);
    bias_control_perYear(ii,2) = sum(mod_optDriv_control - meas_control_currentYear) / numel(meas_control_currentYear);
    bias_control_perYear(ii,3) = sum(mod_enzRig_control - meas_control_currentYear) / numel(meas_control_currentYear);
    
    bias_heated_perYear(ii,1) = sum(mod_noAdapt_heated - meas_heated_currentYear) / numel(meas_heated_currentYear);
    bias_heated_perYear(ii,2) = sum(mod_optDriv_heated - meas_heated_currentYear) / numel(meas_heated_currentYear);
    bias_heated_perYear(ii,3) = sum(mod_enzRig_heated - meas_heated_currentYear) / numel(meas_heated_currentYear);
    
end

% --------
% Plotting
% --------

% ----
% RMSE
% ----

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 25*1.5 15*1.5], 'color', [1 1 1])

subplot(2,1,1)
hold on

scatter(years_perYear, RMSE_control_perYear(:,1), 'MarkerFaceColor', c_noAdapt, 'MarkerEdgeColor', 'none');
scatter(years_perYear, RMSE_control_perYear(:,2), 'MarkerFaceColor', c_optDriv, 'MarkerEdgeColor', 'none');
scatter(years_perYear, RMSE_control_perYear(:,3), 'MarkerFaceColor', c_enzRig, 'MarkerEdgeColor', 'none');

plot(years_perYear, RMSE_control_perYear(:,1), 'color', 'black');

xlabel('Time (years)'); ylabel('RMSE'); title('RMSE - control plots')
set(gca, 'FontSize', 14, 'box', 'on');

subplot(2,1,2)
hold on

p1 = scatter(years_perYear, RMSE_heated_perYear(:,1), 'MarkerFaceColor', c_noAdapt, 'MarkerEdgeColor', 'none');
p2 = scatter(years_perYear, RMSE_heated_perYear(:,2), 'MarkerFaceColor', c_optDriv, 'MarkerEdgeColor', 'none');
p3 = scatter(years_perYear, RMSE_heated_perYear(:,3), 'MarkerFaceColor', c_enzRig, 'MarkerEdgeColor', 'none');

plot(years_perYear, RMSE_heated_perYear(:,1), 'color', 'black');

xlabel('Time (years)'); ylabel('RMSE'); title('RMSE - heated plots')
set(gca, 'FontSize', 14, 'box', 'on');

% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

l = legend([p1, p2, p3], 'No thermal adaptation', 'Optimum driven', 'Enzyme rigidity', 'Orientation', 'Horizontal', 'FontSize', 14);
pos = l.Position;
pos(1) = 0.3; % x
pos(2) = 0.005; % y
l.Position = pos;

% -------------------------
% Nash–Sutcliffe efficiency
% -------------------------

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 25*1.5 15*1.5], 'color', [1 1 1])

subplot(2,1,1)
hold on

scatter(years_perYear, NS_control_perYear(:,1), 'MarkerFaceColor', c_noAdapt, 'MarkerEdgeColor', 'none');
scatter(years_perYear, NS_control_perYear(:,2), 'MarkerFaceColor', c_optDriv, 'MarkerEdgeColor', 'none');
scatter(years_perYear, NS_control_perYear(:,3), 'MarkerFaceColor', c_enzRig, 'MarkerEdgeColor', 'none');

plot(years_perYear, NS_control_perYear(:,1), 'color', 'black');

xlabel('Time (years)'); ylabel('NS'); title('Nash–Sutcliffe efficiency - control plots')
set(gca, 'FontSize', 14, 'box', 'on');

subplot(2,1,2)
hold on

p1 = scatter(years_perYear, NS_heated_perYear(:,1), 'MarkerFaceColor', c_noAdapt, 'MarkerEdgeColor', 'none');
p2 = scatter(years_perYear, NS_heated_perYear(:,2), 'MarkerFaceColor', c_optDriv, 'MarkerEdgeColor', 'none');
p3 = scatter(years_perYear, NS_heated_perYear(:,3), 'MarkerFaceColor', c_enzRig, 'MarkerEdgeColor', 'none');

plot(years_perYear, NS_heated_perYear(:,1), 'color', 'black');

xlabel('Time (years)'); ylabel('NS'); title('Nash–Sutcliffe efficiency - heated plots')
set(gca, 'FontSize', 14, 'box', 'on');

% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

l = legend([p1, p2, p3], 'No thermal adaptation', 'Optimum driven', 'Enzyme rigidity', 'Orientation', 'Horizontal', 'FontSize', 14);
pos = l.Position;
pos(1) = 0.3; % x
pos(2) = 0.005; % y
l.Position = pos;

% ----
% Bias
% ----

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 25*1.5 15*1.5], 'color', [1 1 1])

subplot(2,1,1)
hold on

scatter(years_perYear, bias_control_perYear(:,1), 'MarkerFaceColor', c_noAdapt, 'MarkerEdgeColor', 'none');
scatter(years_perYear, bias_control_perYear(:,2), 'MarkerFaceColor', c_optDriv, 'MarkerEdgeColor', 'none');
scatter(years_perYear, bias_control_perYear(:,3), 'MarkerFaceColor', c_enzRig, 'MarkerEdgeColor', 'none');

plot(years_perYear, bias_control_perYear(:,1), 'color', 'black');

xlabel('Time (years)'); ylabel('NS'); title('Nash–Sutcliffe efficiency - control plots')
set(gca, 'FontSize', 14, 'box', 'on');

subplot(2,1,2)
hold on

p1 = scatter(years_perYear, bias_heated_perYear(:,1), 'MarkerFaceColor', c_noAdapt, 'MarkerEdgeColor', 'none');
p2 = scatter(years_perYear, bias_heated_perYear(:,2), 'MarkerFaceColor', c_optDriv, 'MarkerEdgeColor', 'none');
p3 = scatter(years_perYear, bias_heated_perYear(:,3), 'MarkerFaceColor', c_enzRig, 'MarkerEdgeColor', 'none');

plot(years_perYear, bias_heated_perYear(:,1), 'color', 'black');

xlabel('Time (years)'); ylabel('NS'); title('Bias - heated plots')
set(gca, 'FontSize', 14, 'box', 'on');

% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

l = legend([p1, p2, p3], 'No thermal adaptation', 'Optimum driven', 'Enzyme rigidity', 'Orientation', 'Horizontal', 'FontSize', 14);
pos = l.Position;
pos(1) = 0.3; % x
pos(2) = 0.005; % y
l.Position = pos;

%% Scatterplot of measured versus modelled CO2 per timestep

color_posRes = [228,26,28]./255;
color_negRes = [55,126,184]./255;
color_regr = 'k';%[33,102,172]./255;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 25*1.5 10*1.5], 'color', [1 1 1])

% --------------------------------------
% Forest floor - no adaptation - control
% --------------------------------------

subplot(2,3,1)
hold on

xlim([0 6])
ylim([0 6])

scatter(noAdaptData.CO2flux_negRes_meas_control, noAdaptData.CO2flux_negres_mod_control, 'MarkerFaceColor', color_negRes, 'MarkerEdgeColor', 'none');
scatter(noAdaptData.CO2flux_posRes_meas_control, noAdaptData.CO2flux_posres_mod_control, 'MarkerFaceColor', color_posRes, 'MarkerEdgeColor', 'none'); 

% RMSE
tmpMeas = [noAdaptData.CO2flux_negRes_meas_control; noAdaptData.CO2flux_posRes_meas_control];
tmpMod = [noAdaptData.CO2flux_negres_mod_control; noAdaptData.CO2flux_posres_mod_control];
% RRMSE = sqrt(mean(((tmpMeas - tmpMod)./tmpMeas).^2));
RRMSE = sqrt(mean((tmpMeas - tmpMod).^2));
RRMSE = round(RRMSE*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - .75, ['RRMSE = ' num2str(RRMSE)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

xlabel('Measured CO_{2} (g C per timestep)')
ylabel('Modelled CO_{2} (g C per timestep)')
title('Control treatment')

% Bias
bias = mean(tmpMod - tmpMeas);
bias = round(bias*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - 1.3, ['Bias = ' num2str(bias)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

% The 1:1 line
hline1 = refline(1,0);
hline1.Color = [.5 .5 .5];
hline1.LineStyle = '--';

% The regression line
P = polyfit(tmpMeas,tmpMod,1);
x0 = min(tmpMeas) ; 
x1 = max(tmpMeas) ;
xi = linspace(x0,x1) ;
yi = P(1)*xi+P(2);
hold on
plot(xi, yi, '--', 'color', color_regr, 'linewidth', 2) ;

slope = floor(P(1)*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - 1.8, ['Slope = ' num2str(slope)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

xlabel('Measured CO_{2} (g C per timestep)')
% ylabel('Modelled CO_{2} (g C per timestep)')
title('No adaptation - control')

% --------------------------------------
% Forest floor - Optimum driven - control
% --------------------------------------

subplot(2,3,2)
hold on

xlim([0 6])
ylim([0 6])

scatter(optimumDrivenData.CO2flux_negRes_meas_control, optimumDrivenData.CO2flux_negres_mod_control, 'MarkerFaceColor', color_negRes, 'MarkerEdgeColor', 'none');
scatter(optimumDrivenData.CO2flux_posRes_meas_control, optimumDrivenData.CO2flux_posres_mod_control, 'MarkerFaceColor', color_posRes, 'MarkerEdgeColor', 'none'); 

% RMSE
tmpMeas = [optimumDrivenData.CO2flux_negRes_meas_control; optimumDrivenData.CO2flux_posRes_meas_control];
tmpMod = [optimumDrivenData.CO2flux_negres_mod_control; optimumDrivenData.CO2flux_posres_mod_control];
% RRMSE = sqrt(mean(((tmpMeas - tmpMod)./tmpMeas).^2));
RRMSE = sqrt(mean((tmpMeas - tmpMod).^2));
RRMSE = round(RRMSE*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - .75, ['RRMSE = ' num2str(RRMSE)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

xlabel('Measured CO_{2} (g C per timestep)')
ylabel('Modelled CO_{2} (g C per timestep)')
title('Control treatment')

% Bias
bias = mean(tmpMod - tmpMeas);
bias = round(bias*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - 1.3, ['Bias = ' num2str(bias)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

% The 1:1 line
hline1 = refline(1,0);
hline1.Color = [.5 .5 .5];
hline1.LineStyle = '--';

% The regression line
P = polyfit(tmpMeas,tmpMod,1);
x0 = min(tmpMeas) ; 
x1 = max(tmpMeas) ;
xi = linspace(x0,x1) ;
yi = P(1)*xi+P(2);
hold on
plot(xi, yi, '--', 'color', color_regr, 'linewidth', 2) ;

slope = floor(P(1)*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - 1.8, ['Slope = ' num2str(slope)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

xlabel('Measured CO_{2} (g C per timestep)')
% ylabel('Modelled CO_{2} (g C per timestep)')
title('Optimum driven - control')

% --------------------------------------
% Forest floor - Enzyme rigidity - control
% --------------------------------------

subplot(2,3,3)
hold on

xlim([0 6])
ylim([0 6])

scatter(enzRigData.CO2flux_negRes_meas_control, enzRigData.CO2flux_negres_mod_control, 'MarkerFaceColor', color_negRes, 'MarkerEdgeColor', 'none');
scatter(enzRigData.CO2flux_posRes_meas_control, enzRigData.CO2flux_posres_mod_control, 'MarkerFaceColor', color_posRes, 'MarkerEdgeColor', 'none'); 

% RMSE
tmpMeas = [enzRigData.CO2flux_negRes_meas_control; enzRigData.CO2flux_posRes_meas_control];
tmpMod = [enzRigData.CO2flux_negres_mod_control; enzRigData.CO2flux_posres_mod_control];
% RRMSE = sqrt(mean(((tmpMeas - tmpMod)./tmpMeas).^2));
RRMSE = sqrt(mean((tmpMeas - tmpMod).^2));
RRMSE = round(RRMSE*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - .75, ['RRMSE = ' num2str(RRMSE)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

xlabel('Measured CO_{2} (g C per timestep)')
ylabel('Modelled CO_{2} (g C per timestep)')
% title('Control treatment')

% Bias
bias = mean(tmpMod - tmpMeas);
bias = round(bias*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - 1.3, ['Bias = ' num2str(bias)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

% The 1:1 line
hline1 = refline(1,0);
hline1.Color = [.5 .5 .5];
hline1.LineStyle = '--';

% The regression line
P = polyfit(tmpMeas,tmpMod,1);
x0 = min(tmpMeas) ; 
x1 = max(tmpMeas) ;
xi = linspace(x0,x1) ;
yi = P(1)*xi+P(2);
hold on
plot(xi, yi, '--', 'color', color_regr, 'linewidth', 2) ;

slope = floor(P(1)*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - 1.8, ['Slope = ' num2str(slope)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

xlabel('Measured CO_{2} (g C per timestep)')
% ylabel('Modelled CO_{2} (g C per timestep)')
title('Enzyme rigidity - control')

% --------------------------------------
% Forest floor - no adaptation - warming
% --------------------------------------

subplot(2,3,4)
hold on

xlim([0 6])

ylim([0 6])

scatter(noAdaptData.CO2flux_negRes_meas_warming, noAdaptData.CO2flux_negres_mod_warming, 'MarkerFaceColor', color_negRes, 'MarkerEdgeColor', 'none');
scatter(noAdaptData.CO2flux_posRes_meas_warming, noAdaptData.CO2flux_posres_mod_warming, 'MarkerFaceColor', color_posRes, 'MarkerEdgeColor', 'none'); 

% RMSE
tmpMeas = [noAdaptData.CO2flux_negRes_meas_warming; noAdaptData.CO2flux_posRes_meas_warming];
tmpMod = [noAdaptData.CO2flux_negres_mod_warming; noAdaptData.CO2flux_posres_mod_warming];
% RRMSE = sqrt(mean(((tmpMeas - tmpMod)./tmpMeas).^2));
RRMSE = sqrt(mean((tmpMeas - tmpMod).^2));
RRMSE = round(RRMSE*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - .75, ['RRMSE = ' num2str(RRMSE)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

xlabel('Measured CO_{2} (g C per timestep)')
ylabel('Modelled CO_{2} (g C per timestep)')
title('Control treatment')

% Bias
bias = mean(tmpMod - tmpMeas);
bias = round(bias*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - 1.3, ['Bias = ' num2str(bias)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

% The 1:1 line
hline1 = refline(1,0);
hline1.Color = [.5 .5 .5];
hline1.LineStyle = '--';

% The regression line
P = polyfit(tmpMeas,tmpMod,1);
x0 = min(tmpMeas) ; 
x1 = max(tmpMeas) ;
xi = linspace(x0,x1) ;
yi = P(1)*xi+P(2);
hold on
plot(xi, yi, '--', 'color', color_regr, 'linewidth', 2) ;

slope = floor(P(1)*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - 1.8, ['Slope = ' num2str(slope)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

xlabel('Measured CO_{2} (g C per timestep)')
% ylabel('Modelled CO_{2} (g C per timestep)')
title('No adaptation - warming')

% --------------------------------------
% Forest floor - Optimum driven - warming
% --------------------------------------

subplot(2,3,5)
hold on

xlim([0 6])
ylim([0 6])

scatter(optimumDrivenData.CO2flux_negRes_meas_warming, optimumDrivenData.CO2flux_negres_mod_warming, 'MarkerFaceColor', color_negRes, 'MarkerEdgeColor', 'none');
scatter(optimumDrivenData.CO2flux_posRes_meas_warming, optimumDrivenData.CO2flux_posres_mod_warming, 'MarkerFaceColor', color_posRes, 'MarkerEdgeColor', 'none'); 

% RMSE
tmpMeas = [optimumDrivenData.CO2flux_negRes_meas_warming; optimumDrivenData.CO2flux_posRes_meas_warming];
tmpMod = [optimumDrivenData.CO2flux_negres_mod_warming; optimumDrivenData.CO2flux_posres_mod_warming];
% RRMSE = sqrt(mean(((tmpMeas - tmpMod)./tmpMeas).^2));
RRMSE = sqrt(mean((tmpMeas - tmpMod).^2));
RRMSE = round(RRMSE*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - .75, ['RRMSE = ' num2str(RRMSE)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

xlabel('Measured CO_{2} (g C per timestep)')
ylabel('Modelled CO_{2} (g C per timestep)')
title('Control treatment')

% Bias
bias = mean(tmpMod - tmpMeas);
bias = round(bias*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - 1.3, ['Bias = ' num2str(bias)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

% The 1:1 line
hline1 = refline(1,0);
hline1.Color = [.5 .5 .5];
hline1.LineStyle = '--';

% The regression line
P = polyfit(tmpMeas,tmpMod,1);
x0 = min(tmpMeas) ; 
x1 = max(tmpMeas) ;
xi = linspace(x0,x1) ;
yi = P(1)*xi+P(2);
hold on
plot(xi, yi, '--', 'color', color_regr, 'linewidth', 2) ;

slope = floor(P(1)*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - 1.8, ['Slope = ' num2str(slope)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

xlabel('Measured CO_{2} (g C per timestep)')
% ylabel('Modelled CO_{2} (g C per timestep)')
title('Optimum driven - warming')

% --------------------------------------
% Forest floor - Enzyme rigidity - warming
% --------------------------------------

subplot(2,3,6)
hold on

xlim([0 6])
ylim([0 6])

scatter(enzRigData.CO2flux_negRes_meas_warming, enzRigData.CO2flux_negres_mod_warming, 'MarkerFaceColor', color_negRes, 'MarkerEdgeColor', 'none');
scatter(enzRigData.CO2flux_posRes_meas_warming, enzRigData.CO2flux_posres_mod_warming, 'MarkerFaceColor', color_posRes, 'MarkerEdgeColor', 'none'); 

% RMSE
tmpMeas = [enzRigData.CO2flux_negRes_meas_warming; enzRigData.CO2flux_posRes_meas_warming];
tmpMod = [enzRigData.CO2flux_negres_mod_warming; enzRigData.CO2flux_posres_mod_warming];
% RRMSE = sqrt(mean(((tmpMeas - tmpMod)./tmpMeas).^2));
RRMSE = sqrt(mean((tmpMeas - tmpMod).^2));
RRMSE = round(RRMSE*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - .75, ['RRMSE = ' num2str(RRMSE)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

xlabel('Measured CO_{2} (g C per timestep)')
ylabel('Modelled CO_{2} (g C per timestep)')
% title('Control treatment')

% Bias
bias = mean(tmpMod - tmpMeas);
bias = round(bias*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - 1.3, ['Bias = ' num2str(bias)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

% The 1:1 line
hline1 = refline(1,0);
hline1.Color = [.5 .5 .5];
hline1.LineStyle = '--';

% The regression line
P = polyfit(tmpMeas,tmpMod,1);
x0 = min(tmpMeas) ; 
x1 = max(tmpMeas) ;
xi = linspace(x0,x1) ;
yi = P(1)*xi+P(2);
hold on
plot(xi, yi, '--', 'color', color_regr, 'linewidth', 2) ;

slope = floor(P(1)*100)/100;
T = text(min(get(gca, 'xlim')) + .5, max(get(gca, 'ylim')) - 1.8, ['Slope = ' num2str(slope)]);
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

xlabel('Measured CO_{2} (g C per timestep)')
% ylabel('Modelled CO_{2} (g C per timestep)')
title('Enzyme rigidity - warming')

%% Plotting the amount of available substrate per MBC

% ------------------------------------
% The data is retrieved: no adaptation
% ------------------------------------

% Monomers
noAdapt_monomers_control_litter = noAdaptData.Cpools_litter_control(:,6);
noAdapt_monomers_warming_litter = noAdaptData.Cpools_litter_warming(:,6);

noAdapt_monomers_control_soil = noAdaptData.Cpools_litter_control(:,4);
noAdapt_monomers_warming_soil = noAdaptData.Cpools_litter_warming(:,4);

% Microbial biomass C
noAdapt_mbc_control_litter = sum(noAdaptData.Cpools_litter_control(:,[1 2]),2);
noAdapt_mbc_warming_litter = sum(noAdaptData.Cpools_litter_warming(:,[1 2]),2);

noAdapt_mbc_control_soil = noAdaptData.Cpools_litter_control(:,1);
noAdapt_mbc_warming_soil = noAdaptData.Cpools_litter_warming(:,1);

% % The ratio of substrate to MBC is calculated
% noAdapt_substr_to_mbc_control = (noAdapt_monomers_control_litter + noAdapt_monomers_control_soil)./ (noAdapt_mbc_control_litter + noAdapt_mbc_control_soil);
% noAdapt_substr_to_mbc_warming = (noAdapt_monomers_warming_litter + noAdapt_monomers_warming_soil)./ (noAdapt_mbc_warming_litter + noAdapt_mbc_warming_soil);

% The ratio of MBC to substrate is calculated
noAdapt_mbc_to_substr_control_litter = noAdapt_mbc_control_litter ./ noAdapt_monomers_control_litter;
noAdapt_mbc_to_substr_control_soil = noAdapt_mbc_control_soil ./ noAdapt_monomers_control_soil;
noAdapt_mbc_to_substr_warming_litter = noAdapt_mbc_warming_litter ./ noAdapt_monomers_warming_litter;
noAdapt_mbc_to_substr_warming_soil = noAdapt_mbc_warming_soil ./ noAdapt_monomers_warming_soil;

% -------------------------------------
% The data is retrieved: optimum driven
% -------------------------------------

% Monomers
optDriv_monomers_control_litter = optimumDrivenData.Cpools_litter_control(:,6);
optDriv_monomers_warming_litter = optimumDrivenData.Cpools_litter_warming(:,6);

optDriv_monomers_control_soil = optimumDrivenData.Cpools_litter_control(:,4);
optDriv_monomers_warming_soil = optimumDrivenData.Cpools_litter_warming(:,4);

% Microbial biomass C
optDriv_mbc_control_litter = sum(optimumDrivenData.Cpools_litter_control(:,[1 2]),2);
optDriv_mbc_warming_litter = sum(optimumDrivenData.Cpools_litter_warming(:,[1 2]),2);

optDriv_mbc_control_soil = optimumDrivenData.Cpools_litter_control(:,1);
optDriv_mbc_warming_soil = optimumDrivenData.Cpools_litter_warming(:,1);

% % The ratio of substrate to MBC is calculated
% optDriv_substr_to_mbc_control = (optDriv_monomers_control_litter + optDriv_monomers_control_soil)./ (optDriv_mbc_control_litter + optDriv_mbc_control_soil);
% optDriv_substr_to_mbc_warming = (optDriv_monomers_warming_litter + optDriv_monomers_warming_soil)./ (optDriv_mbc_warming_litter + optDriv_mbc_warming_soil);

% The ratio of MBC to substrate is calculated
optDriv_mbc_to_substr_control_litter = optDriv_mbc_control_litter ./ optDriv_monomers_control_litter;
optDriv_mbc_to_substr_control_soil = optDriv_mbc_control_soil ./ optDriv_monomers_control_soil;
optDriv_mbc_to_substr_warming_litter = optDriv_mbc_warming_litter ./ optDriv_monomers_warming_litter;
optDriv_mbc_to_substr_warming_soil = optDriv_mbc_warming_soil ./ optDriv_monomers_warming_soil;

% --------------------------------------
% The data is retrieved: enzyme rigidity
% --------------------------------------

% Monomers
enzRig_monomers_control_litter = enzRigData.Cpools_litter_control(:,6);
enzRig_monomers_warming_litter = enzRigData.Cpools_litter_warming(:,6);

enzRig_monomers_control_soil = enzRigData.Cpools_litter_control(:,4);
enzRig_monomers_warming_soil = enzRigData.Cpools_litter_warming(:,4);

% Microbial biomass C
enzRig_mbc_control_litter = sum(enzRigData.Cpools_litter_control(:,[1 2]),2);
enzRig_mbc_warming_litter = sum(enzRigData.Cpools_litter_warming(:,[1 2]),2);

enzRig_mbc_control_soil = enzRigData.Cpools_litter_control(:,1);
enzRig_mbc_warming_soil = enzRigData.Cpools_litter_warming(:,1);

% % The ratio of substrate to MBC is calculated
% enzRig_substr_to_mbc_control = (enzRig_monomers_control_litter + enzRig_monomers_control_soil)./ (enzRig_mbc_control_litter + enzRig_mbc_control_soil);
% enzRig_substr_to_mbc_warming = (enzRig_monomers_warming_litter + enzRig_monomers_warming_soil)./ (enzRig_mbc_warming_litter + enzRig_mbc_warming_soil);

% The ratio of MBC to substrate is calculated
enzRig_mbc_to_substr_control_litter = enzRig_mbc_control_litter ./ enzRig_monomers_control_litter;
enzRig_mbc_to_substr_control_soil = enzRig_mbc_control_soil ./ enzRig_monomers_control_soil;
enzRig_mbc_to_substr_warming_litter = enzRig_mbc_warming_litter ./ enzRig_monomers_warming_litter;
enzRig_mbc_to_substr_warming_soil = enzRig_mbc_warming_soil ./ enzRig_monomers_warming_soil;

% ------------------------------------------------------------
% The data for enzyme rigidity and optimum driven are combined
% ------------------------------------------------------------

% Monomers
thermalAdapt_monomers_control_litter = mean([optDriv_monomers_control_litter enzRig_monomers_control_litter],2);
thermalAdapt_monomers_control_soil = mean([optDriv_monomers_control_soil enzRig_monomers_control_soil],2);

thermalAdapt_monomers_warming_litter = mean([optDriv_monomers_warming_litter enzRig_monomers_warming_litter],2);
thermalAdapt_monomers_warming_soil = mean([optDriv_monomers_warming_soil enzRig_monomers_warming_soil],2);

% Microbial biomass C
thermalAdapt_mbc_control_litter = mean([optDriv_mbc_control_litter enzRig_mbc_control_litter],2);
thermalAdapt_mbc_control_soil = mean([optDriv_mbc_control_soil enzRig_mbc_control_soil],2);

thermalAdapt_mbc_warming_litter = mean([optDriv_mbc_warming_litter enzRig_mbc_warming_litter],2);
thermalAdapt_mbc_warming_soil = mean([optDriv_mbc_warming_soil enzRig_mbc_warming_soil],2);

% ---------------------------------------------------
% The average annual amount of microbes is calculated
% ---------------------------------------------------

% The individual yearsare identified
uniqueYears = unique(dates_all.Year);

% Empty array to store the averages in
avgMic_noAdapt_litter_control = NaN(numel(uniqueYears),1);
avgMic_noAdapt_litter_warming = NaN(numel(uniqueYears),1);
avgMic_noAdapt_soil_control = NaN(numel(uniqueYears),1);
avgMic_noAdapt_soil_warming = NaN(numel(uniqueYears),1);

avgMic_optDriv_litter_control = NaN(numel(uniqueYears),1);
avgMic_optDriv_litter_warming = NaN(numel(uniqueYears),1);
avgMic_optDriv_soil_control = NaN(numel(uniqueYears),1);
avgMic_optDriv_soil_warming = NaN(numel(uniqueYears),1);

avgMic_enzRig_litter_control = NaN(numel(uniqueYears),1);
avgMic_enzRig_litter_warming = NaN(numel(uniqueYears),1);
avgMic_enzRig_soil_control = NaN(numel(uniqueYears),1);
avgMic_enzRig_soil_warming = NaN(numel(uniqueYears),1);

% A loop to calculate the annual averages
for ii = 1:numel(uniqueYears)
    
    yr = uniqueYears(ii);
    rowNum = find(dates_all.Year == yr);
    
    % The averages are calculated
    avgMic_noAdapt_litter_control(ii,1) = mean(noAdapt_mbc_control_litter(rowNum));
    avgMic_noAdapt_litter_warming(ii,1) = mean(noAdapt_mbc_warming_litter(rowNum));
    avgMic_noAdapt_soil_control(ii,1) = mean(noAdapt_mbc_control_soil(rowNum));
    avgMic_noAdapt_soil_warming(ii,1) = mean(noAdapt_mbc_warming_soil(rowNum));
    
    avgMic_optDriv_litter_control(ii,1) = mean(optDriv_mbc_control_litter(rowNum));
    avgMic_optDriv_litter_warming(ii,1) = mean(optDriv_mbc_warming_litter(rowNum));
    avgMic_optDriv_soil_control(ii,1) = mean(optDriv_mbc_control_soil(rowNum));
    avgMic_optDriv_soil_warming(ii,1) = mean(optDriv_mbc_warming_soil(rowNum));
    
    avgMic_enzRig_litter_control(ii,1) = mean(enzRig_mbc_control_litter(rowNum));
    avgMic_enzRig_litter_warming(ii,1) = mean(enzRig_mbc_warming_litter(rowNum));
    avgMic_enzRig_soil_control(ii,1) = mean(enzRig_mbc_control_soil(rowNum));
    avgMic_enzRig_soil_warming(ii,1) = mean(enzRig_mbc_warming_soil(rowNum));
    
end

% The average of both thermal adaptation scenarios
avgMic_thermalAdapt_litter_control = mean([avgMic_optDriv_litter_control avgMic_enzRig_litter_control],2);
avgMic_thermalAdapt_litter_warming = mean([avgMic_optDriv_litter_warming avgMic_enzRig_litter_warming],2);
avgMic_thermalAdapt_soil_control = mean([avgMic_optDriv_soil_control avgMic_enzRig_soil_control],2);
avgMic_thermalAdapt_soil_warming = mean([avgMic_optDriv_soil_warming avgMic_enzRig_soil_warming],2);

% ----------------------------------------------------
% The average annual amount of substrate is calculated
% ----------------------------------------------------

% Empty array to store the averages in
avgSub_noAdapt_litter_control = NaN(numel(uniqueYears),1);
avgSub_noAdapt_litter_warming = NaN(numel(uniqueYears),1);
avgSub_noAdapt_soil_control = NaN(numel(uniqueYears),1);
avgSub_noAdapt_soil_warming = NaN(numel(uniqueYears),1);

avgSub_optDriv_litter_control = NaN(numel(uniqueYears),1);
avgSub_optDriv_litter_warming = NaN(numel(uniqueYears),1);
avgSub_optDriv_soil_control = NaN(numel(uniqueYears),1);
avgSub_optDriv_soil_warming = NaN(numel(uniqueYears),1);

avgSub_enzRig_litter_control = NaN(numel(uniqueYears),1);
avgSub_enzRig_litter_warming = NaN(numel(uniqueYears),1);
avgSub_enzRig_soil_control = NaN(numel(uniqueYears),1);
avgSub_enzRig_soil_warming = NaN(numel(uniqueYears),1);

% A loop to calculate the annual averages
for ii = 1:numel(uniqueYears)
    
    yr = uniqueYears(ii);
    rowNum = find(dates_all.Year == yr);
    
    % The averages are calculated
    avgSub_noAdapt_litter_control(ii,1) = mean(noAdapt_monomers_control_litter(rowNum));
    avgSub_noAdapt_litter_warming(ii,1) = mean(noAdapt_monomers_warming_litter(rowNum));
    avgSub_noAdapt_soil_control(ii,1) = mean(noAdapt_monomers_control_soil(rowNum));
    avgSub_noAdapt_soil_warming(ii,1) = mean(noAdapt_monomers_warming_soil(rowNum));
    
    avgSub_optDriv_litter_control(ii,1) = mean(optDriv_monomers_control_litter(rowNum));
    avgSub_optDriv_litter_warming(ii,1) = mean(optDriv_monomers_warming_litter(rowNum));
    avgSub_optDriv_soil_control(ii,1) = mean(optDriv_monomers_control_soil(rowNum));
    avgSub_optDriv_soil_warming(ii,1) = mean(optDriv_monomers_warming_soil(rowNum));
    
    avgSub_enzRig_litter_control(ii,1) = mean(enzRig_monomers_control_litter(rowNum));
    avgSub_enzRig_litter_warming(ii,1) = mean(enzRig_monomers_warming_litter(rowNum));
    avgSub_enzRig_soil_control(ii,1) = mean(enzRig_monomers_control_soil(rowNum));
    avgSub_enzRig_soil_warming(ii,1) = mean(enzRig_monomers_warming_soil(rowNum));
    
end

% The average of both thermal adaptation scenarios
avgSub_thermalAdapt_litter_control = mean([avgSub_optDriv_litter_control avgSub_enzRig_litter_control],2);
avgSub_thermalAdapt_litter_warming = mean([avgSub_optDriv_litter_warming avgSub_enzRig_litter_warming],2);
avgSub_thermalAdapt_soil_control = mean([avgSub_optDriv_soil_control avgSub_enzRig_soil_control],2);
avgSub_thermalAdapt_soil_warming = mean([avgSub_optDriv_soil_warming avgSub_enzRig_soil_warming],2);


% ---------------------------------------------------------------
% The average annual ratio of microbes to substrate is calculated
% ---------------------------------------------------------------

% Empty array to store the averages in
avgRatio_noAdapt_litter_control = NaN(numel(uniqueYears),1);
avgRatio_noAdapt_litter_warming = NaN(numel(uniqueYears),1);
avgRatio_noAdapt_soil_control = NaN(numel(uniqueYears),1);
avgRatio_noAdapt_soil_warming = NaN(numel(uniqueYears),1);

avgRatio_optDriv_litter_control = NaN(numel(uniqueYears),1);
avgRatio_optDriv_litter_warming = NaN(numel(uniqueYears),1);
avgRatio_optDriv_soil_control = NaN(numel(uniqueYears),1);
avgRatio_optDriv_soil_warming = NaN(numel(uniqueYears),1);

avgRatio_enzRig_litter_control = NaN(numel(uniqueYears),1);
avgRatio_enzRig_litter_warming = NaN(numel(uniqueYears),1);
avgRatio_enzRig_soil_control = NaN(numel(uniqueYears),1);
avgRatio_enzRig_soil_warming = NaN(numel(uniqueYears),1);

% A loop to calculate the annual averages
for ii = 1:numel(uniqueYears)
    
    yr = uniqueYears(ii);
    rowNum = find(dates_all.Year == yr);
    
    % The averages are calculated
    avgRatio_noAdapt_litter_control(ii,1) = mean(noAdapt_mbc_to_substr_control_litter(rowNum));
    avgRatio_noAdapt_litter_warming(ii,1) = mean(noAdapt_mbc_to_substr_warming_litter(rowNum));
    avgRatio_noAdapt_soil_control(ii,1) = mean(noAdapt_mbc_to_substr_control_soil(rowNum));
    avgRatio_noAdapt_soil_warming(ii,1) = mean(noAdapt_mbc_to_substr_warming_soil(rowNum));
    
    avgRatio_optDriv_litter_control(ii,1) = mean(optDriv_mbc_to_substr_control_litter(rowNum));
    avgRatio_optDriv_litter_warming(ii,1) = mean(optDriv_mbc_to_substr_warming_litter(rowNum));
    avgRatio_optDriv_soil_control(ii,1) = mean(optDriv_mbc_to_substr_control_soil(rowNum));
    avgRatio_optDriv_soil_warming(ii,1) = mean(optDriv_mbc_to_substr_warming_soil(rowNum));
    
    avgRatio_enzRig_litter_control(ii,1) = mean(enzRig_mbc_to_substr_control_litter(rowNum));
    avgRatio_enzRig_litter_warming(ii,1) = mean(enzRig_mbc_to_substr_warming_litter(rowNum));
    avgRatio_enzRig_soil_control(ii,1) = mean(enzRig_mbc_to_substr_control_soil(rowNum));
    avgRatio_enzRig_soil_warming(ii,1) = mean(enzRig_mbc_to_substr_warming_soil(rowNum));
    
end

% The average of both thermal adaptation scenarios
avgRatio_thermalAdapt_litter_control = mean([avgRatio_optDriv_litter_control avgRatio_enzRig_litter_control],2);
avgRatio_thermalAdapt_litter_warming = mean([avgRatio_optDriv_litter_warming avgRatio_enzRig_litter_warming],2);
avgRatio_thermalAdapt_soil_control = mean([avgRatio_optDriv_soil_control avgRatio_enzRig_soil_control],2);
avgRatio_thermalAdapt_soil_warming = mean([avgRatio_optDriv_soil_warming avgRatio_enzRig_soil_warming],2);

% Dates for plotting the averages
dates_averages = append(string(uniqueYears), '/07/01');
dates_averages = datetime(dates_averages, 'Format', 'yyyy/MM/dd');

% ------------------------------------------------
% The relative decrease in MBC in the heated plots
% ------------------------------------------------

date1 = datetime('01/01/2000', 'format', 'dd/MM/yyyy');
date2 = datetime('31/12/2002', 'format', 'dd/MM/yyyy');

date3 = datetime('01/01/2014', 'format', 'dd/MM/yyyy');
date4 = datetime('31/12/2016', 'format', 'dd/MM/yyyy');

r1 = find(dates_all == date1);
r2 = find(dates_all == date2);

r3 = find(dates_all == date3);
r4 = find(dates_all == date4);

% % No adaptation
% C1 = mean(noAdapt_mbc_warming_litter(r1:r2));
% C2 = mean(noAdapt_mbc_warming_litter(r3:r4));
% dC_ff_noAdapt = (C2 - C1)/C1;
% 
% C3 = mean(noAdapt_mbc_warming_soil(r1:r2));
% C4 = mean(noAdapt_mbc_warming_soil(r3:r4));
% dC_soil_noAdapt = (C4 - C3)/C3;
% 
% C5 = mean(noAdapt_mbc_warming_litter(r1:r2) + noAdapt_mbc_warming_soil(r1:r2));
% C6 = mean(noAdapt_mbc_warming_litter(r3:r4) + noAdapt_mbc_warming_soil(r3:r4));
% dC_tot_noAdapt = (C6 - C5)/C5;
% 
% % Optimum driven
% C1 = mean(optDriv_mbc_warming_litter(r1:r2));
% C2 = mean(optDriv_mbc_warming_litter(r3:r4));
% dC_ff_optDriv = (C2 - C1)/C1;
% 
% C3 = mean(optDriv_mbc_warming_soil(r1:r2));
% C4 = mean(optDriv_mbc_warming_soil(r3:r4));
% dC_soil_optDriv = (C4 - C3)/C3;
% 
% C5 = mean(optDriv_mbc_warming_litter(r1:r2) + optDriv_mbc_warming_soil(r1:r2));
% C6 = mean(optDriv_mbc_warming_litter(r3:r4) + optDriv_mbc_warming_soil(r3:r4));
% dC_tot_optDriv = (C6 - C5)/C5;
% 
% % Enzyme rigidity
% C1 = mean(enzRig_mbc_warming_litter(r1:r2));
% C2 = mean(enzRig_mbc_warming_litter(r3:r4));
% dC_ff_enzRig = (C2 - C1)/C1;
% 
% C3 = mean(enzRig_mbc_warming_soil(r1:r2));
% C4 = mean(enzRig_mbc_warming_soil(r3:r4));
% dC_soil_enzRig = (C4 - C3)/C3;
% 
% C5 = mean(enzRig_mbc_warming_litter(r1:r2) + enzRig_mbc_warming_soil(r1:r2));
% C6 = mean(enzRig_mbc_warming_litter(r3:r4) + enzRig_mbc_warming_soil(r3:r4));
% dC_tot_enzRig = (C6 - C5)/C5;
% 
% dC_tot_thermAdapt = mean([dC_tot_optDriv dC_tot_enzRig]);

% No adaptation
C1 = mean(noAdapt_mbc_control_litter(r3:r4));
C2 = mean(noAdapt_mbc_warming_litter(r3:r4));
dC_ff_noAdapt = (C2 - C1)/C1;

C3 = mean(noAdapt_mbc_control_soil(r3:r4));
C4 = mean(noAdapt_mbc_warming_soil(r3:r4));
dC_soil_noAdapt = (C4 - C3)/C3;

C5 = mean(noAdapt_mbc_control_litter(r3:r4) + noAdapt_mbc_control_soil(r3:r4));
C6 = mean(noAdapt_mbc_warming_litter(r3:r4) + noAdapt_mbc_warming_soil(r3:r4));
dC_tot_noAdapt = (C6 - C5)/C5;

% Optimum driven
C1 = mean(optDriv_mbc_control_litter(r3:r4));
C2 = mean(optDriv_mbc_warming_litter(r3:r4));
dC_ff_optDriv = (C2 - C1)/C1;

C3 = mean(optDriv_mbc_control_soil(r3:r4));
C4 = mean(optDriv_mbc_warming_soil(r3:r4));
dC_soil_optDriv = (C4 - C3)/C3;

C5 = mean(optDriv_mbc_control_litter(r3:r4) + optDriv_mbc_control_soil(r3:r4));
C6 = mean(optDriv_mbc_warming_litter(r3:r4) + optDriv_mbc_warming_soil(r3:r4));
dC_tot_optDriv = (C6 - C5)/C5;

% Enzyme rigidity
C1 = mean(enzRig_mbc_control_litter(r3:r4));
C2 = mean(enzRig_mbc_warming_litter(r3:r4));
dC_ff_enzRig = (C2 - C1)/C1;

C3 = mean(enzRig_mbc_control_soil(r3:r4));
C4 = mean(enzRig_mbc_warming_soil(r3:r4));
dC_soil_enzRig = (C4 - C3)/C3;

C5 = mean(enzRig_mbc_control_litter(r3:r4) + enzRig_mbc_control_soil(r3:r4));
C6 = mean(enzRig_mbc_warming_litter(r3:r4) + enzRig_mbc_warming_soil(r3:r4));
dC_tot_enzRig = (C6 - C5)/C5;

dC_tot_thermAdapt = mean([dC_tot_optDriv dC_tot_enzRig]);

% -----------------------------------------------------------------
% The Change in MBC during the last 3 simulated years is calculated
% -----------------------------------------------------------------

avgMic_lastThreeYears_noAdapt_litter_control = mean(avgMic_noAdapt_litter_control(end-2:end));
avgMic_lastThreeYears_noAdapt_litter_warming = mean(avgMic_noAdapt_litter_warming(end-2:end));
avgMic_lastThreeYears_noAdapt_soil_control = mean(avgMic_noAdapt_soil_control(end-2:end));
avgMic_lastThreeYears_noAdapt_soil_warming = mean(avgMic_noAdapt_soil_warming(end-2:end));

avgMic_lastThreeYears_optDriv_litter_control = mean(avgMic_optDriv_litter_control(end-2:end));
avgMic_lastThreeYears_optDriv_litter_warming = mean(avgMic_optDriv_litter_warming(end-2:end));
avgMic_lastThreeYears_optDriv_soil_control = mean(avgMic_optDriv_soil_control(end-2:end));
avgMic_lastThreeYears_optDriv_soil_warming = mean(avgMic_optDriv_soil_warming(end-2:end));

avgMic_lastThreeYears_enzRig_litter_control = mean(avgMic_enzRig_litter_control(end-2:end));
avgMic_lastThreeYears_enzRig_litter_warming = mean(avgMic_enzRig_litter_warming(end-2:end));
avgMic_lastThreeYears_enzRig_soil_control = mean(avgMic_enzRig_soil_control(end-2:end));
avgMic_lastThreeYears_enzRig_soil_warming = mean(avgMic_enzRig_soil_warming(end-2:end));

% The relative changes (neg is decrease, pos is increase)
relChangeMic_noAdapt_litter = (avgMic_lastThreeYears_noAdapt_litter_warming - avgMic_lastThreeYears_noAdapt_litter_control) / avgMic_lastThreeYears_noAdapt_litter_control;
relChangeMic_noAdapt_soil = (avgMic_lastThreeYears_noAdapt_soil_warming - avgMic_lastThreeYears_noAdapt_soil_control) / avgMic_lastThreeYears_noAdapt_soil_control;

relChangeMic_optDriv_litter = (avgMic_lastThreeYears_optDriv_litter_warming - avgMic_lastThreeYears_optDriv_litter_control) / avgMic_lastThreeYears_optDriv_litter_control;
relChangeMic_optDriv_soil = (avgMic_lastThreeYears_optDriv_soil_warming - avgMic_lastThreeYears_optDriv_soil_control) / avgMic_lastThreeYears_optDriv_soil_control;

relChangeMic_enzRig_litter = (avgMic_lastThreeYears_enzRig_litter_warming - avgMic_lastThreeYears_enzRig_litter_control) / avgMic_lastThreeYears_enzRig_litter_control;
relChangeMic_enzRig_soil = (avgMic_lastThreeYears_enzRig_soil_warming - avgMic_lastThreeYears_enzRig_soil_control) / avgMic_lastThreeYears_enzRig_soil_control;

relChange_Mic_thermalAdapt_litter = mean([relChangeMic_optDriv_litter relChangeMic_enzRig_litter]);
relChange_Mic_thermalAdapt_soil = mean([relChangeMic_optDriv_soil relChangeMic_enzRig_soil]);

% The changes combined for the forest floor and the soil
relChangeMic_noAdapt = mean([relChangeMic_noAdapt_litter relChangeMic_noAdapt_soil]);
relChangeMic_optDriv = mean([relChangeMic_optDriv_litter relChangeMic_optDriv_soil]);
relChangeMic_enzRig = mean([relChangeMic_enzRig_litter relChangeMic_enzRig_soil]);

% ------------------------------------------------
% The relative decrease in substrate in the heated plots
% ------------------------------------------------

date1 = datetime('01/01/2000', 'format', 'dd/MM/yyyy');
date2 = datetime('31/12/2002', 'format', 'dd/MM/yyyy');

date3 = datetime('01/01/2014', 'format', 'dd/MM/yyyy');
date4 = datetime('31/12/2016', 'format', 'dd/MM/yyyy');

r1 = find(dates_all == date1);
r2 = find(dates_all == date2);

r3 = find(dates_all == date3);
r4 = find(dates_all == date4);

% % No adaptation
% C1 = mean(noAdapt_monomers_warming_litter(r1:r2));
% C2 = mean(noAdapt_monomers_warming_litter(r3:r4));
% dC_ff_noAdapt = (C2 - C1)/C1;
% 
% C3 = mean(noAdapt_monomers_warming_soil(r1:r2));
% C4 = mean(noAdapt_monomers_warming_soil(r3:r4));
% dC_soil_noAdapt = (C4 - C3)/C3;
% 
% C5 = mean(noAdapt_monomers_warming_litter(r1:r2) + noAdapt_monomers_warming_soil(r1:r2));
% C6 = mean(noAdapt_monomers_warming_litter(r3:r4) + noAdapt_monomers_warming_soil(r3:r4));
% dC_tot_noAdapt = (C6 - C5)/C5;
% 
% % Optimum driven
% C1 = mean(optDriv_monomers_warming_litter(r1:r2));
% C2 = mean(optDriv_monomers_warming_litter(r3:r4));
% dC_ff_optDriv = (C2 - C1)/C1;
% 
% C3 = mean(optDriv_monomers_warming_soil(r1:r2));
% C4 = mean(optDriv_monomers_warming_soil(r3:r4));
% dC_soil_optDriv = (C4 - C3)/C3;
% 
% C5 = mean(optDriv_monomers_warming_litter(r1:r2) + optDriv_monomers_warming_soil(r1:r2));
% C6 = mean(optDriv_monomers_warming_litter(r3:r4) + optDriv_monomers_warming_soil(r3:r4));
% dC_tot_optDriv = (C6 - C5)/C5;
% 
% % Enzyme rigidity
% C1 = mean(enzRig_monomers_warming_litter(r1:r2));
% C2 = mean(enzRig_monomers_warming_litter(r3:r4));
% dC_ff_enzRig = (C2 - C1)/C1;
% 
% C3 = mean(enzRig_monomers_warming_soil(r1:r2));
% C4 = mean(enzRig_monomers_warming_soil(r3:r4));
% dC_soil_enzRig = (C4 - C3)/C3;
% 
% C5 = mean(enzRig_monomers_warming_litter(r1:r2) + enzRig_monomers_warming_soil(r1:r2));
% C6 = mean(enzRig_monomers_warming_litter(r3:r4) + enzRig_monomers_warming_soil(r3:r4));
% dC_tot_enzRig = (C6 - C5)/C5;
% 
% dC_tot_thermAdapt = mean([dC_tot_optDriv dC_tot_enzRig]);

% No adaptation
C1 = mean(noAdapt_monomers_control_litter(r3:r4));
C2 = mean(noAdapt_monomers_warming_litter(r3:r4));
dC_ff_noAdapt = (C2 - C1)/C1;

C3 = mean(noAdapt_monomers_control_soil(r3:r4));
C4 = mean(noAdapt_monomers_warming_soil(r3:r4));
dC_soil_noAdapt = (C4 - C3)/C3;

C5 = mean(noAdapt_monomers_control_litter(r3:r4) + noAdapt_monomers_control_soil(r3:r4));
C6 = mean(noAdapt_monomers_warming_litter(r3:r4) + noAdapt_monomers_warming_soil(r3:r4));
dC_tot_noAdapt = (C6 - C5)/C5;

% Optimum driven
C1 = mean(optDriv_monomers_control_litter(r3:r4));
C2 = mean(optDriv_monomers_warming_litter(r3:r4));
dC_ff_optDriv = (C2 - C1)/C1;

C3 = mean(optDriv_monomers_control_soil(r3:r4));
C4 = mean(optDriv_monomers_warming_soil(r3:r4));
dC_soil_optDriv = (C4 - C3)/C3;

C5 = mean(optDriv_monomers_control_litter(r3:r4) + optDriv_monomers_control_soil(r3:r4));
C6 = mean(optDriv_monomers_warming_litter(r3:r4) + optDriv_monomers_warming_soil(r3:r4));
dC_tot_optDriv = (C6 - C5)/C5;

% Enzyme rigidity
C1 = mean(enzRig_monomers_control_litter(r3:r4));
C2 = mean(enzRig_monomers_warming_litter(r3:r4));
dC_ff_enzRig = (C2 - C1)/C1;

C3 = mean(enzRig_monomers_control_soil(r3:r4));
C4 = mean(enzRig_monomers_warming_soil(r3:r4));
dC_soil_enzRig = (C4 - C3)/C3;

C5 = mean(enzRig_monomers_control_litter(r3:r4) + enzRig_monomers_control_soil(r3:r4));
C6 = mean(enzRig_monomers_warming_litter(r3:r4) + enzRig_monomers_warming_soil(r3:r4));
dC_tot_enzRig = (C6 - C5)/C5;

dC_tot_thermAdapt = mean([dC_tot_optDriv dC_tot_enzRig]);

% ----------------------------------------------------------------------
% The change in substrate during the last 3 simulated years is calculated
% ----------------------------------------------------------------------

avgSub_lastThreeYears_noAdapt_litter_control = mean(avgSub_noAdapt_litter_control(end-2:end));
avgSub_lastThreeYears_noAdapt_litter_warming = mean(avgSub_noAdapt_litter_warming(end-2:end));
avgSub_lastThreeYears_noAdapt_soil_control = mean(avgSub_noAdapt_soil_control(end-2:end));
avgSub_lastThreeYears_noAdapt_soil_warming = mean(avgSub_noAdapt_soil_warming(end-2:end));

avgSub_lastThreeYears_optDriv_litter_control = mean(avgSub_optDriv_litter_control(end-2:end));
avgSub_lastThreeYears_optDriv_litter_warming = mean(avgSub_optDriv_litter_warming(end-2:end));
avgSub_lastThreeYears_optDriv_soil_control = mean(avgSub_optDriv_soil_control(end-2:end));
avgSub_lastThreeYears_optDriv_soil_warming = mean(avgSub_optDriv_soil_warming(end-2:end));

avgSub_lastThreeYears_enzRig_litter_control = mean(avgSub_enzRig_litter_control(end-2:end));
avgSub_lastThreeYears_enzRig_litter_warming = mean(avgSub_enzRig_litter_warming(end-2:end));
avgSub_lastThreeYears_enzRig_soil_control = mean(avgSub_enzRig_soil_control(end-2:end));
avgSub_lastThreeYears_enzRig_soil_warming = mean(avgSub_enzRig_soil_warming(end-2:end));

% The relative changes (neg is decrease, pos is increase)
relChangeSub_noAdapt_litter = (avgSub_lastThreeYears_noAdapt_litter_warming - avgSub_lastThreeYears_noAdapt_litter_control) / avgSub_lastThreeYears_noAdapt_litter_control;
relChangeSub_noAdapt_soil = (avgSub_lastThreeYears_noAdapt_soil_warming - avgSub_lastThreeYears_noAdapt_soil_control) / avgSub_lastThreeYears_noAdapt_soil_control;

relChangeSub_optDriv_litter = (avgSub_lastThreeYears_optDriv_litter_warming - avgSub_lastThreeYears_optDriv_litter_control) / avgSub_lastThreeYears_optDriv_litter_control;
relChangeSub_optDriv_soil = (avgSub_lastThreeYears_optDriv_soil_warming - avgSub_lastThreeYears_optDriv_soil_control) / avgSub_lastThreeYears_optDriv_soil_control;

relChangeSub_enzRig_litter = (avgSub_lastThreeYears_enzRig_litter_warming - avgSub_lastThreeYears_enzRig_litter_control) / avgSub_lastThreeYears_enzRig_litter_control;
relChangeSub_enzRig_soil = (avgSub_lastThreeYears_enzRig_soil_warming - avgSub_lastThreeYears_enzRig_soil_control) / avgSub_lastThreeYears_enzRig_soil_control;

% The changes combined for the forest floor and the soil
relChangeSub_noAdapt = mean([relChangeSub_noAdapt_litter relChangeSub_noAdapt_soil]);
relChangeSub_optDriv = mean([relChangeSub_optDriv_litter relChangeSub_optDriv_soil]);
relChangeSub_enzRig = mean([relChangeSub_enzRig_litter relChangeSub_enzRig_soil]);

% The average change of both thermal adapation scenarios
reChange_thermalAdapt_litter = mean([relChangeSub_optDriv_litter relChangeSub_enzRig_litter]);
reChange_thermalAdapt_soil = mean([relChangeSub_optDriv_soil relChangeSub_enzRig_soil]);

% ---------------------------------
% The amount of microbes is plotted
% ---------------------------------

noAdapt_mbc_control = noAdapt_mbc_control_litter + noAdapt_mbc_control_soil;
noAdapt_mbc_warming = noAdapt_mbc_warming_litter + noAdapt_mbc_warming_soil;

optDriv_mbc_control = optDriv_mbc_control_litter + optDriv_mbc_control_soil;
optDriv_mbc_warming = optDriv_mbc_warming_litter + optDriv_mbc_warming_soil;

date1 = datetime('01/01/2000', 'Format', 'dd/MM/yyyy');
date2 = datetime('13/12/2016', 'Format', 'dd/MM/yyyy');

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 16*3 8*3], 'color', [1 1 1])

% Forest floor - no adaptation
subplot(2,2,1) 
hold on
p1 = plot(dates_all, noAdapt_mbc_control_litter(2:end), 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, noAdapt_mbc_warming_litter(2:end), 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgMic_noAdapt_litter_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgMic_noAdapt_litter_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 3.7], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 3.7])
xlim([date1 date2])
xlabel('Year')
ylabel('Carbon (g m^{-2})')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'southeast');
legend('box','off')
set(gca, 'FontSize', 14)
title('Microbial carbon in the organic horizon - no adaptation', 'FontSize', 14)

% Soil - no adaptation
subplot(2,2,3) 
hold on
p1 = plot(dates_all, noAdapt_mbc_control_soil(2:end), 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, noAdapt_mbc_warming_soil(2:end), 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgMic_noAdapt_soil_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgMic_noAdapt_soil_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 2], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 2])
xlim([date1 date2])
xlabel('Year')
ylabel('Carbon (g m^{-2})')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'southeast');
legend('box','off')
set(gca, 'FontSize', 14)
title('Microbial carbon in the soil - no adaptation', 'FontSize', 14)

% Forest floor - thermal adaptation
subplot(2,2,2)
hold on
p1 = plot(dates_all, thermalAdapt_mbc_control_litter(2:end), 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, thermalAdapt_mbc_warming_litter(2:end), 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgMic_thermalAdapt_litter_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgMic_thermalAdapt_litter_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 3.7], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 3.7])
xlim([date1 date2])
xlabel('Year')
ylabel('Carbon (g m^{-2})')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'southeast');
legend('box','off')
set(gca, 'FontSize', 14)
title('Microbial carbon in the organic horizon - thermal adaptation', 'FontSize', 14)

% Soil - thermal adaptation
subplot(2,2,4)
hold on
p1 = plot(dates_all, thermalAdapt_mbc_control_soil(2:end), 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, thermalAdapt_mbc_warming_soil(2:end), 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgMic_thermalAdapt_soil_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgMic_thermalAdapt_soil_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 2], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 2])
xlim([date1 date2])
xlabel('Year')
ylabel('Carbon (g m^{-2})')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'southeast');
legend('box','off')
set(gca, 'FontSize', 14)
title('Microbial carbon in the soil - thermal adaptation', 'FontSize', 14)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
% text(0.37,0.98,'Organic carbon stocks','fontweight','bold','fontsize',14)

% Letters
text(0.09, 0.945, '(A)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.52, 0.945, '(B)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.09, 0.47, '(C)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.52, 0.47, '(D)', 'FontSize', 14, 'FontWeight', 'bold')

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Changes in microbial carbon.tiff' -r300 -nocrop -transparent
% =========================================================================

% ----------------------------------
% The amount of substrate is plotted
% ----------------------------------

date1 = datetime('01/01/2000', 'Format', 'dd/MM/yyyy');
date2 = datetime('13/12/2016', 'Format', 'dd/MM/yyyy');

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 16*3 8*3], 'color', [1 1 1])

% Forest floor - no adaptation
subplot(2,2,1) 
hold on
p1 = plot(dates_all, noAdapt_monomers_control_litter(2:end), 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, noAdapt_monomers_warming_litter(2:end), 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgSub_noAdapt_litter_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgSub_noAdapt_litter_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 35], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 10])
xlim([date1 date2])
xlabel('Year')
ylabel('Carbon (g m^{-2})')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'northeast');
legend('box','off')
set(gca, 'FontSize', 14)
title('Available substrate in the organic horizon - no adaptation', 'FontSize', 14)

% Soil - no adaptation
subplot(2,2,3) 
hold on
p1 = plot(dates_all, noAdapt_monomers_control_soil(2:end), 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, noAdapt_monomers_warming_soil(2:end), 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgSub_noAdapt_soil_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgSub_noAdapt_soil_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 7], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 1.2])
xlim([date1 date2])
xlabel('Year')
ylabel('Carbon (g m^{-2})')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'northeast');
legend('box','off')
set(gca, 'FontSize', 14)
title('Available substrate in the soil - no adaptation', 'FontSize', 14)

% Forest floor - thermal adaptation
subplot(2,2,2)
hold on
p1 = plot(dates_all, thermalAdapt_monomers_control_litter(2:end), 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, thermalAdapt_monomers_warming_litter(2:end), 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgSub_thermalAdapt_litter_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgSub_thermalAdapt_litter_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 35], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 10])
xlim([date1 date2])
xlabel('Year')
ylabel('Carbon (g m^{-2})')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'northeast');
legend('box','off')
set(gca, 'FontSize', 14)
title('Available substrate in the organic horizon - thermal adaptation', 'FontSize', 14)

% Soil - thermal adaptation
subplot(2,2,4)
hold on
p1 = plot(dates_all, thermalAdapt_monomers_control_soil(2:end), 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, thermalAdapt_monomers_warming_soil(2:end), 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgSub_thermalAdapt_soil_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgSub_thermalAdapt_soil_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 7], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 1.2])
xlim([date1 date2])
xlabel('Year')
ylabel('Carbon (g m^{-2})')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'northeast');
legend('box','off')
set(gca, 'FontSize', 14)
title('Available substrate in the soil - thermal adaptation', 'FontSize', 14)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
% text(0.37,0.98,'Organic carbon stocks','fontweight','bold','fontsize',14)

% Letters
text(0.09, 0.945, '(A)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.52, 0.945, '(B)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.09, 0.47, '(C)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.52, 0.47, '(D)', 'FontSize', 14, 'FontWeight', 'bold')

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Changes in available substrate.tiff' -r300 -nocrop -transparent
% =========================================================================

% ---------------------------------------
% The ratio of MBC to substrate is plotted
% ---------------------------------------

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 16*3 8*3], 'color', [1 1 1])

% Control - litter
subplot(2,2,1) 
hold on
p1 = plot(dates_all, noAdapt_mbc_to_substr_control_litter(2:end), 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, noAdapt_mbc_to_substr_warming_litter(2:end), 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgRatio_noAdapt_litter_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgRatio_noAdapt_litter_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 0.8], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 0.8])
xlim([date1 date2])
xlabel('Year')
ylabel('Ratio (-)')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'northeast');
legend('box','off')
set(gca, 'FontSize', 14)
title('Ratio of microbes to substrate - litter - no adaptation', 'FontSize', 14)

% Control - soil
subplot(2,2,3) 
hold on
p1 = plot(dates_all, noAdapt_mbc_to_substr_control_soil(2:end), 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, noAdapt_mbc_to_substr_warming_soil(2:end), 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgRatio_noAdapt_soil_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgRatio_noAdapt_soil_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 1.5], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 1.5])
xlim([date1 date2])
xlabel('Year')
ylabel('Ratio (-)')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'northeast');
legend('box','off')
set(gca, 'FontSize', 14)
title('Ratio of microbes to substrate - soil - no adaptation', 'FontSize', 14)

% Warming - litter
subplot(2,2,2)
hold on
p1 = plot(dates_all, optDriv_mbc_to_substr_control_litter(2:end), 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, optDriv_mbc_to_substr_warming_litter(2:end), 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgRatio_optDriv_litter_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgRatio_optDriv_litter_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 0.8], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 0.8])
xlim([date1 date2])
xlabel('Year')
ylabel('Ratio (-)')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'northeast');
legend('box','off')
set(gca, 'FontSize', 14)
title('Ratio of microbes to substrate - litter - optimum driven', 'FontSize', 14)

% Warming - litter
subplot(2,2,4)
hold on
p1 = plot(dates_all, optDriv_mbc_to_substr_control_soil(2:end), 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, optDriv_mbc_to_substr_warming_soil(2:end), 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgRatio_optDriv_soil_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgRatio_optDriv_soil_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 1.5], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 1.5])
xlim([date1 date2])
xlabel('Year')
ylabel('Ratio (-)')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'northeast');
legend('box','off')
set(gca, 'FontSize', 14)
title('Ratio of microbes to substrate - soil - optimum driven', 'FontSize', 14)

%% Mass specific respiration is plotted

% -------------------------
% Repired CO2 is calculated
% -------------------------

noAdapt_CO2_rStrat_litter_control = diff(noAdaptData.Cpools_litter_control(1:end,14))./noAdaptData.dt_array_full;
noAdapt_CO2_kStrat_litter_control = diff(noAdaptData.Cpools_litter_control(1:end,15))./noAdaptData.dt_array_full;
noAdapt_CO2_soil_control = diff(noAdaptData.Cpools_soil_control(1:end,9))./noAdaptData.dt_array_full;
noAdapt_CO2_control_litter = noAdapt_CO2_rStrat_litter_control + noAdapt_CO2_kStrat_litter_control;
noAdapt_CO2_control_soil = noAdapt_CO2_soil_control;

noAdapt_CO2_rStrat_litter_warming = diff(noAdaptData.Cpools_litter_warming(1:end,14))./noAdaptData.dt_array_full;
noAdapt_CO2_kStrat_litter_warming = diff(noAdaptData.Cpools_litter_warming(1:end,15))./noAdaptData.dt_array_full;
noAdapt_CO2_soil_warming = diff(noAdaptData.Cpools_soil_warming(1:end,9))./noAdaptData.dt_array_full;
noAdapt_CO2_warming_litter = noAdapt_CO2_rStrat_litter_warming + noAdapt_CO2_kStrat_litter_warming;
noAdapt_CO2_warming_soil = noAdapt_CO2_soil_warming;

optDriv_CO2_rStrat_litter_control = diff(optimumDrivenData.Cpools_litter_control(1:end,14))./optimumDrivenData.dt_array_full;
optDriv_CO2_kStrat_litter_control = diff(optimumDrivenData.Cpools_litter_control(1:end,15))./optimumDrivenData.dt_array_full;
optDriv_CO2_soil_control = diff(optimumDrivenData.Cpools_soil_control(1:end,9))./optimumDrivenData.dt_array_full;
optDriv_CO2_control_litter = optDriv_CO2_rStrat_litter_control + optDriv_CO2_kStrat_litter_control;
optDriv_CO2_control_soil = optDriv_CO2_soil_control;

optDriv_CO2_rStrat_litter_warming = diff(optimumDrivenData.Cpools_litter_warming(1:end,14))./optimumDrivenData.dt_array_full;
optDriv_CO2_kStrat_litter_warming = diff(optimumDrivenData.Cpools_litter_warming(1:end,15))./optimumDrivenData.dt_array_full;
optDriv_CO2_soil_warming = diff(optimumDrivenData.Cpools_soil_warming(1:end,9))./optimumDrivenData.dt_array_full;
optDriv_CO2_warming_litter = optDriv_CO2_rStrat_litter_warming + optDriv_CO2_kStrat_litter_warming;
optDriv_CO2_warming_soil = optDriv_CO2_soil_warming;

enzRig_CO2_rStrat_litter_control = diff(enzRigData.Cpools_litter_control(1:end,14))./enzRigData.dt_array_full;
enzRig_CO2_kStrat_litter_control = diff(enzRigData.Cpools_litter_control(1:end,15))./enzRigData.dt_array_full;
enzRig_CO2_soil_control = diff(enzRigData.Cpools_soil_control(1:end,9))./enzRigData.dt_array_full;
enzRig_CO2_control_litter = enzRig_CO2_rStrat_litter_control + enzRig_CO2_kStrat_litter_control;
enzRig_CO2_control_soil = enzRig_CO2_soil_control;

enzRig_CO2_rStrat_litter_warming = diff(enzRigData.Cpools_litter_warming(1:end,14))./enzRigData.dt_array_full;
enzRig_CO2_kStrat_litter_warming = diff(enzRigData.Cpools_litter_warming(1:end,15))./enzRigData.dt_array_full;
enzRig_CO2_soil_warming = diff(enzRigData.Cpools_soil_warming(1:end,9))./enzRigData.dt_array_full;
enzRig_CO2_warming_litter = enzRig_CO2_rStrat_litter_warming + enzRig_CO2_kStrat_litter_warming;
enzRig_CO2_warming_soil = enzRig_CO2_soil_warming;

% -------------------------------------------
% The mass specific respiration is calculated
% -------------------------------------------

Rmass_noAdapt_control_litter = noAdapt_CO2_control_litter./noAdapt_mbc_control_litter(2:end);
Rmass_noAdapt_control_soil = noAdapt_CO2_control_soil./noAdapt_mbc_control_soil(2:end);
Rmass_noAdapt_warming_litter = noAdapt_CO2_warming_litter./noAdapt_mbc_warming_litter(2:end);
Rmass_noAdapt_warming_soil = noAdapt_CO2_warming_soil./noAdapt_mbc_warming_soil(2:end);

Rmass_optDriv_control_litter = optDriv_CO2_control_litter./optDriv_mbc_control_litter(2:end);
Rmass_optDriv_control_soil = optDriv_CO2_control_soil./optDriv_mbc_control_soil(2:end);
Rmass_optDriv_warming_litter = optDriv_CO2_warming_litter./optDriv_mbc_warming_litter(2:end);
Rmass_optDriv_warming_soil = optDriv_CO2_warming_soil./optDriv_mbc_warming_soil(2:end);

Rmass_enzRig_control_litter = enzRig_CO2_control_litter./enzRig_mbc_control_litter(2:end);
Rmass_enzRig_control_soil = enzRig_CO2_control_soil./enzRig_mbc_control_soil(2:end);
Rmass_enzRig_warming_litter = enzRig_CO2_warming_litter./enzRig_mbc_warming_litter(2:end);
Rmass_enzRig_warming_soil = enzRig_CO2_warming_soil./enzRig_mbc_warming_soil(2:end);

% The average for both thermal adaptation scenarios is calculated
Rmass_thermalAdapt_control_litter = mean([Rmass_optDriv_control_litter Rmass_enzRig_control_litter], 2);
Rmass_thermalAdapt_control_soil = mean([Rmass_optDriv_control_soil Rmass_enzRig_control_soil], 2);
Rmass_thermalAdapt_warming_litter = mean([Rmass_optDriv_warming_litter Rmass_enzRig_warming_litter], 2);
Rmass_thermalAdapt_warming_soil = mean([Rmass_optDriv_warming_soil Rmass_enzRig_warming_soil], 2);

% --------------------------------------
% The average annual Rmass is calculated
% --------------------------------------

% The individual yearsare identified
uniqueYears = unique(dates_all.Year);

% Empty array to store the averages in
avgRmass_noAdapt_litter_control = NaN(numel(uniqueYears),1);
avgRmass_noAdapt_litter_warming = NaN(numel(uniqueYears),1);
avgRmass_noAdapt_soil_control = NaN(numel(uniqueYears),1);
avgRmass_noAdapt_soil_warming = NaN(numel(uniqueYears),1);

avgRmass_optDriv_litter_control = NaN(numel(uniqueYears),1);
avgRmass_optDriv_litter_warming = NaN(numel(uniqueYears),1);
avgRmass_optDriv_soil_control = NaN(numel(uniqueYears),1);
avgRmass_optDriv_soil_warming = NaN(numel(uniqueYears),1);

avgRmass_enzRig_litter_control = NaN(numel(uniqueYears),1);
avgRmass_enzRig_litter_warming = NaN(numel(uniqueYears),1);
avgRmass_enzRig_soil_control = NaN(numel(uniqueYears),1);
avgRmass_optDriv_soil_warming = NaN(numel(uniqueYears),1);

% A loop to calculate the annual averages
for ii = 1:numel(uniqueYears)
    
    yr = uniqueYears(ii);
    rowNum = find(dates_all.Year == yr);
    
    % The averages are calculated
    avgRmass_noAdapt_litter_control(ii,1) = mean(Rmass_noAdapt_control_litter(rowNum));
    avgRmass_noAdapt_litter_warming(ii,1) = mean(Rmass_noAdapt_warming_litter(rowNum));
    avgRmass_noAdapt_soil_control(ii,1) = mean(Rmass_noAdapt_control_soil(rowNum));
    avgRmass_noAdapt_soil_warming(ii,1) = mean(Rmass_noAdapt_warming_soil(rowNum));
    
    avgRmass_optDriv_litter_control(ii,1) = mean(Rmass_optDriv_control_litter(rowNum));
    avgRmass_optDriv_litter_warming(ii,1) = mean(Rmass_optDriv_warming_litter(rowNum));
    avgRmass_optDriv_soil_control(ii,1) = mean(Rmass_optDriv_control_soil(rowNum));
    avgRmass_optDriv_soil_warming(ii,1) = mean(Rmass_optDriv_warming_soil(rowNum));
    
    avgRmass_enzRig_litter_control(ii,1) = mean(Rmass_enzRig_control_litter(rowNum));
    avgRmass_enzRig_litter_warming(ii,1) = mean(Rmass_enzRig_warming_litter(rowNum));
    avgRmass_enzRig_soil_control(ii,1) = mean(Rmass_enzRig_control_soil(rowNum));
    avgRmass_enzRig_soil_warming(ii,1) = mean(Rmass_enzRig_warming_soil(rowNum));
    
end

% The average for both thermal adaptation scenarios is calculated
avgRmass_thermalAdapt_litter_control = mean([avgRmass_optDriv_litter_control avgRmass_enzRig_litter_control], 2);
avgRmass_thermalAdapt_litter_warming = mean([avgRmass_optDriv_litter_warming avgRmass_enzRig_litter_warming], 2);
avgRmass_thermalAdapt_soil_control = mean([avgRmass_optDriv_soil_control avgRmass_enzRig_soil_control], 2);
avgRmass_thermalAdapt_soil_warming = mean([avgRmass_optDriv_soil_warming avgRmass_enzRig_soil_warming], 2);

% -------------------------------------------------------------------
% The Change in Rmass during the last 3 simulated years is calculated
% -------------------------------------------------------------------

avgRmass_lastThreeYears_noAdapt_litter_control = mean(avgRmass_noAdapt_litter_control(end-2:end));
avgRmass_lastThreeYears_noAdapt_litter_warming = mean(avgRmass_noAdapt_litter_warming(end-2:end));
avgRmass_lastThreeYears_noAdapt_soil_control = mean(avgRmass_noAdapt_soil_control(end-2:end));
avgRmass_lastThreeYears_noAdapt_soil_warming = mean(avgRmass_noAdapt_soil_warming(end-2:end));

avgRmass_lastThreeYears_optDriv_litter_control = mean(avgRmass_optDriv_litter_control(end-2:end));
avgRmass_lastThreeYears_optDriv_litter_warming = mean(avgRmass_optDriv_litter_warming(end-2:end));
avgRmass_lastThreeYears_optDriv_soil_control = mean(avgRmass_optDriv_soil_control(end-2:end));
avgRmass_lastThreeYears_optDriv_soil_warming = mean(avgRmass_optDriv_soil_warming(end-2:end));

avgRmass_lastThreeYears_enzRig_litter_control = mean(avgRmass_enzRig_litter_control(end-2:end));
avgRmass_lastThreeYears_enzRig_litter_warming = mean(avgRmass_enzRig_litter_warming(end-2:end));
avgRmass_lastThreeYears_enzRig_soil_control = mean(avgRmass_enzRig_soil_control(end-2:end));
avgRmass_lastThreeYears_enzRig_soil_warming = mean(avgRmass_enzRig_soil_warming(end-2:end));

% The relative changes (neg is decrease, pos is increase)
relChangeRmass_noAdapt_litter = (avgRmass_lastThreeYears_noAdapt_litter_warming - avgRmass_lastThreeYears_noAdapt_litter_control) / avgRmass_lastThreeYears_noAdapt_litter_control;
relChangeRmass_noAdapt_soil = (avgRmass_lastThreeYears_noAdapt_soil_warming - avgRmass_lastThreeYears_noAdapt_soil_control) / avgRmass_lastThreeYears_noAdapt_soil_control;

relChangeRmass_optDriv_litter = (avgRmass_lastThreeYears_optDriv_litter_warming - avgRmass_lastThreeYears_optDriv_litter_control) / avgRmass_lastThreeYears_optDriv_litter_control;
relChangeRmass_optDriv_soil = (avgRmass_lastThreeYears_optDriv_soil_warming - avgRmass_lastThreeYears_optDriv_soil_control) / avgRmass_lastThreeYears_optDriv_soil_control;

relChangeRmass_enzRig_litter = (avgRmass_lastThreeYears_enzRig_litter_warming - avgRmass_lastThreeYears_enzRig_litter_control) / avgRmass_lastThreeYears_enzRig_litter_control;
relChangeRmass_enzRig_soil = (avgRmass_lastThreeYears_enzRig_soil_warming - avgRmass_lastThreeYears_enzRig_soil_control) / avgRmass_lastThreeYears_enzRig_soil_control;

relChangeRmass_thermalAdapt_litter = mean([relChangeRmass_optDriv_litter relChangeRmass_enzRig_litter]);
relChangeRmass_thermalAdapt_soil = mean([relChangeRmass_optDriv_soil relChangeRmass_enzRig_soil]);

% The changes combined for the forest floor and the soil
relChangeRmass_noAdapt = mean([relChangeRmass_noAdapt_litter relChangeRmass_noAdapt_soil]);
relChangeRmass_optDriv = mean([relChangeRmass_optDriv_litter relChangeRmass_optDriv_soil]);
relChangeRmass_enzRig = mean([relChangeRmass_enzRig_litter relChangeRmass_enzRig_soil]);


% -------------------------------------------------------------------------
% Rmass is plotted
% -------------------------------------------------------------------------

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 25*1.5 10*1.5], 'color', [1 1 1])

% Forest floor - no adaptation
subplot(2,2,1) 
hold on
p1 = plot(dates_all, Rmass_noAdapt_control_litter, 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, Rmass_noAdapt_warming_litter, 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgRmass_noAdapt_litter_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgRmass_noAdapt_litter_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 2.5], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 0.6])
xlim([date1 date2])
xlabel('Year')
ylabel('R_{mass}')
title('R_{mass} - organic horizon - no adaptation')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'northeast');
legend('box','off')
set(gca, 'FontSize', 14)

% Soil - no adaptation
subplot(2,2,3) 
hold on
p1 = plot(dates_all, Rmass_noAdapt_control_soil, 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, Rmass_noAdapt_warming_soil, 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgRmass_noAdapt_soil_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgRmass_noAdapt_soil_warming, 40, [179,0,0]./255, 'filled');
% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 12.5], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 3])
xlim([date1 date2])
xlabel('Year')
ylabel('R_{mass}')
title('R_{mass} - soil - no adaptation')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'northeast');
legend('box','off')
set(gca, 'FontSize', 14)

% Forest floor - thermal adaptation
subplot(2,2,2)
hold on
p1 = plot(dates_all, Rmass_thermalAdapt_control_litter, 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, Rmass_thermalAdapt_warming_litter, 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgRmass_thermalAdapt_litter_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgRmass_thermalAdapt_litter_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 2.5], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 0.6])
xlim([date1 date2])
xlabel('Year')
ylabel('R_{mass}')
title('R_{mass} - organic horizon - thermal adaptation')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'northeast');
legend('box','off')
set(gca, 'FontSize', 14)

% Soil - thermal adaptation
subplot(2,2,4)
hold on
p1 = plot(dates_all, Rmass_thermalAdapt_control_soil, 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, Rmass_thermalAdapt_warming_soil, 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgRmass_thermalAdapt_soil_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgRmass_thermalAdapt_soil_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 12.5], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 3])
xlim([date1 date2])
xlabel('Year')
ylabel('R_{mass}')
title('R_{mass} - soil - thermal adaptation')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'northeast');
legend('box','off')
set(gca, 'FontSize', 14)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
% text(0.37,0.98,'Organic carbon stocks','fontweight','bold','fontsize',14)

% Letters
text(0.09, 0.945, '(A)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.52, 0.945, '(B)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.09, 0.47, '(C)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.52, 0.47, '(D)', 'FontSize', 14, 'FontWeight', 'bold')

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Changes in Rmass.tiff' -r300 -nocrop -transparent
% =========================================================================

%% The carbon use efficiency is plotted

% ------------------------------------
% The data is retrieved: no adaptation
% ------------------------------------

noAdapt_CUE_litter_rStrat_control = noAdaptData.CUE.CUE_litter_rStrat_control;
noAdapt_CUE_litter_rStrat_warming_tmp = noAdaptData.CUE.CUE_litter_rStrat_warming;

noAdapt_CUE_litter_kStrat_control = noAdaptData.CUE.CUE_litter_kStrat_control;
noAdapt_CUE_litter_kStrat_warming_tmp = noAdaptData.CUE.CUE_litter_kStrat_warming;

noAdapt_CUE_soil_control = noAdaptData.CUE.CUE_soil_control;
noAdapt_CUE_soil_warming_tmp = noAdaptData.CUE.CUE_soil_warming;

% The data for the spinup run is added to the data for the warmed run
noAdapt_CUE_litter_rStrat_warming = [noAdapt_CUE_litter_rStrat_control(1:end-numel(noAdapt_CUE_litter_rStrat_warming_tmp)) ;noAdapt_CUE_litter_rStrat_warming_tmp];
noAdapt_CUE_litter_kStrat_warming = [noAdapt_CUE_litter_kStrat_control(1:end-numel(noAdapt_CUE_litter_kStrat_warming_tmp)) ;noAdapt_CUE_litter_kStrat_warming_tmp];
noAdapt_CUE_soil_warming = [noAdapt_CUE_soil_control(1:end-numel(noAdapt_CUE_soil_warming_tmp)) ;noAdapt_CUE_soil_warming_tmp];

% ------------------------------------
% The data is retrieved: optimum driven
% ------------------------------------

optDriv_CUE_litter_rStrat_control = optimumDrivenData.CUE.CUE_litter_rStrat_control;
optDriv_CUE_litter_rStrat_warming_tmp = optimumDrivenData.CUE.CUE_litter_rStrat_warming;

optDriv_CUE_litter_kStrat_control = optimumDrivenData.CUE.CUE_litter_kStrat_control;
optDriv_CUE_litter_kStrat_warming_tmp = optimumDrivenData.CUE.CUE_litter_kStrat_warming;

optDriv_CUE_soil_control = optimumDrivenData.CUE.CUE_soil_control;
optDriv_CUE_soil_warming_tmp = optimumDrivenData.CUE.CUE_soil_warming;

% The data for the spinup run is added to the data for the warmed run
optDriv_CUE_litter_rStrat_warming = [optDriv_CUE_litter_rStrat_control(1:end-numel(optDriv_CUE_litter_rStrat_warming_tmp)) ;optDriv_CUE_litter_rStrat_warming_tmp];
optDriv_CUE_litter_kStrat_warming = [optDriv_CUE_litter_kStrat_control(1:end-numel(optDriv_CUE_litter_kStrat_warming_tmp)) ;optDriv_CUE_litter_kStrat_warming_tmp];
optDriv_CUE_soil_warming = [optDriv_CUE_soil_control(1:end-numel(optDriv_CUE_soil_warming_tmp)); optDriv_CUE_soil_warming_tmp];

% --------------------------------------
% The data is retrieved: enzyme rigidity
% --------------------------------------

enzRig_CUE_litter_rStrat_control = enzRigData.CUE.CUE_litter_rStrat_control;
enzRig_CUE_litter_rStrat_warming_tmp = enzRigData.CUE.CUE_litter_rStrat_warming;

enzRig_CUE_litter_kStrat_control = enzRigData.CUE.CUE_litter_kStrat_control;
enzRig_CUE_litter_kStrat_warming_tmp = enzRigData.CUE.CUE_litter_kStrat_warming;

enzRig_CUE_soil_control = enzRigData.CUE.CUE_soil_control;
enzRig_CUE_soil_warming_tmp = enzRigData.CUE.CUE_soil_warming;

% The data for the spinup run is added to the data for the warmed run
enzRig_CUE_litter_rStrat_warming = [enzRig_CUE_litter_rStrat_control(1:end-numel(enzRig_CUE_litter_rStrat_warming_tmp)) ;enzRig_CUE_litter_rStrat_warming_tmp];
enzRig_CUE_litter_kStrat_warming = [enzRig_CUE_litter_kStrat_control(1:end-numel(enzRig_CUE_litter_kStrat_warming_tmp)) ;enzRig_CUE_litter_kStrat_warming_tmp];
enzRig_CUE_soil_warming = [enzRig_CUE_soil_control(1:end-numel(enzRig_CUE_soil_warming_tmp)); enzRig_CUE_soil_warming_tmp];

% ------------------------------------------------------------
% The data for enzyme rigidity and optimum driven are combined
% ------------------------------------------------------------

thermalAdapt_CUE_litter_rStrat_control = mean([optDriv_CUE_litter_rStrat_control enzRig_CUE_litter_rStrat_control],2);
thermalAdapt_CUE_litter_rStrat_warming = mean([optDriv_CUE_litter_rStrat_warming enzRig_CUE_litter_rStrat_warming],2);

thermalAdapt_CUE_litter_kStrat_control = mean([optDriv_CUE_litter_kStrat_control enzRig_CUE_litter_kStrat_control],2);
thermalAdapt_CUE_litter_kStrat_warming = mean([optDriv_CUE_litter_kStrat_warming enzRig_CUE_litter_kStrat_warming],2);

thermalAdapt_CUE_soil_control = mean([optDriv_CUE_soil_control enzRig_CUE_soil_control],2);
thermalAdapt_CUE_soil_warming = mean([optDriv_CUE_soil_warming enzRig_CUE_soil_warming],2);

% ------------------------------------------------------
% The average annual carbon use efficiency is calculated
% ------------------------------------------------------

% The individual yearsare identified
uniqueYears = unique(dates_all.Year);

% Empty array to store the averages in
avgCUE_noAdapt_litter_rStrat_control = NaN(numel(uniqueYears),1);
avgCUE_noAdapt_litter_kStrat_control = NaN(numel(uniqueYears),1);
avgCUE_noAdapt_soil_control = NaN(numel(uniqueYears),1);
avgCUE_noAdapt_litter_rStrat_warming = NaN(numel(uniqueYears),1);
avgCUE_noAdapt_litter_kStrat_warming = NaN(numel(uniqueYears),1);
avgCUE_noAdapt_soil_warming = NaN(numel(uniqueYears),1);

avgCUE_optDriv_litter_rStrat_control = NaN(numel(uniqueYears),1);
avgCUE_optDriv_litter_kStrat_control = NaN(numel(uniqueYears),1);
avgCUE_optDriv_soil_control = NaN(numel(uniqueYears),1);
avgCUE_optDriv_litter_rStrat_warming = NaN(numel(uniqueYears),1);
avgCUE_optDriv_litter_kStrat_warming = NaN(numel(uniqueYears),1);
avgCUE_optDriv_soil_warming = NaN(numel(uniqueYears),1);

avgCUE_enzRig_litter_rStrat_control = NaN(numel(uniqueYears),1);
avgCUE_enzRig_litter_kStrat_control = NaN(numel(uniqueYears),1);
avgCUE_enzRig_soil_control = NaN(numel(uniqueYears),1);
avgCUE_enzRig_litter_rStrat_warming = NaN(numel(uniqueYears),1);
avgCUE_enzRig_litter_kStrat_warming = NaN(numel(uniqueYears),1);
avgCUE_enzRig_soil_warming = NaN(numel(uniqueYears),1);

% A loop to calculate the annual averages
for ii = 1:numel(uniqueYears)
    
    yr = uniqueYears(ii);
    rowNum = find(dates_all.Year == yr);
    
    % The averages are calculated
    avgCUE_noAdapt_litter_rStrat_control(ii,1) = mean(noAdapt_CUE_litter_rStrat_control(rowNum));
    avgCUE_noAdapt_litter_rStrat_warming(ii,1) = mean(noAdapt_CUE_litter_rStrat_warming(rowNum));
    
    avgCUE_noAdapt_litter_kStrat_control(ii,1) = mean(noAdapt_CUE_litter_kStrat_control(rowNum));
    avgCUE_noAdapt_litter_kStrat_warming(ii,1) = mean(noAdapt_CUE_litter_kStrat_warming(rowNum));
    
    avgCUE_noAdapt_soil_control(ii,1) = mean(noAdapt_CUE_soil_control(rowNum));
    avgCUE_noAdapt_soil_warming(ii,1) = mean(noAdapt_CUE_soil_warming(rowNum));
    
    avgCUE_optDriv_litter_rStrat_control(ii,1) = mean(optDriv_CUE_litter_rStrat_control(rowNum));
    avgCUE_optDriv_litter_rStrat_warming(ii,1) = mean(optDriv_CUE_litter_rStrat_warming(rowNum));
    
    avgCUE_optDriv_litter_kStrat_control(ii,1) = mean(optDriv_CUE_litter_kStrat_control(rowNum));
    avgCUE_optDriv_litter_kStrat_warming(ii,1) = mean(optDriv_CUE_litter_kStrat_warming(rowNum));
    
    avgCUE_optDriv_soil_control(ii,1) = mean(optDriv_CUE_soil_control(rowNum));
    avgCUE_optDriv_soil_warming(ii,1) = mean(optDriv_CUE_soil_warming(rowNum));
    
    avgCUE_enzRig_litter_rStrat_control(ii,1) = mean(enzRig_CUE_litter_rStrat_control(rowNum));
    avgCUE_enzRig_litter_rStrat_warming(ii,1) = mean(enzRig_CUE_litter_rStrat_warming(rowNum));
    
    avgCUE_enzRig_litter_kStrat_control(ii,1) = mean(enzRig_CUE_litter_kStrat_control(rowNum));
    avgCUE_enzRig_litter_kStrat_warming(ii,1) = mean(enzRig_CUE_litter_kStrat_warming(rowNum));
    
    avgCUE_enzRig_soil_control(ii,1) = mean(enzRig_CUE_soil_control(rowNum));
    avgCUE_enzRig_soil_warming(ii,1) = mean(enzRig_CUE_soil_warming(rowNum));
    
end

% The average of both thermal adaptation scenarios
avgCUE_thermalAdapt_litter_rStrat_control = mean([avgCUE_optDriv_litter_rStrat_control avgCUE_enzRig_litter_rStrat_control],2);
avgCUE_thermalAdapt_litter_rStrat_warming = mean([avgCUE_optDriv_litter_rStrat_warming avgCUE_enzRig_litter_rStrat_warming],2);

avgCUE_thermalAdapt_litter_kStrat_control = mean([avgCUE_optDriv_litter_kStrat_control avgCUE_enzRig_litter_kStrat_control],2);
avgCUE_thermalAdapt_litter_kStrat_warming = mean([avgCUE_optDriv_litter_kStrat_warming avgCUE_enzRig_litter_kStrat_warming],2);

avgCUE_thermalAdapt_soil_control = mean([avgCUE_optDriv_soil_control avgCUE_enzRig_soil_control],2);
avgCUE_thermalAdapt_soil_warming = mean([avgCUE_optDriv_soil_warming avgCUE_enzRig_soil_warming],2);

% -----------------------------------------------------------------
% The Change in CUE during the last 3 simulated years is calculated
% -----------------------------------------------------------------

avgCUE_lastThreeYears_noAdapt_litter_rStrat_control = mean(avgCUE_noAdapt_litter_rStrat_control(end-2:end));
avgCUE_lastThreeYears_noAdapt_litter_kStrat_control = mean(avgCUE_noAdapt_litter_kStrat_control(end-2:end));
avgCUE_lastThreeYears_noAdapt_litter_rStrat_warming = mean(avgCUE_noAdapt_litter_rStrat_warming(end-2:end));
avgCUE_lastThreeYears_noAdapt_litter_kStrat_warming = mean(avgCUE_noAdapt_litter_kStrat_warming(end-2:end));
avgCUE_lastThreeYears_noAdapt_soil_control = mean(avgCUE_noAdapt_soil_control(end-2:end));
avgCUE_lastThreeYears_noAdapt_soil_warming = mean(avgCUE_noAdapt_soil_warming(end-2:end));

avgCUE_lastThreeYears_optDriv_litter_rStrat_control = mean(avgCUE_optDriv_litter_rStrat_control(end-2:end));
avgCUE_lastThreeYears_optDriv_litter_kStrat_control = mean(avgCUE_optDriv_litter_kStrat_control(end-2:end));
avgCUE_lastThreeYears_optDriv_litter_rStrat_warming = mean(avgCUE_optDriv_litter_rStrat_warming(end-2:end));
avgCUE_lastThreeYears_optDriv_litter_kStrat_warming = mean(avgCUE_optDriv_litter_kStrat_warming(end-2:end));
avgCUE_lastThreeYears_optDriv_soil_control = mean(avgCUE_optDriv_soil_control(end-2:end));
avgCUE_lastThreeYears_optDriv_soil_warming = mean(avgCUE_optDriv_soil_warming(end-2:end));

avgCUE_lastThreeYears_enzRig_litter_rStrat_control = mean(avgCUE_enzRig_litter_rStrat_control(end-2:end));
avgCUE_lastThreeYears_enzRig_litter_kStrat_control = mean(avgCUE_enzRig_litter_kStrat_control(end-2:end));
avgCUE_lastThreeYears_enzRig_litter_rStrat_warming = mean(avgCUE_enzRig_litter_rStrat_warming(end-2:end));
avgCUE_lastThreeYears_enzRig_litter_kStrat_warming = mean(avgCUE_enzRig_litter_kStrat_warming(end-2:end));
avgCUE_lastThreeYears_enzRig_soil_control = mean(avgCUE_enzRig_soil_control(end-2:end));
avgCUE_lastThreeYears_enzRig_soil_warming = mean(avgCUE_enzRig_soil_warming(end-2:end));

avgCUE_lastThreeYears_noAdapt_litter_control = mean([avgCUE_lastThreeYears_noAdapt_litter_rStrat_control avgCUE_lastThreeYears_noAdapt_litter_kStrat_control]);
avgCUE_lastThreeYears_noAdapt_litter_warming = mean([avgCUE_lastThreeYears_noAdapt_litter_rStrat_warming avgCUE_lastThreeYears_noAdapt_litter_kStrat_warming]);
avgCUE_lastThreeYears_optDriv_litter_control = mean([avgCUE_lastThreeYears_optDriv_litter_rStrat_control avgCUE_lastThreeYears_optDriv_litter_kStrat_control]);
avgCUE_lastThreeYears_optDriv_litter_warming = mean([avgCUE_lastThreeYears_optDriv_litter_rStrat_warming avgCUE_lastThreeYears_optDriv_litter_kStrat_warming]);
avgCUE_lastThreeYears_enzRig_litter_control = mean([avgCUE_lastThreeYears_enzRig_litter_rStrat_control avgCUE_lastThreeYears_enzRig_litter_kStrat_control]);
avgCUE_lastThreeYears_enzRig_litter_warming = mean([avgCUE_lastThreeYears_enzRig_litter_rStrat_warming avgCUE_lastThreeYears_enzRig_litter_kStrat_warming]);

avgCUE_lastThreeyears_thermalAdapt_soil_control = mean([avgCUE_lastThreeYears_optDriv_soil_control avgCUE_lastThreeYears_enzRig_soil_control]);
avgCUE_lastThreeyears_thermalAdapt_soil_warming = mean([avgCUE_lastThreeYears_optDriv_soil_warming avgCUE_lastThreeYears_enzRig_soil_warming]);

% The relative changes (neg is decrease, pos is increase)
relChangeCUE_noAdapt_litter_rStrat = (avgCUE_lastThreeYears_noAdapt_litter_rStrat_warming - avgCUE_lastThreeYears_noAdapt_litter_rStrat_control) / avgCUE_lastThreeYears_noAdapt_litter_rStrat_control;
relChangeCUE_noAdapt_litter_kStrat = (avgCUE_lastThreeYears_noAdapt_litter_kStrat_warming - avgCUE_lastThreeYears_noAdapt_litter_kStrat_control) / avgCUE_lastThreeYears_noAdapt_litter_kStrat_control;
relChangeCUE_noAdapt_soil = (avgCUE_lastThreeYears_noAdapt_soil_warming - avgCUE_lastThreeYears_noAdapt_soil_control) / avgCUE_lastThreeYears_noAdapt_soil_control;

relChangeCUE_optDriv_litter_rStrat = (avgCUE_lastThreeYears_optDriv_litter_rStrat_warming - avgCUE_lastThreeYears_optDriv_litter_rStrat_control) / avgCUE_lastThreeYears_optDriv_litter_rStrat_control;
relChangeCUE_optDriv_litter_kStrat = (avgCUE_lastThreeYears_optDriv_litter_kStrat_warming - avgCUE_lastThreeYears_optDriv_litter_kStrat_control) / avgCUE_lastThreeYears_optDriv_litter_kStrat_control;
relChangeCUE_optDriv_soil = (avgCUE_lastThreeYears_optDriv_soil_warming - avgCUE_lastThreeYears_optDriv_soil_control) / avgCUE_lastThreeYears_optDriv_soil_control;

relChangeCUE_enzRig_litter_rStrat = (avgCUE_lastThreeYears_enzRig_litter_rStrat_warming - avgCUE_lastThreeYears_enzRig_litter_rStrat_control) / avgCUE_lastThreeYears_enzRig_litter_rStrat_control;
relChangeCUE_enzRig_litter_kStrat = (avgCUE_lastThreeYears_enzRig_litter_kStrat_warming - avgCUE_lastThreeYears_enzRig_litter_kStrat_control) / avgCUE_lastThreeYears_enzRig_litter_kStrat_control;
relChangeCUE_enzRig_soil = (avgCUE_lastThreeYears_enzRig_soil_warming - avgCUE_lastThreeYears_enzRig_soil_control) / avgCUE_lastThreeYears_enzRig_soil_control;

% The mean of the optimum driven and enzyme rigidity scenarios is calculated
relChangeCUE_thermalAdapt_litter_rStrat = mean([relChangeCUE_optDriv_litter_rStrat relChangeCUE_enzRig_litter_rStrat]);
relChangeCUE_thermalAdapt_litter_kStrat = mean([relChangeCUE_optDriv_litter_kStrat relChangeCUE_enzRig_litter_kStrat]);
relChangeCUE_thermalAdapt_soil = mean([relChangeCUE_optDriv_soil relChangeCUE_enzRig_soil]);

% The changes combined for the forest floor and the soil
% relChangeCUE_noAdapt_litter = mean([relChangeCUE_noAdapt_litter_rStrat relChangeCUE_noAdapt_litter_kStrat]);
% relChangeCUE_optDriv_litter = mean([relChangeCUE_optDriv_litter_rStrat relChangeCUE_optDriv_litter_kStrat]);

% ---------------------------------------
% The carbon use efficiencies are plotted
% ---------------------------------------

% Dates for plotting the averages
dates_averages = append(string(uniqueYears), '/07/01');
dates_averages = datetime(dates_averages, 'Format', 'yyyy/MM/dd');

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 25*1.5 25*1.5], 'color', [1 1 1])

% No adaptation - Control - litter - rStrat
subplot(3,2,1) 
hold on
p1 = plot(dates_all, noAdapt_CUE_litter_rStrat_control, 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, noAdapt_CUE_litter_rStrat_warming, 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgCUE_noAdapt_litter_rStrat_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgCUE_noAdapt_litter_rStrat_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 0.45], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 0.45])
xlim([date1 date2])
xlabel('Year')
ylabel('CUE')
title('CUE - organic horizon, r-strategists - no adaptation')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'southeast');
legend('box','off')
set(gca, 'FontSize', 14)

% No adaptation - Control - litter - kStrat
subplot(3,2,3) 
hold on
p1 = plot(dates_all, noAdapt_CUE_litter_kStrat_control, 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, noAdapt_CUE_litter_kStrat_warming, 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgCUE_noAdapt_litter_kStrat_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgCUE_noAdapt_litter_kStrat_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 0.45], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 0.45])
xlim([date1 date2])
xlabel('Year')
ylabel('CUE')
title('CUE - organic horizon, k-strategists - no adaptation')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'southeast');
legend('box','off')
set(gca, 'FontSize', 14)

% No adaptation - Control - soil
subplot(3,2,5) 
hold on
p1 = plot(dates_all, noAdapt_CUE_soil_control, 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, noAdapt_CUE_soil_warming, 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgCUE_noAdapt_soil_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgCUE_noAdapt_soil_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 0.45], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 0.45])
xlim([date1 date2])
xlabel('Year')
ylabel('CUE')
title('CUE - soil - no adaptation')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'southeast');
legend('box','off')
set(gca, 'FontSize', 14)

% Thermal adaptation - Control - litter - rStrat
subplot(3,2,2) 
hold on
p1 = plot(dates_all, thermalAdapt_CUE_litter_rStrat_control, 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, thermalAdapt_CUE_litter_rStrat_warming, 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgCUE_thermalAdapt_litter_rStrat_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgCUE_thermalAdapt_litter_rStrat_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 0.45], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 0.45])
xlim([date1 date2])
xlabel('Year')
ylabel('CUE')
title('CUE - organic horizon, r-strategists - thermal adaptation')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'southeast');
legend('box','off')
set(gca, 'FontSize', 14)

% Thermal adaptation - Control - litter - kStrat
subplot(3,2,4) 
hold on
p1 = plot(dates_all, thermalAdapt_CUE_litter_kStrat_control, 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, thermalAdapt_CUE_litter_kStrat_warming, 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgCUE_thermalAdapt_litter_kStrat_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgCUE_thermalAdapt_litter_kStrat_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 0.45], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 0.45])
xlim([date1 date2])
xlabel('Year')
ylabel('CUE')
title('CUE - forest floor, k-strategists - thermal adaptation')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'southeast');
legend('box','off')
set(gca, 'FontSize', 14)

% Thermal adaptation - Control - soil
subplot(3,2,6) 
hold on
p1 = plot(dates_all, thermalAdapt_CUE_soil_control, 'color', [54,144,192]./255, 'LineWidth', 1.5);
p2 = plot(dates_all, thermalAdapt_CUE_soil_warming, 'color', [239,101,72]./255, 'LineWidth', 1.5);
p3 = scatter(dates_averages, avgCUE_thermalAdapt_soil_control, 40, [2,56,88]./255, 'filled');
p4 = scatter(dates_averages, avgCUE_thermalAdapt_soil_warming, 40, [179,0,0]./255, 'filled');

% The start of warming is indicated
d1 = date_startWarming;
d2 = date_startWarming;
plot([d1 d2],[0 0.45], ':', 'color', 'k', 'lineWidth', 2);

ylim([0 0.45])
xlim([date1 date2])
xlabel('Year')
ylabel('CUE')
title('CUE - soil - thermal adaptation')
legend([p1, p2], 'Control', 'Warming', 'Orientation', 'Horizontal', 'FontSize', 14, 'Location', 'southeast');
legend('box','off')
set(gca, 'FontSize', 14)

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
% text(0.37,0.98,'Organic carbon stocks','fontweight','bold','fontsize',14)

% Letters
text(0.07, 0.945, '(A)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.51, 0.945, '(B)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.07, 0.64, '(C)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.51, 0.64, '(D)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.07, 0.34, '(E)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.51, 0.34, '(F)', 'FontSize', 14, 'FontWeight', 'bold')

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Changes in carbon use efficiency.tiff' -r300 -nocrop -transparent
% =========================================================================

% ---------------------------------------------------
% The relation between CUE and temperature is plotted
% ---------------------------------------------------

% The temperatures are loaded
soilTemperature_spinupAndControl = noAdaptData.soilTemperature_spinupAndControl;
soilTemperature_spinupAndWarming = noAdaptData.soilTemperature_spinupAndWarming;

% The data for the duration of the treatment is isolated
rowNum = find(dates_all == date_startWarming);
% soilTemperature_control = soilTemperature_spinupAndControl(rowNum:end);
% soilTemperature_warming = soilTemperature_spinupAndWarming(rowNum:end);

c_control = [54,144,192]./255;
c_warming = [239,101,72]./255;
alphaValue = .3;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 25*1.5 25*1.5], 'color', [1 1 1])

% No adaptation - litter - rStrat
subplot(3,2,1) 
hold on
p1 = scatter(soilTemperature_spinupAndControl(rowNum:end), noAdapt_CUE_litter_rStrat_control(rowNum:end), 20, c_control, 'filled');
p2 = scatter(soilTemperature_spinupAndWarming(rowNum:end), noAdapt_CUE_litter_rStrat_warming(rowNum:end), 20, c_warming, 'filled', 'MarkerFaceAlpha', alphaValue);
ylim([0.1 0.5])
xlim([-5 30])
xlabel('Temperature (°C)')
ylabel('CUE')
title('Organic horizon, r-strategists - control - no adaptation')
set(gca, 'FontSize', 14)
legend([p1, p2], 'Control', 'Warming' ,'FontSize', 14, 'Location', 'southwest')
legend('box', 'off')

% No adaptation - litter - kStrat
subplot(3,2,3) 
hold on
p1 = scatter(soilTemperature_spinupAndControl(rowNum:end), noAdapt_CUE_litter_kStrat_control(rowNum:end), 20, c_control, 'filled');
p2 = scatter(soilTemperature_spinupAndWarming(rowNum:end), noAdapt_CUE_litter_kStrat_warming(rowNum:end), 20, c_warming, 'filled', 'MarkerFaceAlpha', alphaValue);
ylim([0 0.5])
xlim([-5 30])
xlabel('Temperature (°C)')
ylabel('CUE')
title('Organic horizon, K-strategists - control - no adaptation')
set(gca, 'FontSize', 14)
legend([p1, p2], 'Control', 'Warming' ,'FontSize', 14, 'Location', 'southwest')
legend('box', 'off')

% No adaptation - soil
subplot(3,2,5) 
hold on
p1 = scatter(soilTemperature_spinupAndControl(rowNum:end), noAdapt_CUE_soil_control(rowNum:end), 20, c_control, 'filled');
p2 = scatter(soilTemperature_spinupAndWarming(rowNum:end), noAdapt_CUE_soil_warming(rowNum:end), 20, c_warming, 'filled', 'MarkerFaceAlpha', alphaValue);
ylim([0 0.5])
xlim([-5 30])
xlabel('Temperature (°C)')
ylabel('CUE')
title('Soil - control - no adaptation')
set(gca, 'FontSize', 14)
legend([p1, p2], 'Control', 'Warming' ,'FontSize', 14, 'Location', 'southwest')
legend('box', 'off')

% Thermal adaptation - r-strategists
subplot(3,2,2) 
hold on
p1 = scatter(soilTemperature_spinupAndControl(rowNum:end), thermalAdapt_CUE_litter_rStrat_control(rowNum:end), 20, c_control, 'filled');
p2 = scatter(soilTemperature_spinupAndWarming(rowNum:end), thermalAdapt_CUE_litter_rStrat_warming(rowNum:end), 20, c_warming, 'filled', 'MarkerFaceAlpha', alphaValue);
ylim([0.1 0.5])
xlim([-5 30])
xlabel('Temperature (°C)')
ylabel('CUE')
title('Organic horizon, r-strategists - warming - no adaptation')
set(gca, 'FontSize', 14)
legend([p1, p2], 'Control', 'Warming' ,'FontSize', 14, 'Location', 'southwest')
legend('box', 'off')

% Thermal adaptation - K-strategists
subplot(3,2,4) 
hold on
p1 = scatter(soilTemperature_spinupAndControl(rowNum:end), thermalAdapt_CUE_litter_kStrat_control(rowNum:end), 20, c_control, 'filled');
p2 = scatter(soilTemperature_spinupAndWarming(rowNum:end), thermalAdapt_CUE_litter_kStrat_warming(rowNum:end), 20, c_warming, 'filled', 'MarkerFaceAlpha', alphaValue);
ylim([0 0.5])
xlim([-5 30])
xlabel('Temperature (°C)')
ylabel('CUE')
title('Organic horizon, K-strategists - warming - no adaptation')
set(gca, 'FontSize', 14)
legend([p1, p2], 'Control', 'Warming' ,'FontSize', 14, 'Location', 'southwest')
legend('box', 'off')

% Thermal adaptation - soil
subplot(3,2,6) 
hold on
p1 = scatter(soilTemperature_spinupAndControl(rowNum:end), thermalAdapt_CUE_soil_control(rowNum:end), 20, c_control, 'filled');
p2 = scatter(soilTemperature_spinupAndWarming(rowNum:end), thermalAdapt_CUE_soil_warming(rowNum:end), 20, c_warming, 'filled', 'MarkerFaceAlpha', alphaValue);
ylim([0 0.5])
xlim([-5 30])
xlabel('Temperature (°C)')
ylabel('CUE')
title('Soil - warming - no adaptation')
set(gca, 'FontSize', 14)
legend([p1, p2], 'Control', 'Warming' ,'FontSize', 14, 'Location', 'southwest')
legend('box', 'off')

% =========================================================================
% A new invisible axis is created
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
xlim([0 1])
ylim([0 1])

% The title of the graph
% text(0.37,0.98,'Organic carbon stocks','fontweight','bold','fontsize',14)

% Letters
text(0.07, 0.945, '(A)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.51, 0.945, '(B)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.07, 0.64, '(C)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.51, 0.64, '(D)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.07, 0.34, '(E)', 'FontSize', 14, 'FontWeight', 'bold')
text(0.51, 0.34, '(F)', 'FontSize', 14, 'FontWeight', 'bold')

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'Relation temperature - CUE.tiff' -r300 -nocrop -transparent
% =========================================================================

%% The optimum temperature for CO2 respiration for every year is calculated

% ----------------------------------------------------------------
% For every year, the optimum temperature is calculated using MMRT
% ----------------------------------------------------------------

% The MMRT parameters are defined
kb = 1.380649e-23; % Boltzmann's constant (J K-1)
h = 6.62607015e-34; % Planck's constant (J s)
R = 8.314; % Universal gas constant (J K-1 mol-1)
dS = -0.105;

% Only data from the year 2000 onwards is retained
rowNum = find(noAdaptData.dates_all.Year >= 2000);
uniqueYears = unique(noAdaptData.dates_all(rowNum).Year);
nYears = numel(uniqueYears);

% The Topt values are obtained in a loop

% No adapt - control
T_tmp = noAdaptData.soilTemperature_spinupAndControl + 273;
CO2_tmp = noAdaptData.CO2_soil_control;
allDates = noAdaptData.dates_all;
Topt_allYears_noAdapt_control = get_Topt(nYears, uniqueYears, T_tmp, CO2_tmp, allDates, kb, h, R, dS);

% No adapt - warming
T_tmp = noAdaptData.soilTemperature_spinupAndWarming + 273;
CO2_tmp = noAdaptData.CO2_soil_warmed;
allDates = noAdaptData.dates_all;
Topt_allYears_noAdapt_warming = get_Topt(nYears, uniqueYears, T_tmp, CO2_tmp, allDates, kb, h, R, dS);

% Optimum driven - control
T_tmp = optimumDrivenData.soilTemperature_spinupAndControl + 273;
CO2_tmp = optimumDrivenData.CO2_soil_control;
allDates = optimumDrivenData.dates_all;
Topt_allYears_optDriv_control = get_Topt(nYears, uniqueYears, T_tmp, CO2_tmp, allDates, kb, h, R, dS);

% Optimum driven - warming
T_tmp = optimumDrivenData.soilTemperature_spinupAndWarming + 273;
CO2_tmp = optimumDrivenData.CO2_soil_warmed;
allDates = optimumDrivenData.dates_all;
Topt_allYears_optDriv_warming = get_Topt(nYears, uniqueYears, T_tmp, CO2_tmp, allDates, kb, h, R, dS);

% Enzyme rigidity - control
T_tmp = enzRigData.soilTemperature_spinupAndControl + 273;
CO2_tmp = enzRigData.CO2_soil_control;
allDates = enzRigData.dates_all;
Topt_allYears_enzRig_control = get_Topt(nYears, uniqueYears, T_tmp, CO2_tmp, allDates, kb, h, R, dS);

% Enzyme rigidity - warming
T_tmp = enzRigData.soilTemperature_spinupAndWarming + 273;
CO2_tmp = enzRigData.CO2_soil_warmed;
allDates = enzRigData.dates_all;
Topt_allYears_enzRig_warming = get_Topt(nYears, uniqueYears, T_tmp, CO2_tmp, allDates, kb, h, R, dS);

% ---------------------------
% The averages are calculated
% ---------------------------

avg_noAdapt_control = mean(Topt_allYears_noAdapt_control(3:end));
avg_noAdapt_warming = mean(Topt_allYears_noAdapt_warming(3:end));
changePerDegree_noAdapt = (avg_noAdapt_warming - avg_noAdapt_control)/5;
changePerDegree_noAdapt = round(changePerDegree_noAdapt*100)/100;

avg_optDriv_control = mean(Topt_allYears_optDriv_control(3:end));
avg_optDriv_warming = mean(Topt_allYears_optDriv_warming(3:end));
changePerDegree_optDriv = (avg_optDriv_warming - avg_optDriv_control)/5;
changePerDegree_optDriv = round(changePerDegree_optDriv*100)/100;

avg_enzRig_control = mean(Topt_allYears_enzRig_control(3:end));
avg_enzRig_warming = mean(Topt_allYears_enzRig_warming(3:end));
changePerDegree_enzRig = (avg_enzRig_warming - avg_enzRig_control)/5;
changePerDegree_enzRig = round(changePerDegree_enzRig*100)/100;

% --------
% Plotting
% --------

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 25*1.5 25*1.5], 'color', [1 1 1])

subplot(3,1,1)
hold on
scatter(uniqueYears,Topt_allYears_noAdapt_control, 'filled');
scatter(uniqueYears,Topt_allYears_noAdapt_warming, 'filled');
ylim([285 305])
title('No thermal adaptation')
text(2000.5, 302, ['Change per deg. K = +' num2str(changePerDegree_noAdapt)], 'FontSize', 16)

subplot(3,1,2)
hold on
scatter(uniqueYears,Topt_allYears_optDriv_control, 'filled');
scatter(uniqueYears,Topt_allYears_optDriv_warming, 'filled');
ylim([285 305])
title('Optimum driven')
text(2000.5, 302, ['Change per deg. K = +' num2str(changePerDegree_optDriv)], 'FontSize', 16)

subplot(3,1,3)
hold on
scatter(uniqueYears,Topt_allYears_enzRig_control, 'filled');
scatter(uniqueYears,Topt_allYears_enzRig_warming, 'filled');
ylim([285 305])
title('Enzyme rigidity')
text(2000.5, 302, ['Change per deg. K = +' num2str(changePerDegree_enzRig)], 'FontSize', 16)



figure






