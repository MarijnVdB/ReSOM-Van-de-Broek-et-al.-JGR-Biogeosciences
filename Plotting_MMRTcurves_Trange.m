%% This script is to plot the MMRT curves and temperature ranges

%% MMRT curve for the no adaptation scenario

clc; close all

% !!! First run the main script untill the script 'Initialize_temperature_moisture' is run

tmp = load('GA output/Barre Woods - noTau - 5 June 2024/xga_noAdapt_16-May-2024.mat');
xga = tmp.xga;

% Colors are defined
c_curve = [0,0,0]./255;
c_range1 = [116,169,207]./255;
c_range2 = [239,101,72]./255;

% The year for this the curves are plotted are defined
yrRange1 = 1998:2000;
yrRange2 = 2010:2012;

% The min and max temperature for these years is calculated, and all
% temperatures are isolated
[r c] = find (dates_all.Year >= yrRange1(1) & dates_all.Year <= yrRange1(end));
allTemp_range1 = soilTemperature_spinupAndWarming(c);
[r c] = find (dates_all.Year >= yrRange2(1) & dates_all.Year <= yrRange2(end));
allTemp_range2 = soilTemperature_spinupAndWarming(c);

% The MMRT curves are calculated
par_MMRT = [xga(15);...   % Cp
            NaN; ...  % dH
            NaN];     % dS

Topt = xga(16);
Tref_MMRT = 299.9;

R = 8.314;

par_MMRT(2) = (Topt*(-par_MMRT(1)*1e3-R) + par_MMRT(1)*1e3*Tref_MMRT)/1e3; 
par_MMRT(3) = (0.0033 * (par_MMRT(2)*1e3) - 204.01) / 1e3;

VmaxRef = 1;
temp = 273.15 + (-10:1:50);

% The temperature-dependent parameters are calculated
Topt_yr1 = Topt;
Tref_yr1 = Topt_yr1 - 6;
par_MMRT_yr1 = par_MMRT;
par_MMRT_yr1(2) = (Topt_yr1*(-par_MMRT_yr1(1)*1e3-R) + par_MMRT_yr1(1)*1e3*Tref_yr1)/1e3;
% par_MMRT_yr1(2) = par_MMRT_yr1(2);

Topt_yr2 = Topt;
Tref_yr2 = Topt_yr2 - 6;
par_MMRT_yr2 = par_MMRT;
par_MMRT_yr2(2) = (Topt_yr1*(-par_MMRT_yr2(1)*1e3-R) + par_MMRT_yr2(1)*1e3*Tref_yr1)/1e3;
% par_MMRT_yr2(2) = par_MMRT_yr2(2);

% The MMRT factors are calculated
MMRTfactor_yr1 = MMRT(par_MMRT_yr1, Tref_yr1, temp, Topt);
MMRTfactor_yr2 = MMRT(par_MMRT_yr2, Tref_yr2, temp, Topt);

Vmax_yr1 = VmaxRef.*MMRTfactor_yr1;
Vmax_yr2 = VmaxRef.*MMRTfactor_yr2;

% The MMRT rate modifier is calculatd for all temperatures in the ranges
MMRT_rateMod_range1 = MMRT(par_MMRT_yr1, Tref_yr1, 273.15 + allTemp_range1, Topt);
MMRT_rateMod_range2 = MMRT(par_MMRT_yr2, Tref_yr2, 273.15 + allTemp_range2, Topt);

% These data are formatted for plotting
yrCell = cell(numel(MMRT_rateMod_range1) + numel(MMRT_rateMod_range2),1);
yrCell(1:numel(MMRT_rateMod_range1)) = {'1998:2000'};
yrCell(numel(MMRT_rateMod_range1) + 1 : end) = {'2010:2012'};

% --------------------------------------
% Plot all temperatures versus all rates
% --------------------------------------

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [15 15 30 30], 'color', [1 1 1])

lineWidth = 2;

% The temperatures and MMRT values
s1 = scatterhist([allTemp_range1; allTemp_range2], [MMRT_rateMod_range1; MMRT_rateMod_range2], 'Group', yrCell, 'Kernel', 'off','MarkerSize', .1,...
    'color', [c_range1;c_range2],'LineStyle',{'-','-'}, 'Location', 'northeast', 'Direction','out', 'NBins', 25, 'Style', 'bar');

set(gca, 'XAxisLocation', 'bottom', 'YAxisLocation', 'left')

hold(s1(1), 'on')

% The MMRT curves
p1 = plot(temp-273.15, Vmax_yr1, 'color', c_curve, 'linewidth', lineWidth);
s2 = scatter(allTemp_range1, MMRT_rateMod_range1, 400, 'markerfacecolor', c_range1, 'markeredgecolor', c_range1);
s3 = scatter(allTemp_range2, MMRT_rateMod_range2, 60, 'markerfacecolor', c_range2, 'markeredgecolor', c_range2);

set(gca, 'fontsize',30, 'ytick', 0:0.2:1)
xlabel('Temperature (°C)');
ylabel('f_{MMRT}');
ylim([0 1])

lgd = legend([p1, s2, s3], 'MMRT curve', '1998 - 2000', '2010 - 2012');
set(legend,'FontSize',26);
legend('boxoff')
pos = lgd.Position;
pos(1) = 0.15; % horizontal
pos(2) = 0.61; % vertical
set(lgd, 'Position', pos)

text(-3,0.95, 'No thermal adaptation', 'FontWeight', 'bold', 'FontSize', 30, 'Parent', s1(1))
text(-7.5,1.28,'(A)', 'FontWeight', 'bold', 'FontSize', 24)

% Position of the main graph
s1(1).Position(1) = 0.13; % horizontal
s1(1).Position(2) = 0.12; % vertical
s1(1).Position(3) = .67; % Width
s1(1).Position(4) = .67; % height

% Position of the upper graph
drawnow

s1(2).Position(2) = 0.80; % vertical
s1(2).Position(4) = 0.19; % height

% Position of the right graph
s1(3).Position(1) = 0.82; % horizontal
s1(3).Position(3) = 0.17; % width

% =========================================================================
% Export the figure
% drawnow
% addpath('export_fig')
% export_fig 'MMRT_A.tiff' -r300 -nocrop -transparent
% =========================================================================

%% MMRT curve for the optimum driven scenario

% !!! First run the main script untill the script 'Initialize_temperature_moisture' is run

tmp = load('GA output/Barre Woods - noTau - 5 June 2024/xga_optDriv_04-May-2024.mat');
xga = tmp.xga;

% Colors are defined
c_range1 = [116,169,207]./255;
c_range2 = [239,101,72]./255;

% The year for this the curves are plotted are defined
yrRange1 = 1998:2000;
yrRange2 = 2010:2012;

% The mean annual, min and max temperature for these years is calculated
[r c] = find (dates_all.Year >= yrRange1(1) & dates_all.Year <= yrRange1(end));
allTemp_range1 = soilTemperature_spinupAndWarming(c);
meanTemp_yr1 = mean(soilTemperature_spinupAndWarming(c));
[r c] = find (dates_all.Year >= yrRange2(1) & dates_all.Year <= yrRange2(end));
allTemp_range2 = soilTemperature_spinupAndWarming(c);
meanTmp_yr2 = mean(soilTemperature_spinupAndWarming(c));

% The MMRT curves are calculated
par_MMRT = [NaN;...   % Cp
            NaN; ...  % dH
            NaN];     % dS

alphaT = xga(16);
betaT = xga(17);
nYears_thermalAdapt = xga(18);

R = 8.314;
VmaxRef = 1;
temp = 273.15 + (-10:1:50);

% The temperature-dependent parameters are calculated
Topt_yr1 = alphaT*meanTemp_yr1 + betaT + 273.15;
Tref_yr1 = 299.9;
% Topt_yr1 = Tref_yr1;
par_MMRT_yr1 = par_MMRT;
par_MMRT_yr1(1) = xga(15);
par_MMRT_yr1(2) = (Topt_yr1*(-par_MMRT_yr1(1)*1e3-R) + par_MMRT_yr1(1)*1e3*Tref_yr1)/1e3; % dH is calculated based on Topt (which is optimized)
par_MMRT_yr1(3) = (0.0033 * (par_MMRT_yr1(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9))

Topt_yr2 = alphaT*meanTmp_yr2 + betaT + 273.15;
Tref_yr2 = 299.9;
% Topt_yr2 = Tref_yr2;
par_MMRT_yr2 = par_MMRT;
par_MMRT_yr2(1) = xga(15);
par_MMRT_yr2(2) = (Topt_yr2*(-par_MMRT_yr2(1)*1e3-R) + par_MMRT_yr2(1)*1e3*Tref_yr2)/1e3; % dH is calculated based on Topt (which is optimized)
par_MMRT_yr2(3) = (0.0033 * (par_MMRT_yr2(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9))

% The MMRT factors are calculated
MMRTfactor_yr1 = MMRT(par_MMRT_yr1, Tref_yr1, temp, Topt_yr1);
MMRTfactor_yr2 = MMRT(par_MMRT_yr2, Tref_yr2, temp, Topt_yr2);

Vmax_yr1 = VmaxRef.*MMRTfactor_yr1;
Vmax_yr2 = VmaxRef.*MMRTfactor_yr2;

% The MMRT rate modifier is calculatd for all temperatures in the ranges
MMRT_rateMod_range1 = MMRT(par_MMRT_yr1, Tref_yr1, 273.15 + allTemp_range1, Topt_yr1);
MMRT_rateMod_range2 = MMRT(par_MMRT_yr2, Tref_yr2, 273.15 + allTemp_range2, Topt_yr2);

% These data are formatted for plotting
yrCell = cell(numel(MMRT_rateMod_range1) + numel(MMRT_rateMod_range2),1);
yrCell(1:numel(MMRT_rateMod_range1)) = {'1998:2000'};
yrCell(numel(MMRT_rateMod_range1) + 1 : end) = {'2010:2012'};

% --------------------------------------
% Plot all temperatures versus all rates
% --------------------------------------

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [5 5 30 30], 'color', [1 1 1])

% The temperatures and MMRT values
s1 = scatterhist([allTemp_range1; allTemp_range2], [MMRT_rateMod_range1; MMRT_rateMod_range2], 'Group', yrCell, 'Kernel', 'off','MarkerSize', .1,...
    'color', [c_range1;c_range2],'LineStyle',{'-','-'}, 'Location', 'northeast', 'Direction','out', 'NBins', 25, 'Style', 'bar');

set(gca, 'XAxisLocation', 'bottom', 'YAxisLocation', 'left')
hold(s1(1), 'on')

lineWidth = 2;

% The MMRT curves
p1 = plot(temp-273.15, Vmax_yr1, 'color', c_range1, 'linewidth', lineWidth);
p2 = plot(temp-273.15, Vmax_yr2, 'color', c_range2, 'linewidth', lineWidth);
% s2 = scatter(allTemp_range1, MMRT_rateMod_range1, 200, 'markerfacecolor', c_range1, 'markeredgecolor', c_range1);
% s3 = scatter(allTemp_range2, MMRT_rateMod_range2, 200, 'markerfacecolor', c_range2, 'markeredgecolor', c_range2);

s2 = scatter(allTemp_range1, MMRT_rateMod_range1, 200, 'markeredgecolor', c_range1);
s3 = scatter(allTemp_range2, MMRT_rateMod_range2, 200, 'markeredgecolor', c_range2);

set(gca, 'fontsize',30, 'ytick', 0:0.2:1)
xlabel('Temperature (°C)');
ylabel('f_{MMRT}');
ylim([0 1])

lgd = legend([p1, p2, s2, s3], 'MMRT curve 1998 - 2000', 'MMRT curve 2010 - 2012', '1998 - 2000', '2010 - 2012');
set(legend,'FontSize',26);
legend('boxoff')
pos = lgd.Position;
pos(1) = 0.15; % horizontal
pos(2) = 0.58; % vertical
set(lgd, 'Position', pos)

text(-3,0.95, 'Optimum driven', 'FontWeight', 'bold', 'FontSize', 30, 'Parent', s1(1))
text(-7.5,1.28,'(B)', 'FontWeight', 'bold', 'FontSize', 24)

% Position of the main graph
s1(1).Position(1) = 0.13; % horizontal
s1(1).Position(2) = 0.12; % vertical
s1(1).Position(3) = .67; % Width
s1(1).Position(4) = .67; % height

% Position of the upper graph
drawnow

s1(2).Position(2) = 0.80; % vertical
s1(2).Position(4) = 0.19; % height

% Position of the right graph
s1(3).Position(1) = 0.82; % horizontal
s1(3).Position(3) = 0.17; % width

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'MMRT_B.tiff' -r300 -nocrop -transparent
% =========================================================================

%% MMRT curve for the enzyme rigidity scenario

% !!! First run the main script untill the script 'Initialize_temperature_moisture' is run

clc; close all

tmp = load('GA output/Barre Woods - noTau - 5 June 2024/xga_enzRig_07-May-2024.mat');
xga = tmp.xga;

% Colors are defined
c_range1 = [116,169,207]./255;
c_range2 = [239,101,72]./255;

% The year for this the curves are plotted are defined
yrRange1 = 1998:2000;
yrRange2 = 2010:2012;

% The mean annual, min and max temperature for these years is calculated
[r c] = find (dates_all.Year >= yrRange1(1) & dates_all.Year <= yrRange1(end));
allTemp_range1 = soilTemperature_spinupAndWarming(c);
meanTemp_yr1 = mean(soilTemperature_spinupAndWarming(c));
[r c] = find (dates_all.Year >= yrRange2(1) & dates_all.Year <= yrRange2(end));
allTemp_yr2 = soilTemperature_spinupAndWarming(c);
meanTemp_yr2 = mean(soilTemperature_spinupAndWarming(c));

% The MMRT curves are calculated
par_MMRT = [NaN;...   % Cp
            NaN; ...  % dH
            NaN];     % dS
        
alphaT = xga(15);
betaT = xga(16);
alphaC = xga(17);
betaC = xga(18);

R = 8.314;
VmaxRef = 1;
temp = 273.15 + (-10:1:50);

% The temperature-dependent parameters are calculated
Topt_yr1 = alphaT*meanTemp_yr1 + betaT + 273.15;
Tref_yr1 = 299.9;
par_MMRT_yr1 = par_MMRT;
par_MMRT_yr1(1) = alphaC*meanTemp_yr1 + betaC;
par_MMRT_yr1(2) = (Topt_yr1*(-par_MMRT_yr1(1)*1e3-R) + par_MMRT_yr1(1)*1e3*Tref_yr1)/1e3; % dH is calculated based on Topt (which is optimized)
par_MMRT_yr1(3) = (0.0033 * (par_MMRT_yr1(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9))

Topt_yr2 = alphaT*meanTemp_yr2 + betaT + 273.15;
Tref_yr2 = 299.9;
par_MMRT_yr2 = par_MMRT;
par_MMRT_yr2(1) = alphaC*meanTemp_yr2 + betaC;
par_MMRT_yr2(2) = (Topt_yr2*(-par_MMRT_yr2(1)*1e3-R) + par_MMRT_yr2(1)*1e3*Tref_yr2)/1e3; % dH is calculated based on Topt (which is optimized)
par_MMRT_yr2(3) = (0.0033 * (par_MMRT_yr2(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9))

% The MMRT factors are calculated
MMRTfactor_yr1 = MMRT(par_MMRT_yr1, Tref_yr1, temp, Topt_yr1);
MMRTfactor_yr2 = MMRT(par_MMRT_yr2, Tref_yr2, temp, Topt_yr2);

Vmax_yr1 = VmaxRef.*MMRTfactor_yr1;
Vmax_yr2 = VmaxRef.*MMRTfactor_yr2;

% The MMRT rate modifier is calculatd for all temperatures in the ranges
MMRT_rateMod_range1 = MMRT(par_MMRT_yr1, Tref_yr1, 273.15 + allTemp_range1, Topt_yr1);
MMRT_rateMod_range2 = MMRT(par_MMRT_yr2, Tref_yr2, 273.15 + allTemp_range2, Topt_yr2);

% These data are formatted for plotting
yrCell = cell(numel(MMRT_rateMod_range1) + numel(MMRT_rateMod_range2),1);
yrCell(1:numel(MMRT_rateMod_range1)) = {'1998:2000'};
yrCell(numel(MMRT_rateMod_range1) + 1 : end) = {'2010:2012'};

% --------------------------------------
% Plot all temperatures versus all rates
% --------------------------------------

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [5 5 30 30], 'color', [1 1 1])

% The temperatures and MMRT values
s1 = scatterhist([allTemp_range1; allTemp_range2], [MMRT_rateMod_range1; MMRT_rateMod_range2], 'Group', yrCell, 'Kernel', 'off','MarkerSize', .1,...
    'color', [c_range1;c_range2],'LineStyle',{'-','-'}, 'Location', 'northeast', 'Direction','out', 'NBins', 25, 'Style', 'bar');

set(gca, 'XAxisLocation', 'bottom', 'YAxisLocation', 'left')

hold(s1(1), 'on')

% The MMRT curves
p1 = plot(temp-273.15, Vmax_yr1, 'color', c_range1, 'linewidth', lineWidth);
p2 = plot(temp-273.15, Vmax_yr2, 'color', c_range2, 'linewidth', lineWidth);
s2 = scatter(allTemp_range1, MMRT_rateMod_range1, 200, 'markerfacecolor', c_range1, 'markeredgecolor', c_range1);
s3 = scatter(allTemp_range2, MMRT_rateMod_range2, 200, 'markerfacecolor', c_range2, 'markeredgecolor', c_range2);

set(gca, 'fontsize',30, 'ytick', 0:0.2:1)
xlabel('Temperature (°C)');
ylabel('f_{MMRT}');
ylim([0 1])

lgd = legend([p1, p2, s2, s3], 'MMRT curve 1998 - 2000', 'MMRT curve 2010 - 2012', '1998 - 2000', '2010 - 2012');
set(legend,'FontSize',26);
legend('boxoff')
pos = lgd.Position;
pos(1) = 0.15; % horizontal
pos(2) = 0.58; % vertical
set(lgd, 'Position', pos)

text(-3,0.95, 'Enzyme rigidity', 'FontWeight', 'bold', 'FontSize', 30, 'Parent', s1(1))
text(-7.5,1.28,'(C)', 'FontWeight', 'bold', 'FontSize', 24)

% Position of the main graph
s1(1).Position(1) = 0.13; % horizontal
s1(1).Position(2) = 0.12; % vertical
s1(1).Position(3) = .67; % Width
s1(1).Position(4) = .67; % height

% Position of the upper graph
drawnow

s1(2).Position(2) = 0.80; % vertical
s1(2).Position(4) = 0.19; % height

% Position of the right graph
s1(3).Position(1) = 0.82; % horizontal
s1(3).Position(3) = 0.17; % width

% =========================================================================
% Export the figure
% addpath('export_fig')
% export_fig 'MMRT_C.tiff' -r300 -nocrop -transparent
% =========================================================================






