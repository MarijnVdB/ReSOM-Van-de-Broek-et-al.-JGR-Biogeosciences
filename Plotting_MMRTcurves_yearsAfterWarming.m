%% This script is to plot the MMRT curves and temperature ranges for the first years after warming

clc; close all

%% MMRT curves for the no adaptation scenario

tmp = load('GA output/Barre Woods - noTau - 5 June 2024/xga_noAdapt_16-May-2024.mat');
xga = tmp.xga;

% The colors
c_yr1 = [253,187,132]./255;
c_yr2 = [252,141,89]./255;
c_yr3 = [239,101,72]./255;
c_yr4 = [215,48,31]./255;
c_yr5 = [179,0,0]./255;
c_curve = 'k';

% The year for this the curves are plotted are defined
dates1 = datetime(2002,05,21):datetime(2003,05,20);
dates2 = datetime(2003,05,21):datetime(2004,05,20);
dates3 = datetime(2004,05,21):datetime(2005,05,20);
dates4 = datetime(2005,05,21):datetime(2006,05,20);
dates5 = datetime(2006,05,21):datetime(2007,05,20);

% All temperatures of the evaluated years are isolated
[r c] = ismember(dates1,dates_all);
allTemp_range1 = soilTemperature_spinupAndWarming(c);

[r c] = ismember(dates2,dates_all);
allTemp_range2 = soilTemperature_spinupAndWarming(c);

[r c] = ismember(dates3,dates_all);
allTemp_range3 = soilTemperature_spinupAndWarming(c);

[r c] = ismember(dates4,dates_all);
allTemp_range4 = soilTemperature_spinupAndWarming(c);

[r c] = ismember(dates5,dates_all);
allTemp_range5 = soilTemperature_spinupAndWarming(c);

% The mean temperature for every based on which the MMRT curve is
% calculated is obtained, to calculate Tref_MMRT
% (not necessary for the no adaptation scenario)

% The MMRT curves are calculated
par_MMRT = [xga(15);...   % Cp
            NaN; ...  % dH
            NaN];     % dS

Topt = xga(16);
Tref_MMRT = 299.9;

R = 8.314;
VmaxRef = 1;
temp = 273.15 + (-10:1:50);

% The temperature-dependent parameters are calculated
Tref_yr1 = Tref_MMRT;
Tref_yr2 = Tref_MMRT;
Tref_yr3 = Tref_MMRT;
Tref_yr4 = Tref_MMRT;
Tref_yr5 = Tref_MMRT;

% Tref = Tref_MMRT;
Topt_yr1 = Topt;
Topt_yr2 = Topt;
Topt_yr3 = Topt;
Topt_yr4 = Topt;
Topt_yr5 = Topt;

par_MMRT_yr1 = par_MMRT;
par_MMRT_yr1(2) = (Topt_yr1*(-par_MMRT_yr1(1)*1e3-R) + par_MMRT_yr1(1)*1e3*Tref_MMRT)/1e3; 
par_MMRT_yr1(3) = (0.0033 * (par_MMRT_yr1(2)*1e3) - 204.01) / 1e3;

par_MMRT_yr2 = par_MMRT;
par_MMRT_yr2(2) = (Topt_yr2*(-par_MMRT_yr2(1)*1e3-R) + par_MMRT_yr2(1)*1e3*Tref_MMRT)/1e3; 
par_MMRT_yr2(3) = (0.0033 * (par_MMRT_yr2(2)*1e3) - 204.01) / 1e3;

par_MMRT_yr3 = par_MMRT;
par_MMRT_yr3(2) = (Topt_yr3*(-par_MMRT_yr3(1)*1e3-R) + par_MMRT_yr3(1)*1e3*Tref_MMRT)/1e3; 
par_MMRT_yr3(3) = (0.0033 * (par_MMRT_yr3(2)*1e3) - 204.01) / 1e3;

par_MMRT_yr4 = par_MMRT;
par_MMRT_yr4(2) = (Topt_yr4*(-par_MMRT_yr4(1)*1e3-R) + par_MMRT_yr4(1)*1e3*Tref_MMRT)/1e3; 
par_MMRT_yr4(3) = (0.0033 * (par_MMRT_yr4(2)*1e3) - 204.01) / 1e3;

par_MMRT_yr5 = par_MMRT;
par_MMRT_yr5(2) = (Topt_yr5*(-par_MMRT_yr5(1)*1e3-R) + par_MMRT_yr5(1)*1e3*Tref_MMRT)/1e3; 
par_MMRT_yr5(3) = (0.0033 * (par_MMRT_yr5(2)*1e3) - 204.01) / 1e3;

% The MMRT factors are calculated
MMRTfactor_yr1 = MMRT(par_MMRT_yr1, Tref_yr1, temp, Topt_yr1);
MMRTfactor_yr2 = MMRT(par_MMRT_yr2, Tref_yr2, temp, Topt_yr2);
MMRTfactor_yr3 = MMRT(par_MMRT_yr3, Tref_yr3, temp, Topt_yr3);
MMRTfactor_yr4 = MMRT(par_MMRT_yr4, Tref_yr4, temp, Topt_yr4);
MMRTfactor_yr5 = MMRT(par_MMRT_yr5, Tref_yr5, temp, Topt_yr5);

Vmax_yr1 = VmaxRef.*MMRTfactor_yr1;
Vmax_yr2 = VmaxRef.*MMRTfactor_yr2;
Vmax_yr3 = VmaxRef.*MMRTfactor_yr3;
Vmax_yr4 = VmaxRef.*MMRTfactor_yr4;
Vmax_yr5 = VmaxRef.*MMRTfactor_yr5;

% The MMRT rate modifier is calculatd for all temperatures in the ranges
MMRT_rateMod_range1 = MMRT(par_MMRT_yr1, Tref_yr1, 273.15 + allTemp_range1, Topt_yr1);
MMRT_rateMod_range2 = MMRT(par_MMRT_yr2, Tref_yr2, 273.15 + allTemp_range2, Topt_yr2);
MMRT_rateMod_range3 = MMRT(par_MMRT_yr3, Tref_yr3, 273.15 + allTemp_range3, Topt_yr3);
MMRT_rateMod_range4 = MMRT(par_MMRT_yr4, Tref_yr4, 273.15 + allTemp_range4, Topt_yr4);
MMRT_rateMod_range5 = MMRT(par_MMRT_yr5, Tref_yr5, 273.15 + allTemp_range5, Topt_yr5);

% These data are formatted for plotting
yrCell = cell(numel(dates1) + numel(dates2) + numel(dates3) + numel(dates4) + numel(dates5),1);
yrCell(1:numel(dates1)) = {'Year 0'};
yrCell(numel(dates1) + 1 : numel(dates1) + numel(dates2)) = {'Year 1'};
yrCell(numel(dates1) + numel(dates2) + 1 : numel(dates1) + numel(dates2) + numel(dates3)) = {'Year 2'};
yrCell(numel(dates1) + numel(dates2) + numel(dates3) + 1 : numel(dates1) + numel(dates2) + numel(dates3) + numel(dates4)) = {'Year 3'};
yrCell(numel(dates1) + numel(dates2) + numel(dates3) + numel(dates4) + 1 : numel(dates1) + numel(dates2) + numel(dates3) + numel(dates4) + numel(dates5)) = {'Year 4'};

% --------------------------------------
% Plot all temperatures versus all rates
% --------------------------------------

lineWidth = 2;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [5 5 30 30], 'color', [1 1 1])

bw1 = .5;
bw2 = .02;

% The temperatures and MMRT values
s1 = scatterhist([allTemp_range1; allTemp_range2; allTemp_range3; allTemp_range4; allTemp_range5], ...
                        [MMRT_rateMod_range1; MMRT_rateMod_range2; MMRT_rateMod_range3; MMRT_rateMod_range4; MMRT_rateMod_range5], 'Group', yrCell, 'Kernel', 'on','MarkerSize', .5,...
    'color', [c_yr1;c_yr2; c_yr3; c_yr4; c_yr5],'LineStyle',{'-','-','-','-','-'}, 'Location', 'northeast', 'Direction','out', 'Bandwidth', [bw1, bw1, bw1, bw1, bw1; bw2, bw2, bw2, bw2,bw2]);

set(gca, 'XAxisLocation', 'bottom', 'YAxisLocation', 'left')

hold(s1(1), 'on')

% The MMRT curves
p1 = plot(temp-273.15, Vmax_yr1, 'color', c_curve, 'linewidth', lineWidth);
s2 = scatter(allTemp_range1, MMRT_rateMod_range1, 30, 'markerfacecolor', c_yr1, 'markeredgecolor', c_yr1);
s3 = scatter(allTemp_range2, MMRT_rateMod_range2, 30, 'markerfacecolor', c_yr2, 'markeredgecolor', c_yr2);
s4 = scatter(allTemp_range3, MMRT_rateMod_range3, 30, 'markerfacecolor', c_yr3, 'markeredgecolor', c_yr3);
s5 = scatter(allTemp_range4, MMRT_rateMod_range4, 30, 'markerfacecolor', c_yr4, 'markeredgecolor', c_yr4);
s6 = scatter(allTemp_range5, MMRT_rateMod_range5, 30, 'markerfacecolor', c_yr5, 'markeredgecolor', c_yr5);

set(gca, 'fontsize', 24, 'YTick',[0:0.2:1])
xlabel('Temperature (°C)');
ylabel('MMRT rate modifier');
ylim([0 1])

lgd = legend([p1, s2, s3, s4, s5, s6], 'MMRT curve', 'Last year before soil warming', '1^{st} year of soil warming', '2^{nd} year of soil warming', '3^{th} year of soil warming', '4^{th} year of soil warming');
set(legend,'FontSize',24);
legend('boxoff')
pos = lgd.Position;
pos(1) = 0.1; % horizontal
pos(2) = 0.46; % vertical
set(lgd, 'Position', pos)

text(-3,0.95, 'No thermal adaptation', 'FontWeight', 'bold', 'FontSize', 28, 'Parent', s1(1))
text(-7.5,1.28,'(A)', 'FontWeight', 'bold', 'FontSize', 28)

% Position of the main graph
s1(1).Position(1) = 0.1; % horizontal
s1(1).Position(2) = 0.09; % vertical
s1(1).Position(3) = .67; % Width
s1(1).Position(4) = .67; % height

% Position of the upper graph
drawnow

s1(2).Position(2) = 0.78; % vertical
s1(2).Position(4) = 0.2; % height

% Position of the right graph
s1(3).Position(1) = 0.8; % horizontal
s1(3).Position(3) = 0.18; % width

%% MMRT curve for the optimum driven scenario

tmp = load('GA output/Barre Woods - noTau - 5 June 2024/xga_optDriv_04-May-2024.mat');
xga = tmp.xga;

% The colors
c_yr1 = [253,187,132]./255;
c_yr2 = [252,141,89]./255;
c_yr3 = [239,101,72]./255;
c_yr4 = [215,48,31]./255;
c_yr5 = [179,0,0]./255;
c_curve = 'k';

% The year for this the curves are plotted are defined
dates1 = datetime(2002,05,21):datetime(2003,05,20);
dates2 = datetime(2003,05,21):datetime(2004,05,20);
dates3 = datetime(2004,05,21):datetime(2005,05,20);
dates4 = datetime(2005,05,21):datetime(2006,05,20);
dates5 = datetime(2006,05,21):datetime(2007,05,20);

% All temperatures of the evaluated years are isolated
[r c] = ismember(dates1,dates_all);
allTemp_range1 = soilTemperature_spinupAndWarming(c);

[r c] = ismember(dates2,dates_all);
allTemp_range2 = soilTemperature_spinupAndWarming(c);

[r c] = ismember(dates3,dates_all);
allTemp_range3 = soilTemperature_spinupAndWarming(c);

[r c] = ismember(dates4,dates_all);
allTemp_range4 = soilTemperature_spinupAndWarming(c);

[r c] = ismember(dates5,dates_all);
allTemp_range5 = soilTemperature_spinupAndWarming(c);

% The mean temperature for every based on which the MMRT curve is
% calculated is obtained, to calculate Tref_MMRT

nYears = xga(18);
% The number of time steps equal to nYears
nTimeSteps = ceil((365*nYears));

% The average temperature from every day to that day minus nYears is calculated
allAvgTemp_yr1 = NaN(numel(dates1),1);
for ii = 1:numel(dates1)
    % The index of the current day in temperature_spinupAndControl
    [r, ind] = find(dates_all == dates1(ii));
    % The average temperature is stored
    allAvgTemp_yr1(ii,1) = mean(soilTemperature_spinupAndWarming(ind-nTimeSteps:ind));
end
avgTempForMMRT_yr1 = mean(allAvgTemp_yr1); 

allAvgTemp_yr2 = NaN(numel(dates2),1);
for ii = 1:numel(dates2)
    % The index of the current day in temperature_spinupAndControl
    [r, ind] = find(dates_all == dates2(ii));
    % The average temperature is stored
    allAvgTemp_yr2(ii,1) = mean(soilTemperature_spinupAndWarming(ind-nTimeSteps:ind));
end
avgTempForMMRT_yr2 = mean(allAvgTemp_yr2); 

allAvgTemp_yr3 = NaN(numel(dates3),1);
for ii = 1:numel(dates3) 
    % The index of the current day in temperature_spinupAndControl
    [r, ind] = find(dates_all == dates3(ii));
    % The average temperature is stored
    allAvgTemp_yr3(ii,1) = mean(soilTemperature_spinupAndWarming(ind-nTimeSteps:ind));
end
avgTempForMMRT_yr3 = mean(allAvgTemp_yr3); 

allAvgTemp_yr4 = NaN(numel(dates4),1);
for ii = 1:numel(dates4) 
    % The index of the current day in temperature_spinupAndControl
    [r, ind] = find(dates_all == dates4(ii));
    % The average temperature is stored
    allAvgTemp_yr4(ii,1) = mean(soilTemperature_spinupAndWarming(ind-nTimeSteps:ind));
end
avgTempForMMRT_yr4 = mean(allAvgTemp_yr4);

allAvgTemp_yr5 = NaN(numel(dates5),1);
for ii = 1:numel(dates5) 
    % The index of the current day in temperature_spinupAndControl
    [r, ind] = find(dates_all == dates5(ii));
    % The average temperature is stored
    allAvgTemp_yr5(ii,1) = mean(soilTemperature_spinupAndWarming(ind-nTimeSteps:ind));
end
avgTempForMMRT_yr5 = mean(allAvgTemp_yr5);

% The MMRT curves are calculated
par_MMRT = [NaN;...   % Cp
            NaN; ...  % dH
            NaN];     % dS

% The calibrated parameter values
alphaT = xga(16);
betaT = xga(17);

R = 8.314;
VmaxRef = 1;
temp = 273.15 + (-10:1:50);

% The temperature-dependent parameters are calculated
Topt_yr1 = alphaT*avgTempForMMRT_yr1 + betaT + 273.15;
Tref_yr1 = 299.9;
par_MMRT_yr1 = par_MMRT;
par_MMRT_yr1(1) = xga(15);
par_MMRT_yr1(2) = (Topt_yr1*(-par_MMRT_yr1(1)*1e3-R) + par_MMRT_yr1(1)*1e3*Tref_yr1)/1e3; % dH is calculated based on Topt (which is optimized)
par_MMRT_yr1(3) = (0.0033 * (par_MMRT_yr1(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9))

Topt_yr2 = alphaT*avgTempForMMRT_yr2 + betaT + 273.15;
Tref_yr2 = 299.9;
par_MMRT_yr2 = par_MMRT;
par_MMRT_yr2(1) = xga(15);
par_MMRT_yr2(2) = (Topt_yr2*(-par_MMRT_yr2(1)*1e3-R) + par_MMRT_yr2(1)*1e3*Tref_yr2)/1e3; % dH is calculated based on Topt (which is optimized)
par_MMRT_yr2(3) = (0.0033 * (par_MMRT_yr2(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9))

Topt_yr3 = alphaT*avgTempForMMRT_yr3 + betaT + 273.15;
Tref_yr3 = 299.9;
par_MMRT_yr3 = par_MMRT;
par_MMRT_yr3(1) = xga(15);
par_MMRT_yr3(2) = (Topt_yr3*(-par_MMRT_yr3(1)*1e3-R) + par_MMRT_yr3(1)*1e3*Tref_yr3)/1e3; % dH is calculated based on Topt (which is optimized)
par_MMRT_yr3(3) = (0.0033 * (par_MMRT_yr3(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9))

Topt_yr4 = alphaT*avgTempForMMRT_yr4 + betaT + 273.15;
Tref_yr4 = 299.9;
par_MMRT_yr4 = par_MMRT;
par_MMRT_yr4(1) = xga(15);
par_MMRT_yr4(2) = (Topt_yr4*(-par_MMRT_yr4(1)*1e3-R) + par_MMRT_yr4(1)*1e3*Tref_yr4)/1e3; % dH is calculated based on Topt (which is optimized)
par_MMRT_yr4(3) = (0.0033 * (par_MMRT_yr4(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9))

Topt_yr5 = alphaT*avgTempForMMRT_yr5 + betaT + 273.15;
Tref_yr5 = 299.9;
par_MMRT_yr5 = par_MMRT;
par_MMRT_yr5(1) = xga(15);
par_MMRT_yr5(2) = (Topt_yr5*(-par_MMRT_yr5(1)*1e3-R) + par_MMRT_yr5(1)*1e3*Tref_yr5)/1e3; % dH is calculated based on Topt (which is optimized)
par_MMRT_yr5(3) = (0.0033 * (par_MMRT_yr5(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9))

% The MMRT factors are calculated
MMRTfactor_yr1 = MMRT(par_MMRT_yr1, Tref_yr1, temp, Topt_yr1);
MMRTfactor_yr2 = MMRT(par_MMRT_yr2, Tref_yr2, temp, Topt_yr2);
MMRTfactor_yr3 = MMRT(par_MMRT_yr3, Tref_yr3, temp, Topt_yr3);
MMRTfactor_yr4 = MMRT(par_MMRT_yr4, Tref_yr4, temp, Topt_yr4);
MMRTfactor_yr5 = MMRT(par_MMRT_yr5, Tref_yr5, temp, Topt_yr5);

Vmax_yr1 = VmaxRef.*MMRTfactor_yr1;
Vmax_yr2 = VmaxRef.*MMRTfactor_yr2;
Vmax_yr3 = VmaxRef.*MMRTfactor_yr3;
Vmax_yr4 = VmaxRef.*MMRTfactor_yr4;
Vmax_yr5 = VmaxRef.*MMRTfactor_yr5;

% The MMRT rate modifier is calculatd for all temperatures in the ranges
MMRT_rateMod_range1 = MMRT(par_MMRT_yr1, Tref_yr1, 273.15 + allTemp_range1, Topt_yr1);
MMRT_rateMod_range2 = MMRT(par_MMRT_yr2, Tref_yr2, 273.15 + allTemp_range2, Topt_yr2);
MMRT_rateMod_range3 = MMRT(par_MMRT_yr3, Tref_yr3, 273.15 + allTemp_range3, Topt_yr3);
MMRT_rateMod_range4 = MMRT(par_MMRT_yr4, Tref_yr4, 273.15 + allTemp_range4, Topt_yr4);
MMRT_rateMod_range5 = MMRT(par_MMRT_yr5, Tref_yr5, 273.15 + allTemp_range5, Topt_yr5);

% These data are formatted for plotting
yrCell = cell(numel(dates1) + numel(dates2) + numel(dates3) + numel(dates4) + numel(dates5),1);
yrCell(1:numel(dates1)) = {'Year 0'};
yrCell(numel(dates1) + 1 : numel(dates1) + numel(dates2)) = {'Year 1'};
yrCell(numel(dates1) + numel(dates2) + 1 : numel(dates1) + numel(dates2) + numel(dates3)) = {'Year 2'};
yrCell(numel(dates1) + numel(dates2) + numel(dates3) + 1 : numel(dates1) + numel(dates2) + numel(dates3) + numel(dates4)) = {'Year 3'};
yrCell(numel(dates1) + numel(dates2) + numel(dates3) + numel(dates4) + 1 : numel(dates1) + numel(dates2) + numel(dates3) + numel(dates4) + numel(dates5)) = {'Year 4'};

% --------------------------------------
% Plot all temperatures versus all rates
% --------------------------------------

lineWidth = 2;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [5 5 30 30], 'color', [1 1 1])

bw1 = .5;
bw2 = .02;

% The temperatures and MMRT values
s1 = scatterhist([allTemp_range1; allTemp_range2; allTemp_range3; allTemp_range4; allTemp_range5], ...
                        [MMRT_rateMod_range1; MMRT_rateMod_range2; MMRT_rateMod_range3; MMRT_rateMod_range4; MMRT_rateMod_range5], 'Group', yrCell, 'Kernel', 'on','MarkerSize', .5,...
    'color', [c_yr1;c_yr2; c_yr3; c_yr4; c_yr5],'LineStyle',{'-','-','-','-','-'}, 'Location', 'northeast', 'Direction','out', 'Bandwidth', [bw1, bw1, bw1, bw1, bw1; bw2, bw2, bw2, bw2,bw2]);

set(gca, 'XAxisLocation', 'bottom', 'YAxisLocation', 'left')

hold(s1(1), 'on')

% The MMRT curves
p1 = plot(temp-273.15, Vmax_yr1, 'color', c_curve, 'linewidth', lineWidth);
p2 = plot(temp-273.15, Vmax_yr2, 'color', c_curve, 'linewidth', lineWidth);
p3 = plot(temp-273.15, Vmax_yr3, 'color', c_curve, 'linewidth', lineWidth);
p4 = plot(temp-273.15, Vmax_yr4, 'color', c_curve, 'linewidth', lineWidth);
p5 = plot(temp-273.15, Vmax_yr5, 'color', c_curve, 'linewidth', lineWidth);

s2 = scatter(allTemp_range1, MMRT_rateMod_range1, 30, 'markerfacecolor', c_yr1, 'markeredgecolor', c_yr1);
s3 = scatter(allTemp_range2, MMRT_rateMod_range2, 30, 'markerfacecolor', c_yr2, 'markeredgecolor', c_yr2);
s4 = scatter(allTemp_range3, MMRT_rateMod_range3, 30, 'markerfacecolor', c_yr3, 'markeredgecolor', c_yr3);
s5 = scatter(allTemp_range4, MMRT_rateMod_range4, 30, 'markerfacecolor', c_yr4, 'markeredgecolor', c_yr4);
s6 = scatter(allTemp_range5, MMRT_rateMod_range5, 30, 'markerfacecolor', c_yr5, 'markeredgecolor', c_yr5);

set(gca, 'fontsize',24, 'YTick', [0:0.2:1])
xlabel('Temperature (°C)');
ylabel('MMRT rate modifier');
ylim([0 1])

lgd = legend([p1, s2, s3, s4, s5, s6], 'MMRT curve', 'Last year before soil warming', '1^{st} year of soil warming', '2^{nd} year of soil warming', '3^{th} year of soil warming', '4^{th} year of soil warming');
set(legend,'FontSize',24);
legend('boxoff')
pos = lgd.Position;
pos(1) = 0.1; % horizontal
pos(2) = 0.46; % vertical
set(lgd, 'Position', pos)

text(-3,0.95, 'Optimum driven', 'FontWeight', 'bold', 'FontSize', 28, 'Parent', s1(1))
text(-7.5,1.28,'(B)', 'FontWeight', 'bold', 'FontSize', 28)

% Position of the main graph
s1(1).Position(1) = 0.1; % horizontal
s1(1).Position(2) = 0.09; % vertical
s1(1).Position(3) = .67; % Width
s1(1).Position(4) = .67; % height

% Position of the upper graph
drawnow

s1(2).Position(2) = 0.78; % vertical
s1(2).Position(4) = 0.2; % height

% Position of the right graph
s1(3).Position(1) = 0.8; % horizontal
s1(3).Position(3) = 0.18; % width

%% MMRT curve for the enzyme rigidity scenario

tmp = load('GA output/Barre Woods - noTau - 5 June 2024/xga_enzRig_07-May-2024.mat');
xga = tmp.xga;

% The colors
c_yr1 = [253,187,132]./255;
c_yr2 = [252,141,89]./255;
c_yr3 = [239,101,72]./255;
c_yr4 = [215,48,31]./255;
c_yr5 = [179,0,0]./255;
c_curve = 'k';

% The year for this the curves are plotted are defined
dates1 = datetime(2002,05,21):datetime(2003,05,20);
dates2 = datetime(2003,05,21):datetime(2004,05,20);
dates3 = datetime(2004,05,21):datetime(2005,05,20);
dates4 = datetime(2005,05,21):datetime(2006,05,20);
dates5 = datetime(2006,05,21):datetime(2007,05,20);

% All temperatures of the evaluated years are isolated
[r c] = ismember(dates1,dates_all);
allTemp_range1 = soilTemperature_spinupAndWarming(c);

[r c] = ismember(dates2,dates_all);
allTemp_range2 = soilTemperature_spinupAndWarming(c);

[r c] = ismember(dates3,dates_all);
allTemp_range3 = soilTemperature_spinupAndWarming(c);

[r c] = ismember(dates4,dates_all);
allTemp_range4 = soilTemperature_spinupAndWarming(c);

[r c] = ismember(dates5,dates_all);
allTemp_range5 = soilTemperature_spinupAndWarming(c);

% The mean temperature for every based on which the MMRT curve is
% calculated is obtained, to calculate Tref_MMRT

nYears = xga(19);
% The number of time steps equal to nYears
nTimeSteps = ceil((365*nYears));

% The average temperature from every day to that day minus nYears is calculated
allAvgTemp_yr1 = NaN(numel(dates1),1);
for ii = 1:numel(dates1)
    % The index of the current day in temperature_spinupAndControl
    [r, ind] = find(dates_all == dates1(ii));
    % The average temperature is stored
    allAvgTemp_yr1(ii,1) = mean(soilTemperature_spinupAndWarming(ind-nTimeSteps:ind));
end
avgTempForMMRT_yr1 = mean(allAvgTemp_yr1); 

allAvgTemp_yr2 = NaN(numel(dates2),1);
for ii = 1:numel(dates2)
    % The index of the current day in temperature_spinupAndControl
    [r, ind] = find(dates_all == dates2(ii));
    % The average temperature is stored
    allAvgTemp_yr2(ii,1) = mean(soilTemperature_spinupAndWarming(ind-nTimeSteps:ind));
end
avgTempForMMRT_yr2 = mean(allAvgTemp_yr2); 

allAvgTemp_yr3 = NaN(numel(dates3),1);
for ii = 1:numel(dates3) 
    % The index of the current day in temperature_spinupAndControl
    [r, ind] = find(dates_all == dates3(ii));
    % The average temperature is stored
    allAvgTemp_yr3(ii,1) = mean(soilTemperature_spinupAndWarming(ind-nTimeSteps:ind));
end
avgTempForMMRT_yr3 = mean(allAvgTemp_yr3); 

allAvgTemp_yr4 = NaN(numel(dates4),1);
for ii = 1:numel(dates4) 
    % The index of the current day in temperature_spinupAndControl
    [r, ind] = find(dates_all == dates4(ii));
    % The average temperature is stored
    allAvgTemp_yr4(ii,1) = mean(soilTemperature_spinupAndWarming(ind-nTimeSteps:ind));
end
avgTempForMMRT_yr4 = mean(allAvgTemp_yr4);

allAvgTemp_yr5 = NaN(numel(dates5),1);
for ii = 1:numel(dates5) 
    % The index of the current day in temperature_spinupAndControl
    [r, ind] = find(dates_all == dates5(ii));
    % The average temperature is stored
    allAvgTemp_yr5(ii,1) = mean(soilTemperature_spinupAndWarming(ind-nTimeSteps:ind));
end
avgTempForMMRT_yr5 = mean(allAvgTemp_yr5);

% The MMRT curves are calculated
par_MMRT = [NaN;...   % Cp
            NaN; ...  % dH
            NaN];     % dS

% The calibrated parameter values
alphaT = xga(15);
betaT = xga(16);
alphaC = xga(17);
betaC = xga(18);

R = 8.314;
VmaxRef = 1;
temp = 273.15 + (-10:1:50);

% The temperature-dependent parameters are calculated
Topt_yr1 = alphaT*avgTempForMMRT_yr1 + betaT + 273.15;
Tref_yr1 = 299.9;
par_MMRT_yr1 = par_MMRT;
par_MMRT_yr1(1) = alphaC*avgTempForMMRT_yr1 + betaC;
par_MMRT_yr1(2) = (Topt_yr1*(-par_MMRT_yr1(1)*1e3-R) + par_MMRT_yr1(1)*1e3*Tref_yr1)/1e3; % dH is calculated based on Topt (which is optimized)
par_MMRT_yr1(3) = (0.0033 * (par_MMRT_yr1(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9))

Topt_yr2 = alphaT*avgTempForMMRT_yr2 + betaT + 273.15;
Tref_yr2 = 299.9;
par_MMRT_yr2 = par_MMRT;
par_MMRT_yr2(1) = alphaC*avgTempForMMRT_yr2 + betaC;
par_MMRT_yr2(2) = (Topt_yr2*(-par_MMRT_yr2(1)*1e3-R) + par_MMRT_yr2(1)*1e3*Tref_yr2)/1e3; % dH is calculated based on Topt (which is optimized)
par_MMRT_yr2(3) = (0.0033 * (par_MMRT_yr2(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9))

Topt_yr3 = alphaT*avgTempForMMRT_yr3 + betaT + 273.15;
Tref_yr3 = 299.9;
par_MMRT_yr3 = par_MMRT;
par_MMRT_yr3(1) = alphaC*avgTempForMMRT_yr3 + betaC;
par_MMRT_yr3(2) = (Topt_yr3*(-par_MMRT_yr3(1)*1e3-R) + par_MMRT_yr3(1)*1e3*Tref_yr3)/1e3; % dH is calculated based on Topt (which is optimized)
par_MMRT_yr3(3) = (0.0033 * (par_MMRT_yr3(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9))

Topt_yr4 = alphaT*avgTempForMMRT_yr4 + betaT + 273.15;
Tref_yr4 = 299.9;
par_MMRT_yr4 = par_MMRT;
par_MMRT_yr4(1) = alphaC*avgTempForMMRT_yr4 + betaC;
par_MMRT_yr4(2) = (Topt_yr4*(-par_MMRT_yr4(1)*1e3-R) + par_MMRT_yr4(1)*1e3*Tref_yr4)/1e3; % dH is calculated based on Topt (which is optimized)
par_MMRT_yr4(3) = (0.0033 * (par_MMRT_yr4(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9))

Topt_yr5 = alphaT*avgTempForMMRT_yr5 + betaT + 273.15;
Tref_yr5 = 299.9;
par_MMRT_yr5 = par_MMRT;
par_MMRT_yr5(1) = alphaC*avgTempForMMRT_yr5 + betaC;
par_MMRT_yr5(2) = (Topt_yr5*(-par_MMRT_yr5(1)*1e3-R) + par_MMRT_yr5(1)*1e3*Tref_yr5)/1e3; % dH is calculated based on Topt (which is optimized)
par_MMRT_yr5(3) = (0.0033 * (par_MMRT_yr5(2)*1e3) - 204.01) / 1e3; % dS is calculated based on dH (these are strongly correlated (R2 = 0.9))


% Topt_yr1 = alphaT*avgTempForMMRT_yr1 + betaT + 273.15;
% par_MMRT_yr1 = par_MMRT;
% Tref_yr1 = Topt_yr1 - 6;
% par_MMRT_yr1(1) = alphaC*avgTempForMMRT_yr1 + betaC;
% par_MMRT_yr1(2) = (Topt_yr1*(-par_MMRT_yr1(1)*1e3-R) + par_MMRT_yr1(1)*1e3*Tref_yr1)/1e3;
% 
% Topt_yr2 = alphaT*avgTempForMMRT_yr2 + betaT + 273.15;
% par_MMRT_yr2 = par_MMRT;
% Tref_yr2 = Topt_yr2 - 6;
% par_MMRT_yr2(1) = alphaC*avgTempForMMRT_yr2 + betaC;
% par_MMRT_yr2(2) = (Topt_yr2*(-par_MMRT_yr2(1)*1e3-R) + par_MMRT_yr2(1)*1e3*Tref_yr2)/1e3;
% 
% Topt_yr3 = alphaT*avgTempForMMRT_yr3 + betaT + 273.15;
% par_MMRT_yr3 = par_MMRT;
% Tref_yr3 = Topt_yr3 - 6;
% par_MMRT_yr3(1) = alphaC*avgTempForMMRT_yr3 + betaC;
% par_MMRT_yr3(2) = (Topt_yr3*(-par_MMRT_yr3(1)*1e3-R) + par_MMRT_yr3(1)*1e3*Tref_yr3)/1e3;
% 
% Topt_yr4 = alphaT*avgTempForMMRT_yr4 + betaT + 273.15;
% par_MMRT_yr4 = par_MMRT;
% Tref_yr4 = Topt_yr4 - 6;
% par_MMRT_yr4(1) = alphaC*avgTempForMMRT_yr4 + betaC;
% par_MMRT_yr4(2) = (Topt_yr4*(-par_MMRT_yr4(1)*1e3-R) + par_MMRT_yr4(1)*1e3*Tref_yr4)/1e3;
% 
% Topt_yr5 = alphaT*avgTempForMMRT_yr5 + betaT + 273.15;
% par_MMRT_yr5 = par_MMRT;
% Tref_yr5 = Topt_yr5 - 6;
% par_MMRT_yr5(1) = alphaC*avgTempForMMRT_yr5 + betaC;
% par_MMRT_yr5(2) = (Topt_yr5*(-par_MMRT_yr5(1)*1e3-R) + par_MMRT_yr5(1)*1e3*Tref_yr5)/1e3;

% The MMRT factors are calculated
MMRTfactor_yr1 = MMRT(par_MMRT_yr1, Tref_yr1, temp, Topt_yr1);
MMRTfactor_yr2 = MMRT(par_MMRT_yr2, Tref_yr2, temp, Topt_yr2);
MMRTfactor_yr3 = MMRT(par_MMRT_yr3, Tref_yr3, temp, Topt_yr3);
MMRTfactor_yr4 = MMRT(par_MMRT_yr4, Tref_yr4, temp, Topt_yr4);
MMRTfactor_yr5 = MMRT(par_MMRT_yr5, Tref_yr5, temp, Topt_yr5);

Vmax_yr1 = VmaxRef.*MMRTfactor_yr1;
Vmax_yr2 = VmaxRef.*MMRTfactor_yr2;
Vmax_yr3 = VmaxRef.*MMRTfactor_yr3;
Vmax_yr4 = VmaxRef.*MMRTfactor_yr4;
Vmax_yr5 = VmaxRef.*MMRTfactor_yr5;

% The MMRT rate modifier is calculatd for all temperatures in the ranges
MMRT_rateMod_range1 = MMRT(par_MMRT_yr1, Tref_yr1, 273.15 + allTemp_range1, Topt_yr1);
MMRT_rateMod_range2 = MMRT(par_MMRT_yr2, Tref_yr2, 273.15 + allTemp_range2, Topt_yr2);
MMRT_rateMod_range3 = MMRT(par_MMRT_yr3, Tref_yr3, 273.15 + allTemp_range3, Topt_yr3);
MMRT_rateMod_range4 = MMRT(par_MMRT_yr4, Tref_yr4, 273.15 + allTemp_range4, Topt_yr4);
MMRT_rateMod_range5 = MMRT(par_MMRT_yr5, Tref_yr5, 273.15 + allTemp_range5, Topt_yr5);

% These data are formatted for plotting
yrCell = cell(numel(dates1) + numel(dates2) + numel(dates3) + numel(dates4) + numel(dates5),1);
yrCell(1:numel(dates1)) = {'Year 0'};
yrCell(numel(dates1) + 1 : numel(dates1) + numel(dates2)) = {'Year 1'};
yrCell(numel(dates1) + numel(dates2) + 1 : numel(dates1) + numel(dates2) + numel(dates3)) = {'Year 2'};
yrCell(numel(dates1) + numel(dates2) + numel(dates3) + 1 : numel(dates1) + numel(dates2) + numel(dates3) + numel(dates4)) = {'Year 3'};
yrCell(numel(dates1) + numel(dates2) + numel(dates3) + numel(dates4) + 1 : numel(dates1) + numel(dates2) + numel(dates3) + numel(dates4) + numel(dates5)) = {'Year 4'};

% --------------------------------------
% Plot all temperatures versus all rates
% --------------------------------------

lineWidth = 2;

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [5 5 30 30], 'color', [1 1 1])

bw1 = .5;
bw2 = .02;

% The temperatures and MMRT values
s1 = scatterhist([allTemp_range1; allTemp_range2; allTemp_range3; allTemp_range4; allTemp_range5], ...
                        [MMRT_rateMod_range1; MMRT_rateMod_range2; MMRT_rateMod_range3; MMRT_rateMod_range4; MMRT_rateMod_range5], 'Group', yrCell, 'Kernel', 'on','MarkerSize', .5,...
    'color', [c_yr1;c_yr2; c_yr3; c_yr4; c_yr5],'LineStyle',{'-','-','-','-','-'}, 'Location', 'northeast', 'Direction','out', 'Bandwidth', [bw1, bw1, bw1, bw1, bw1; bw2, bw2, bw2, bw2,bw2]);

set(gca, 'XAxisLocation', 'bottom', 'YAxisLocation', 'left')

hold(s1(1), 'on')

% The MMRT curves
p1 = plot(temp-273.15, Vmax_yr1, 'color', c_curve, 'linewidth', lineWidth);
p2 = plot(temp-273.15, Vmax_yr2, 'color', c_curve, 'linewidth', lineWidth);
p3 = plot(temp-273.15, Vmax_yr3, 'color', c_curve, 'linewidth', lineWidth);
p4 = plot(temp-273.15, Vmax_yr4, 'color', c_curve, 'linewidth', lineWidth);
p5 = plot(temp-273.15, Vmax_yr5, 'color', c_curve, 'linewidth', lineWidth);

s2 = scatter(allTemp_range1, MMRT_rateMod_range1, 30, 'markerfacecolor', c_yr1, 'markeredgecolor', c_yr1);
s3 = scatter(allTemp_range2, MMRT_rateMod_range2, 30, 'markerfacecolor', c_yr2, 'markeredgecolor', c_yr2);
s4 = scatter(allTemp_range3, MMRT_rateMod_range3, 30, 'markerfacecolor', c_yr3, 'markeredgecolor', c_yr3);
s5 = scatter(allTemp_range4, MMRT_rateMod_range4, 30, 'markerfacecolor', c_yr4, 'markeredgecolor', c_yr4);
s6 = scatter(allTemp_range5, MMRT_rateMod_range5, 30, 'markerfacecolor', c_yr5, 'markeredgecolor', c_yr5);

set(gca, 'fontsize',24,'YTick',[0:0.2:1])
xlabel('Temperature (°C)');
ylabel('MMRT rate modifier');
ylim([0 1])

lgd = legend([p1, s2, s3, s4, s5, s6], 'MMRT curve', 'Last year before soil warming', '1^{st} year of soil warming', '2^{nd} year of soil warming', '3^{th} year of soil warming', '4^{th} year of soil warming');
set(legend,'FontSize',23);
legend('boxoff')
pos = lgd.Position;
pos(1) = 0.1; % horizontal
pos(2) = 0.46; % vertical
set(lgd, 'Position', pos)

text(-3,0.95, 'Enzyme rigidity', 'FontWeight', 'bold', 'FontSize', 28, 'Parent', s1(1))
text(-7.5,1.28,'(C)', 'FontWeight', 'bold', 'FontSize', 28)

% Position of the main graph
s1(1).Position(1) = 0.1; % horizontal
s1(1).Position(2) = 0.09; % vertical
s1(1).Position(3) = .67; % Width
s1(1).Position(4) = .67; % height

% Position of the upper graph
drawnow

s1(2).Position(2) = 0.78; % vertical
s1(2).Position(4) = 0.2; % height

% Position of the right graph
s1(3).Position(1) = 0.8; % horizontal
s1(3).Position(3) = 0.18; % width








































