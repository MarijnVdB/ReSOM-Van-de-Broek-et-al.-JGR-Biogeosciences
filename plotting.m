%%  Plotting everything together

close all

% Decide which plots are displayed
plot_temperature = 0;
plot_litter = 1;
plot_soil = 1;
plot_CO2 = 1;
plot_annual_CO2 = 1;
plot_CO2_per_microbe = 0;
plot_massSpec_resp = 0;

if variableTimeStep == 0
   dt_array_full = ones(numel(dates_all),1).*dt;
end

%% The measured data are formatted for plotting

% In case the treatment run is longer than the period for which data is
% available, the variable 'dates_treatmentRun' is adjusted

% If the end date is before the last years with available data
if max(dates_treatmentRun) > max(dates_CO2_measurements_control)
    dates_treatmentRun = datesWithData; 
end

%% The temperature data is plotted

if plot_temperature == 1
    
    hfig = figure;
    set(hfig, 'Units', 'centimeter', 'Position', [3 3 25 10], 'color', [1 1 1])
    hold on
    s1 = scatter(dates_all', soilTemperature_spinupAndControl);
    s2 = scatter(dates_all', soilTemperature_spinupAndWarming);

    legend([s1, s2], 'Temperature spinup and contol run', 'Temperature warming run', ...
        'location', 'northwest')
    legend('boxoff')
    
end

%% Plotting litter pools for the control run

if plot_litter == 1
    
ymin_litter = 0;
ymax_litter = 3300;

lineWidth = 1.5;

hfig = figure;
set(hfig,'Units','centimeters', 'Position',[3 3 16 20],'color',[1 1 1]);

% The change in C stocks
subplot(3,3,[1 2])
hold on

% The colors are defined
color_micr_rStrat_litter = [252,141,89]./255;
color_micr_kStrat_litter = 'r';%[252,141,89]./255;
color_micReserve_rStrat_litter = [223,101,176]./255;
color_micReserve_kStrat_litter = [223,101,176]./255;
color_monomers_litter = [140,107,177]./255;
color_adsMonomers_litter = [140,107,177]./255;
color_metabolic_litter = [254,196,79]./255;
color_structural_litter = [217,95,14]./255;
color_enz_rStrat_litter = [0,109,44]./255;
color_enz_kStrat_litter = [8,81,156]./255;
color_adsEnz_rStrat_litter = [0,109,44]./255;
color_adsEnz_kStrat_litter = [8,81,156]./255;
color_CO2_rStrat_litter = [166,189,219]./255;
color_CO2_kStrat_litter = [8,81,156]./255;
color_totalCarbon_litter = 'k';
color_measuredPOC = [153,52,4]./255;

% The C pools are defined
micr_rStrat_litter = Cpools_litter_control(2:end,1);
micr_kStrat_litter = Cpools_litter_control(2:end,2);
micReserve_rStrat_litter = Cpools_litter_control(2:end,3);
micReserve_kStrat_litter = Cpools_litter_control(2:end,4);
surfaces_litter = Cpools_litter_control(2:end,5);
monomers_litter = Cpools_litter_control(2:end,6);
adsMonomers_litter = Cpools_litter_control(2:end,7);
metabolic_litter = Cpools_litter_control(2:end,8);
structural_litter = Cpools_litter_control(2:end,9);
enz_rStrat_litter = Cpools_litter_control(2:end,10);
enz_kStrat_litter = Cpools_litter_control(2:end,11);
adsEnz_rStrat_litter = Cpools_litter_control(2:end,12);
adsEnz_kStrat_litter = Cpools_litter_control(2:end,13);
CO2_rStrat_litter_control = diff(Cpools_litter_control(2:end,14))./dt_array_full(2:end);
CO2_kStrat_litter_control = diff(Cpools_litter_control(2:end,15))./dt_array_full(2:end);

totalCarbon_litter_control = micr_rStrat_litter + micr_kStrat_litter + micReserve_rStrat_litter + ...
    micReserve_kStrat_litter + monomers_litter + adsMonomers_litter + metabolic_litter + structural_litter + ...
    enz_rStrat_litter + enz_kStrat_litter + adsEnz_rStrat_litter + adsEnz_kStrat_litter;

% The measured data is plotted
s1 = scatter(date_SOC_measurement, Cstock_forestFloor_control,'k');
s2 = scatter(date_SOC_measurement, Cstock_forestFloor_control*0.058, 'markeredgecolor',color_monomers_litter);           % MAOC
s4 = scatter(date_SOC_measurement, Cstock_forestFloor_control*0.942*0.66, 'markeredgecolor',color_metabolic_litter);     % Metabolic
s5 = scatter(date_SOC_measurement, Cstock_forestFloor_control*0.942*0.34, 'markeredgecolor',color_structural_litter);     % Structural

% Data are plotted
p1 = plot(dates_all,micr_rStrat_litter,'color', color_micr_rStrat_litter, 'LineWidth', lineWidth); % Microbial biomass
p2 = plot(dates_all,micr_kStrat_litter, '--','color', color_micr_kStrat_litter, 'LineWidth', lineWidth); % Reserve biomass
p3 = plot(dates_all,micReserve_rStrat_litter, '--','color', color_micReserve_rStrat_litter, 'LineWidth', lineWidth); % Monomers
p4 = plot(dates_all,micReserve_kStrat_litter,'color', color_micReserve_kStrat_litter, 'LineWidth', lineWidth); % Adsorbed monomers
p5 = plot(dates_all,monomers_litter,'color', color_monomers_litter, 'LineWidth', lineWidth); % Polymers
p6 = plot(dates_all,adsMonomers_litter, '--','color', color_monomers_litter, 'LineWidth', lineWidth); % Polymers
p7 = plot(dates_all,metabolic_litter, '--','color', color_metabolic_litter, 'LineWidth', lineWidth); % Enzymes
p8 = plot(dates_all,structural_litter, 'color', color_structural_litter, 'LineWidth', lineWidth); % Adsorbed enzymes
p9 = plot(dates_all,enz_rStrat_litter, 'color', color_enz_rStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes
p10 = plot(dates_all,enz_kStrat_litter, 'color', color_enz_kStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes
p11 = plot(dates_all,adsEnz_rStrat_litter, '--', 'color', color_adsEnz_rStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes
p12 = plot(dates_all,adsEnz_kStrat_litter, '--', 'color', color_adsEnz_kStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes
p13 = plot(dates_all,totalCarbon_litter_control, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

ylim([ymin_litter ymax_litter])

title('Forest floor C stock (g m^{-2}) - Control run')
xlabel('Years')
ylabel('g C m^{-2}')

% CO2 emissions
subplot(3,3,[4 5 6])
hold on

p14 = plot(dates_all(2:end),CO2_rStrat_litter_control, 'color', color_CO2_rStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes
p15 = plot(dates_all(2:end),CO2_kStrat_litter_control, '--', 'color', color_CO2_kStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

set(gca,'box','off')

title('CO_{2} emissions (g C m^{-2} d^{-1})')
xlabel('Years')
ylabel('(g C m^{-2}) per timestep')

% Carbon use efficiency
subplot(3,3,[7 8 9])
hold on

cue_rStrat = diff(Cpools_litter_control(:,16))./dt_array_full;
cue_kStrat = diff(Cpools_litter_control(:,17))./dt_array_full;

p16 = plot(dates_all,cue_kStrat, 'color', color_CO2_kStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes
p17 = plot(dates_all,cue_rStrat, '--', 'color', color_CO2_rStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

legend([p16, p17],'CUE K-strategists','CUE r-strategists', 'location', 'southwest')

title('Carbon use efficiency')
xlabel('Years')
ylabel('CUE')

% -------------------------------------------------------------------------
% A new invisible axis is created
% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% axes(ax1)
% xlim([0 1])
% ylim([0 1])

Legend1 = legend([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12], ...
    'r-strategists','K-strategists','Reserve r-strategists',...
    'Reserve K-strategists','Monomers', 'Adsorbed monomers','Metabolic litter',...
    'Structural litter','Enzymes r-strategists','Enzymes K-strategists', ...
    'Adsorbed Enz r-strategists','Adsorbed Enz K-strategists');

set(Legend1,'FontSize',9, 'color','w','edgecolor','w');
legend('boxoff')

pos = Legend1.Position;
pos(1) = 0.65;
pos(2) = 0.73;
set(Legend1, 'Position', pos)

%% Plotting litter pools for the warmed run

lineWidth = 1.5;

hfig = figure;
set(hfig,'Units','centimeters', 'Position',[3 3 16 20],'color',[1 1 1]);

% The change in C stocks
subplot(3,3,[1 2])
hold on

% The colors are defined
color_micr_rStrat_litter = [252,141,89]./255;
color_micr_kStrat_litter = 'r';%[252,141,89]./255;
color_micReserve_rStrat_litter = [223,101,176]./255;
color_micReserve_kStrat_litter = [223,101,176]./255;
color_monomers_litter = [140,107,177]./255;
color_adsMonomers_litter = [140,107,177]./255;
color_metabolic_litter = [254,196,79]./255;
color_structural_litter = [217,95,14]./255;
color_enz_rStrat_litter = [0,109,44]./255;
color_enz_kStrat_litter = [8,81,156]./255;
color_adsEnz_rStrat_litter = [0,109,44]./255;
color_adsEnz_kStrat_litter = [8,81,156]./255;
color_CO2_rStrat_litter = [166,189,219]./255;
color_CO2_kStrat_litter = [8,81,156]./255;
color_totalCarbon_litter = 'k';

% The C pools are defined
micr_rStrat_litter = Cpools_litter_warming(2:end,1);
micr_kStrat_litter = Cpools_litter_warming(2:end,2);
micReserve_rStrat_litter = Cpools_litter_warming(2:end,3);
micReserve_kStrat_litter = Cpools_litter_warming(2:end,4);
surfaces_litter = Cpools_litter_warming(2:end,5);
monomers_litter = Cpools_litter_warming(2:end,6);
adsMonomers_litter = Cpools_litter_warming(2:end,7);
metabolic_litter = Cpools_litter_warming(2:end,8);
structural_litter = Cpools_litter_warming(2:end,9);
enz_rStrat_litter = Cpools_litter_warming(2:end,10);
enz_kStrat_litter = Cpools_litter_warming(2:end,11);
adsEnz_rStrat_litter = Cpools_litter_warming(2:end,12);
adsEnz_kStrat_litter = Cpools_litter_warming(2:end,13);
CO2_rStrat_litter_control = diff(Cpools_litter_warming(2:end,14))./dt_array_full(2:end);
CO2_kStrat_litter_control = diff(Cpools_litter_warming(2:end,15))./dt_array_full(2:end);

totalCarbon_litter_warming = micr_rStrat_litter + micr_kStrat_litter + micReserve_rStrat_litter + ...
    micReserve_kStrat_litter + monomers_litter + adsMonomers_litter + metabolic_litter + structural_litter + ...
    enz_rStrat_litter + enz_kStrat_litter + adsEnz_rStrat_litter + adsEnz_kStrat_litter;

% The measures data is plotted
s1 = scatter(date_SOC_measurement, Cstock_forestFloor_warming,'k');
s2 = scatter(date_SOC_measurement, Cstock_forestFloor_warming*0.058, 'markeredgecolor',color_monomers_litter);           % MAOC
s4 = scatter(date_SOC_measurement, Cstock_forestFloor_warming*0.942*0.66, 'markeredgecolor',color_metabolic_litter);     % Metabolic
s5 = scatter(date_SOC_measurement, Cstock_forestFloor_warming*0.942*0.34, 'markeredgecolor',color_structural_litter);     % Structural

% Data are plotted
p1 = plot(dates_all,micr_rStrat_litter,'color', color_micr_rStrat_litter, 'LineWidth', lineWidth); % Microbial biomass
p2 = plot(dates_all,micr_kStrat_litter, '--','color', color_micr_kStrat_litter, 'LineWidth', lineWidth); % Reserve biomass
p3 = plot(dates_all,micReserve_rStrat_litter, '--','color', color_micReserve_rStrat_litter, 'LineWidth', lineWidth); % Monomers
p4 = plot(dates_all,micReserve_kStrat_litter,'color', color_micReserve_kStrat_litter, 'LineWidth', lineWidth); % Adsorbed monomers
p5 = plot(dates_all,monomers_litter,'color', color_monomers_litter, 'LineWidth', lineWidth); % Polymers
p6 = plot(dates_all,adsMonomers_litter, '--','color', color_monomers_litter, 'LineWidth', lineWidth); % Polymers
p7 = plot(dates_all,metabolic_litter, '--','color', color_metabolic_litter, 'LineWidth', lineWidth); % Enzymes
p8 = plot(dates_all,structural_litter, 'color', color_structural_litter, 'LineWidth', lineWidth); % Adsorbed enzymes
p9 = plot(dates_all,enz_rStrat_litter, 'color', color_enz_rStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes
p10 = plot(dates_all,enz_kStrat_litter, 'color', color_enz_kStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes
p11 = plot(dates_all,adsEnz_rStrat_litter, '--', 'color', color_adsEnz_rStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes
p12 = plot(dates_all,adsEnz_kStrat_litter, '--', 'color', color_adsEnz_kStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes
p13 = plot(dates_all,totalCarbon_litter_warming, 'color', color_totalCarbon_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

% A line indicating when warming started
initYear = date_startWarming;%datetime('1990', 'InputFormat', 'yyyy');
plot([initYear initYear], [ymin_litter ymax_litter], '--', 'Color', [.3 .3 .3])

ylim([ymin_litter ymax_litter])

title('Forest floor C stock (g m^{-2}) - warming run')
xlabel('Years')
ylabel('g m^{-2}')

% CO2 emissions
subplot(3,3,[4 5 6])
hold on

p14 = plot(dates_all(2:end),CO2_rStrat_litter_control, 'color', color_CO2_rStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes
p15 = plot(dates_all(2:end),CO2_kStrat_litter_control, '--', 'color', color_CO2_kStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

% A line indicating when warming started
initYear = datetime('1990', 'InputFormat', 'yyyy');
plot([initYear initYear], [0 40], '--', 'Color', [.3 .3 .3])

ylim([0 2])

set(gca,'box','off')

title('CO_{2} emissions (g C m^{-2} d^{-1})')
xlabel('Years')
ylabel('(g C m^{-2}) per timestep')

% Carbon use efficiency
subplot(3,3,[7 8 9])
hold on

cue_rStrat = diff(Cpools_litter_warming(:,16))./dt_array_full;
cue_kStrat = diff(Cpools_litter_warming(:,17))./dt_array_full;

p16 = plot(dates_all,cue_kStrat, 'color', color_CO2_kStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes
p17 = plot(dates_all,cue_rStrat, '--', 'color', color_CO2_rStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

% A line indicating when warming started
initYear = datetime('1990', 'InputFormat', 'yyyy');
plot([initYear initYear], [0 40], '--', 'Color', [.3 .3 .3])

ylim([0 1])

legend([p16, p17],'CUE K-strategists','CUE r-strategists', 'location', 'southwest')

title('Carbon use efficiency')
xlabel('Years')
ylabel('CUE')

% -------------------------------------------------------------------------
% A new invisible axis is created
% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% axes(ax1)
% xlim([0 1])
% ylim([0 1])

Legend1 = legend([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12], ...
    'r-strategists','K-strategists','Reserve r-strategists',...
    'Reserve K-strategists','Monomers', 'Adsorbed monomers','Metabolic litter',...
    'Structural litter','Enzymes r-strategists','Enzymes K-strategists', ...
    'Adsorbed Enz r-strategists','Adsorbed Enz K-strategists');

set(Legend1,'FontSize',9, 'color','w','edgecolor','w');
legend('boxoff')

pos = Legend1.Position;
pos(1) = 0.65;
pos(2) = 0.73;
set(Legend1, 'Position', pos)

end

%% Plotting soil pools for the control run

if plot_soil == 1
    
ymin_soil = 0;
ymax_soil = 10000;

% Line colors
color_microbes_soil = [252,141,89]./255;
color_reserve_soil = [223,101,176]./255;
color_freeEnzymes_soil = [78,179,211]./255;
color_adsEnzymes_soil = [78,179,211]./255;
color_freeMonomers_soil = [140,107,177]./255;
color_adsMonomers_soil = [140,107,177]./255;
color_polymers_soil = [153,52,4]./255;
color_freeSurfaces_soil = [2,129,138]./255;
color_CO2_soil = [54,144,192]./255;
color_totalCarbon_soil = 'k';

% The data is converted to OC%
micrBiomass_conc_soil = Cpools_soil_control(2:end,1);
micrReserve_conc_soil = Cpools_soil_control(2:end,2);
monomers_conc_soil = Cpools_soil_control(2:end,4);
adsMonomers_conc_soil = Cpools_soil_control(2:end,5);
polymers_conc_soil = Cpools_soil_control(2:end,6);
enzymes_conc_soil = Cpools_soil_control(2:end,7);
adsEnzymes_conc_soil = Cpools_soil_control(2:end,8);
totalCarbon_conc_soil_control = micrBiomass_conc_soil + micrReserve_conc_soil + monomers_conc_soil + ...
    adsMonomers_conc_soil + polymers_conc_soil + enzymes_conc_soil + adsEnzymes_conc_soil;

lineWidth = 1.5;

hfig = figure;
set(hfig,'Units','centimeters', 'Position',[3 3 16 20],'color',[1 1 1]);

% The change in C stocks
subplot(3,3,[1 2])
hold on

% The measures data is plotted
s1 = scatter(date_SOC_measurement, Cstock_soil_control,'k');
s1 = scatter(date_SOC_measurement, Cstock_soil_control.*0.58,'MarkerEdgeColor', color_adsMonomers_soil);
s1 = scatter(date_SOC_measurement, Cstock_soil_control.*0.42,'MarkerEdgeColor', color_polymers_soil);

% Data are plotted as %OC
p1 = plot(dates_all,micrBiomass_conc_soil,'color', color_microbes_soil, 'LineWidth', lineWidth); % Microbial biomass
p2 = plot(dates_all,micrReserve_conc_soil,'color', color_reserve_soil, 'LineWidth', lineWidth); % Reserve biomass
p3 = plot(dates_all,monomers_conc_soil,'color', color_freeMonomers_soil, 'LineWidth', lineWidth); % Monomers
p4 = plot(dates_all,adsMonomers_conc_soil, '--','color', color_adsMonomers_soil, 'LineWidth', lineWidth); % Adsorbed monomers
p5 = plot(dates_all,polymers_conc_soil,'color', color_polymers_soil, 'LineWidth', lineWidth); % Polymers
p6 = plot(dates_all,enzymes_conc_soil,'color', color_freeEnzymes_soil, 'LineWidth', lineWidth); % Enzymes
p7 = plot(dates_all,adsEnzymes_conc_soil, '--','color', color_adsEnzymes_soil, 'LineWidth', lineWidth); % Adsorbed enzymes
p8 = plot(dates_all,totalCarbon_conc_soil_control, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

title('SOC stock (g m^{-1}) - control run')
xlabel('Years')
ylabel('g m^{-2}')
ylim([ymin_soil ymax_soil])

% The carbon use efficiency

subplot(3,3,[4 5 6])

cue_soil = diff(Cpools_soil_control(:,10))./dt_array_full;

plot(dates_all,cue_soil,'color', color_freeSurfaces_soil, 'LineWidth', lineWidth); % Free surfaces

set(gca,'box','off')
xlabel('Years')
title('Carbon use efficiency')

% CO2 emissions
subplot(3,3,[7 8 9])

plot(dates_all(2:end),diff(Cpools_soil_control(2:end,9))./dt_array_full(2:end),'color', color_CO2_soil, 'LineWidth', lineWidth);

% Legend
set(gca,'box','off')

title('CO2 emissions (g C m^{-2})')
xlabel('Years')
ylabel('(g C m^{-2})')

Legend1 = legend([p1, p2, p3, p4, p5, p6, p7, p8], ...
    'Microbial biomass','Reserve biomass','Monomers',...
    'Adsorbed monomers','Polymers', 'Enzymes','Adsorbed enzymes',...
    'Total SOC');

set(Legend1,'FontSize',9, 'color','w','edgecolor','w');

pos = Legend1.Position;
pos(1) = 0.65;
pos(2) = 0.73;
set(Legend1, 'Position', pos)

%% Plotting soil pools for the warmed run

% Line colors
color_microbes_soil = [252,141,89]./255;
color_reserve_soil = [223,101,176]./255;
color_freeEnzymes_soil = [78,179,211]./255;
color_adsEnzymes_soil = [78,179,211]./255;
color_freeMonomers_soil = [140,107,177]./255;
color_adsMonomers_soil = [140,107,177]./255;
color_polymers_soil = [153,52,4]./255;
color_freeSurfaces_soil = [2,129,138]./255;
color_CO2_soil = [54,144,192]./255;
color_totalCarbon_soil = 'k';

% The data is converted to OC%
micrBiomass_conc_soil = Cpools_soil_warming(2:end,1);
micrReserve_conc_soil = Cpools_soil_warming(2:end,2);
monomers_conc_soil = Cpools_soil_warming(2:end,4);
adsMonomers_conc_soil = Cpools_soil_warming(2:end,5);
polymers_conc_soil = Cpools_soil_warming(2:end,6);
enzymes_conc_soil = Cpools_soil_warming(2:end,7);
adsEnzymes_conc_soil = Cpools_soil_warming(2:end,8);
totalCarbon_conc_soil_warming = micrBiomass_conc_soil + micrReserve_conc_soil + monomers_conc_soil + adsMonomers_conc_soil + polymers_conc_soil + enzymes_conc_soil + adsEnzymes_conc_soil;

lineWidth = 1.5;

hfig = figure;
set(hfig,'Units','centimeters', 'Position',[3 3 16 20],'color',[1 1 1]);

% The change in C stocks
subplot(3,3,[1 2])
hold on

s1 = scatter(date_SOC_measurement, Cstock_soil_warming,'k');
s1 = scatter(date_SOC_measurement, Cstock_soil_warming.*0.58,'MarkerEdgeColor', color_adsMonomers_soil);
s1 = scatter(date_SOC_measurement, Cstock_soil_warming.*0.42,'MarkerEdgeColor', color_polymers_soil);

% Data are plotted as %OC
p1 = plot(dates_all,micrBiomass_conc_soil,'color', color_microbes_soil, 'LineWidth', lineWidth); % Microbial biomass
p2 = plot(dates_all,micrReserve_conc_soil,'color', color_reserve_soil, 'LineWidth', lineWidth); % Reserve biomass
p3 = plot(dates_all,monomers_conc_soil,'color', color_freeMonomers_soil, 'LineWidth', lineWidth); % Monomers
p4 = plot(dates_all,adsMonomers_conc_soil, '--','color', color_adsMonomers_soil, 'LineWidth', lineWidth); % Adsorbed monomers
p5 = plot(dates_all,polymers_conc_soil,'color', color_polymers_soil, 'LineWidth', lineWidth); % Polymers
p6 = plot(dates_all,enzymes_conc_soil,'color', color_freeEnzymes_soil, 'LineWidth', lineWidth); % Enzymes
p7 = plot(dates_all,adsEnzymes_conc_soil, '--','color', color_adsEnzymes_soil, 'LineWidth', lineWidth); % Adsorbed enzymes
p8 = plot(dates_all,totalCarbon_conc_soil_warming, 'color', color_totalCarbon_soil, 'LineWidth', lineWidth); % Adsorbed enzymes

% A line indicating when warming started
initYear = date_startWarming;%datetime('1990', 'InputFormat', 'yyyy');
plot([initYear initYear], [ymin_soil ymax_soil], '--', 'Color', [.3 .3 .3])

title('SOC stock (g m^{-1}) - warmed run')
xlabel('Years')
ylabel('g m^{-2}')
ylim([ymin_soil ymax_soil])

% The carbon use efficiency

subplot(3,3,[4 5 6])

cue_soil = diff(Cpools_soil_warming(:,10))./dt_array_full;

plot(dates_all,cue_soil,'color', color_freeSurfaces_soil, 'LineWidth', lineWidth); % Free surfaces

set(gca,'box','off')
xlabel('Years')
title('Carbon use efficiency')

% CO2 emissions
subplot(3,3,[7 8 9])
hold on

plot(dates_all(2:end),diff(Cpools_soil_warming(2:end,9))./dt_array_full(2:end),'color', color_CO2_soil, 'LineWidth', lineWidth);

% A line indicating when warming started
initYear = datetime('1990', 'InputFormat', 'yyyy');
plot([initYear initYear], [0 3], '--', 'Color', [.3 .3 .3])

set(gca,'box','off')

title('CO2 emissions (g C m^{-2})')
xlabel('Years')
ylabel('(g C m^{-2})')

% Legend
set(gca,'box','off')

title('CO2 emissions (g C m^{-2})')
xlabel('Years')
ylabel('(g C m^{-2})')

Legend1 = legend([p1, p2, p3, p4, p5, p6, p7, p8], ...
    'Microbial biomass','Reserve biomass','Monomers',...
    'Adsorbed monomers','Polymers', 'Enzymes','Adsorbed enzymes',...
    'Total SOC');

set(Legend1,'FontSize',9, 'color','w','edgecolor','w');

pos = Legend1.Position;
pos(1) = 0.65;
pos(2) = 0.73;
set(Legend1, 'Position', pos)

end

%% The total CO2 flux is plotted for the control run

if plot_CO2 == 1

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [5 5 20 15], 'color', [1 1 1]);
hold on

% The colors are defined
color_CO2_rStrat_litter = [166,189,219]./255;
color_CO2_kStrat_litter = [43,140,190]./255;
color_CO2_soil = [35,132,67]./255;
markerColor = [215,48,31]./255;%[.5 .5 .5];

% The data are formatted
CO2_rStrat_litter_control = diff(Cpools_litter_control(1:end,14))./dt_array_full;
CO2_kStrat_litter_control = diff(Cpools_litter_control(1:end,15))./dt_array_full;
CO2_soil_control = diff(Cpools_soil_control(1:end,9))./dt_array_full;

% CO2 flux from forest floor
p1 = plot(dates_all(1:end),CO2_rStrat_litter_control, 'color', color_CO2_rStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes
p2 = plot(dates_all(1:end),CO2_kStrat_litter_control, '--', 'color', color_CO2_kStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

% CO2 flux from the mineral soil
p3 = plot(dates_all(1:end),CO2_soil_control,'color', color_CO2_soil, 'LineWidth', lineWidth);

% The total CO2 flux
totalCO2_control = CO2_rStrat_litter_control + CO2_kStrat_litter_control + CO2_soil_control;
p4 = plot(dates_all(1:end),totalCO2_control,'color', 'k', 'LineWidth', lineWidth);

% The measured CO2 fluxes are plotted
s1 = scatter(dates_CO2_measurements_control, Measured_CO2_flux_control, 10, 'filled', 'markerfacecolor', markerColor, 'markeredgecolor', markerColor);

xlim([datetime(1985, 01, 01), dates_all(end)])
ylim([0 6])

title('CO_{2} fluxes for the control run')

Legend1 = legend([p1, p2, p3], ...
    'r-strat forest floor','K-strat forest floor','Soil');

set(Legend1,'FontSize',9, 'color','w','edgecolor','w');
legend('boxoff')

pos = Legend1.Position;
pos(1) = 0.65;
pos(2) = 0.73;
set(Legend1, 'Position', pos)

%% The total CO2 flux is plotted for the warmed run

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [5 5 20 15], 'color', [1 1 1]);
hold on

% The colors are defined
color_CO2_rStrat_litter = [166,189,219]./255;
color_CO2_kStrat_litter = [43,140,190]./255;
color_CO2_soil = [35,132,67]./255;
markerColor = [215,48,31]./255;%[.5 .5 .5];

% The data are formatted
CO2_rStrat_litter_warmed = diff(Cpools_litter_warming(1:end,14))./dt_array_full;
CO2_kStrat_litter_warmed = diff(Cpools_litter_warming(1:end,15))./dt_array_full;
CO2_soil_warmed = diff(Cpools_soil_warming(1:end,9))./dt_array_full;

% CO2 flux from forest floor
p1 = plot(dates_all(1:end),CO2_rStrat_litter_warmed, 'color', color_CO2_rStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes
p2 = plot(dates_all(1:end),CO2_kStrat_litter_warmed, '--', 'color', color_CO2_kStrat_litter, 'LineWidth', lineWidth); % Adsorbed enzymes

% CO2 flux from the mineral soil
p3 = plot(dates_all(1:end),CO2_soil_warmed,'color', color_CO2_soil, 'LineWidth', lineWidth);

% The total CO2 flux
totalCO2_warmed = CO2_rStrat_litter_warmed + CO2_kStrat_litter_warmed + CO2_soil_warmed;
p4 = plot(dates_all(1:end),totalCO2_warmed,'color', 'k', 'LineWidth', lineWidth);

% The measured CO2 fluxes are plotted
if ~strcmp(site,'artificialRun') % In case of the artificial run, no measured data is available
    s1 = scatter(dates_CO2_measurements_warming, Measured_CO2_flux_warming, 10, 'filled', 'markerfacecolor', markerColor, 'markeredgecolor', markerColor);
end

xlim([datetime(1985, 01, 01), dates_all(end)])
ylim([0 6])

title('CO_{2} fluxes for the warmed run')

%% Scatterplot of measured versus modelled CO2 per timestep

% This is only possible if the timestep is 1 day, since the measured CO2
% data has this unit

if dt == 1 || variableTimeStep == 1
    hfig = figure;
    set(hfig, 'Units', 'centimeter', 'Position', [5 5 25 15], 'color', [1 1 1]);

    % The indices of the modelled timeseries for which measured data is
    % available are isolated
    [r c] = ismember(dates_all, allMeasurements.dates_CO2_measurements_control);
    
    dotColor = [54,144,192]./255;
    color_posRes = [228,26,28]./255;
    color_negRes = [55,126,184]./255;
    
    % -------------
    % Control plots
    % -------------
    
    subplot(1,2,1)
    
    % The measurement dates that are modelled are isolated
    tmpMeas = allMeasurements.Measured_CO2_flux_control(allMeasurements.dates_CO2_measurements_control< endDate);
    tmpMeas_control = tmpMeas;
    
    % The modeled fluxes are isolated
    tmpMod = totalCO2_control(r,1);
    tmpMod_control = tmpMod;
    
    % Dots with positive and negative residuals are given a different color
    res = tmpMod - tmpMeas;
    [r1 c1] = find(res < 0);
    [r2 c2] = find(res > 0);
    
    % The data is saved for export
    CO2flux_negRes_meas_control = tmpMeas(r1);
    CO2flux_negres_mod_control = tmpMod(r1);
    CO2flux_posRes_meas_control = tmpMeas(r2);
    CO2flux_posres_mod_control = tmpMod(r2);
    
    hold on
    xlim([0 6])
    ylim([0 6])
    scatter(tmpMeas(r1), tmpMod(r1), 'MarkerFaceColor', color_negRes, 'MarkerEdgeColor', 'none');
    scatter(tmpMeas(r2), tmpMod(r2), 'MarkerFaceColor', color_posRes, 'MarkerEdgeColor', 'none');    
    hline1 = refline(1,0);
    hline1.Color = [.5 .5 .5];
    hline1.LineStyle = '--';

    % RMSE
    RRMSE = sqrt(mean(((tmpMeas - totalCO2_control(r,1))./tmpMeas).^2));
    RRMSE = round(RRMSE*100)/100;
    T = text(min(get(gca, 'xlim')) + 1, max(get(gca, 'ylim')) - 1, ['RRMSE = ' num2str(RRMSE)]);
    set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

    xlabel('Measured CO_{2} (g C per timestep)')
    ylabel('Modelled CO_{2} (g C per timestep)')
    title('Control treatment')
    
    % Bias
    bias = mean(totalCO2_control(r,1) - tmpMeas);
    bias = round(bias*100)/100;
    T = text(min(get(gca, 'xlim')) + 1, max(get(gca, 'ylim')) - 1.3, ['Bias = ' num2str(bias)]);
    set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

    % ------------
    % Warmed plots
    % -------------
    
    subplot(1,2,2)
    
    % The measurement dates that are modelled are isolated
    tmpMeas = allMeasurements.Measured_CO2_flux_warming;
    tmpMeas_heated = tmpMeas;

    % The modeled fluxes are isolated
    tmpMod = totalCO2_warmed(r,1);
    tmpMod_heated = tmpMod;
    
    % Dots with positive and negative residuals are given a different color
    res = tmpMod - tmpMeas;
    [r1 c1] = find(res < 0);
    [r2 c2] = find(res > 0);
    
     % The data is saved for export
    CO2flux_negRes_meas_warming = tmpMeas(r1);
    CO2flux_negres_mod_warming = tmpMod(r1);
    CO2flux_posRes_meas_warming = tmpMeas(r2);
    CO2flux_posres_mod_warming = tmpMod(r2);
    
    hold on
    xlim([0 6])
    ylim([0 6])
    scatter(tmpMeas(r1), tmpMod(r1), 'MarkerFaceColor', color_negRes, 'MarkerEdgeColor', 'none');
    scatter(tmpMeas(r2), tmpMod(r2), 'MarkerFaceColor', color_posRes, 'MarkerEdgeColor', 'none');
    hline2 = refline(1,0);
    hline2.Color = [.5 .5 .5];
    hline2.LineStyle = '--';

    % The RMSE is calculated
    RRMSE = sqrt(mean(((tmpMeas - totalCO2_warmed(r,1))./tmpMeas).^2));
    RRMSE = round(RRMSE*100)/100;
    T = text(min(get(gca, 'xlim')) + 1, max(get(gca, 'ylim')) - 1, ['RRMSE = ' num2str(RRMSE)]); 
    set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    % Bias
    bias = mean(totalCO2_warmed(r,1) - tmpMeas);
    bias = round(bias*100)/100;
    T = text(min(get(gca, 'xlim')) + 1, max(get(gca, 'ylim')) - 1.3, ['Bias = ' num2str(bias)]);
    set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');

    xlabel('Measured CO_{2} (g C per timestep)')
    ylabel('Modelled CO_{2} (g C per timestep)')
    title('Heated treatment (+ 5 °C)')
end

end

%% The annual CO2 fluxes are calculated and plotted

if plot_annual_CO2 == 1
    
    % -----------------------------------------
    % The annual fluxes are calculated
    % -----------------------------------------
    
    % The model results are isolated
    [r, c] = find(dates_all > date_startWarming);
    dates_after1991 = dates_all(1,c(1:end-1));
    all_CO2_control = CO2_rStrat_litter_control + CO2_kStrat_litter_control + CO2_soil_control;
    CO2_after1991_control = all_CO2_control(c(1:end-1),1);
    all_CO2_warming = CO2_rStrat_litter_warmed + CO2_kStrat_litter_warmed + CO2_soil_warmed;
    CO2_after1991_warming = all_CO2_warming(c(1:end-1),1);
    
    % The daily fluxes are calculated
    % The flux per timestep (dt) is converted to the daily flux
    if variableTimeStep == 1
        dt = d2;
    end
    dailyCO2Flux_modDates_control = CO2_after1991_control./dt;
    dailyCO2Flux_modDates_warming = CO2_after1991_warming./dt;
    
    % These daily fluxes are interpolated to a daily timestep
    days_CO2_forPlotting = dates_after1991(1):caldays(1):dates_after1991(end);
    dailyCO2Flux_control = interp1(dates_after1991,dailyCO2Flux_modDates_control,days_CO2_forPlotting);
    dailyCO2Flux_warming = interp1(dates_after1991,dailyCO2Flux_modDates_warming,days_CO2_forPlotting);
    
    % The total annual flux is calculated
    allYears = days_CO2_forPlotting.Year;
    uniqueYears_CO2diff = unique(allYears);
    uniqueYears_CO2diff(uniqueYears_CO2diff<1991) = [];
    
    annualCO2_control = NaN(size(uniqueYears_CO2diff,2),1);
    annualCO2_warming = NaN(size(uniqueYears_CO2diff,2),1);
    
    % A loop through all years of warming
    for i = 1:size(uniqueYears_CO2diff,2)
        
        % Find the rownumbers for the specific year
        yy = uniqueYears_CO2diff(1,i);
        [r, c] = find(dates_after1991.Year == yy);
        annualCO2_control(i,1) = sum(dailyCO2Flux_control(1,c));
        annualCO2_warming(i,1) = sum(dailyCO2Flux_warming(1,c));
        
    end
    
    % -----------------------------------------
    % The annual fluxes are plotted
    % -----------------------------------------
    
    barColor_control = [54,144,192]./255;
    barColor_warmed = [239,101,72]./255;
    
    markerColor_control = [2,56,88]./255;
    markerColor_warmed = [127,0,0]./255;
    
    hfig = figure;
    set(hfig, 'Units', 'centimeter', 'Position', [3 3 20 20], 'color', [1 1 1])
    
    subplot(2,1,1)
    hold on
    
    % The model results are plotted
    histData = [annualCO2_control annualCO2_warming];
    b1 = bar(uniqueYears_CO2diff,histData);
    set(b1(1),'FaceColor',barColor_control);
    set(b1(2),'FaceColor',barColor_warmed);

    legend(b1, {'Control' 'Warming'})
    legend('boxoff')
    
    ylim([0 1100])
    xlabel('Years')
    ylabel('Annual CO_{2} emissions (g CO_{2}-C yr^{-1})')
    
    % -----------------------------------------
    % The difference in annual fluxes is plotted
    % -----------------------------------------
    
    subplot(2,1,2)
    hold on
    
    color_diff = [54,144,192]./255;
    
    diffAnnualCO2 = annualCO2_warming - annualCO2_control;
    
    if ~strcmp(site,'artificialRun') % In case of the artificial run, no measured data is available
        measuredDiff = measuredAnnualCO2_warming - measuredAnnualCO2_control;
        % The standard deviations is calculated
        stDev_diff = sqrt((allMeasurements.measuredAnnualCO2_stDev_control.^2) + (allMeasurements.measuredAnnualCO2_stDev_warming.^2));
    end
    
    % The modelled data is plotted
    b1 = bar(uniqueYears_CO2diff,diffAnnualCO2);
    set(b1(1), 'FaceColor', color_diff)

    xlabel('Years')
    ylabel('CO_{2} warming - control (g CO_{2}-C yr^{-1})')
    
end

%% The difference in cumulative CO2 emissions between control and warmed plots

color_control = [55,126,184]./255;
color_warming = [228,26,28]./255;

% The model results are isolated
[r c] = find(dates_all > date_startWarming);
dates_after1991 = dates_all(1,c(1:end-1));

% The cumulative CO2 emisions per time step are isolated
all_CO2_control_cumul = Cpools_litter_control(1:end,14) + Cpools_litter_control(1:end,15) + Cpools_soil_control(1:end,9);
CO2_after1991_control_cumul = all_CO2_control_cumul(c(1:end-1),1) - all_CO2_control_cumul(c(1),1);
all_CO2_warming_cumul = Cpools_litter_warming(1:end,14) + Cpools_litter_warming(1:end,15) + Cpools_soil_warming(1:end,9);
CO2_after1991_warming_cumul = all_CO2_warming_cumul(c(1:end-1),1) - all_CO2_warming_cumul(c(1),1);

% The difference between the cumulative emissions from control and warmed plots are calculated
diffCumul = CO2_after1991_warming_cumul - CO2_after1991_control_cumul;

% The cumulative fluxes and their difference is plotted
hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 20 20], 'color', [1 1 1])

% The cumulative fluxes
subplot(2,1,1)
hold on
p1 = plot(dates_after1991, CO2_after1991_control_cumul, 'color', color_control, 'linewidth', 2);
p2 = plot(dates_after1991, CO2_after1991_warming_cumul, 'color', color_warming, 'linewidth', 2);

L1 = legend([p1, p2], 'Control', 'Warming', 'location', 'southeast');
legend('boxoff')

title('Cumulative CO2 fluxes')
xlabel('Year')
ylabel('g CO_{2}')

% Their difference
subplot(2,1,2)
hold on

p1 = plot(dates_after1991, diffCumul, 'color', color_control, 'linewidth', 2);

title('Difference in cumulative fluxes (warmed - control)')
xlabel('Year')
ylabel('g CO_{2}')

%% Plotting CO2 production per microbe

if plot_CO2_per_microbe == 1

% The colors are defined
color_CO2_rStrat_litter = [166,189,219]./255;
color_CO2_kStrat_litter = [43,140,190]./255;
color_CO2_soil = [35,132,67]./255;
color_maint = [255,127,0]./255;
color_growth = [152,78,163]./255;
color_enzyme = [77,175,74]./255;
color_uptake = [55,126,184]./255;
markerColor = [.5 .5 .5];

% The CO2 flux per timestep is calculated
CO2_rStrat_litter_maint_control = diff(Cpools_litter_control(2:end,20));
CO2_rStrat_litter_growth_control = diff(Cpools_litter_control(2:end,22));
CO2_rStrat_litter_enzyme_control = diff(Cpools_litter_control(2:end,24));
CO2_rStrat_litter_uptake_control = diff(Cpools_litter_control(2:end,26));

CO2_rStrat_litter_control_split = [CO2_rStrat_litter_maint_control CO2_rStrat_litter_growth_control ...
    CO2_rStrat_litter_enzyme_control CO2_rStrat_litter_uptake_control];
clear CO2_rStrat_litter_maint_control CO2_rStrat_litter_growth_control ...
    CO2_rStrat_litter_enzyme_control CO2_rStrat_litter_uptake_control;

CO2_kStrat_litter_maint_control = diff(Cpools_litter_control(2:end,21));
CO2_kStrat_litter_growth_control = diff(Cpools_litter_control(2:end,23));
CO2_kStrat_litter_enzyme_control = diff(Cpools_litter_control(2:end,25));
CO2_kStrat_litter_uptake_control = diff(Cpools_litter_control(2:end,27));

CO2_kStrat_litter_control_split = [CO2_kStrat_litter_maint_control CO2_kStrat_litter_growth_control ...
    CO2_kStrat_litter_enzyme_control CO2_kStrat_litter_uptake_control];
clear CO2_kStrat_litter_maint_control CO2_kStrat_litter_growth_control ...
    CO2_kStrat_litter_enzyme_control CO2_kStrat_litter_uptake_control;

CO2_soil_maint_control = diff(Cpools_soil_control(2:end,12));
CO2_soil_growth_control = diff(Cpools_soil_control(2:end,13));
CO2_soil_enzyme_control = diff(Cpools_soil_control(2:end,14));
CO2_soil_uptake_control = diff(Cpools_soil_control(2:end,15));

CO2_soil_control_split = [CO2_soil_maint_control CO2_soil_growth_control ...
    CO2_soil_enzyme_control CO2_soil_uptake_control];
clear CO2_soil_maint_control CO2_soil_growth_control ...
    CO2_soil_enzyme_control CO2_soil_uptake_control;

CO2_rStrat_litter_maint_warming = diff(Cpools_litter_warming(2:end,20));
CO2_rStrat_litter_growth_warming = diff(Cpools_litter_warming(2:end,22));
CO2_rStrat_litter_enzyme_warming = diff(Cpools_litter_warming(2:end,24));
CO2_rStrat_litter_uptake_warming = diff(Cpools_litter_warming(2:end,26));

CO2_rStrat_litter_warming_split = [CO2_rStrat_litter_maint_warming CO2_rStrat_litter_growth_warming ...
    CO2_rStrat_litter_enzyme_warming CO2_rStrat_litter_uptake_warming];
clear CO2_rStrat_litter_maint_warming CO2_rStrat_litter_growth_warming ...
    CO2_rStrat_litter_enzyme_warming CO2_rStrat_litter_uptake_warming;

CO2_kStrat_litter_maint_warming = diff(Cpools_litter_warming(2:end,21));
CO2_kStrat_litter_growth_warming = diff(Cpools_litter_warming(2:end,23));
CO2_kStrat_litter_enzyme_warming = diff(Cpools_litter_warming(2:end,25));
CO2_kStrat_litter_uptake_warming = diff(Cpools_litter_warming(2:end,27));

CO2_kStrat_litter_warming_split = [CO2_kStrat_litter_maint_warming CO2_kStrat_litter_growth_warming ...
    CO2_kStrat_litter_enzyme_warming CO2_kStrat_litter_uptake_warming];
clear CO2_kStrat_litter_maint_warming CO2_kStrat_litter_growth_warming ...
    CO2_kStrat_litter_enzyme_warming CO2_kStrat_litter_uptake_warming;

CO2_soil_maint_warming = diff(Cpools_soil_warming(2:end,12));
CO2_soil_growth_warming = diff(Cpools_soil_warming(2:end,13));
CO2_soil_enzyme_warming = diff(Cpools_soil_warming(2:end,14));
CO2_soil_uptake_warming = diff(Cpools_soil_warming(2:end,15));

CO2_soil_warming_split = [CO2_soil_maint_warming CO2_soil_growth_warming ...
    CO2_soil_enzyme_warming CO2_soil_uptake_warming];
clear CO2_soil_maint_warming CO2_soil_growth_warming ...
    CO2_soil_enzyme_warming CO2_soil_uptake_warming;

% The model results after 1991 are isolated
[r c] = find(dates_all > datetime(1990,1,1));
dates_after1991 = dates_all(1,c(1:end-1));

CO2_rStrat_litter_control_plot = CO2_rStrat_litter_control(c(1:end-1),1);
CO2_rStrat_litter_control_split_plot = CO2_rStrat_litter_control_split(c(1:end-1),:);

CO2_kStrat_litter_control_plot = CO2_kStrat_litter_control(c(1:end-1),1);
CO2_kStrat_litter_control_split_plot = CO2_kStrat_litter_control_split(c(1:end-1),:);

CO2_soil_control_plot = CO2_soil_control(c(1:end-1),1);
CO2_soil_control_split_plot = CO2_soil_control_split(c(1:end-1),:);

CO2_rStrat_litter_warming = CO2_rStrat_litter_warmed(c(1:end-1),1);
CO2_rStrat_litter_warming_split_plot = CO2_rStrat_litter_warming_split(c(1:end-1),:);

CO2_kStrat_litter_warming = CO2_kStrat_litter_warmed(c(1:end-1),1);
CO2_kStrat_litter_warming_split_plot = CO2_kStrat_litter_warming_split(c(1:end-1),:);

CO2_soil_warming = CO2_soil_warmed(c(1:end-1),1);
CO2_soil_warming_split_plot = CO2_soil_warming_split(c(1:end-1),:);

% The daily fluxes are calculated
% The flux per timestep (dt) is converted to the daily flux
daily_CO2_rStrat_litter_control = CO2_rStrat_litter_control_plot./dt;
daily_CO2_rStrat_litter_control_split = CO2_rStrat_litter_control_split_plot./dt;

daily_CO2_kStrat_litter_control = CO2_kStrat_litter_control_plot./dt;
daily_CO2_kStrat_litter_control_split = CO2_kStrat_litter_control_split_plot./dt;

daily_CO2_soil_control = CO2_soil_control_plot./dt;
daily_CO2_soil_control_split = CO2_soil_control_split_plot./dt;

daily_CO2_rStrat_litter_warming = CO2_rStrat_litter_warming./dt;
daily_CO2_rStrat_litter_warming_split = CO2_rStrat_litter_warming_split_plot./dt;

daily_CO2_kStrat_litter_warming = CO2_kStrat_litter_warming./dt;
daily_CO2_kStrat_litter_warming_split = CO2_kStrat_litter_warming_split_plot./dt;

daily_CO2_soil_warming = CO2_soil_warming./dt;
daily_CO2_soil_warming_split = CO2_soil_warming_split_plot./dt;

% These daily fluxes are interpolated to a daily timestep
days_CO2_forPlotting = dates_after1991(1):caldays(1):dates_after1991(end);

dailyCO2Flux_rStrat_control = interp1(dates_after1991,daily_CO2_rStrat_litter_control,days_CO2_forPlotting);
dailyCO2Flux_rStrat_control_split = interp1(dates_after1991,daily_CO2_rStrat_litter_control_split,days_CO2_forPlotting)';

dailyCO2Flux_kStrat_control = interp1(dates_after1991,daily_CO2_kStrat_litter_control,days_CO2_forPlotting);
dailyCO2Flux_kStrat_control_split = interp1(dates_after1991,daily_CO2_kStrat_litter_control_split,days_CO2_forPlotting)';

dailyCO2Flux_soil_control = interp1(dates_after1991,daily_CO2_soil_control,days_CO2_forPlotting);
dailyCO2Flux_soil_control_split = interp1(dates_after1991,daily_CO2_soil_control_split,days_CO2_forPlotting)';

dailyCO2Flux_rStrat_warming = interp1(dates_after1991,daily_CO2_rStrat_litter_warming,days_CO2_forPlotting);
dailyCO2Flux_rStrat_warming_split = interp1(dates_after1991,daily_CO2_rStrat_litter_warming_split,days_CO2_forPlotting)';

dailyCO2Flux_kStrat_warming = interp1(dates_after1991,daily_CO2_kStrat_litter_warming,days_CO2_forPlotting);
dailyCO2Flux_kStrat_warming_split = interp1(dates_after1991,daily_CO2_kStrat_litter_warming_split,days_CO2_forPlotting)';

dailyCO2Flux_soil_warming = interp1(dates_after1991,daily_CO2_soil_warming,days_CO2_forPlotting);
dailyCO2Flux_soil_warming_split = interp1(dates_after1991,daily_CO2_soil_warming_split,days_CO2_forPlotting)';

% The total annual flux is calculated
allYears = days_CO2_forPlotting.Year;
uniqueYears = unique(allYears);
uniqueYears(uniqueYears<1991) = [];

annualCO2_rStrat_control = NaN(size(uniqueYears,2),1);
annualCO2_rStrat_litter_control_split = NaN(size(uniqueYears,2),4);

annualCO2_kStrat_control = NaN(size(uniqueYears,2),1);
annualCO2_kStrat_litter_control_split = NaN(size(uniqueYears,2),4);

annualCO2_soil_control = NaN(size(uniqueYears,2),1);
annualCO2_soil_control_split = NaN(size(uniqueYears,2),4);

annualCO2_rStrat_warming = NaN(size(uniqueYears,2),1);
annualCO2_rStrat_litter_warming_split = NaN(size(uniqueYears,2),4);

annualCO2_kStrat_warming = NaN(size(uniqueYears,2),1);
annualCO2_kStrat_litter_warming_split = NaN(size(uniqueYears,2),4);

annualCO2_soil_warming = NaN(size(uniqueYears,2),1);
annualCO2_soil_warming_split = NaN(size(uniqueYears,2),4);

% A loop through all years of warming
for i = 1:size(uniqueYears,2)

    % Find the rownumbers for the specific year
    yy = uniqueYears(1,i);
    [r c] = find(allYears == yy);
    
    annualCO2_rStrat_control(i,1) = sum(dailyCO2Flux_rStrat_control(1,c));
    annualCO2_rStrat_litter_control_split(i,:) = sum(dailyCO2Flux_rStrat_control_split(:,c),2);
    
    annualCO2_kStrat_control(i,1) = sum(dailyCO2Flux_kStrat_control(1,c));
    annualCO2_kStrat_litter_control_split(i,:) = sum(dailyCO2Flux_kStrat_control_split(:,c),2);
    
    annualCO2_soil_control(i,1) = sum(dailyCO2Flux_soil_control(1,c));
    annualCO2_soil_control_split(i,:) = sum(dailyCO2Flux_soil_control_split(:,c),2);
    
    annualCO2_rStrat_warming(i,1) = sum(dailyCO2Flux_rStrat_warming(1,c));
    annualCO2_rStrat_litter_warming_split(i,:) = sum(dailyCO2Flux_rStrat_warming_split(:,c),2);
    
    annualCO2_kStrat_warming(i,1) = sum(dailyCO2Flux_kStrat_warming(1,c));
    annualCO2_kStrat_litter_warming_split(i,:) = sum(dailyCO2Flux_kStrat_warming_split(:,c),2);
    
    annualCO2_soil_warming(i,1) = sum(dailyCO2Flux_soil_warming(1,c));
    annualCO2_soil_warming_split(i,:) = sum(dailyCO2Flux_soil_warming_split(:,c),2);

end

% -----------------------------------------
% The annual fluxes are plotted - grouped
% -----------------------------------------

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 25 20], 'color', [1 1 1])

% CO2 r-strategists litter - control
subplot(3,2,1)
hold on
b1 = bar(uniqueYears,annualCO2_rStrat_control);
set(b1(1),'FaceColor',color_CO2_rStrat_litter);
title('CO_{2} flux: r-stategists litter - control')
ylim([0 300])

% CO2 K-strategists litter - control
subplot(3,2,3)
hold on
b1 = bar(uniqueYears,annualCO2_kStrat_control);
set(b1(1),'FaceColor',color_CO2_kStrat_litter);
title('CO_{2} flux: K-stategists litter - control')
ylim([0 150])

% CO2 soil - control
subplot(3,2,5)
hold on
b1 = bar(uniqueYears,annualCO2_soil_control);
set(b1(1),'FaceColor',color_CO2_soil);
title('CO_{2} flux: Soil - control')
ylim([0 400])

% CO2 r-strategists litter - warming
subplot(3,2,2)
hold on
b1 = bar(uniqueYears,annualCO2_rStrat_warming);
b2 = bar(uniqueYears,annualCO2_rStrat_control, .2);
set(b1(1),'FaceColor',color_CO2_rStrat_litter);
title('CO_{2} flux: r-stategists litter - warming')
ylim([0 300])

% CO2 K-strategists litter - warming
subplot(3,2,4)
hold on
b1 = bar(uniqueYears,annualCO2_kStrat_warming);
b2 = bar(uniqueYears,annualCO2_kStrat_control, .2);
set(b1(1),'FaceColor',color_CO2_kStrat_litter);
title('CO_{2} flux: K-stategists litter - warming')
ylim([0 150])

% CO2 soil - warming
subplot(3,2,6)
hold on
b1 = bar(uniqueYears,annualCO2_soil_warming);
b2 = bar(uniqueYears,annualCO2_soil_control, .2);
set(b1(1),'FaceColor',color_CO2_soil);
title('CO_{2} flux: Soil - warming')
ylim([0 420])

% -----------------------------------------
% The annual fluxes are plotted - subdivided per flux
% -----------------------------------------

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 25 20], 'color', [1 1 1])

% CO2 r-strategists litter - control
subplot(3,2,1)
hold on
b1 = bar(uniqueYears,annualCO2_rStrat_litter_control_split, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);
title('CO_{2} flux: r-stategists litter - control')
ylim([0 300])

% CO2 K-strategists litter - control
subplot(3,2,3)
hold on
b1 = bar(uniqueYears,annualCO2_kStrat_litter_control_split, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);title('CO_{2} flux: K-stategists litter - control')
ylim([0 150])

% CO2 soil - control
subplot(3,2,5)
hold on
b1 = bar(uniqueYears,annualCO2_soil_control_split, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);title('CO_{2} flux: Soil - control')
ylim([0 400])

% CO2 r-strategists litter - warming
subplot(3,2,2)
hold on
b1 = bar(uniqueYears,annualCO2_rStrat_litter_warming_split, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);title('CO_{2} flux: r-stategists litter - warming')
ylim([0 300])

% CO2 K-strategists litter - warming
subplot(3,2,4)
hold on
b1 = bar(uniqueYears,annualCO2_kStrat_litter_warming_split, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);title('CO_{2} flux: K-stategists litter - warming')
ylim([0 150])

% CO2 soil - warming
subplot(3,2,6)
hold on
b1 = bar(uniqueYears,annualCO2_soil_warming_split, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);title('CO_{2} flux: Soil - warming')
ylim([0 400])

l1 = legend('Maintenance','Growth','Enzyme', 'Uptake');
l1.FontSize = 8;

% -----------------------------------------
% The annual fluxes are plotted - subdivided per flux; normalized fractions
% -----------------------------------------

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 25 20], 'color', [1 1 1])

% CO2 r-strategists litter - control
subplot(3,2,1)
hold on
b1 = bar(uniqueYears,annualCO2_rStrat_litter_control_split./sum(annualCO2_rStrat_litter_control_split,2)*100, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);
title('CO_{2} flux: r-stategists litter - control')
ylim([0 100])

% CO2 K-strategists litter - control
subplot(3,2,3)
hold on
b1 = bar(uniqueYears,annualCO2_kStrat_litter_control_split./sum(annualCO2_kStrat_litter_control_split,2)*100, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);title('CO_{2} flux: K-stategists litter - control')
ylim([0 100])

% CO2 soil - control
subplot(3,2,5)
hold on
b1 = bar(uniqueYears,annualCO2_soil_control_split./sum(annualCO2_soil_control_split,2)*100, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);title('CO_{2} flux: Soil - control')
ylim([0 100])

% CO2 r-strategists litter - warming
subplot(3,2,2)
hold on
b1 = bar(uniqueYears,annualCO2_rStrat_litter_warming_split./sum(annualCO2_rStrat_litter_warming_split,2)*100, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);title('CO_{2} flux: r-stategists litter - warming')
ylim([0 100])

% CO2 K-strategists litter - warming
subplot(3,2,4)
hold on
b1 = bar(uniqueYears,annualCO2_kStrat_litter_warming_split./sum(annualCO2_kStrat_litter_warming_split,2)*100, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);title('CO_{2} flux: K-stategists litter - warming')
ylim([0 100])

% CO2 soil - warming
subplot(3,2,6)
hold on
b1 = bar(uniqueYears,annualCO2_soil_warming_split./sum(annualCO2_soil_warming_split,2)*100, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);title('CO_{2} flux: Soil - warming')
ylim([0 100])

l1 = legend('Maintenance','Growth','Enzyme', 'Uptake');
l1.FontSize = 8;

end

%% Mass specific respiration
% The average amount of microbes per year is calculated

if plot_massSpec_resp == 1

% The model results after 1991 are isolated
[r c] = find(dates_all > datetime(1990,1,1));
dates_after1991 = dates_all(1,c(1:end-1));

rStrat_litter_control_plot = Cpools_litter_control(c(1:end-1),1);
kStrat_litter_control_plot = Cpools_litter_control(c(1:end-1),1);
soilMicrobes_control_plot = Cpools_soil_control(c(1:end-1),1);

rStrat_litter_warming = Cpools_litter_warming(c(1:end-1),1);
kStrat_litter_warming = Cpools_litter_warming(c(1:end-1),1);
soilMicrobes_warming = Cpools_soil_warming(c(1:end-1),1);

% These are interpolated to a daily timestep
days_CO2_forPlotting = dates_after1991(1):caldays(1):dates_after1991(end);

daily_rStrat_control = interp1(dates_after1991,rStrat_litter_control_plot,days_CO2_forPlotting);
daily_kStrat_control = interp1(dates_after1991,kStrat_litter_control_plot,days_CO2_forPlotting);
daily_soilMicrobes_control = interp1(dates_after1991,soilMicrobes_control_plot,days_CO2_forPlotting);

daily_rStrat_warming = interp1(dates_after1991,rStrat_litter_warming,days_CO2_forPlotting);
daily_kStrat_warming = interp1(dates_after1991,kStrat_litter_warming,days_CO2_forPlotting);
daily_soilMicrobes_warming = interp1(dates_after1991,soilMicrobes_warming,days_CO2_forPlotting);

% The average amount of microbes per year is calculated
allYears = days_CO2_forPlotting.Year;
uniqueYears = unique(allYears);
uniqueYears(uniqueYears<1991) = [];

rStrat_control = NaN(size(uniqueYears,2),1);
kStrat_control = NaN(size(uniqueYears,2),1);
soil_control = NaN(size(uniqueYears,2),1);

rStrat_warming = NaN(size(uniqueYears,2),1);
kStrat_warming = NaN(size(uniqueYears,2),1);
soil_warming = NaN(size(uniqueYears,2),1);

% A loop through all years of warming
for i = 1:size(uniqueYears,2)

    % Find the rownumbers for the specific year
    yy = uniqueYears(1,i);
    [r c] = find(allYears == yy);
    
    rStrat_control(i,1) = mean(daily_rStrat_control(1,c));
    kStrat_control(i,1) = mean(daily_kStrat_control(1,c));
    soilMicrobe_control(i,1) = mean(daily_soilMicrobes_control(1,c));
    
    rStrat_warming(i,1) = mean(daily_rStrat_warming(1,c));
    kStrat_warming(i,1) = mean(daily_kStrat_warming(1,c));
    soilMicrobe_warming(i,1) = mean(daily_soilMicrobes_warming(1,c));

end

% -----------------------------------------
% The average annual Rmass is plotted
% -----------------------------------------

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 25 20], 'color', [1 1 1])

% CO2 r-strategists litter - control
subplot(3,2,1)
hold on
b1 = bar(uniqueYears,annualCO2_rStrat_control./rStrat_control);
set(b1(1),'FaceColor',color_CO2_rStrat_litter);
avg_rStrat = mean(annualCO2_rStrat_control./rStrat_control);
plot([uniqueYears(1) uniqueYears(end)],[avg_rStrat avg_rStrat], '--r')
title('Rmass: r-stategists litter - control')
ylim([0 300])

% CO2 K-strategists litter - control
subplot(3,2,3)
hold on
b1 = bar(uniqueYears,annualCO2_kStrat_control./kStrat_control);
set(b1(1),'FaceColor',color_CO2_kStrat_litter);
avg_kStrat = mean(annualCO2_kStrat_control./kStrat_control);
plot([uniqueYears(1) uniqueYears(end)],[avg_kStrat avg_kStrat], '--r')
title('Rmass: K-stategists litter - control')
ylim([0 180])

% CO2 soil - control
subplot(3,2,5)
hold on
b1 = bar(uniqueYears,annualCO2_soil_control./soilMicrobe_control);
set(b1(1),'FaceColor',color_CO2_soil);
avg_soil = mean(annualCO2_soil_control./soilMicrobe_control);
plot([uniqueYears(1) uniqueYears(end)],[avg_soil avg_soil], '--r')
title('Rmass: Soil - control')
ylim([0 400])

% CO2 r-strategists litter - warming
subplot(3,2,2)
hold on
b1 = bar(uniqueYears,annualCO2_rStrat_warming./rStrat_warming);
set(b1(1),'FaceColor',color_CO2_rStrat_litter);
plot([uniqueYears(1) uniqueYears(end)],[avg_rStrat avg_rStrat], '--r')
title('Rmass: r-stategists litter - warming')
ylim([0 300])

% CO2 K-strategists litter - warming
subplot(3,2,4)
hold on
b1 = bar(uniqueYears,annualCO2_kStrat_warming./kStrat_warming);
set(b1(1),'FaceColor',color_CO2_kStrat_litter);
plot([uniqueYears(1) uniqueYears(end)],[avg_kStrat avg_kStrat], '--r')
title('Rmass: K-stategists litter - warming')
ylim([0 180])

% CO2 soil - warming
subplot(3,2,6)
hold on
b1 = bar(uniqueYears,annualCO2_soil_warming./soilMicrobe_warming);
set(b1(1),'FaceColor',color_CO2_soil);
plot([uniqueYears(1) uniqueYears(end)],[avg_soil avg_soil], '--r')
title('Rmass: Soil - warming')
ylim([0 400])

% -----------------------------------------
% The average annual Rmass is plotted: subdivided in different sources
% -----------------------------------------

hfig = figure;
set(hfig, 'Units', 'centimeter', 'Position', [3 3 25 20], 'color', [1 1 1])

% CO2 r-strategists litter - control
subplot(3,2,1)
hold on
b1 = bar(uniqueYears,annualCO2_rStrat_litter_control_split./rStrat_control, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);
avg_rStrat = mean(annualCO2_rStrat_control./rStrat_control);
plot([uniqueYears(1) uniqueYears(end)],[avg_rStrat avg_rStrat], '--r')
title('Rmass: r-stategists litter - control')
ylim([0 300])

% CO2 K-strategists litter - control
subplot(3,2,3)
hold on
b1 = bar(uniqueYears,annualCO2_kStrat_litter_control_split./kStrat_control, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);avg_kStrat = mean(annualCO2_kStrat_control./kStrat_control);
plot([uniqueYears(1) uniqueYears(end)],[avg_kStrat avg_kStrat], '--r')
title('Rmass: K-stategists litter - control')
ylim([0 180])

% CO2 soil - control
subplot(3,2,5)
hold on
b1 = bar(uniqueYears,annualCO2_soil_control_split./soilMicrobe_control, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);avg_soil = mean(annualCO2_soil_control./soilMicrobe_control);
plot([uniqueYears(1) uniqueYears(end)],[avg_soil avg_soil], '--r')
title('Rmass: Soil - control')
ylim([0 100])

% CO2 r-strategists litter - warming
subplot(3,2,2)
hold on
b1 = bar(uniqueYears,annualCO2_rStrat_litter_warming_split./rStrat_warming, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);plot([uniqueYears(1) uniqueYears(end)],[avg_rStrat avg_rStrat], '--r')
title('Rmass: r-stategists litter - warming')
ylim([0 300])

% CO2 K-strategists litter - warming
subplot(3,2,4)
hold on
b1 = bar(uniqueYears,annualCO2_kStrat_litter_warming_split./kStrat_warming, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);plot([uniqueYears(1) uniqueYears(end)],[avg_kStrat avg_kStrat], '--r')
title('Rmass: K-stategists litter - warming')
ylim([0 180])

% CO2 soil - warming
subplot(3,2,6)
hold on
b1 = bar(uniqueYears,annualCO2_soil_warming_split./soilMicrobe_warming, 'stacked');
set(b1(1),'FaceColor', color_maint);
set(b1(2),'FaceColor',color_growth);
set(b1(3),'FaceColor',color_enzyme);
set(b1(4),'FaceColor',color_uptake);plot([uniqueYears(1) uniqueYears(end)],[avg_soil avg_soil], '--r')
title('Rmass: Soil - warming')
ylim([0 100])

l1 = legend('Maintenance','Growth','Enzyme', 'Uptake');
l1.FontSize = 8;

end














