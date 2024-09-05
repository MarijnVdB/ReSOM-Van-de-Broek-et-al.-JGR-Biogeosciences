function dxdt = ReSOM_litter(t,x,param)

% DESCRIPTION
% Core file for microbial model with many bugs attacking many substrates
% Author: Jinyun Tang, jinyuntang@lbl.gov
% Created Feb, 9, 2014.
% currently, I only deal with carbon substrate, so the first step carbon
% yield is well defined for each substrate by each microbe

% Adapted by Marijn Van de Broek to simulate an organic horizon

vid_litter = param.vid_litter;                             % id index
par_mic_litter_rStrat = param.par_mic_litter_rStrat;       % r-strategists
par_mic_litter_kStrat = param.par_mic_litter_kStrat;       % K-strategists
par_enz_litter_rStrat = param.par_enz_litter_rStrat;       % Enzymes from r-strategists
par_enz_litter_kStrat = param.par_enz_litter_kStrat;       % Enzymes from K-strategists
par_ss_litter = param.par_ss_litter;                       % parameter structure for substrate
input_litter = param.input_litter;                         % substrate input_litter structure
par_surface_litter = param.par_surface_litter;             % each par represents one mineral surface

% Obtain number of pools
nmicrobes = length(vid_litter.microbep);
nenzymes = length(vid_litter.enzymes);
nsurfaces = length(vid_litter.surfaces);
nmonomers = length(vid_litter.monomers);
npolymers = length(vid_litter.polymers);

% Initialize the trend vector
dxdt = zeros(size(x));

% The substrate-microbe relationship is represented as a matrix with each
% column representing substrates, and each row representing microbes (or
% mineral surface) the matrices are K parameters, Vmax parameters and 
% substrate yielding rates polymer degradation needs extracellular enzyme,
% which is produced by producers. In this model, mineral surfaces only 
% adsorb enzymes and monomers but do not actively protect them through
% aggregation.

% Obtain polymer degradation flux for metabolic litter
[polymer_degrad_flux_metabolic, enzyme_adsorp_flux_rStrat] = polymer_degradation(x(vid_litter.polymers(1)), x(vid_litter.enzymes(1)), ...
    x(vid_litter.surfaces), par_enz_litter_rStrat);

% Obtain polymer degradation flux for structural litter
[polymer_degrad_flux_structural, enzyme_adsorp_flux_kStrat] = polymer_degradation(x(vid_litter.polymers(2)), x(vid_litter.enzymes(2)), ...
    x(vid_litter.surfaces), par_enz_litter_kStrat);

% The parameters for r- and K-strategists are temporarily merged
par_allMicrobes = cell2struct(cellfun(@horzcat,struct2cell(par_mic_litter_rStrat),struct2cell(par_mic_litter_kStrat),'uni',0),fieldnames(par_mic_litter_rStrat),1);

% Obtain monomer uptake flux (and monomer adsorption flux) by r-strategists
[monomer_uptake_matrix, monomer_adsorp_flux, enz_pro_stress_tot] = monomer_uptake(x(vid_litter.monomers), ...
    x(vid_litter.microbep([1 2])), x(vid_litter.surfaces), par_allMicrobes, par_surface_litter);

% The monomer uptake is assigned to the r- and K-strategists
monomer_uptake_matrix_rStrat = monomer_uptake_matrix(1,1);
monomer_uptake_matrix_kStrat = monomer_uptake_matrix(1,2);

enz_pro_stress_tot_rStrat = enz_pro_stress_tot(1,1);
enz_pro_stress_tot_kStrat = enz_pro_stress_tot(1,2);

% Turnover (desorption) of adsorbed enzyme pool for r-strategists
adsenz_decay_flux_rStrat = zeros(1,nsurfaces);
for ii = 1:1
    for jj = 1:nsurfaces
        adsenz_decay_flux_rStrat(ii,jj) = par_enz_litter_rStrat(ii).adsenz_decay(jj)*x(vid_litter.enzymes_ads(1));
    end
end

% Turnover (desorption) of adsorbed enzyme pool for r-strategists
adsenz_decay_flux_kStrat = zeros(1,nsurfaces);
for ii = 1:1
    for jj = 1:nsurfaces
        adsenz_decay_flux_kStrat(ii,jj) = par_enz_litter_kStrat(ii).adsenz_decay(jj)*x(vid_litter.enzymes_ads(2));
    end
end

% turnover of adsorbed monomer pool
adsmon_decay_flux = zeros(nmonomers,nsurfaces);
for ii = 1:nmonomers
    for jj = 1:nsurfaces
        adsmon_decay_flux(ii,jj) = par_surface_litter(jj).adsmon_decay(ii)*x(vid_litter.monomers_ads((ii-1)*nsurfaces+jj));
    end
end

% Now assign the effective carbon yield flux to each microbes: r-strategists
[monomer_to_microbe_flux_rStrat, assim_co2_flux_rStrat] = monomer_microbe_yield(monomer_uptake_matrix_rStrat,par_mic_litter_rStrat,x);

% Now assign the effective carbon yield flux to each microbes: K-strategists
[monomer_to_microbe_flux_kStrat, assim_co2_flux_kStrat] = monomer_microbe_yield(monomer_uptake_matrix_kStrat,par_mic_litter_kStrat,x);

% Now do cell level metabolism to find rates of maintenance, population
% growth and enzyme production for r-strategists
[maintenance_rStrat, growth_rStrat, pro_enz_rStrat, cell_death_rStrat, co2_metab_flux_rStrat, uptake_rStrat] = microbe_metabolism(monomer_to_microbe_flux_rStrat,...
    enz_pro_stress_tot_rStrat, x, vid_litter, par_mic_litter_rStrat, 1);

% Now do cell level metabolism to find rates of maintenance, population
% growth and enzyme production for K-strategists
[maintenance_kStrat, growth_kStrat, pro_enz_kStrat, cell_death_kStrat, co2_metab_flux_kStrat, uptake_kStrat] = microbe_metabolism(monomer_to_microbe_flux_kStrat,...
    enz_pro_stress_tot_kStrat, x, vid_litter, par_mic_litter_kStrat, 2);

% Calculate CUE for r-strategists
cue_rStrat = 1 - ((co2_metab_flux_rStrat + assim_co2_flux_rStrat)./(monomer_to_microbe_flux_rStrat/par_mic_litter_rStrat.Yld_micx_monomer));
global value_CUE_litter_rStrat
value_CUE_litter_rStrat = cue_rStrat(1);

% Calculate CUE for K-strategists
cue_kStrat = 1 - ((co2_metab_flux_kStrat + assim_co2_flux_kStrat)./(monomer_to_microbe_flux_kStrat/par_mic_litter_kStrat.Yld_micx_monomer));
global value_CUE_litter_kStrat
value_CUE_litter_kStrat = cue_kStrat(1);

jj = 1;
dxdt(vid_litter.microbep(jj)) = (growth_rStrat-cell_death_rStrat)*x(vid_litter.microbep(jj));     
dxdt(vid_litter.micc(jj)) = sum(monomer_to_microbe_flux_rStrat,1)-(par_mic_litter_rStrat.kappa_micb...
            -growth_rStrat(jj)+cell_death_rStrat(jj))*x(vid_litter.micc(jj));

% K-strategists
jj = 2;
dxdt(vid_litter.microbep(jj)) = (growth_kStrat-cell_death_kStrat)*x(vid_litter.microbep(jj));     
dxdt(vid_litter.micc(jj)) = sum(monomer_to_microbe_flux_kStrat,1)-(par_mic_litter_kStrat.kappa_micb...
            -growth_kStrat+cell_death_kStrat)*x(vid_litter.micc(jj));

% Enzyme production matrix for enzymes by r-strategists
jj = 1;
pro_enz_bulk_rStrat = pro_enz_rStrat.*(x(vid_litter.microbep(jj)))'; % Bulk enzyme production
pro_enz_matrix_rStrat = pro_enz_bulk_rStrat*par_mic_litter_rStrat.pro_enz_dist(1)';

% Enzyme production matrix for enzymes by K-strategists
jj = 2;
pro_enz_bulk_kStrat = pro_enz_kStrat.*(x(vid_litter.microbep(jj)))'; % Bulk enzyme production
pro_enz_matrix_kStrat = pro_enz_bulk_kStrat*par_mic_litter_kStrat.pro_enz_dist(2)';

% Bulk microbial mortality: r-strategists
cell_death_bulk_microbep_rStrat = cell_death_rStrat.*(x(vid_litter.microbep(1)))';
cell_death_bulk_micc_rStrat = cell_death_rStrat.*(x(vid_litter.micc(1)))';

% Bulk microbial mortality: K-strategists
cell_death_bulk_microbep_kStrat = cell_death_kStrat.*(x(vid_litter.microbep(2)))';
cell_death_bulk_micc_kStrat = cell_death_kStrat.*(x(vid_litter.micc(2)))';

enz_decay_rStrat = par_enz_litter_rStrat.decay_ee*x(vid_litter.enzymes(1));
% the net enzyme flux
enz_flux_rStrat = sum(pro_enz_matrix_rStrat,2) - enz_decay_rStrat - sum(enzyme_adsorp_flux_rStrat,1)' + sum(adsenz_decay_flux_rStrat,2); % nenzymes x 1
dxdt(vid_litter.enzymes(1)) = enz_flux_rStrat;

% update extracellular enzyme: produced by K-strategists: NOT calculated in
% a function like in the soil
nenzymes = 1;

enz_decay_kStrat = par_enz_litter_kStrat.decay_ee*x(vid_litter.enzymes(2));
% the net enzyme flux
enz_flux_kStrat = sum(pro_enz_matrix_kStrat,2) - enz_decay_kStrat - sum(enzyme_adsorp_flux_kStrat,1)' + sum(adsenz_decay_flux_kStrat,2); % nenzymes x 1
dxdt(vid_litter.enzymes(2)) = enz_flux_kStrat;

% update adsorbed enzymes by mineral surfaces: r-strategists
dxdt(vid_litter.enzymes_ads(1)) = reshape(enzyme_adsorp_flux_rStrat - adsenz_decay_flux_rStrat',nsurfaces*nenzymes,1);

% update adsorbed enzymes by mineral surfaces: K-strategists
dxdt(vid_litter.enzymes_ads(2)) = reshape(enzyme_adsorp_flux_kStrat - adsenz_decay_flux_kStrat',nsurfaces*nenzymes,1);

% update polymers: metabolic litter
dxdt(vid_litter.polymers(1)) = input_litter.polymers(1) ...
    - sum(polymer_degrad_flux_metabolic,2)' ...
    + (par_ss_litter.deadmicrobep2polymers(1)*cell_death_bulk_microbep_rStrat)' ...
    + (par_ss_litter.deadmicrobep2polymers(1)*cell_death_bulk_microbep_kStrat)' ...
    + (par_ss_litter.deadmicc2polymers(1)*cell_death_bulk_micc_rStrat)' ...
    + (par_ss_litter.deadmicc2polymers(1)*cell_death_bulk_micc_kStrat)' ...
    + (par_ss_litter.deadenz2polymers(1)*enz_decay_rStrat)' + ...
    + (par_ss_litter.deadenz2polymers(1)*enz_decay_kStrat)';

% update polymers: structural litter
dxdt(vid_litter.polymers(2)) = input_litter.polymers(2) ...
    - sum(polymer_degrad_flux_structural,2)' ...
    + (par_ss_litter.deadmicrobep2polymers(2)*cell_death_bulk_microbep_rStrat)' ...
    + (par_ss_litter.deadmicrobep2polymers(2)*cell_death_bulk_microbep_kStrat)' ...
    + (par_ss_litter.deadmicc2polymers(2)*cell_death_bulk_micc_rStrat)' ...
    + (par_ss_litter.deadmicc2polymers(2)*cell_death_bulk_micc_kStrat)' ...
    + (par_ss_litter.deadenz2polymers(2)*enz_decay_rStrat)' + ...
    + (par_ss_litter.deadenz2polymers(2)*enz_decay_kStrat)';

% update monomers
dxdt(vid_litter.monomers) = input_litter.monomers ...
    + (par_ss_litter.polymer2monomer(1)*sum(polymer_degrad_flux_metabolic,2))' ...
    + (par_ss_litter.polymer2monomer(2)*sum(polymer_degrad_flux_structural,2))' ...
    - sum(monomer_uptake_matrix_rStrat,2)' ...
    - sum(monomer_uptake_matrix_kStrat,2)' ...
    + (par_ss_litter.deadmicrobep2monomers*cell_death_bulk_microbep_rStrat)' ...
    + (par_ss_litter.deadmicrobep2monomers*cell_death_bulk_microbep_kStrat)' ...
    + (par_ss_litter.deadmicc2monomers*cell_death_bulk_micc_rStrat)' ...
    + (par_ss_litter.deadmicc2monomers*cell_death_bulk_micc_kStrat)' ...
    + (par_ss_litter.deadenz2monomers*enz_decay_rStrat)' ...
    + (par_ss_litter.deadenz2monomers*enz_decay_kStrat)' ...
    - sum(monomer_adsorp_flux,2)' ...
    + sum(adsmon_decay_flux,2)';

% update adsorped monomers
dxdt(vid_litter.monomers_ads) = reshape(monomer_adsorp_flux' - adsmon_decay_flux',nsurfaces*nmonomers,1);

% update mineral surfaces (effective binding surfaces)
dxdt(vid_litter.surfaces) =  sum(adsmon_decay_flux,1) ...    % Increase in free surfaces through decaying monomers
    + sum(adsenz_decay_flux_rStrat,1) ...                  % Increase in free surfaces through decaying enzymes
    + sum(adsenz_decay_flux_kStrat,1) ...                  % Increase in free surfaces through decaying enzymes
    - sum(enzyme_adsorp_flux_rStrat,2)' ...                 % Decrease in free surfaces through adsorbed enzymes
    - sum(enzyme_adsorp_flux_kStrat,2)' ...                 % Decrease in free surfaces through adsorbed enzymes
    - sum(monomer_adsorp_flux,1);                    % Decrease in free surfaces through adsorbed monomers

% update co2 for r-strategists
jj = 1;
dxdt(vid_litter.co2(1)) = sum(assim_co2_flux_rStrat(:)) + sum(co2_metab_flux_rStrat,2);
dxdt(vid_litter.maintCO2(1)) = maintenance_rStrat*x(vid_litter.microbep(jj));
dxdt(vid_litter.growthCO2(1)) = growth_rStrat*x(vid_litter.microbep(jj))*((1-par_mic_litter_rStrat(jj).Yld_micb)/par_mic_litter_rStrat(jj).Yld_micb);
dxdt(vid_litter.enzymeCO2(1)) = pro_enz_rStrat*x(vid_litter.microbep(jj))*((1-par_mic_litter_rStrat(jj).Yld_ee)/par_mic_litter_rStrat(jj).Yld_ee);
dxdt(vid_litter.uptakeCO2(1)) = assim_co2_flux_rStrat;

% update co2 for r-strategists
jj = 2;
dxdt(vid_litter.co2(2)) = sum(assim_co2_flux_kStrat(:)) + sum(co2_metab_flux_kStrat,2);
dxdt(vid_litter.maintCO2(2)) = maintenance_kStrat*x(vid_litter.microbep(jj));
dxdt(vid_litter.growthCO2(2)) = growth_kStrat*x(vid_litter.microbep(jj))*((1-par_mic_litter_kStrat.Yld_micb)/par_mic_litter_kStrat.Yld_micb);
dxdt(vid_litter.enzymeCO2(2)) = pro_enz_kStrat*x(vid_litter.microbep(jj))*((1-par_mic_litter_kStrat.Yld_ee)/par_mic_litter_kStrat.Yld_ee);
dxdt(vid_litter.uptakeCO2(2)) = assim_co2_flux_kStrat;

%  for jj = 1 : nmicrobes
% update cue

% r-strategists
dxdt(vid_litter.cue(1)) = cue_rStrat;

% r-strategists
dxdt(vid_litter.cue(2)) = cue_kStrat;

%  end

% do mass balance check
[residual] = mass_balance_check_litter(dxdt,x, param);

if(abs(residual)>1d-10)
    fprintf('residual = %f\n',residual);
    error('mass balance error');
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [residual] = mass_balance_check_litter(dxdt,x, par)
%do mass balance check
%     global input_litter vid_litter

    vid_litter = par.vid_litter;
    input_litter = par.input_litter;

    idx = [vid_litter.microbep vid_litter.micc vid_litter.monomers vid_litter.monomers_ads ...
        vid_litter.polymers vid_litter.enzymes vid_litter.enzymes_ads vid_litter.co2];
    
    residual = sum(dxdt(idx)) - (sum(input_litter.polymers) + sum(input_litter.monomers));
    
end





