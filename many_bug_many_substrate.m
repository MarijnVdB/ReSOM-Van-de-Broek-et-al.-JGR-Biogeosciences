function dxdt=many_bug_many_substrate(t,x,param, tmpDummy)

% DESCRIPTION
% Core file for microbial model with many bugs attacking many substrates
% Author: Jinyun Tang, jinyuntang@lbl.gov
% Created Feb, 9, 2014.
% currently, I only deal with carbon substrate, so the first step carbon
% yield is well defined for each substrate by each microbe

% Adapted by Marijn Van de Broek

vid_soil = param.vid_soil;                  % id index
par_mic_soil = param.par_mic_soil;          % each par represents a microbe and its enzyme
par_enz_soil = param.par_enz_soil;          % each par represents a enzyme
par_surface_soil = param.par_surface_soil;  % each par represents one mineral surface
par_ss_soil = param.par_ss_soil;            % parameter structure for substrate
input_soil = param.input_soil;              % substrate input_soil structure

% Obtain number of pools
nmicrobes = length(vid_soil.microbep);
nenzymes = length(vid_soil.enzymes);
nsurfaces = length(vid_soil.surfaces);
nmonomers = length(vid_soil.monomers);
npolymers = length(vid_soil.polymers);

% Initialize the trend vector
dxdt = zeros(size(x));

% The substrate-microbe relationship is represented as a matrix with each
% column representing substrates, and each row representing microbes (or
% mineral surface) the matrices are K parameters, Vmax parameters and 
% substrate yielding rates polymer degradation needs extracellular enzyme,
% which is produced by producers. In this model, mineral surfaces only 
% adsorb enzymes and monomers but do not actively protect them through
% aggregation.

% Obtain polymer degradation flux (by enzyme, npolymers x nenzymes) and enzyme adsorption flux (by mineral, nsurfaces x nenzymes)
[polymer_degrad_flux, enzyme_adsorp_flux] = polymer_degradation(x(vid_soil.polymers), x(vid_soil.enzymes), ...
    x(vid_soil.surfaces), par_enz_soil);

% Obtain monomer uptake flux (by microbe, nmonomers x nmicrobes) and monomer adsorption flux (by mineral, nmonomers x nsurfaces)
[monomer_uptake_matrix, monomer_adsorp_flux, enz_pro_stress_tot] = monomer_uptake(x(vid_soil.monomers), ...
    x(vid_soil.microbep), x(vid_soil.surfaces), par_mic_soil, par_surface_soil);

% turnover of adsorbed enzyme pool
adsenz_decay_flux = zeros(nenzymes,nsurfaces);
for ii = 1:nenzymes
    for jj = 1:nsurfaces
        adsenz_decay_flux(ii,jj) = par_enz_soil(ii).adsenz_decay(jj)*x(vid_soil.enzymes_ads((ii-1)*nsurfaces+jj));
    end
end

% turnover of adsorbed monomer pool
adsmon_decay_flux = zeros(nmonomers,nsurfaces);
for ii = 1:nmonomers
    for jj = 1:nsurfaces
        adsmon_decay_flux(ii,jj) = par_surface_soil(jj).adsmon_decay(ii)*x(vid_soil.monomers_ads((ii-1)*nsurfaces+jj));
    end
end

%now assign the effective carbon yield flux to each microbes, nmonomers x nmicrobes
[monomer_to_microbe_flux, assim_co2_flux] = monomer_microbe_yield(monomer_uptake_matrix(:,1:nmicrobes),par_mic_soil,x);

% Now do cell level metabolism to find rates of maintenance, population
% growth and enzyme production
[maintenance, growth, pro_enz, cell_death, co2_metab_flux, uptake] = microbe_metabolism(monomer_to_microbe_flux,...
    enz_pro_stress_tot, x, vid_soil, par_mic_soil, 1);

%calculate CUE
for jj = 1 : nmicrobes
    cue(jj) = 1 - ((co2_metab_flux(jj) + assim_co2_flux(jj))./(monomer_to_microbe_flux/par_mic_soil.Yld_micx_monomer));
end

global value_CUE_soil
value_CUE_soil = cue(1);

% do mass update

% microbial population
for jj = 1 : nmicrobes
    if(plasticity_test(par_mic_soil(jj))) % need change
        %rigid microbe
        dxdt(vid_soil.microbep(jj))=(growth(jj)-...
            cell_death(jj))*x(vid_soil.microbep(jj)); % !! add 13C and 14C stuff for rigid microbe
    else
        dxdt(vid_soil.microbep(jj))=(growth(jj)-cell_death(jj))*x(vid_soil.microbep(jj));     
        dxdt(vid_soil.micc(jj))=sum(monomer_to_microbe_flux(:,jj),1)-(par_mic_soil(jj).kappa_micb...
            -growth(jj)+cell_death(jj))*x(vid_soil.micc(jj));        
    end
end

% Enzyme production matrix
pro_enz_bulk=pro_enz.*(x(vid_soil.microbep))'; % Bulk enzyme production
pro_enz_matrix=zeros(nenzymes,nmicrobes);
for jj=1:nmicrobes
    pro_enz_matrix(:,jj)=pro_enz_bulk(jj)*par_mic_soil(jj).pro_enz_dist';
end

% Bulk microbial mortality
cell_death_bulk_microbep=cell_death.*(x(vid_soil.microbep))';
cell_death_bulk_micc=cell_death.*(x(vid_soil.micc))';

% update extracellular enzyme
[dxdt(vid_soil.enzymes),enz_decay]=enzyme_update(x,pro_enz_matrix,vid_soil, par_enz_soil,enzyme_adsorp_flux,adsenz_decay_flux);

% update adsorbed enzymes by mineral surfaces
dxdt(vid_soil.enzymes_ads)=reshape(enzyme_adsorp_flux-adsenz_decay_flux',nsurfaces*nenzymes,1);

% update polymers
dxdt(vid_soil.polymers)=input_soil.polymers ...
    -sum(polymer_degrad_flux,2)' ...
    +(par_ss_soil.deadmicrobep2polymers*cell_death_bulk_microbep)' ...
    +(par_ss_soil.deadmicc2polymers*cell_death_bulk_micc)' ...
    +(par_ss_soil.deadenz2polymers*enz_decay)';

% update monomers
dxdt(vid_soil.monomers)=input_soil.monomers ...
    +(par_ss_soil.polymer2monomer*sum(polymer_degrad_flux,2))' ...
    -sum(monomer_uptake_matrix,2)' ...
    +(par_ss_soil.deadmicrobep2monomers*cell_death_bulk_microbep)' ...
    +(par_ss_soil.deadmicc2monomers*cell_death_bulk_micc)' ...
    +(par_ss_soil.deadenz2monomers*enz_decay)' ...
    -sum(monomer_adsorp_flux,2)' ...
    +sum(adsmon_decay_flux,2)';

% update adsorped monomers
dxdt(vid_soil.monomers_ads)=reshape(monomer_adsorp_flux'-adsmon_decay_flux',nsurfaces*nmonomers,1);

% update mineral surfaces (effective binding surfaces)
dxdt(vid_soil.surfaces)= sum(adsmon_decay_flux,1) ...   % Increase in free surfaces through decaying monomers
    + sum(adsenz_decay_flux,1) ...                      % Increase in free surfaces through decaying enzymes
    -sum(enzyme_adsorp_flux,2)' ...                     % Decrease in free surfaces through adsorbed enzymes
    -sum(monomer_adsorp_flux,1);                        % Decrease in free surfaces through adsorbed monomers

% update co2
dxdt(vid_soil.co2)=sum(assim_co2_flux(:))+sum(co2_metab_flux,2);

dxdt(vid_soil.maintCO2) = maintenance*x(vid_soil.microbep);
dxdt(vid_soil.growthCO2) = growth*x(vid_soil.microbep)*((1-par_mic_soil.Yld_micb)/par_mic_soil.Yld_micb);
dxdt(vid_soil.enzymeCO2) = pro_enz*x(vid_soil.microbep)*((1-par_mic_soil.Yld_ee)/par_mic_soil.Yld_ee);
dxdt(vid_soil.uptakeCO2) = assim_co2_flux;

 for jj = 1 : nmicrobes
% update cue
         dxdt(vid_soil.cue(jj))=cue(jj);
 end

% do mass balance check
[residual]=mass_balance_check(dxdt,x, param);

if(abs(residual)>1d-10)
    fprintf('residual=%f\n',residual);
    error('mass balance error');
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [residual]=mass_balance_check(dxdt,x, par)
%do mass balance check
%     global input_soil vid_soil
    
    input_soil = par.input_soil;
    vid_soil = par.vid_soil;

    idx=[vid_soil.microbep vid_soil.micc vid_soil.monomers vid_soil.monomers_ads vid_soil.polymers vid_soil.enzymes vid_soil.enzymes_ads vid_soil.co2];
    residual=sum(dxdt(idx))-(sum(input_soil.polymers)+sum(input_soil.monomers));
    
end





