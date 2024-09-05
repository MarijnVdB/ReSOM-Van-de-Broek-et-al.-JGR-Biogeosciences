function par = set_msurface_par_default_litter(par)
%
%Description
%set up parameters for monomer-mineral surface interaction
%author: Jinyun Tang

cpar_litter = par.cpar_litter;
vid_litter = par.vid_litter;

nmicrobes = length(vid_litter.microbep);
nenzymes = length(vid_litter.enzymes);
nsurfaces = length(vid_litter.surfaces);
nmonomers = length(vid_litter.monomers);
npolymers = length(vid_litter.polymers);

par.Kaff_monomer = repmat(25d0,1,nmonomers);    % 1 x nmonomers, adsorption surface doc affinity           (g C/m3)
par.Vmax_ads_monomer = cpar_litter{13};         % 1 x nmonomers, maximum monomer adsorption rate (1/day)
par.adsmon_decay = cpar_litter{14};             % 1 x nmonomers, arbitary, need change
par.adspoly_decay = cpar_litter{14}/100;        % 1 x npolymers, arbitary, need change
end