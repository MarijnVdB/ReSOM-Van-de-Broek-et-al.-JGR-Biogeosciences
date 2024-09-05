function par = set_enzyme_par_default_kStrat(scal, opt, par)

% Description:
%   Setup parameters for enzymes
%   Added by XZhu

cpar_litter = par.cpar_litter;
vid_litter = par.vid_litter;

nmicrobes = length(vid_litter.microbep);
nenzymes = length(vid_litter.enzymes);
nsurfaces=length(vid_litter.surfaces);
nmonomers = length(vid_litter.monomers);
npolymers = length(vid_litter.polymers);

if(nargin == 0)
    scal = 1;
end
if(nargin <= 1)
    opt = 1;
end

% - scalars

par.decay_ee   = cpar_litter{9}(2).*scal;       % enzyme decay rate (1/day)

% - vectors

par.Vmax_ee    = cpar_litter{10}(2).*scal;      % 1 x npolymers, maximum som degradation rate (V_Emax) (1/day)


par.Kaff_ee_polymer = cpar_litter{16}(2);       % 1 x npolymers, enzyme affinity to polymer (K_ES) (g ee c/m3)
par.Kaff_ee_msurf = repmat(5d1,1,nsurfaces);    % 1 x nsurfaces, enzyme affinity for adsorptive surface (K_ME) (g ee c/m3)

par.Vmax_ads_enzyme = cpar_litter{11}(1).*scal; % 1 x nsurfaces, maximum enzyme adsorption rate (1/day)
par.adsenz_decay = cpar_litter{12}(1).*scal;    % 1 x nsurfaces, decay rate of adsorbed enzyme (1/day) (ganmma_E)

end