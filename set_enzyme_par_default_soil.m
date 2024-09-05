function par = set_enzyme_par_default_soil(scal, opt, par)

% Description:
%   Setup parameters for enzymes
%   Added by XZhu

cpar_soil = par.cpar_soil;
vid_soil = par.vid_soil;

nmicrobes=length(vid_soil.microbep);
nenzymes=length(vid_soil.enzymes);
nsurfaces=length(vid_soil.surfaces);
nmonomers=length(vid_soil.monomers);
npolymers=length(vid_soil.polymers);

if(nargin==0)
    scal=1;
end
if(nargin<=1)
    opt=1;
end

% - scalars

par.decay_ee   = cpar_soil{9}.*scal;            % enzyme decay rate          (1/day)

% - vectors

par.Vmax_ee    = cpar_soil{10}.*scal;           % 1 x npolymers, maximum som degradation rate (V_Emax)             (1/day)

par.Kaff_ee_polymer = 5;                        % 1 x npolymers, enzyme affinity to polymer (K_ES)                (g ee c/m3)
par.Kaff_ee_msurf = 0.5;                        % 1 x nsurfaces, enzyme affinity for adsorptive surface (K_ME)    (g ee c/m3)

par.Vmax_ads_enzyme = cpar_soil{11}.*scal;      % 1 x nsurfaces, maximum enzyme adsorption rate (1/day) (??)
par.adsenz_decay = cpar_soil{12}.*scal;         % 1 x nsurfaces, decay rate of adsorbed enzyme (1/day) (ganmma_E)

end