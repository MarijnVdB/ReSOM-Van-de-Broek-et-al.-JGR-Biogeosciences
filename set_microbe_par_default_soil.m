function par=set_microbe_par_default_soil(scal, opt, par)

% Description
%   Setup parameters for microbes
%   Author: Jinyun Tang

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

par.decay_micb0 = cpar_soil{6}.*scal;   % Microbial death rate (1/day)
par.decay_micb1 = 0d0;                  % Scaling parameter to account for mortality change from respiration stress
par.gmax_micb  = cpar_soil{8}.*scal;    % Growth rate (1/day)
par.pro_ee0    = cpar_soil{7}.*scal;    % Maximum enzyme production rate (1/day)
par.pro_ee1    = 0d-6.*scal;            % Inductive enzyme production rate (1/day)
par.Yld_micb   = 0.8;                   % Yield rate of structural biomass from reserve metabolites (unitless)
par.Yld_ee     = 0.8;                   % Yield rate of enzyme from reserve metabolites (unitless)
par.mr_micb    = cpar_soil{3}.*scal;    % Microbial reference maintenance rate (1/day)
par.kappa_micb = cpar_soil{4}.*scal;    % Metabolic turnover rate for plastic microbe (1/day)
par.zb         = 0.05;                  % Scaling factor between transporter and microbial cell biomass
par.micb0      = 1d-4;                  % Half saturating microbial biomass for mortality (gC/m3), I have used an alternative value 1d-3 for sensitivity test

% - vectors

par.Kaff_monomer = repmat(1d0,1,nmonomers);         % 1 x nmonomers, microbial doc affinity                    (g C/m3)
par.Vmax_micb  = cpar_soil{5}.*scal;                % 1 x nmonomers, maximum doc uptake rate                   (1/day)
par.Yld_micx_monomer = repmat(0.5,1,nmonomers);     % 1 x nmonomers, assimilation efficiency from doc/monomer uptake (yield rate of reserve metabolite from DOM)   (g res C/g DOC C)
par.pro_enz_dist = repmat(1/nenzymes,1,nenzymes);   % 1 x nenzymes, distribution of bulk enzyme biomass to each enzyme

% the two parameters below are only for analytic solution test
par.pro_ee = cpar_soil{15}.*scal;
par.decay_mic = par.decay_micb0; 
end