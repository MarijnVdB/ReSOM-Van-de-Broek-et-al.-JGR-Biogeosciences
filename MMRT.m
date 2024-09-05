function MMRTfactor = MMRT(par_MMRT, T0, T, Topt)

% The parameters are converted to J mol-1 K-1
Cp = par_MMRT(1).*1e+3;     % Change in heat capacity
dH = par_MMRT(2).*1e+3;     % Enthalpy 
dS = par_MMRT(3).*1e+3;     % Entropy

% Standard parameters
kb = 1.380649e-23;      % Boltzmann's constant (J K-1)
h = 6.62607015e-34;     % Planck's constant (J s)
R = 8.314;              % Universal gas constant (J K-1 mol-1)

% The log rates at T and Topt
logk_T = log((kb.*T)./h) - ((dH + Cp.*(T-T0))./(R.*T)) + ((dS+Cp.*(log(T)-log(T0)))./R);
logk_Topt = log((kb.*Topt)./h) - ((dH + Cp.*(Topt-T0))./(R.*Topt)) + ((dS+Cp.*(log(Topt)-log(T0)))./R);

MMRTfactor = exp(logk_T)./exp(logk_Topt);

end

