function theta_op = calculate_theta_op

% calculate the optimum water content (theta_op) using equation (11) in Methods

%% parameter initialization (parameter values change with soil characteristics)

    phi = 0.76;                     % Soil porosity (-), calculated for HF
    c = 0.1;  %  g/g                % Clay content (g/g)
    m1 = 1.5;                       % Cementation exponent for solute diffusion (-)
    n1 = 2;                         % Saturation exponent for solute diffusion (-)
    m2 = 1.5;                       % Cementation exponent for gas diffusion (-)
    n2 = 2.5;                       % Saturation exponent for gas diffusion (-)
    nu_DO = 2.45; % g/g             % Stoichiometric coefficient of DO (g/g)
    K_theta = 0.1; % m^3/m^3        % Moisture constant (m3/m3)
    alpha = 2e-8; % 1/s             % Mass transfer coefficient (s-1)
    rho_s = 2.65e3;  % kg/m^3       % Density of soil mineral (kg/m3)
    H = 0.2; % m                    % Thickness of layer?
    C_soc = 0.05; % g/g             % OC concentration (g/g; %/100)
    D_GO = 2.1e-5; % m^2/s          % Diffusion coefficient of gaseous O2 (m2/s)

%% calculate SOC-microorganism collocation factor (a) according to clay content (c)

% for heterogenous soils
    if c<0.036 
        a = 0;
    elseif c > 0.34
        a = 1;
    else
        a = 3.31*c - 0.12;
    end

% for homogeneous soils
%     a = 0;
    
%% parameter b is assumed constant
% b = 1.7 for homogeneous soils
b = 0.75; % for heterogenous soils
%     b = 1.7;

%% calculate total mass of SOC per unit surface area (m_soc)
% m_soc can also be given directly
m_soc = rho_s*(1-phi)*H*C_soc;


%% calculate optimal water content (theta_op)
    function f=symfun(x)
        f = nu_DO*x/(K_theta+x)*alpha*m_soc*phi^(a*(m1-n1))*x^(a*n1)-(1.72*C_soc+0.065)*phi^(m2-n2)*(phi-x)^b*D_GO;
    end

    theta_op = bisection(@symfun,0,phi,1e-9);
        
end

%% sub-function 

function p = bisection(f,a,b,eps)

if f(a)*f(b)>0 
    disp('Wrong choice bro')
else
    p = (a + b)/2;
    err = abs(f(p));
    while err > eps
   if f(a)*f(p)<0 
       b = p;
   else
       a = p;          
   end
    p = (a + b)/2; 
   err = abs(f(p));
    end
end
    
end