function isgrw=deb_microbe_init(m0,g0,cnp,yld,je,ev)
%
% DESCRIPTION
% initialize the DEB module
% author: Jinyun Tang
% first the metabolic reserve turnover is used to support maintenance
% respiration. If there is excessive carbon flux then it is used to support
% growth and other activity. Were there any stresses, that should be
% effectively change the maintenance requirement
% On input it is assued the maintenance requirement has included cost to
% deal with stress
% input variables:
% m0: maintenance demand
% g0: a list of (nx1) maximum production rates
% cnp: a matrix of (nx2) elemental stoichiometry for the process indicated
% in g0
% yld: a list of (nx1) yielding rate for given processes, in the unit of carbon.
% je: the export metabolic flux to cell structure, before dilution
% correction
% the returning variables should in the order of [netgrowth, activity
% investment]
% determine the number of elements
% return variable, 
nelm=length(je);  
% Compute the excessive carbon 
% je is the amount of C out of reserve, relative to the amount of structural biomass
% m0 is the maintenance respiration, relative to the amount of structural biomass
% dc is the amount of C available from growth and enzyme production,
%   relative to the amount of structural biomass
dc = je(1) - m0;

isgrw=-1;
if(dc>0)
    % There is carbon to support growth activity
    % compute the actual carbon flux to support growth, the yield rate here
    % is less than 1, indicating the fraction of carbon being turned into
    % the required structure after taking off the overhead
    
    jc = dc.*yld; % The actual amount of C available for growth, after accounting for CO2 losses
    
    switch nelm
        case 1               
        %carbon only
           % maximum carbon export
           % ev = ratio of reserve C to MBC
           % g0 contains gmax and pemax
           % dc1 > 1 if maximum growth is achieved, dc < 1 otherwise
           dc1 = dc - g0(1)*ev(1);

           scal_c = dc1/(sum(g0./yld)); 
           if(scal_c>=1)
               % Maximum growth
               isgrw=1;
           else
               % Less than maximum growth
               isgrw=0;
           end
        case 2
        %carbon and nitrogen
           %carbon export
           dc1=dc-g0(1)*ev(1);
           %nitrogen export
           dn1=je(2)-g0(1)*ev(2);
           %compute c-based down-regulation factor
           scal_c=dc1/sum(g0./yld);
           %compute nitrogen based down-regulation factor
           scal_n=dn1/sum(g0./cnp(:,1));
           %compute the actual growth rate        
           scal=min([scal_c,scal_n]);       
           if(scal>=1)
               %maximum growth
               isgrw=1;
           else
               isgrw=0;
           end
        case 3
           %carbon and nitrogen and phosphorus
           %carbon export
           dc1=dc-g0(1)*ev(1);
           %nitrogen export
           dn1=je(2)-g0(1)*ev(2);
           %phosphorus export
           dp1=je(3)-g0(1)*ev(3);

            %compute c-based down-regulation factor
            scal_c=dc1/sum(gp./yld);
            %compute nitrogen based down-regulation factor
            scal_n=dn1/sum(gp./cnp(:,1));
            %compute phosphorus based down-regulation factor
            scal_p=dp1/sum(gp./cnp(:,2));
            %compute the actual growth rate
            scal=min([scal_c,scal_n,scal_p]);       
            
            if(scal>=1)
                %maximum growth
                isgrw=1;
            else
                isgrw=0;
            end
    end
        
end
   
end