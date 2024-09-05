function siej=ECAcomplex(kd,ss,ee)
% ECAcomplex(kd,ss,ee)
% ECA kinetics
% reference: Tang and Riley, BG, 2013. 
% I, number of different substrates (e.g. polymeric C and minerals)
% J, number of different binding sites (e.g. enzymes)

[I,J]=size(kd);     % kd is the K Matrix
siej=zeros(I,J);    % To store the output Ci,j
dnrm=zeros(I,J);    % The denominator(s) of the ECA equation

for i = 1 : I                           % For every substrate
    dnm1 = 0.0;                         % The initial denominator value is 0
    for k = 1 : J                       % For every binding site
        if(kd(i,k)>0)                   % If the K value > 0
            dnm1=dnm1+ee(k)/kd(i,k);    % The denoninator associated with binding site(s)
        end
    end
    
    for j = 1 : J                       % For every binding site
        dnm2=0.0;                       % The initial denominator value is 0
        if(kd(i,j)>0)                   % If the K value > 0
            for k = 1 : I               % For every substrate
                if(kd(k,j)>0.0)                
                    dnm2=dnm2+ss(k)/kd(k,j); % The denominator associated with the substrate(s)
                end
            end    
            dnrm(i,j)=kd(i,j)*(1+dnm1+dnm2); % The total denominator is calculated
            siej(i,j)=ss(i)*ee(j)/dnrm(i,j); % The amount of complex is calculated
        end
    end
end

end


