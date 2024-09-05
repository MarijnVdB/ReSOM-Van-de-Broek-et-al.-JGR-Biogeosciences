function y=mbbks1(odefun,y0,t,dt,param)

% First order mbbks ode integrator
dydt=feval(odefun,t,y0,param);

% Identify pmax
nJ=0;
neq=length(dydt);
pmax=0.0;
aj=zeros(size(dydt));
for jj = 1 : neq
    if(dydt(jj)<0.)
        nJ=nJ+1;
        pm=-y0(jj)/(dydt(jj)*dt);
        aj(nJ)=-1./pm;
        if(nJ==1)
            pmax=pm;
        else
            pmax=min(pm,pmax);
        end
    end
end

if(nJ>0)
    pmax=min(1.0,pmax.^nJ);
    %solve
    p=GetGdtScalar(aj,nJ,pmax);
    p=p.^(1./nJ);
    y=y0+dydt.*dt.*p;
else
    p=1.0;
    y=y0+dydt.*dt;
end

end