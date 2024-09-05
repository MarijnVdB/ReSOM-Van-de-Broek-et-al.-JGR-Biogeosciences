function y = adptmbbks1(odefun,y0,t,dt,param)

% first order adaptive mbbks ode integrator
% the step controller is implemented as
% (1) find coarse solution, uc
% (2) find fine solution, uf
% (3) take the difference and evaluate rerr=abs(uf-uc).^2/abs(uf).^2

y = y0;             % The pools at the start of the iteration
dt2 = dt;           % The time step
dtmin = dt/2.^5;    % ???
dtr = dt;           % The time step
tt = t;             % The current time

while(1)
    
    if(dt2<=dtmin)
        y = mbbks1(odefun,y,tt,dt2,param);
        dtr=dtr-dt2;
        tt=tt+dt2;

    else

        %find coarse solution
        yc=mbbks1(odefun,y,tt,dt2,param);
        
        %find fine solution
        dt05=dt2/2;
        tt2=tt+dt05;
        yf=mbbks1(odefun,y,tt,dt05,param);
        yf=mbbks1(odefun,yf,tt2,dt05,param);
        
        %determine relative error
        rerr=max(abs(yc-yf)./(abs(yf)+eps));
        [dt_scal,acc]=get_tscal(rerr);
    
        if(acc)
            dtr=dtr-dt2;
            tt=tt+dt2;
            dt2=dt2*dt_scal;
            y=yf;               
        else
            dt2=dt2*dt_scal;
        end
        dt2=min(dt2,dtr);
    end

    if(abs(dtr/dt)<1.d-3)
        break;
    end
end

end

function [dt_scal,acc]=get_tscal(rerr)
%step controller
%if(rerr<0.5rdif)
%dt=dt*2
%elseif(rerr<rdif)
%dt=dt
%elseif(rerr<2*rdif)
%dt=dt/2
%else
%retry
%endif
rerrdif=1.d-4;

if(rerr<0.5*rerrdif)
    dt_scal=2;
    acc=1;
elseif(rerr<rerrdif)
    dt_scal=1;
    acc=1;
elseif(rerr<2*rerrdif)
    dt_scal=0.5;
    acc=1;
else
    dt_scal=0.5;
    acc=0;
end
end