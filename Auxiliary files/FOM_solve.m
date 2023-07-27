function [yFOM,xFOM,FOM_time,delta,omega]=FOM_solve(f,Jack,tspan,x0,solver,n_oc)
xdotNL= @(t,x) f(x); % ode
optionsFOM = odeset('RelTol',1e-6,'AbsTol',1e-9);
optionsNL.solver = 'fsolve';   %'fsolve','NewtonRaphson'
optionsNL.funJacobian=Jack;
optionsNL.RelTol = 1e-6;
optionsNL.AbsTol = 1e-9;
tic
switch(solver)
    case 'ode15s'
       [tFOM,xFOM]=ode15s(xdotNL,tspan,x0,optionsFOM);
    case 'Euler'
        [tFOM,xFOM]=implicitEuler(xdotNL,tspan,x0,optionsNL);
end
FOM_time=toc
delta=xFOM(:,1:n_oc);
omega=xFOM(:,n_oc+1:end);
%if n_oc ==300
%[v,val]=max(sum(delta));
%delta=delta-delta(:,val);
%avg_del=sum(delta,2)/n_oc;   %1st state-variable  
%avg_omega=sum(omega,2)/n_oc; %2st state-variable 
%else
avg_del=sum(delta,2)/n_oc;   %1st state-variable 
avg_omega=sum(omega,2)/n_oc; %2st state-variable 
%end
yFOM(:,1)=avg_del;
yFOM(:,2)=avg_omega;
end