function [yROM,ROM_time,deltar,omegar]=ROM_solve_DMD(f,Jack,tspan,x0,solver,V,n_oc)
xdotr= @(t,x) f(t);
x0r=V'*x0;
optionsROM = odeset('RelTol',1e-6,'AbsTol',1e-9);
optionsNL.solver = 'fsolve';   %'fsolve','NewtonRaphson'
optionsNL.funJacobian=Jack;
optionsNL.RelTol = 1e-6;
optionsNL.AbsTol = 1e-9;
tic
switch(solver)
    case 'ode15s'
       [tFOM,xROM]=ode15s(xdotr,tspan,x0r,optionsROM);
    case 'Euler'
        [tFOM,xROM]=implicitEuler(xdotr,tspan,x0r,optionsNL);
end
ROM_time=toc
xxROM=(V*xROM')';
deltar=xxROM(:,1:n_oc); 
omegar=xxROM(:,n_oc+1:end);
        avg_delr=sum(deltar,2)/n_oc; % 1st state-variable
        avg_omegar=sum(omegar,2)/n_oc; %2st state-variable 
 %end
yROM(:,1)=avg_delr;
yROM(:,2)=avg_omegar;
end