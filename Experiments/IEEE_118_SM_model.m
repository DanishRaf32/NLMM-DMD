% This file is part of "Synergestic use of intrusive and non-intrusive Model 
% Order Reducion techniques for Dynamic Power Grids"
%
%------------------------------------------------------------------
% Authors:      Danish Rafiq, Junaid Farooq, M. A. Bazaz 
% Email:        <a href="mailto:danish_pha2007@nitsri.ac.in">
% Last Change:  15 June 2021
% Copyright (c) Department of Electrical Engineering, NIT Srinagar
%------------------------------------------------------------------
% This scipt simulates:
% 1) the full-order SM models of the IEEE-118 bus system,
% 2) Reduced-Order Models via NLMM-DMD and POD for comparison.
%   
% Selected-Refs:
%   [1] T. Nishikawa and A. E. Motter, Comparative analysis of existing
%       models for power-grid synchronization, New J. Phys. 17, 015012 (2015).
%
%   [2] A. E. Motter, S. A. Myers, M. Anghel, and T. Nishikawa, Spontaneous
%       synchrony in power-grid networks, Nat. Phys. 9, 191-197 (2013).
%
%   [3] Tobias, Sara Grundel, Nonlinear model reduction of dynamical power 
%       grid models using quadratization and balanced truncation (2021)
%
%   [4] Rafiq, Bazaz, Nonlinear Model Order Reduction via nonlinear moment-
%       matching with Dynamic Mode Decomposition (2020).
%------------------------------------------------------------------
% Pre-reqiuisites to run the script:
% 1. NLMM code along-with solvers from zenodo: (https://zenodo.org/record/3542641#.YMhNLDaA79E)
% 2. MATPOWER toolbox: (https://matpower.org/download/)
% 3. pg-sync-models from sourceforge: (https://sourceforge.net/projects/pg-sync-models/)
% 4. sss & sssMOR MATLAB toolbox https://github.com/MORLab

%% Initialize the IEEE Power System and model configuration
clear;clc;close all
 mpc=case118;      %(see MATPOWER data file)
 mpc.ref_freq=60;  %reference frequency
 global data
 [data,details]=SM_model(mpc);  
 n_oc=length(data.H);   % No of oscillators (FOM Size = 2*n_oc)
 u=@(t) 1;
 tspan=0:0.01:5; % Initial simulation time
 solver='Euler'; % (Select solver: Euler, ode15s) 
 rdefl=20; % Rank for POD and NLMM 
 rdmd=25; % Reduced rank of Koopman operator 
% FOM:
 disp('Simulating Full Order Model...')
 f= @(x) power_func(x); %nonlinear function
 Jack=@(x) Power_Jack(x); %Jacobian function
 x0=zeros(2*n_oc,1);    %initial conditions
 [yFOM,xFOM,FOM_time,delta_pre,omega_pre]=FOM_solve(f,Jack,tspan,x0,solver,n_oc);
% POD-DEIM:
 [Vpod]=basisRed(xFOM',rdefl);
 Opts.NoSV = 30; %DMD-Rank
 [U_trunc,P,~,~] = funcSnapBasisDEIM(f,xFOM',Opts);
 pre = Vpod'*U_trunc/(P'*U_trunc)*P';
 disp('Simulating POD-DMD-ROM..')
 [yROM_POD,ROM_POD_time,delta_POD_pre,omega_POD_pre]=ROM_solve_DEIM(pre,f,Jack,tspan,x0,solver,Vpod,n_oc);
% NLMM: signal generator + approximated nonlinear Sylvester eqn  
 fun=@(x,u) f(x);         
 Jackfun= @(x,u) Jack(x); 
 x_eq0=x0;                
 A_lin=Jack(x_eq0);
 B_lin=zeros(2*n_oc,1);
 [xr, t, xr0, Sv, R, v0] = SigGenPowerMOR(x_eq0,A_lin,B_lin,fun);
 sigGen_time=toc;
 optionsMM.rdefl = rdefl;  %deflation of V_NLMM
 optionsMM.RelTol = 1e-3; % relative tolerance for residual error: norm(f(xcurr))
 optionsMM.AbsTol = 1e-6;
 optionsMM.Display = 'none'; % 'none', 'iter'
 optionsMM.real = 1;
 optionsMM.solver='fsolve';  %NewtonRaphson, fsolve
 tic
 [Vnlmm, nlmmNLSE, v0raw] = nlmm(fun,Sv,R,xr,v0,Jackfun,optionsMM);
 NLMM_time=toc
 % Hyper-reduction via DMD: nonlinear snapshots + Koopman operator + rank truncation 
 X=xFOM';
 dt=0.01;
 F = zeros(size(X));
 for i = 1:size(X,2)
      F(:,i) = f(X(:,i));  % Snapshots of nonlinear term
 end
 X1=F(:,1:end-1);
 X2=F(:,2:end);
 tic
    [U,S,VV]=svd(X1,'econ');
    U=U(:,1:rdmd);
    S=S(1:rdmd,1:rdmd);
    VV=VV(:,1:rdmd);
    atilde=U'*X2*VV/(S);
    [W,D]=eig(atilde);
    Phi=X2*VV/S*W;         %DMD modes
    Lambda=abs(diag(D)) ;  %discrete-time eigenvalues
    omega= log(D)./dt;     %continuous-time eigenvalues
    b=pinv(Phi)*X1(:,1);
    base=Vnlmm'*Phi;
    f_DMD =@(t) base*diag(exp(omega*t).*b);
 DMD_time=toc;
 disp('Simulating NLMM-DMD ROM...')
 [yROM_NLMM,ROM_NLMM_time,delta_DMD_pre,omega_DMD_pre]=ROM_solve_DMD(f_DMD,Jack,tspan,x0,solver,Vnlmm,n_oc);
% Errors
 [Error,quality]=error_cal(yFOM,yROM_POD,yROM_NLMM);
% Plotting Outputs
 plotOut(tspan,yFOM,yROM_POD,yROM_NLMM)
% Plotting Errors
 error_Plot(tspan,Error)
% Qualitative Table
 results_Table(Error,FOM_time,ROM_POD_time,ROM_NLMM_time,FOM_time,NLMM_time,yFOM,n_oc,rdefl,rdmd)
%% CPU Time Comparison for different Ranks of: POD, DEIM, NLMM, DMD
% Note: CPU times can vary depending upon the system-processor and core-speed
disp('Starting Iterations..')
rank_no=5:5:30;
for i =1:length(rank_no)
    disp(['Iteration:',num2str(i),' of ',num2str(length(rank_no))])
% POD-DEIM
tic
 [Vpod]=basisRed(xFOM',rank_no(i));
POD_offline(i)=toc;
POD_offline(i)=POD_offline(i)+FOM_time;
Opts.NoSV = rank_no(i); 
[U_trunc,P,~,~] = funcSnapBasisDEIM(f,xFOM',Opts);
pre = Vpod'*U_trunc/(P'*U_trunc)*P';
[yROM_POD,ROM_POD_time]=ROM_solve_DEIM(pre,f,Jack,tspan,x0,solver,Vpod,n_oc);
POD_online(i)=ROM_POD_time;
POD_total(i)=POD_offline(i) + POD_online(i);
error_POD(i)=norm(abs(yFOM(:,1)-yROM_POD(:,1)),2);
% NLMM-DMD
fun=@(x,u) f(x);
Jackfun= @(x,u) Jack(x);
x_eq0=x0;
A_lin=ones(2*n_oc,1);
B_lin=ones(2*n_oc,1);
[xr, t, xr0, Sv, R, v0] = SigGenPowerMOR(x_eq0,A_lin,B_lin,fun);
sigGen_time=toc;
optionsMM.rdefl = rank_no(i);  %deflation of V_NLMM
optionsMM.RelTol = 1e-3; % relative tolerance for residual error: norm(f(xcurr))
optionsMM.AbsTol = 1e-6;
optionsMM.Display = 'none'; % 'none', 'iter'
optionsMM.real = 1;
optionsMM.solver='fsolve';  %NewtonRaphson, fsolve
tic
 [Vnlmm, nlmmNLSE, v0raw] = nlmm(fun,Sv,R,xr,v0,Jackfun,optionsMM);
 NLMM_time=toc;
NLMM_offline(i)=NLMM_time;
X=xFOM';
dt=0.01;
F = zeros(size(X));
for j = 1:size(X,2)
     F(:,j) = f(X(:,j));  % Snapshots of nonlinear term
end
X1=F(:,1:end-1);
X2=F(:,2:end);
rdmd=rank_no(i);
tic
   [U,S,VV]=svd(X1,'econ');
   U=U(:,1:rdmd);
   S=S(1:rdmd,1:rdmd);
   VV=VV(:,1:rdmd);
   atilde=U'*X2*VV/(S);
   [W,D]=eig(atilde);
   Phi=X2*VV/S*W;         %DMD modes
   Lambda=abs(diag(D)) ;  %discrete-time eigenvalues
   omega= log(D)./dt;     %continuous-time eigenvalues
   b=pinv(Phi)*X1(:,1);
   base=Vnlmm'*Phi;
   f_DMD =@(t) base*diag(exp(omega*t).*b);
DMD_time=toc;
[yROM_NLMM,ROM_NLMM_time]=ROM_solve_DMD(f_DMD,Jack,tspan,x0,solver,Vnlmm,n_oc);
NLMM_online(i)=ROM_NLMM_time;
NLMM_total(i)= NLMM_offline(i) + NLMM_online(i);
error_NLMM(i)=norm(abs(yFOM(:,1)-yROM_NLMM(:,1)),2);
end
% plot CPU times
figure
subplot(1,3,1)
semilogy(rank_no,POD_offline,'g-*',rank_no,NLMM_offline,'r-o','linewidth',2)
title('Offline time')
ylabel('CPU time(sec)','Interpreter','LaTeX')
xlabel('Rank: $r$','Interpreter','LaTeX')
grid on 
legend('POD','NLMM','Interpreter','LaTeX')
set(gca,'FontSize',20,'TicklabelInterpreter','LaTeX')
subplot(1,3,2)
semilogy(rank_no,POD_online,'g-*',rank_no,NLMM_online,'r-o','linewidth',2)
title('Online time')
ylabel('CPU time(sec)','Interpreter','LaTeX')
xlabel('r','Interpreter','LaTeX')
grid on 
legend('DEIM','DMD','Interpreter','LaTeX')
set(gca,'FontSize',20,'TicklabelInterpreter','LaTeX')
subplot(1,3,3)
semilogy(rank_no,POD_total,'g-*',rank_no,NLMM_total,'r-o','linewidth',2)
grid on
title('Total CPU time')
ylabel('CPU time(sec)','Interpreter','LaTeX')
xlabel('Rank $r$','Interpreter','LaTeX')
legend('POD-DEIM','NLMM-DMD','Interpreter','LaTeX')
set(gca,'FontSize',20,'TicklabelInterpreter','LaTeX')
% plot error
figure
semilogy(rank_no,error_POD,'g-*',rank_no,error_NLMM,'r-o','linewidth',2)
title('Error')
grid on
ylabel('$e(t)$','Interpreter','LaTeX')
xlabel('Rank: $r$','Interpreter','LaTeX')
legend('POD-DEIM','NLMM-DMD','Interpreter','LaTeX')
set(gca,'FontSize',20,'TicklabelInterpreter','LaTeX')
disp('Code finished successfully!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%