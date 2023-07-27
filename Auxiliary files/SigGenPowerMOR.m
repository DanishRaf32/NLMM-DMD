function [xr, t, xr0, Sv, R, v0] = SigGenPowerMOR(x_eq0,A_lin,B_lin,fun)
% SIGNALGENERATORTEMPLATE - Template for Signal Generator script
%
% Syntax:
%		[xr, t, xr0, Sv, R, v0] = SIGNALGENERATORTEMPLATE()
%
% Description:
%       This function computes the inputs for NLMM according to the signal generators 
%       chosen by the user. This function distinguishes 3 different types of 
%       signal generators: 
%
%       # nonlinear signal generators, where solution is simulated (case 1)
%       # linear signal generators, where solution is simulated (case 2)
%       # signal generators with explicitly given analytical solution (case 3)
% 
%       In case 1 and 2 the solution is calculated by MATLAB's built-in solver
%       ODE45.
% 
% //Note: 
%       Cell arrays xrdotFun and xrsdotFun HAVE to contain the correct 
%       number of columns. 
%
%       The number of columns in xrdotFun is the sum of the 
%       number of signal generators in case 1 and 2.
%
% Input Arguments:
%
% Output Arguments:
%       -xr:
%       -t:
%       -xr0:
%       -Sv:
%       -R:
%       -v0:

%% preliminaries
m = size(B_lin,2);

%% ========================================================================
% =========================== CASES 1 AND 2 ===============================
% =========================================================================
%% nonlinear SG-ODE, where solution is simulated (case 1)


%% linear SG-ODE, where solution is simulated (case 2)
% xrdot1 = @(t,xr) tanh(xr) + 0.0005; %Signal Generator 1
% xr01 = 0;
% u1 = @(xr) 1 .* xr; 
% xrdot2 = @(t,xr) -exp(-t); %Signal Generator 1
%  xr02 = 1;
%  u2 = @(xr) 1 .* xr;

% xrdot1 = @(t,xr) 10*2*t.*exp(-t/5) + 10*t.^2.*exp(-t/5)*(-1/5); %Signal Generator 1
% xr01 = 0;
% u1 = @(xr) 1 .* xr;

% xrdot1 = @(t,xr) -0.5*sin((2*pi/10)*t)*(2*pi/10); %Signal Generator 2
% xr01 = 1;
% u1 = @(xr) 1 .* xr;

%% cell arrays containing every signal generator in case 1 and 2
%---Default:
xrdotFun = {};
xr0Fun = {};
uFun = {};

%---Customized:
%  xrdotFun = {xrdot1}; % containing every (!) xrdot, which has to be simulated 
%  xr0Fun = {xr01}; % containing every (!) xr0
%  uFun = {u1}; % containing every (!) u

% xrdotFun = {xrdot1,xrdot2}; % containing every (!) xrdot, which has to be simulated 
% xr0Fun = {xr01,xr02}; % containing every (!) xr0
% uFun = {u1,u2}; % containing every (!) u

%% simulating solution xr(t) for signal generators from case 1 and 2 and construction of xr, t, xr0, Sv, R, v0
t0 = 0; 
tEnd = 10;

%---Default:
%--- EACH signal generator has different K snapshots
t = cell(1,length(xrdotFun));
xr = cell(1,length(xrdotFun));
xr0 = cell(1,length(xrdotFun));
Sv = cell(1,length(xrdotFun));
R = cell(1,length(xrdotFun));
v0 = cell(1,length(xrdotFun));

%--- ALL signal generators are simulated with K snapshots
K = 50; tSimSigGen = linspace(t0,tEnd,K); % time span for simulation

%---Customized:
for i1 = 1 : length(xrdotFun)
    %--- EACH signal generator has different K snapshots
    [t{1,i1}, xr{1,i1}] = ode45(xrdotFun{i1},[t0 tEnd],xr0Fun{i1});
    %--- ALL signal generators are simulated with K snapshots
%     [t{1,i1}, xr{1,i1}] = ode45(xrdotFun{i1},tSimSigGen,xr0Fun{i1});
    
    Sv{1,i1} = xrdotFun{i1}(t{1,i1},xr{1,i1}); % containing every (!) xrdot{#}(t{1,#},xr{1,#})
    R{1,i1} = uFun{i1}(xr{1,i1}).';
    xr0{1,i1} = xr0Fun{i1};
    v0{1,i1} = fsolve(@(v) fun(v*xr0{i1}, uFun{i1}(xr0{i1}).'), x_eq0);
end

%% ========================================================================
% ============================== CASE 3 ===================================
% =========================================================================
%% calculate xr(t) for explicitly given analytical solutions of xr (case 3)
K = 150; thelp = linspace(t0,tEnd,K)';

%-------------------------------------------------------------------------%
% EXPONENTIAL FUNCTION
%-------------------------------------------------------------------------%
% Generators with xrDot = s0*xr with analytical solution xr(t) = exp(s0*t)*xr0.
% xrsdot# = @(t) s0(#)*exp(s0(#)*t)*xrs0#; 
% xrs# = @(t) exp(s0(#)*t)*xrs0#;
% us# = @(xrs) Rt(:,#).' .* xrs;

%---Default:
s01 = [];
Rt1 = [];
thelp1 = {};

%---Customized:
%---Either:
% s01 = [-(2*pi/10)*1i (2*pi/10)*1i]; % when using complex shifts, use the complex-conjugated pair
% s01 = [-0.3i 0.3i];
% xrs01 = [0.5 0.5];
% r = length(s01);
% Rt1 = ones(m,r);
% thelp1 = {tEnd, thelp};

%---Or:
% r = 15;
% s01 = logspace(-3,1,r);
% Rt1 = ones(m,r);
% xrs01 = randn(1,r);
% thelp1 = {tEnd, tEnd, tEnd, tEnd, tEnd, tEnd, tEnd, tEnd, tEnd, tEnd, tEnd, tEnd, tEnd, tEnd, tEnd};

if isempty(xrdotFun)
    i1 = 0;
end
for iShift1 = 1 : length(s01)
    xrsdotFun = @(t) s01(iShift1)*exp(s01(iShift1)*t)*xrs01(iShift1); 
    xrsFun = @(t) exp(s01(iShift1)*t)*xrs01(iShift1);
    usFun = @(xrs) Rt1(:,iShift1).' .* xrs;
    
    t{1,i1+iShift1} = thelp1{iShift1};
    xr{1,i1+iShift1} = xrsFun(thelp1{iShift1}); % apply t to xrsFun
    Sv{1,i1+iShift1} = xrsdotFun(thelp1{iShift1});
    R{1,i1+iShift1} = usFun(xrsFun(thelp1{iShift1})).';
    xr0{1,i1+iShift1} = xrs01(iShift1);
    v0{1,i1+iShift1} = (s01(iShift1)*speye(size(x_eq0,1))-A_lin)\(B_lin * Rt1(:,iShift1));
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% SHIFTED EXPONENTIAL FUNCTION
%-------------------------------------------------------------------------%
% Generators with xrDot = s0*(xr + alpha) with analytical solution xr(t) = exp(s0*t)*(xr0 + alpha) - alpha.
% xrsdotFun# = @(t) s0(#)*(exp(s0(#)*t)*(xrs0(#) + alpha(#)) - alpha(#) + alpha(#)); 
% xrsFun# = @(t) exp(s0(#)*t)*(xrs0(#) + alpha(#)) - alpha(#);
% usFun# = @(xrs) Rt(:,#).' .* xrs;

%---Default:
s02 = [];
Rt2 = [];
thelp2 = {};

%---Customized:
%  s02 = [];
%  alpha = [];
%  xrs02 = [];
%  Rt2 = ;
%  thelp2 = {thelp};

if isempty(s01)
    iShift1 = 0;
end
for iShift2 = 1 : length(s02)
    xrsdotFun = @(t) s02(iShift2)*(exp(s02(iShift2)*t)*(xrs02(iShift2) + alpha(iShift2)) - alpha(iShift2) + alpha(iShift2));
    xrsFun = @(t) exp(s02(iShift2)*t)*(xrs02(iShift2) + alpha(iShift2)) - alpha(iShift2);
    usFun = @(xrs) Rt2(:,iShift2).' .* xrs;
    
    t{1,i1+iShift1+iShift2} = thelp2{iShift2};
    xr{1,i1+iShift1+iShift2} = xrsFun(thelp2{iShift2}); % apply t to xrsFun
    Sv{1,i1+iShift1+iShift2} = xrsdotFun(thelp2{iShift2});
    R{1,i1+iShift1+iShift2} = usFun(xrsFun(thelp2{iShift2})).';
    xr0{1,i1+iShift1+iShift2} = xrs02(iShift2);
    v0{1,i1+iShift1+iShift2} = (s02(iShift2)*speye(size(x_eq0,1))-A_lin)\(B_lin * Rt2(:,iShift2));
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% OTHERS 
%-------------------------------------------------------------------------%
% Other generators with analytical solution xr(t)
% xrsdotFun = @(t,xrs) 
% xrsFun = @(t)
% usFun = @(xrs)

%------
t0=0;
tEnd=100;
xrsdotFun1 = @(t,xrs) -pi*exp(-pi*t);
xrs0Fun1 = 0;
xrsFun1 = @(t) exp(-pi*t);
usFun1 = @(xrs) 1 .* xrs;
K = 28; thelp3{1} = linspace(t0,tEnd,K)';

%------
xrsdotFun2 = @(t,xrs) 10*2*t.*exp(-t/5) + 10*t.^2.*exp(-t/5)*(-1/5);
xrs0Fun2 = 1;
xrsFun2 = @(t) 10*(t.^2.*exp(-t/5));
usFun2 = @(xrs) 1 .* xrs;
K = 101; thelp3{1} = linspace(t0,tEnd,K)';

%------
A = 0.5;
T = 10;
xrsdotFun3 = @(t,xrs) -A*sin((2*pi/T)*t)*(2*pi/T);
xrs0Fun3 = 2*A;
xrsFun3 = @(t) A*(cos((2*pi/T)*t) + 1);
usFun3 = @(xrs) 1 .* xrs;
K = 28; thelp3{1} = linspace(t0,tEnd,K)';

%------
A = 5;
T = 10;
xrsdotFun4 = @(t,xrs) A*cos((2*pi/T)*t)*(2*pi/T);
xrs0Fun4 = A;
xrsFun4 = @(t) A*(sin((2*pi/T)*t) + 1);
usFun4 = @(xrs) 1 .* xrs;
K = 28; thelp3{1} = linspace(t0,tEnd,K)';

%--
t0=0;
tEnd=100;
xrsdotFun5 = @(t,xrs) 5*pi*sin(pi*t);
xrs0Fun5 = 1;
xrsFun5 = @(t) -5*cos(pi*t);
usFun5 = @(xrs) 1 .* xrs;
K = 28; thelp3{1} = linspace(t0,tEnd,K)';

%--
t0=0;
tEnd=20; %20-EN, 10,2-SM, (IEEE300=8)
amp=5; %5-EN, 
freq= 2*pi*60; %2*pi*60-EN: (IEEE300=20*pi*600)
phase=0; %0
xrsdotFun6 = @(t,xrs) freq*amp*sin(freq*t + phase);
xrs0Fun6 = 0; %0-EN
xrsFun6 = @(t) - amp*cos(freq*t + phase);
usFun6 = @(xrs) 1* xrs;
K = 40; %40-EN (IEEE300=80)
thelp3{1} = linspace(t0,tEnd,K)';
% 
% xrsdotFun7 = @(t,xrs) 0
% xrs0Fun7 = 0.95;
% xrsFun7 = @(t) 0.95;
% usFun7 = @(xrs) 1 .* xrs;
% K = 101; thelp3{4} = linspace(t0,tEnd,K)';

%---Default:
% xrsdotFun = {};
% xrsFun = {};
% xrs0Fun = {};
% usFun = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Customized:
xrsdotFun = {xrsdotFun6}; % containing every (!) xrsdotFun 
xrsFun = {xrsFun6}; % containing every (!) xrsdotFun
xrs0Fun = {xrs0Fun6}; % containing every (!) xrs0Fun
usFun = {usFun6}; % containing every (!) usFun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xrsdotFun = {xrsdotFun1,xrsdotFun2,xrsdotFun3}; % containing every (!) xrsdotFun 
% xrsFun = {xrsFun1,xrsFun2,xrsFun3}; % containing every (!) xrsdotFun
% xrs0Fun = {xrs0Fun1,xrs0Fun2,xrs0Fun3}; % containing every (!) xrs0Fun
% usFun = {usFun1,usFun2,usFun3}; % containing every (!) usFun

if isempty(s02)
    iShift2 = 0;
end
for iShift3 = 1 : length(xrsdotFun)
    t{1,i1+iShift1+iShift2+iShift3} = thelp3{iShift3};
    xr{1,i1+iShift1+iShift2+iShift3} = xrsFun{iShift3}(thelp3{iShift3}); % apply t to xrsFun
    Sv{1,i1+iShift1+iShift2+iShift3} = xrsdotFun{iShift3}(thelp3{iShift3},xrsFun{iShift3}(thelp3{iShift3})); % containing every (!) xrsdotFun{#}
    R{1,i1+iShift1+iShift2+iShift3} = usFun{iShift3}(xrsFun{iShift3}(thelp3{iShift3})).';
    xr0{1,i1+iShift1+iShift2+iShift3} = xrs0Fun{iShift3};
    v0{1,i1+iShift1+iShift2+iShift3} = fsolve(@(v) fun(v*xrs0Fun{iShift3}, usFun{iShift3}(xrs0Fun{iShift3}).'), x_eq0);
end
%-------------------------------------------------------------------------%
end

%% Symbolic computation of scalar, first-order ODE
% syms xr(t)
% 
% %-------
% ode = diff(xr,t) == 10*2*t*exp(-t/5) + 10*t.^2*(-1/5)*exp(-t/5)
% xrSol(t) = dsolve(ode)
% 
% cond = xr(0) == 0;
% xrSol(t) = dsolve(ode,cond)
% 
% %-------
% ode = diff(xr,t) == -exp(-t)
% xrSol(t) = dsolve(ode)
% 
% cond = xr(0) == 1;
% xrSol(t) = dsolve(ode,cond)
% 
% %-------
% ode = diff(xr,t) == -0.5*sin((2*pi/10)*t)*(2*pi/10)
% xrSol(t) = dsolve(ode)
% 
% cond = xr(0) == 1;
% xrSol(t) = dsolve(ode,cond)
% 
% %-------
% ode = diff(xr,t) == 5*cos((2*pi/10)*t)*(2*pi/10)
% xrSol(t) = dsolve(ode)
% 
% cond = xr(0) == 5;
% xrSol(t) = dsolve(ode,cond)

