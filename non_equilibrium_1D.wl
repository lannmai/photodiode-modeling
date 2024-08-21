(* ::Package:: *)

q = 1.6*10^(-19);
k = 1.38*10^(-23);
T = 300;
L = 1*10^(-6);
Dp = 12*10^(-4);
tauP = 10^(-6);
alpha = 1*10^5;
modFreq = 10^6;
opticalPower = 10^(-3);
Ephoton = 1.37*1.6*10^(-19);
xMax = L;
tMax = tauP;

G[x_,t_]:=(-alpha)*Exp[-alpha*x]*(opticalPower*(1+Sin[modFreq*t]))/Ephoton
pde = D[delp[x,t],{t,1}] == (Dp)*D[delp[x,t],{x,2}] - delp[x,t]/tauP + G[x,t];

sol = DSolve[{pde, delp[0,t]==c}, delp[x,t], x, t]

excessP[x_,t_] = delp[x,t]/.sol;
currentDen[x_,t_] = -q*Dp*D[excessP[x,t],{x,1}];

Plot[excessP[xMax,t],{t,0,tMax}, PlotRange->All, AxesLabel->{"Time (s)", "Excess hole concentration (1/m^3)"}]
Plot[currentDen[xMax,t],{t,0,tMax}, PlotRange->All, AxesLabel->{"Time (s)", "Current density (A/m^2)"}]
