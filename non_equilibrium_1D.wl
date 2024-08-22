(* ::Package:: *)

SetOptions[Plot, BaseStyle -> FontSize -> 14];

q = 1.6*10^(-19);
k = 1.38*10^(-23);
T = 300;
W = 1*10^(-6);
Dp = 12*10^(-4);
tauP = 10^(-6);
alpha = 1*10^5;
modFreq = 10^9;
opticalPower = 10^(-3);
Ephoton = 1.37*1.6*10^(-19);
xMax = W;
tMax = 0.1tauP;

G[x_,t_] := (-alpha)*Exp[-alpha*x] * (opticalPower/Ephoton)*(1+Sin[modFreq*t]);
pde = D[delp[x,t], {t,1}] == (Dp)*D[delp[x,t], {x,2}] - delp[x,t]/tauP + G[x,t];

sol = DSolve[{pde}, delp[x,t], x, t];

excessP[x_,t_] = delp[x,t]/.sol;
currentDen[x_,t_] = -q*(Dp)*D[excessP[x,t], {x,1}];

maxCurrentDen = FindMaxValue[currentDen[xMax,t], {t,0}];
minG = FindMinValue[G[xMax,t], {t,0}];

Plot[excessP[xMax,t], {t,0,tMax}, 
PlotRange->All, AxesLabel->{"Time (s)", "Excess hole concentration (1/m^3)"}]

Plot[Abs[currentDen[xMax,t]]*10^3, {t,0,tMax}, 
PlotRange->All, AxesLabel->{"Time (s)", "Current density (mA/m^2)"}]

Plot[{Abs[currentDen[xMax,t]/maxCurrentDen], Abs[G[xMax,t]/minG]}, {t,0,tMax}, 
PlotRange->All, AxesLabel->{"Time (s)", "Normalized value"}]



