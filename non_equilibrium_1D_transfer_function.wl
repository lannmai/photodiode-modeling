(* ::Package:: *)

ClearAll;
SetOptions[Plot, BaseStyle -> FontSize -> 14];
q = 1.6*10^(-19);
k = 1.38*10^(-23);
T = 300;
activeArea = 5*10^(-8);
W = 10^(-5);
Wdep = 10^(-6);
Dp = 12*10^(-4);
Dn = 36*10^(-4);
tauP = 10^(-3);
tauN = 10^(-3);
alpha = 2*10^5;
powerIn = 10*10^(-6);
opticalPower = powerIn/activeArea;
Ephoton = 1.65*1.6*10^(-19);
xP = W;
xN = W+Wdep;
tMax = tauP;
impedance = 50;

G[x_,t_] := (-alpha)*Exp[-alpha*x] * (opticalPower/Ephoton)*(1+Sin[modFreq*t]);


(* ::Text:: *)
(*Solve 1D diffusion PDE analytically for holes in n-type layer. delp[x,t] is excess hole concentration.*)


pde1 = D[delp[x,t], {t,1}] == (Dp)*D[delp[x,t], {x,2}] - delp[x,t]/tauP + G[x,t];
sol1 = DSolve[{pde1}, delp[x,t], x, t];


(* ::Text:: *)
(*Assign delp[x,t] to a new function excessP[x,t] for easier handling. Calculate hole current density.*)


excessP[x_,t_] = delp[x,t]/.sol1;
currentDenP[x_,t_] = -q*(Dp)*D[excessP[x,t], {x,1}];


(* ::Text:: *)
(*Solve 1D diffusion PDE analytically for electrons in p-type layer. deln[x,t] is excess electron concentration.*)


pde2 = D[deln[x,t], {t,1}] == (Dn)*D[deln[x,t], {x,2}] - deln[x,t]/tauN + G[x,t];
sol2 = DSolve[{pde2}, deln[x,t], x, t];


(* ::Text:: *)
(*Assign delp[x,t] to a new function excessP[x,t] for easier handling. Calculate hole current density.*)


excessN[x_,t_] = deln[x,t]/.sol2;
currentDenN[x_,t_] = q*(Dn)*D[excessN[x,t], {x,1}];


outPower[t_] = ((currentDenP[xP,t] + currentDenN[xN,t])*activeArea)^2*impedance;
outPowerLP[s_] = LaplaceTransform[outPower[t], t, s];
inPower[t_] = powerIn;
inPowerLP[s_] = LaplaceTransform[inPower[t], t, s];


H[s_] = FullSimplify[outPowerLP[s]/inPowerLP[s]];


freqRes[modFreq_] = FullSimplify[Sqrt[H[I*2*Pi*modFreq]*H[-I*2*Pi*modFreq]]];


LogLinearPlot[Log10[freqRes[modFreq]], {modFreq, 10^6, 4*10^9}, 
PlotRange->All]
