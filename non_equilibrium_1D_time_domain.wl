(* ::Package:: *)

SetOptions[Plot, BaseStyle -> FontSize -> 14];
q = 1.6*10^(-19);
T = 300;
activeArea = 5*10^(-8);
W = 10^(-5);
Wdep = 10^(-6);
Dp = 12*10^(-4);
Dn = 36*10^(-4);
tauP = 10^(-6);
tauN = 10^(-6);
absorptionCoeff = 2*10^5;
powerIn = 10*10^(-6);
opticalPower = powerIn/activeArea;
Ephoton = 1.65*1.6*10^(-19);
modFreq = 1*10^8;
xP = W;
xN = W+Wdep;
tMax = 0.5*tauP;
impedance = 50;

G[x_, t_] := (-absorptionCoeff)*Exp[-absorptionCoeff*x]*(opticalPower/Ephoton)*(1+Sin[modFreq*t]);


(* ::Text:: *)
(*Solve 1D diffusion PDE analytically for holes in n-type layer. delp[x,t] is excess hole concentration.*)


pde1 = D[delp[x, t], {t, 1}] == (Dp)*D[delp[x, t], {x, 2}]-delp[x, t]/tauP+G[x, t];
sol1 = DSolve[{pde1}, delp[x, t], x, t];


(* ::Text:: *)
(*Assign delp[x,t] to a new function excessP[x,t] for easier handling. Calculate hole current density.*)


excessP[x_, t_] = delp[x, t]/.sol1;
currentDenP[x_, t_] = -q*(Dp)*D[excessP[x, t], {x, 1}];


Plot[excessP[xP, t], {t, 0, tMax}, 
PlotRange->All, AxesLabel->{"Time (s)", "Excess hole concentration (1/m^3)"}, PlotStyle->Red]

Plot[currentDenP[xP, t], {t, 0, tMax}, 
PlotRange->All, AxesLabel->{"Time (s)", "Hole current density (A/m^2)"}]


(* ::Text:: *)
(*Solve 1D diffusion PDE analytically for electrons in p-type layer. deln[x,t] is excess electron concentration.*)


pde2 = D[deln[x, t], {t, 1}] == (Dn)*D[deln[x, t], {x, 2}]-deln[x, t]/tauN+G[x, t];
sol2 = DSolve[{pde2}, deln[x, t], x, t];


(* ::Text:: *)
(*Assign delp[x,t] to a new function excessP[x,t] for easier handling. Calculate hole current density.*)


excessN[x_, t_] = deln[x, t]/.sol2;
currentDenN[x_, t_] = q*(Dn)*D[excessN[x, t], {x, 1}];


Plot[excessN[xN, t], {t, 0, tMax}, 
PlotRange->All, AxesLabel->{"Time (s)", "Excess electron concentration (1/m^3)"}, PlotStyle->Red]

Plot[currentDenN[xN, t], {t, 0, tMax}, 
PlotRange->All, AxesLabel->{"Time (s)", "Electron current density (A/m^2)"}]


Plot[(currentDenN[xN, t]+currentDenP[xP, t]), {t, 0, tMax}, 
PlotRange->All, AxesLabel->{"Time (s)", "Current density (A/m^2)"}]


outPower[t_] = ((currentDenP[xP, t]+currentDenN[xN, t])*activeArea)^2*impedance;
inPower[t_] = powerIn*(1+Sin[modFreq*t]);

maxInPower = FindMaxValue[inPower[t], {t, 0}];
maxOutPower = FindMaxValue[outPower[t], {t, 0}];

Plot[{inPower[t]/maxInPower, outPower[t]/maxOutPower}, {t, 0, tMax}, 
PlotRange->All, AxesLabel->{"Time (s)", "Normalized power"}, PlotLegends->{"Input (optical) power", "Output (photocurrent) power"}]
