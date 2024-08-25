(* ::Package:: *)

SetOptions[{Plot, LogLinearPlot}, BaseStyle -> FontSize -> 14];
k = 1.38*10^(-23);
q = 1.60*10^(-19);
h = 6.63*10^(-34);
elecmass = 9.11*10^(-31);
me = 0.043*elecmass;
mh = 0.46*elecmass;
EG = 1.12*q;
ND = 10^21;
NA = 10^25;
eps = 11.7;
A = 10^(-12);
pF = 10^12;

EionD = 13.6*(me/elecmass)*q/eps^2;
EionA = 13.6*(mh/elecmass)*q/eps^2;

NDion[T_] := ND/(Exp[(EionD/2)/(k*T)]+1);
NAion[T_] := NA/(Exp[(EionA/2)/(k*T)]+1);

V0[T_] := (k*T/q)*Log[(NDion[T]*NAion[T]/4)*(h^6/(2*Pi*k*T)^3)*(me*mh)^(-3/2)]+EG/q;
Cj[V_, T_] := (A/2)*Sqrt[2*q*eps*NDion[T]/(V0[T]-V)];


LogLinearPlot[Evaluate@Table[Cj[V, T]*pF, {V, {-4, -3, -2, -1, 0}}],{T, 6, 300}, 
PlotRange->All, PlotLegends->LineLegend[Table[V, {V, {-4, -3, -2, -1, 0}}], LegendLabel->"Reverse bias (V)"], 
AxesLabel->{"log T (K)", "Junction capacitance (pF)"}]


Plot[Evaluate@Table[Cj[V, T]*pF, {T, {4.7, 6, 10, 30, 50, 100}}], {V, -5, 0},
PlotRange->All, PlotLegends->LineLegend[Table[T, {T, {4.7, 6, 10, 30, 50, 100}}], LegendLabel->"T (K)"], 
AxesLabel->{"V (V)", "Junction capacitance (pF)"}]


(* ::Text:: *)
(*Capacitance sensitivity to reverse bias and temperature*)


CjV'[V_, T_]=D[Cj[V,T], V];
CjT'[V_, T_]=D[Cj[V,T], T];


Plot[Evaluate@Table[CjV'[V, T]*pF,{T, {6, 10, 20, 50, 77, 100, 200, 300}}], {V, -5, 0}, 
PlotLegends->LineLegend[Table[T, {T, {6, 10, 20, 50, 77, 100, 200, 300}}], LegendLabel->"T (K)"], AxesLabel->{"Reverse bias (V)", "Sensitivity (pF/V)"}, PlotRange->All]


LogLinearPlot[Evaluate@Table[CjT'[V, T]*pF, {V, {-4, -3, -2, -1, 0}}], {T, 6, 300},
PlotLegends->LineLegend[Table[V, {V, {-4, -3, -2, -1, 0}}], LegendLabel->"Reverse bias (V)"], AxesLabel->{"T (K)", "Sensitivity (pF/T)"}, PlotRange->All]
