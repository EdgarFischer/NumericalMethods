(* ::Package:: *)

BeginPackage["Functions`"]


(* First define conversion factors between various units*)
CPrec        =500; (* numerical precission for later calculations*)
constants=SetPrecision[{hbar->1.05457*^-34,c0->2.997925*^8,ec->1.602177*^-19},5 CPrec];
(* Note that ec is the unitless(!!) value of the electric charge, which is used to convert eV to J. ev=ec*J *)

m2invMeV  =1/(hbar c0/(1*^6*ec))/.constants; (*to be read as convert meter to MeV^-1*)
kg2MeV       =1/(ec 10^6/(c0^2))/.constants; 
s2invMeV  =1/(hbar/(1*^6 ec))/.constants; 
N2invMeV2=kg2MeV m2invMeV/s2invMeV^2; (*Newton to inverse MeV^2*)
kgm32MeV4=kg2MeV/m2invMeV^3; (*to be read as kg/m^3 converted to MeV^4 / density to natural units*)

(*Next define various constants / parameters that are use in the LLR analysis*)

\[Rho]M2= SetPrecision[2514.0`*kgm32MeV4,5 CPrec]; (*mirror density cannex*)

ddr=SetPrecision[5.0677398985414185`*^7,1000]; (*10 \[Mu]m in MeV*)


mpl=SetPrecision[2.4353635268*^21,5 CPrec](*planck mass in MeV*)

H0 := SetPrecision[197.3269788/(1.374*10^41),5CPrec];
\[CapitalOmega] := SetPrecision[0.73,5CPrec];
VC= 3 \[CapitalOmega] (mpl^2) (H0^2)         (* SetPrecision[6.73304*^-34,5 CPrec]; *)(*in MeV^4, this is the effective potential for the vacuum density and the vaccuum expactation value of \[Phi]*)





EndPackage[];
