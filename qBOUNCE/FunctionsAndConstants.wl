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


Clear[b,A]
B[a_,b_,c_,rho_,\[Phi]Old_,H_,DVeff_,DDVeff_]:= Block[{Output},
Output={};

AppendTo[Output,DVeff[\[Phi]Old[[1]],a,b,c,rho[[1]]]  -DDVeff[\[Phi]Old[[1]],a,b,c,rho[[1]]]*\[Phi]Old[[1]] -Block[{$MaxExtraPrecision=1000},\[Phi]\[Rho][a,b,c,\[Rho]M]]/(H[[1]]^2*m2invMeV^2)](*From first boundary condition*);
For[i=2, i<=Length[\[Phi]Old]-1,i++, AppendTo[Output, DVeff[\[Phi]Old[[i]],a,b,c,rho[[i]]]  -DDVeff[\[Phi]Old[[i]],a,b,c,rho[[i]]]*\[Phi]Old[[i]]]];
AppendTo[Output,DVeff[\[Phi]Old[[-1]],a,b,c,rho[[-1]]]  -DDVeff[\[Phi]Old[[-1]],a,b,c,rho[[-1]]]*\[Phi]Old[[-1]] -Block[{$MaxExtraPrecision=1000},\[Phi]\[Rho][a,b,c,\[Rho]V]]/(H[[-1]]^2*m2invMeV^2)](*from 2nd bountady condition*);
Output]
(*Subscript[\[Rho], i] computes Subscript[\[Rho], i] prior to the iterations, such that it is not done repeadedly, same for Subscript[h, i]*)
\[Rho]i[points_]:=Block[{Output},
Output={};
For[i=1, i<=Length[points],i++, AppendTo[Output,\[Rho][points[[i]]]]];
Output]

h[points_]:=Block[{Output},
Output={};
For[i=1, i<=Length[points]-1,i++, AppendTo[Output,(points[[i+1]]-points[[i]])]];
AppendTo[Output,Output[[-1]]]; (*Subscript[h, N] would point outside the degrees of freedom and is therefore not part of points, I just repeat the last entry*);
Output]

(*H refers to the externally computed h, rho to Subscript[\[Rho], i]*)
A[a_,b_,c_,H_,rho_,\[Phi]Old_,DDVeff_]:=Block[{Output},
Output=ConstantArray[0,{Length[rho],Length[rho]}](*tridiagonal matrix with mostly 0 elements*); 
Output[[1]][[1]]=-(2/(H[[1]]^2*m2invMeV^2))-DDVeff[\[Phi]Old[[1]],a,b,c,rho[[1]]] (*from boundary*);
Output[[1]][[2]]=1/(H[[1]]^2*m2invMeV^2);

For[i=2, i<=Length[H]-1,i++, 
Output[[i]][[i+1]]=2/(H[[i]](H[[i]]+H[[i-1]])*m2invMeV^2);
Output[[i]][[i-1]]=2/(H[[i-1]](H[[i]]+H[[i-1]])*m2invMeV^2);
Output[[i]][[i]]= -(2/(H[[i]](H[[i]]+H[[i-1]])*m2invMeV^2))-2/(H[[i-1]](H[[i]]+H[[i-1]])*m2invMeV^2)-DDVeff[\[Phi]Old[[i]],a,b,c,rho[[i]]](*from boundary*)];
Output[[-1]][[-2]]=2/(H[[-2]](H[[-1]]+H[[-2]])*m2invMeV^2);
Output[[-1]][[-1]]=-(2/(H[[-1]]^2*m2invMeV^2))-DDVeff[\[Phi]Old[[-1]],a,b,c,rho[[-1]]](*from boundary*);
Output

];

(*This function converts the continous function f into a discrete version that can be used as seed for newtons method*)
CreateSeed[f_,points_]:=Block[{Seed},
Seed={};
For[i=1, i<=Length[points],i++, AppendTo[Seed,f[points[[i]]]]];
Seed]

(*The iterate function takes \[Phi]^(n-1) and returns \[Phi]^n*)

Clear[NORM]
NORM[Y1_]:=Sqrt[Sum[Y1[[i]]^2,{i,Length[points]}]]

Iterate[a_,b_,c_,H_,rho_,\[Phi]Old_,DVeff_,DDVeff_]:=Block[{Output,B2,C2},
B2=B[a,b,c,rho,\[Phi]Old,H,DVeff,DDVeff];
C2=A[a,b,c,H,rho,\[Phi]Old,DDVeff];
Output=LinearSolve[C2,B2];
Output]

(*careful, Seed needs to be converted to a vector prior to using SOLVE*)
SOLVEFD[a_,b_,c_,points_,Prec_,Max_,Seed_,DVeff_,DDVeff_]:=Block[{YOLD,YNEW,DIFF,M,H,RHO,Func,Output},
M=2;
H=h[points];
RHO=\[Rho]i[points];
YOLD=If[StringQ[Seed],If[Seed=="Default",TrivialSeed2[points,a,b,c],Seed],Seed];
YNEW=Iterate[a,b,c,H,RHO,YOLD,DVeff,DDVeff];
DIFF=YNEW-YOLD;
While[(NORM[DIFF]/NORM[YOLD]>10^-Prec&&M<=Max),
M=M+1;
YOLD=YNEW;
YNEW=Iterate[a,b,c,H,RHO,YOLD,DVeff,DDVeff];
DIFF=YNEW-YOLD;
];
Func={};
For[i=1,i<=Length[points],i++,
AppendTo[Func,{points[[i]],YNEW[[i]]}]];
Output=Interpolation[Func,InterpolationOrder->3];
M=M-1;
(*Print["iterations"M]*);
{YNEW,Output}]

TrivialSeed2[points_,a_,b_,c_]:=Block[{Seed,RHO},
RHO=\[Rho]i[points];
Seed={};
For[i=1, i<=Length[points],i++, AppendTo[Seed,Block[{$MaxExtraPrecision=1000},\[Phi]\[Rho][a,b,c,RHO[[i]]]]]];
Seed]







EndPackage[];
