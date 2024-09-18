#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*Notebook tested on Mathematica 14.1*)
(*This notebook computes \!\(\*SubsuperscriptBox[\(C\), \(lm\), \(e/o\)]\)(\[Omega]) on a grid of \[Omega] values for a given set of (l,m) modes and exports the result to a file. *)
(**)
(*Written by Geoffrey Comp\[EGrave]re and Lorenzo K\[UDoubleDot]chler*)


(* ::Chapter::Closed:: *)
(*I. Parameters of this computation *)


NPoints= 2000; (* Number of points on the grid of [-1,1] frequency values *)
Maxw=4;(* Range of \[Omega] frequency: \[Omega] \[Epsilon] [-Maxw,Maxw]] *)
thiscutoff=5.999; (* Cutoff rc close to the ISCO *)
PrecisionAin=100;(* Precision required for the computation of Ain by the Black Hole Perturbation Toolkit procedure *)


 (* List of (l,m) modes to be computed *)
 ListModes={{2,2},{2,1},{2,-1},{2,-2}};
 (* ListModes={{3,3},{3,2},{3,1},{3,-1},{3,-2},{3,-3}}; *)
 (* ListModes={{4,4},{4,3},{4,2},{4,1},{4,-1},{4,-2},{4,-3},{4,-4}}; *)
 (* ListModes={{5,5},{5,4},{5,3},{5,2},{5,1},{5,-1},{5,-2},{5,-3},{5,-4},{5,-5}};*)


 (* Filename where the results will be exported: *)
FilenameWhereToExportFinalData="Modes2.dat";


(* ::Chapter::Closed:: *)
(*II. Setup parallel computing*)


(* ::Text:: *)
(*Setup of parallel computing on High Performance Computers with job management system based on Slurm*)
(*If Slurm is not available, e.g., on a local machine, launch available kernels on the local machine. *)


(*Get association of resources,name of local host,and remove local host from available resources*)
nodeList=Import["!scontrol show hostnames $SLURM_JOB_NODELIST","Table"];
hosts=Counts[nodeList];
localNodeID=ToExpression[Environment["SLURM_NODEID"]];
If[localNodeID== $Failed,localHost={};,localHost=nodeList[[localNodeID+1]];]
hosts[localHost]--;

(*Launch subkernels and connect them to the controlling Wolfram Kernel*)
(* Needs["SubKernels`RemoteKernels`"];*)
Map[If[hosts[#]>0,LaunchKernels[RemoteMachine[#,hosts[#]]]]&,Keys[hosts]];



Print[{$MachineName,$KernelID}]
Print["kernel count is  ",$KernelCount];
Print[" machine name is  ",ParallelEvaluate[$MachineName]];
Print[" kernel id is  ",ParallelEvaluate[$KernelID]];
Print["processor count is  ",$ProcessorCount];


(* ::Chapter::Closed:: *)
(*III. Computation*)


Print["Main computation started"];


<<ReggeWheeler`
(*Ain*)
Ain$BHPT[l_,w_]:=Module[{output},output=ReggeWheelerRadial[2,l,w]["In"]["Amplitudes"]["Incidence"];If[Not[$MessageList== {}],Print["Evaluation of ReggeWheelerRadial mode l=",l," w=",w,": Encountered messages:",$MessageList]];output]


(* ::Section::Closed:: *)
(*1. Analytic derivation of the integrand of \!\(\*SubsuperscriptBox[\(C\), \(lm\), \(e/o\)]\)(\[Omega])*)


(* ::Text:: *)
(*Derivation of the integrand \!\(\*SubsuperscriptBox[\(C\), \(lm\), \(e/o\)]\)(\[Omega])*)


(* ::Subsection::Closed:: *)
(*Source: plunging geodesic with conserved charges at the ISCO values*)


MVal=1;t0InitVal=0;\[Phi]0InitVal=0;
\[Gamma][l_]:=((l+2)(l-1))/2
\[CapitalLambda][l_,rG_]:=\[Gamma][l]+(3M)/rG
rGdot[rG_]:=-(1/2^(3/2)) ((6M)/rG-1)^(3/2) (1-(2M)/rG)
to[rG_]:=(2Sqrt[2](rG-24M))/((6M)/rG-1)^(1/2)-44M Sqrt[2]ArcTan[((6M)/rG-1)^(1/2)]+4M ArcTanh[1/Sqrt[2] ((6M)/rG-1)^(1/2)]+t0InitVal
\[Phi]o[rG_]:=-((2Sqrt[3])/((6M)/rG-1)^(1/2))+\[Phi]0InitVal
\[CapitalOmega]o[rG_]:=Sqrt[3/2] (3M)/rG^2 (1-(2M)/rG)
SubsSphTerms={Ylm->Function[{l}, SphericalHarmonicY[l,m,\[Pi]/2,0]],D\[Theta]Ylm->Function[{l}, D[SphericalHarmonicY[l,m,\[Theta],0],\[Theta]]/.{\[Theta]-> \[Pi]/2}],lfac-> Function[{l},(l-2)!/(l+2)!]}//FunctionExpand//FullSimplify;
(* Source terms A and B *)
Aodd[l_,m_,rG_]:=-(1 /(l (1+l) \[Gamma][l] rG^7))4  M \[Pi] (-2 M+rG)^2 (756 Sqrt[3] M^3-18 I m M (6 M-rG)^(3/2) Sqrt[rG]-216 Sqrt[3] M^2 rG+45 Sqrt[3] M rG^2+7 Sqrt[3] rG^3)  D\[Theta]Ylm[l]
Bodd[l_,m_,rG_]:=-((36 Sqrt[3] M \[Pi] (2 M-rG)^3 (12 M^2+rG^2) D\[Theta]Ylm[l])/(l (1+l) \[Gamma][l] rG^6))
Aeven[l_,m_,rG_]:=-(1/((1+\[Gamma][l]) rG^5 (3 M+\[Gamma][l] rG)^2))2 Sqrt[2] \[Pi] (-2 M+rG)^2 (324 M^4+36 M^3 (3+8 \[Gamma][l]) rG+9 M^2 (1+4 \[Gamma][l]^2) rG^2+2 (M+9 M \[Gamma][l]) rG^3+3 \[Gamma][l] (1+\[Gamma][l]) rG^4)Ylm[l]-(16 Sqrt[6]  M \[Pi] (6 M-rG)^(3/2) (-2 M+rG)^2)/(l (1+l) rG^(9/2) (3 M+\[Gamma][l] rG)) (-I m Ylm[l])-(288 Sqrt[2]  lfac[l] M^2 \[Pi] (-2 M+rG)^2)/rG^5 (-m^2+l (l+1)/2)Ylm[l]
Beven[l_,m_,rG_]:=(6 Sqrt[2]  \[Pi] (-2 M+rG)^3 (12 M^2+rG^2) Ylm[l])/((1+\[Gamma][l]) rG^4 (3 M+\[Gamma][l] rG))


(* ::Subsection::Closed:: *)
(*Homogeneous data*)


(* In solution *)
\[Psi]oIn[l_,w_,r_,M_]:=1/r 2 E^(-I (4 M-r) w) M (-1+r/(2 M))^(-2 I M w) HeunC[l+l^2-2 (1+4 I M w+8 M^2 w^2),-4 M w (-I+4 M w),1-4 I M w,-3,-4 I M w,1-r/(2 M)]
\[Psi]eIn[l_,w_,r_,M_]:=1/(r^2 ((-1+l) l (1+l) (2+l)-12 I M w)) 2 E^(-I (4 M-r) w) M (-1+r/(2 M))^(-2 I M w) (1/(6 M+(-2+l+l^2) r) ((-2+l+l^2) (24 M^2+6 (-2+l+l^2) M r+(-1+l) l (1+l) (2+l) r^2)-12 I M (4 M-r) (6 M+(-2+l+l^2) r) w) HeunC[-2+l+l^2-8 M w (I+2 M w),4 M w (I-4 M w),1-4 I M w,-3,-4 I M w,1-r/(2 M)]+6 (2 M-r) HeunCPrime[-2+l+l^2-8 M w (I+2 M w),4 M w (I-4 M w),1-4 I M w,-3,-4 I M w,1-r/(2 M)])
RInHodd[l_,w_,r_]:=\[Psi]oIn[l,w,r,MVal]
dRInHodd[l_,w_,r_]:=1/r^2 2^(2 I w) E^(I (-4+r) w) (-2+r)^(-1-2 I w) ((4-2 r+2 I (-4+r) r w) HeunC[l+l^2-2 (1+4 I w+8 w^2),4 (I-4 w) w,1-4 I w,-3,-4 I w,1-r/2]-(-2+r) r HeunCPrime[l+l^2-2 (1+4 I w+8 w^2),4 (I-4 w) w,1-4 I w,-3,-4 I w,1-r/2])
RInHeven[l_,w_,r_]:=\[Psi]eIn[l,w,r,MVal]
dRInHeven[l_,w_,r_]:=(2^(1+2 I w) E^(I (-4+r) w) (-2+r)^(-2 I w) (((-((-1+l) (2+l) (-2+r) (-144+(-2+l+l^2) r (-36+l (1+l) (-2+l+l^2) r^2)))+I (-4+r) (6+(-2+l+l^2) r) (-144+r (72+(-1+l) l (1+l) (2+l) r (6+(-2+l+l^2) r))) w-12 r^3 (6+(-2+l+l^2) r)^2 w^2) HeunC[-2+l+l^2-8 w (I+2 w),4 (I-4 w) w,1-4 I w,-3,-4 I w,1-r/2])/((-2+r) (6+(-2+l+l^2) r)^2)+1/(12+2 (-2+l+l^2) r) (144+r (-72+(-1+l) l (1+l) (2+l) r (-6-(-2+l+l^2) r))) HeunCPrime[-2+l+l^2-8 w (I+2 w),4 (I-4 w) w,1-4 I w,-3,-4 I w,1-r/2]))/(r^3 ((-1+l) l (1+l) (2+l)-12 I w))
f[r_]:=1-(2 M)/r


(* \!\(
\*SubsuperscriptBox[\(K\), \(lm\), \(e/o\)]\ ; \ \(Cutoff\ close\ to\ the\ \(ISCO : \ 3 M\  < \ rc\  < \ 6 M\)\)\) *)


KOdd[l_,m_,w_,rG_]:=((1-(2M)/rG)^-1 RInHodd[l,w,rG] Aodd[l,m,rG]-(-(1-(2M)/rG)^-2 (2M)/rG^2 RInHodd[l,w,rG]+(1-(2M)/rG)^-1 dRInHodd[l,w,rG])Bodd[l,m,rG])/.M->MVal
KEven[l_,m_,w_,rG_]:=((1-(2M)/rG)^-1 RInHeven[l,w,rG] Aeven[l,m,rG]-(-(1-(2M)/rG)^-2 (2M)/rG^2 RInHeven[l,w,rG]+(1-(2M)/rG)^-1 dRInHeven[l,w,rG])Beven[l,m,rG])/.M->MVal


$Assumptions={r>2,r<6,Element[m,Integers],Element[w,Reals],Element[l,Integers],l>1};


(* ::Text:: *)
(*Integrand scales as (r - 6 M)^(-3/2) close to the ISCO*)


(* for 2M<r<rc *) 
IntCo1=1/rGdot[r] Exp[I w to[r]-I m \[Phi]o[r]]KOdd[l,m,w,r]/.M->MVal/.SubsSphTerms ;
IntCe1=1/rGdot[r] Exp[I w to[r]-I m \[Phi]o[r]]KEven[l,m,w,r]/.M->MVal/.SubsSphTerms ;


(* ::Section::Closed:: *)
(*2. Numerical integration scheme of \!\(\*SubsuperscriptBox[\(C\), \(lm\), \(e/o\)]\)(\[Omega]) *)


(* ::Input:: *)
(*(* Integration scheme tested for cutoff at rc=5.999, frequencies up to \[Omega] = 4 and l mode up to l=5 *)*)


(* Integration strategy for rcutoff=5.999M, l=4 *)
CVale[thisl_,thism_,thisw_,thisrc_]:=Block[{l=thisl,m=thism,w=thisw},
NIntegrate[IntCe1,{r,2,3},Method->{ "LocalAdaptive",Method->{ "GaussKronrodRule","Points"-> 50}}]+
NIntegrate[IntCe1,{r,3,4},Method->{ "LocalAdaptive",Method->{ "GaussKronrodRule","Points"-> 50}}]+
NIntegrate[IntCe1,{r,4,5},Method->{ "LocalAdaptive",Method->{ "LevinRule","Points"-> 50}}]+
If[thisl < 4,
	NIntegrate[IntCe1,{r,5,59/10},Method->{ "LocalAdaptive",Method->{ "LevinRule","Points"-> 50}}]+ 
	If[Abs[thisw-thism/6^(3/2)]<10^(-1),
		NIntegrate[IntCe1,{r,59/10,thisrc},Method->{ "LocalAdaptive",Method->{ "LevinRule","Points"-> 50}}],
		NIntegrate[IntCe1,{r,59/10,thisrc},Method->{ "GlobalAdaptive",Method->{ "LevinRule","Points"-> 50}}]
	]
	,
	NIntegrate[IntCe1,{r,5,thisrc},Method->{ "LocalAdaptive",Method->{ "LevinRule","Points"-> 50}}]
	]
]

CValo[thisl_,thism_,thisw_,thisrc_]:=Block[{l=thisl,m=thism,w=thisw},
NIntegrate[IntCo1,{r,2,3},Method->{ "LocalAdaptive",Method->{ "GaussKronrodRule","Points"-> 50}}]+
NIntegrate[IntCo1,{r,3,4},Method->{ "LocalAdaptive",Method->{ "GaussKronrodRule","Points"-> 50}}]+
NIntegrate[IntCo1,{r,4,5},Method->{ "LocalAdaptive",Method->{ "LevinRule","Points"-> 50}}]+
If[thisl < 4,
	NIntegrate[IntCo1,{r,5,59/10},Method->{ "LocalAdaptive",Method->{ "LevinRule","Points"-> 50}}]+ 
	If[Abs[thisw-thism/6^(3/2)]<10^(-1),
		NIntegrate[IntCo1,{r,59/10,thisrc},Method->{ "LocalAdaptive",Method->{ "LevinRule","Points"-> 50}}],
		NIntegrate[IntCo1,{r,59/10,thisrc},Method->{ "GlobalAdaptive",Method->{ "LevinRule","Points"-> 50}}]
	]
	,
	NIntegrate[IntCo1,{r,5,thisrc},Method->{ "LocalAdaptive",Method->{ "LevinRule","Points"-> 50}}]
	]
]


(* ::Section::Closed:: *)
(*3. Parallel computation of \!\(\*SubsuperscriptBox[\(C\), \(lm\), \(e/o\)]\)(\[Omega])*)


NModes=Length[ListModes];


(*Init*)
MaxN= Maxw NPoints;
\[Omega]Table=Table[Maxw (i-MaxN/2)/(MaxN/2),{i,0,MaxN}];
data=Table[0,{i,1,NModes},{j,0,MaxN}];
SetSharedVariable[data];


ParallelDo[
thisl=ListModes[[nn,1]];
thism=ListModes[[nn,2]];
thisw=\[Omega]Table[[compt]];
data[[nn,compt]]= If[thisw== 0,0,Module[{output},output=If[EvenQ[thisl+thism],CVale[thisl,thism,thisw,thiscutoff],CValo[thisl,thism,thisw,thiscutoff]];If[Not[$MessageList == {}],Print["Evaluation of CVal l=",thisl,"m=",thism," w=",thisw,": Encountered messages:",$MessageList]];output]/(2 I thisw Ain$BHPT[thisl,N[thisw,PrecisionAin]])];
,{compt,1,MaxN+1},{nn,1,NModes}];



Result=Table[Association@{"l"-> ListModes[[nn,1]],"m"->ListModes[[nn,2]],"rCutoff"-> thiscutoff,"wCutoff"-> Maxw,"FourierModes"->Table[{\[Omega]Table[[compt]], data[[nn,compt]]},{compt,1,MaxN+1}]},{nn,1,NModes}];


Print["Main computation ended"];


(* ::Section::Closed:: *)
(*4. Export result*)


Print["Exporting results"];


Export[FilenameWhereToExportFinalData,Result,"Table"];


Print["Results exported"];


(* ::Chapter::Closed:: *)
(*IV. Closing  up*)


CloseKernels[];
