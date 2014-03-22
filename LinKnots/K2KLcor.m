(* Mathematica Version: 4.0; 4.2 *)

(* Name: `K2KL` *)

(* Title: K2KL *)

(* Authors: M.Ochiai and N.Imafuji are the authors of Knot 2000 (K2K), 
and S. Jablan and R. Sazdanovic of LinKnot (K2KC) *)

(* Copyright: Copyright 2003 *)

(* History: current Version 1.1 at 1.05.2003 *)

(* Summary: The complete package K2K.m works with two packages: 
knotbycomp.m by M.Ochiai and N.Imafuji, and LinKnot.m  
by S. Jablan and R. Sazdanovic. The package LinKnot.m, 
contained in K2K.m, works with knots and links in Conway 
notation, with basic polyhedra with N <= 20 vertices, 
producing their pdata, Dowker codes and Millet codes, 
and calculating linking numbers and estimated unknotting 
and unlinking numbers. *)
  
(* Context: `LinKnots` *)

(* Keywords: Knot Theory, Knots, Links, Conway Notation, Basic Polyhedra *)

(* Requirements: DiscreteMath`Combinatorica, knotbycomp.m ,
"LinearAlgebra`MatrixManipulation`",
"DiscreteMath`ComputationalGeometry`"*)

BeginPackage["LinKnots`", "DiscreteMath`Combinatorica`",
"DiscreteMath`ComputationalGeometry`",
"LinearAlgebra`MatrixManipulation`","PolyBase`",
"KnotLinkBase`","NalgBase`"]; 

 Get["knotbycomp.m"];
 Get["RecBase.txt"];
 Get["Classic.txt"];
(* USAGE MESSAGES *)
            
fCreatePData::usage = 
"fCreatePData[Conway_String] calculates pdata for 
given knot or link given in Conway notation."

Dow::usage = 
"Dow[Conway_String] produces Dowker code of a knot or link K.
The first data of result are lengths of components, and the second Dowker 
code with the signs of crossing points."

fGaussExt::usage="fFormDL[Ul_] produces extended Gauss code with signs 
of knot or link given by Conway symbol or Dowker code." 

LinkingNo::usage="LinkingNo[Ul_] computes the linking number of a 
link given by Conway symbol or Dowker code."

fGenSign::usage =
"fGenSign[Conway_String] returns the signlist of crossing points 
of a KL given by Conway symbol."

fDowfromPD::usage = 
"fDowfromPD[pdata_List] computes Dowker code from pdata."

fMillett::usage="fMillett[Conway_String] writes in the file 
izlaz11.txt the Millet code of a knot or link and calculates
polynomial invariants, including Homfly polynomial."

UnKnotLink::usage="UnKnotLink[Ul_] returns an (estimated)
unknotting or unlinking number of a given knot or link given 
by its Conway symbol, Dowker code or pdata. For the 
reduction of knots or links it uses ReductionKnotLink."

RatReduce::usage="RatReduce[Conway_String] reduces 
Conway symbol of a rational KL."

UnR::usage="UnR[Conway_String] returns an estimated unknotting 
or unlinking number of a given rational knot or link. It is
computed by using continued fractions."

fCreateGraphics::usage="fCreateGraphics[Conway_String] writes in 
the file c:\\Program Files\\KnotPlot\\graphics.txt the coordinates 
of a knot or link. The file obtained can be loaded in the program 
KnotPlot by writing load graphics.txt in the KnotPlot command line."
        
GetKnotLink::usage="GetKnotLink[Str_String, No_Integer].
Str is the name of a list of alternating or non-alternating 
knots and links with a specified number of crossings (for example, 
a10 or n10, where a stands for alternating and n for non-alternating 
knots and links) and No is the number of a particullar desired knot ot 
link from that list. As a result, GetKnotLink returns its Conway symbol."

NumberOfKL::usage="NumberOfKL[Str_String] Str is the name of a list
of knots and links with a specified number of crossings 
(for example, a10 or n10, where a stands for alternating and n 
for non-alternating knots and links, written as a string). The 
function NumberOfKL gives the number of alternating or 
non-alternating knots and links with a specified number of 
crossings."

Dowker::usage ="Interior function"
fConvert::usage="Interior function"
fConvertPoly::usage="Interior function"
fBlock::usage="Interior function"
fComma::usage="Interior function"
fPlus::usage="Interior function"
fInv::usage="Interior function"
fSpace::usage="Interior function"
fRef::usage="Interior function"
fGen::usage="Interior function"
fMoja::usage="Interior function"
fDowkerCode::usage="Interior function"

fDToD::usage="fDtoD[Con_String] calculates Dowker code 
of a KL given in Conway notation."

fProjections::usage="fProjections[Con_String] calculates 
Conway symbols of all projections of a KL given by its Conway symbol."

fFixP::usage=""
fChangePart::usage=""


fDiffProjectionsAltKL::usage="fDiffProjectionsAltKL[Con_String] calculates 
Conway symbols and minimal Dowker codes of all non-isomorphic 
projections of a KL given by its Conway symbol."

MinDowProjAltKL ::usage="MinDowProjAltKL[Ul_] calculates minimal 
Dowker code without signs (in Knotscape form) for a given 
alternating KL projection. An input is its Conway symbol, Dowker code, 
or pdata."

MinDowAltKL::usage="MinDowAltKL[Con_String] calculates 
minimal Dowker code of any alternating KL given by its Conway symbol."

SameAltProjKL::usage="SameAltProjKL[Con1_String,Con2_String] compares two 
alternating KL projections given by Conway symbols, 
Dowker codes, or pdata. The result is 1 
for equal, and 0 for non-equal projections."

SameAltConKL::usage="SameAltConKL[Con1_String,Con2_String] compares 
two alternating KLs given by Conway symbols. 
The result is 1 for equal, and 0 for non-equal KLs."

fFindCon::usage="fFindCon[Ul_] finds Conway symbol of 
any alternating KL with at most 12 vertices 
given by its Dowker code with signs, or by pdata."

fAlexPoly::usage="fAlexPoly[Conway_String] calculates multi-variable 
Alexander polynomial of a KL from its Conway symbol 
by using Writinger presentation (generators of KL)."

fSeifert::usage="fSeifert[Ul_] calculates Seifert matrix of a 
KL given by its Conway symbol, Dowker code, 
or pdata. Its basic function SeifertMatrix  
is written by Stepan Orevkov."

SeifertMatrix::usage="Interior function"

fSignat::usage="fSignat[Ul_] calculates signature of a KL 
given by its Conway symbol, Dowker code, or pdata. 
Its basic function ssmW is written by Stepan Orevkov."

fOrientedLink::usage="fOrientedLink[Ul_] calculates Gauss 
codes for a link projection given by Conway symbol, 
Dowker code, or pdata. The Gauss codes obtained correspond 
to different orientations of components (the first part of data 
obtained). The second part gives the orientations of components 
(where + is denoted by 1, and - by 0), and the third part is the 
(signed) linking number of the corrsponding oriented link."

fGaussExtSigns::usage="fGaussExtSigns[Ulaz_] calculates 
Gauss code with signs for a link projection 
given by its Conway symbol, Dowker code, or pdata."

fGraphInc::usage="fGraphInc[Ul_] calculates from a 
Conway symbol, Dowker code, or pdata of a KL, the 
corresponding graph of KL given by edges (unordered pairs), 
and by the list of vertex signs."

fPrimeGraph::usage="fPrimeGraph[GG_List] tests that a graph given 
by a list of unordered pairs, after transforming it 
into the corresponding alternating KL, will result in a  
prime or composite alternating KL. The output is 1 for a prime, 
and 0 for a composite KL."

fPlanarEmbKL::usage="fPlanarEmb[Ul_] calculates the planar 
embedding of a prime KL given by Conway symbol, pdata or 
Dowker code. An output is the list that consists of the graph of 
input KL, its planar embedding, and faces of planar 
embedded graph. As the basis of 
this program it is used the external program 
planarity.exe written by J.M.Boyer."

fPlanarEmbGraph::usage="fPlanarEmbGraph[LP_] gives the planar embedding 
of a 3-connected planar graph given by a list of unordered pairs. 
An output is the list that consists of the graph 
given by a list of unordered pairs, its
planar embedding given by vertex cycles, 
and faces of planar embedded graph. As the basis of 
this program it is used the external program planarity.exe written by \
J.M.Boyer."

DrawPlanarEmbKL::usage="DrawPlanarEmbKL[Ul_] draws the planar embedding 
of a prime KL given by Conway symbol, Dowker code, or pdata, without 
showing digonal faces. As the basis of this program it is used the 
program 3-Dimensional Convex Drawings of 3-Connected Planar Graphs 
by M.Ochiai, N.Imafuji and N.Morimura."

DrawPlanarEmbGraph::usage="DrawPlanarEmbGraph[Ul_] draws the planar 
embedding of a 3-connected graph given by the list of unordered pairs.  
As the basis of this program it is used the program 3-Dimensional 
Convex Drawings of 3-Connected Planar Graphs by M.Ochiai, N.Imafuji 
and N.Morimura."

fMidEdgeGraph::usage="fMidEdgeGraph[L_List] gives the graph 
defined by mid-edge points of a given polyhedral graph G 
given by a list of unordered pairs of vertices."

fKLfromGraph::usage="fKLfromGraph[L_List] gives the KL
defined by a graph G given as a list of unordered pairs."

fKLinGraph::usage="fKLinGraph[UOPair_List] gives the graph 
defined by mid-edge points of a polyhedral graph G 
given by a list of unordered pairs of vertices."

fAddDig::usage="fAddDig[GL_List] produces from a given graph 
all 4-regular non-isomorphic graphs by replacing single edges by  
double (digonal) edges."

fSignsKL::usage="fSignsKL[Dow_List] calculates Dowker code with signs 
of a KL given by its  Dowker code in Knotscape form. 
For an alternating knot, an input is Dowker code without signs, 
and for a nonalternatng KL, an input is Dowker code containing 
only signs of points with signs changed with regard to the 
corresponding alternating KL."

fGraphKL::usage="fGraphKL[Ulaz_] calculates and draws graph  
of a KL given by Conway symbol, Dowker code, or pdata.  
The result is the graph given by the list of unordered pairs."

fGenerators::usage="fGenerators[Ul_] calculates 
generators of a KL and the list of the signs corresponding 
to crossing points for a KL given by Conway symbol,  
Dowker code, or pdata. The first generator incommes to  
a vertex, the seccond is a the outgoing generator, and the  
third is the incomming generator (IOP=Incomming-Outgoing-Passing).  
The result is divided according to the components of KL."

fColTest::usage="fColTest[Ul_,cNo_Integer] calculates from a given 
Conway symbol, Dowker code,  pdata, and from a number of colors 
(greater then 2), a coloring of KL. The result is the list of colored 
generators according to the generator list given by  fGenerators. 
The resuct can be a perfect coloring, where in every crossing all 
generators have different color, standard coloring where all generators 
have different colors or the same color, and Kauffman coloring where 
the only forbiden case is that incomming and passing generator, or 
outgoing and passing generator have the same color, different from 
the third generator."

fPDataFromDow::usage="fPDataFromDow[Dow_List] calculates pdata 
from Dowker code in Knotscape format (without signs for an 
alternating prime KL, or from Dowker code with signs of changed 
crossings for nonalternating prime KL)."

fPDataFromDowker::usage="fPDataFromDowker[Dow_List] calculates pdata 
from Dowker code with signs."

fKnotscapeDow::usage="fKnotscapeDow[Con_String] calculates from a 
Conway symbol of a KL its Dowker code in a the Knotscape format: 
Dowker code without signs for alternating KL, or Dowker code with 
signs of changed crossings for nonalternating KL."

fPrimeKL::usage="fPrimeKL[Ul_] checks that KL given by Dowker code, 
or pdata is a prime or composite KL, i.e. a direct product of some 
prime KLs. The result is 1 for a prime KL, and 0 for a composite KL."

fTorusKL::usage="fTorusKL[a_Integer,b_Integer] calculates for a 
torus KL [a,b] its braid word, min. number of crossings, unknotting 
number or number of components, bridge number, Alexander polynomial 
and (Murasugi) signature (see: Murasugi K., Knot Theory and its 
Applications, Birkhauser, Boston, Basel, Berlin, 1996)."

fComponentNo::usage="fComponentNo[Ul_] calculates the number of 
components of KL given by Conway symbol, Dowker code, or pdata."

fBreakComp::usage="fBreakComp[Ul_,k_Integer] in a link given by 
Conway symbol, Dowker code, or pdata, and by ordering number 
k of a cutted component calculates the pdata of a link with 
k-th component cutted."

BreakCoAll::usage="BreakCoAll[Ul_] in a link given by its Conway 
symbol, Dowker code, or pdata, cutts all components. The results 
are all different KLs obtained by cutting all components, where 
for unknot (unlink) the result is {0}."

CuttNo::usage="CuttNo[Ul_] calculates the cutting number of a  
link given by its Conway symbol, Dowker code, or pdata."

NoSelfCrossNo::usage="U-infinity number is the minimum number of 
cuts in self-crossing points of a KL necessary to obtain unknot, 
or a link without self-crossings. The function NoSelfCrossNo[Ul_] 
calculates the U-infinity number (see: Jablan S.: Unknotting number 
and Infty-Unknotting Number of a Knot, Filomat (Nis), 12:1 (1998), 
113-120) of a  KL given by its Conway symbol, Dowker code, or pdata."

fCuttRealKL::usage="fCuttRealKL[UL_] calculates the number of 
real cuttings of a given projection of KL given by its Conway 
symbol, Dowker code, or pdata, i.e., the number of cuttings 
with a different cutting point in the projection. The result is 
the number of different real cutting classes with preserved 
signs, or with preserved or reversed signs."

SplittNo::usage="SplittNo[PData_] calculates splitting number 
of a link given by its Conway symbol, Dowker code, or pdata, i.e., 
the minimum number of crossing changes necessary in order to obtain a \
splitted 
link. It is illustrated by the example from Adams. (see: Adams C: 
Splitting Versus Unlinking, Journal of Knot Theory and Ramifications, 
Vol.5, No. 3 (1996) 295-299)."

AmphiProjAltKL::usage="AmphiProAltjKL[Ul_] tests the amphicheirality of 
a given projection of an alternating KL given by its Conway symbol, 
Dowker code, or pdata."

PeriodProjAltKL::usage="PeriodProjAltKL[Ul_] calculates the period 
of a given projection of an alternating knot given by its 
Conway symbol, Dowker code, or pdata."

PeriodAltKL::usage="PeriodAltKL[Conway_String] calculates the 
period of a given alternating KL given by its Conway symbol."

AmphiAltKL::usage="AmphiAltKL[Conway_String] tests the 
amphicheirality of an alternating KL given by its Conway symbol."

Symm::usage="Symm[Ul_] calculates all automorphisms of an 
alternating KL given by its Conway symbol, Dowker code, 
or pdata. The first part of the result is the list of 
automorphisms given by permutations, and the other the 
list of the corresponding cycles."

MaxSymmProjAltKL::usage="MaxSymmProjAltKL[Con_String] finds maximum 
symmetrical projection of an alternating KL given by its Conway symbol."

fBasicPoly::usage="fBasicPoly[Ul_List] finds the basic polyhedron 
for a KL given by Dowker code, or pdata."

fCompositePoly::usage="fCompositePoly[Conway_String] checks that basic 
polyhedron given by its Conway symbol is elementary or composite.  
The result is 1 for composite, and 0 for elementary basic polyhedron."

fPolyFlype::usage="fPolyFlype[Conway_String] checks that basic 
polyhedron given by its Conway symbol permits flypes or not.  
The result is 1 for basic polyhedra permiting flypes, and 0 for the others."

fAllStatesProj::usage="fAllStatesProj[UL_] calculates all states of 
a given alternating KL, i.e., all different variations of signs in 
a given KL projection. Projection can be given by its Conway symbol, 
Dowker code, or pdata. The result is the list of reduced KL obtained 
from it by all sign-changes. The first datum is  of crossing changes, 
and the second the list of all KL obtained."

RationalKL::usage="RationalKL[n_Integer] calculates the number and 
Conway symbols of all rational KLs for a given number of crossings n."

RationalAmphiK::usage="RationalAmphiK[n_Integer] calculates the number 
and Conway symbols of all rational amphicheiral knots for a given 
number of crossings n."

RationalAmphiL::usage="RationalAmphiL[n_Integer] calculates the number 
and Conway symbols of all rational amphicheiral links for a given 
number of crossings n."

RatGenSourKL::usage="RatGenSourKL[n_Integer, m_Integer] calculates the 
number and Conway symbols of all rational generating KLs (for m=3) 
and rational source KLs (for m=2) with n crossings."

RatSourceKLNo::usage="RatSourceKLNo[n_Integer] calculates the 
number of rational source KLs (for k=2) with n crossings 
according to the general recursion formula: b[0]=1, b[1]=1, 
b[2n-2]+b[2n-1]=b[2n], b[2n]+b[2n-1]-f[n-1]=b[2n+1], where 
f is the Fibonacci sequence given by recursion 
f[0]=1, f[1]=1, f[n-2]+f[n-1]=f[n]."

RatLinkU1::usage="RatLinkU1[n_Integer] calculates the number 
and Conway symbols of all rational links with the unlinking 
number 1 with n crossings."

RatLinkU0::usage="RatLinkU0[n_Integer] calculates the number 
and Conway symbols of all rational unlinks (U=0) with n crossings."

RatKnotGenU0::usage="RatKnotGenU0[n_Integer] gives Conway symbols 
of all rational links with the unlinking number 0 with n crossings."

RatKnotGenU1::usage="RatKnotGenU1[n_Integer] gives Conway symbols 
of all rational knots with the unknotting number 1 with n crossings."

MSigRat::usage="MSigRat[Conway_String] calculates the fraction 
and Murasugi signature of a rational KL given by its Conway symbol."

AllStatesRational::usage="AllStatesRational[Conway_String] calculates 
all states of a given rational KL given by its Conway symbol, i.e., 
all different variations of signs in a given projection of a rational 
KL. The first datum in the result is the unknotting (unlinking) 
number obtained from the fixed projection of a given rational KL,  
followed by  the list of minimum numbers of crossing changes necessary 
to obtain the corresponding rational KL given by its Conway symbol, 
i.e., the list of the KL distances from a given rational KL projection."

AllStatesRatFast::usage="Interior function"

fDirectProduct::usage="Interior function"

fDToDDirect::usage="fDToDDirect[ConDir_String] computes Dowker code 
of a direct product of KLs given in Conway notation."

fGenSignDirProd::usage="fGenSignDirProd[Conway_String] computes 
signs of a direct product of KLs given in Conway notation."

fConvertDirect::usage="Interior function"

fClassicToCon::usage="fClassicToCon[Ulaz_String] gives Conway 
symbol of a KL given in classical notation (according to D. Rolfsen)."

fGetBraidRepresent::usage="fGetBraidRepresent[Ul_] calls the 
external program Braids-9.0 by A. Bartholomew. It returns a B-word of 
a braid representation for the 
given Conway symbol or P-data."

fGaussExtSignsBraid::usage=""

BraidWReduced::usage="BraidWReduced[PD_List] gives a 
reduced braid word of a KL given by its pdata."

fEdmonds::usage="fEdmonds[UlL_List] gives list of polygons, 
characteristic and genus for a graph given as a list of 
unordered pairs."

fStellarBasic::usage="fStellarBasic[n_Integer] produces the list 
of basic stellar KLs with n crossings."

fStellar::usage="fStellar[n_Integer] produces the list 
of stellar KLs with n crossings."

fCompl::usage="Interior function"
fStelRot::usage="Interior function"
fGener::usage="Interior function"
fMakeStel::usage="Interior function"
fPlanarEmb::usage="fPlanarEmb[neu_List] calculates 
the planar embedding of a prime KL given by 
Conway symbol, pdata or  Dowker code. An output is 
the list that consists of the graph of  input KL, 
its planar embedding given  by vertex cycles, and 
faces of planar embedded graph. As the basis of  
this program it is used the external MS DOS 
program planarity.exe written by J.M.Boyer."

fStellarPlus::usage="fStellarPlus[n_Integer] produces the list 
of distinct stellar KLs with n crossings and with pluses."

fStellarNalt::usage="fStellarNalt[n_Integer] produces the list 
of nonalternating stellar KLs with n crossings."

NMoveRat::usage="NMoveRat[n_Integer,Conway_String] calculates 
for a given integer n and Conway symbol of a KL the minimum 
number of n-moves neccessary to unknott a minimal projection 
of a rational KL."

fJablanPoly::usage="fJablanPoly[Conway_String] calculates 
multivariable Jablan polynomial for a projection of KL 
given by its Conway symbol."

fGener::usage="Interior function"

fGapRat::usage="fGapRat[Conway_String] finds rational KLs with 
an unlinking gap. The result is the unlinking gap, followed by 
the unlinking number of KL and the unlinking number of its 
fixed minimal projection."

fUnKLFixed::usage="fUnKLFixed[UL_,n_Integer] checks that a 
fixed projection of KL can be unknotted (unlinked) by n crossing 
changes. The result is {} for a KLs that cannot be unknotted 
(unlinked) by n crossing changes, and n for the others."

fGap::usage="fGap[Conway_String] finds KLs with the unlinking gap. 
The first part of result is 1 for KLs with an unlinking gap, 
and 0 otherwise, and second part is unlinking number."

fUnRFixProj::usage="fUnRFixProj[Conway_String,k_Integer] checks that 
a fixed projection of a rational KL can be unknotted (unlinked) 
by k crossing changes. The result is 1 for a KLs that can be 
unknotted (unlinked) by k crossing changes, and 0 for the others."

fAutoSigInp::usage="fAutoSigInp[Gr_List] calculates stable 
states of an automaton, given by a list of outgoing edges, 
signs of vertices and inputs in vertices. In the vertices 
with a sign 1 is used the operation NOR, and in the 
vertices -1 the operation OR. The result is the list of 
edge colorings corresponding to stable states and the 
list of stable states according to vertices (see Kauffman). 
If the signlist is empty, it is computed as {1,...,1}, 
i.e. with NOR in all vertices."

fAutoKL::usage="fAutoKL[Ul_List] calculates the stable 
states of an automaton obtained from a KL given in Conway 
notation, followed by a list of signs of vertices, and 
inputs in vertices. In the vertices with a sign 1 is 
used the operation NOR, and in the vertices -1 the 
operation OR. The result is the list of edge colorings 
corresponding to stable states and the list of stable states 
according to vertices (see Kauffmann). If the signlist 
is empty, it is computed as the signlist of the given KL."

LiangPoly::usage="LiangPoly[Conway_String] calculates Liang 
polynomial of a KL projection given by its Conway symbol 
and distinguishes amphicheiral 
KL projections with a Liang polynomial satisfying the 
relationship L(t)=L(1/t)."

fBoolean::usage="fBoolean[LL_] calculates the square-free 
polynomial for a given logical term."

fDiffSeq::usage="fDiffSeq[n_Integer,p_Integer] calculates 
all different periodic sequences with n variables, 
of period p. As the result, every sequence is given 
by its code, Kauffman code, and square-free polynomials 
corresponding to it."

fBalanced::usage="fBalanced[Diff_List,redBr_Integer] 
calculates stable (balanced) states for 
a given periodic sequence."

fDowCodes::usage="fDowCodes[n_Integer] 
generates all different knot and link projections
with n crossings, including non-prime KLs. It uses
external program plantri.exe by B.McKay and G.Brinkmann."

fGenKL::usage="fGenKL[b_Integer, pp_String] generates 
all different alternating KLs with n crossings 
from a given source KL."

fMinDowKnotPr::usage=""

fWrGraph::usage=""

fPlanEmbedding::usage=""

fForknotFind::usage=""

fKnotFind::usage=""

fGapKnotsc::usage=""

fVarNewP::usage=""

RedKauffmanPolynomial::usage=""

fProdTangles::usage="fProdTangles[s1_String, s2_String, tt_String] 
gives product of non-algebraic tangles m5*,m7*,m81*-m82*,m91*-m96*,
m101*-m1011*,m111*-m1138* with an algebraic tangle tt placed in it."

fSumTangles::usage="fSumTangles[s1_String, s2_String, tt_String] 
gives sum of non-algebraic tangles m5*,m7*,m81*-m82*,m91*-m96*,
m101*-m1011*,m111*-m1138* with an algebraic tangle tt placed in it."

ListOfOneFactors::usage:="ListOfOneFactors[n] gives the list of all 
possible one-factors in a cubic graph with a given Hamilton cycle 
on n vertices. Some combinations can be isomorphic."

fDiffViae::usage="fDiffViae[n_Integer] gives all different
viae with n mirrors"

fViaToKL::usage="fViaToKL[LL_List] from given via 
computes Dowker code of corresponding KL."

AdmissibleEdge::usage=""

ConnectedComponents::usage =""

DFS::usage =""

Begin["`Private`"]

(*######### Implementation of +,",","_",- #############*)

fRef[{L1_List,{{a_,b_},{a1_,b1_}},L3_List}]:={L1,{{a,a1},{b,b1}},-L3}
fBlock[x_,n_Integer]:={{},{{4n-3,4n},{4n-2,4n-1}},{1}}
fInv[{L1_List,{{a_,b_},{a1_,b1_}},L3_List}]:={L1,{{a,b},{a1,b1}},-L3}
fPlus[L_List]:=
  Module[{s,i,l,p,d,p1,d1},l=Length[L];p=L[[1]][[2]][[1]];d=L[[1]][[2]][[2]];
    p1=L[[2]][[2]][[1]];d1=L[[2]][[2]][[2]];
    If[l>2,s=L[[1]];  Do[s=fPlus[{s,L[[i]]}],{i,2,l}],
      s={Join[Join[L[[1]][[1]],
              L[[2]][[1]]],{{p[[2]],p1[[1]]},{d[[2]],d1[[1]]}}],{{p[[1]],
              p1[[2]]},{d[[1]],d1[[2]]}}, Join[L[[1]][[3]],L[[2]][[3]]]}];
              s
    ] 
fComma[L_List]:=fPlus[Map[fRef,L]]
fSpace[L_List]:=Module[{s,i,l},l=Length[L];
    If[l>2,s=L[[1]];  Do[s=fSpace[{s,L[[i]]}],{i,2,l}],
      s=fPlus[{fRef[L[[1]]],L[[2]]}]];
    s
    ] 
    
 (* ######## Conversion from Conway Notation ############ *) 
   
fConvert[ sInput_String, opts___Rule ] :=
  Module[ {sOut,sCh,iI,iJ,iLevel,iPos,iPrev,iNext} ,
     sCh = Characters[ sInput ] ;
    lComm = {};     
    For[ iI=1, iI<=Length[sCh], iI++,  
      If[ sCh[[iI]]//DigitQ ,        
        While[ iI+1<=Length[sCh] && DigitQ[sCh[[iI+1]]] ,
           sCh[[iI]] = sCh[[iI]] <> sCh[[iI+1]];
           sCh       = Drop[ sCh , {iI+1} ];
        ];
        AppendTo[ lComm , 4 ];  
      , AppendTo[ lComm ,
          Switch[ sCh[[iI]],
               "(" , 2  ,          
               ")" , -2 ,          
               "+" , 5  ,         
               " " , 6  ,         
               "," , 7  ,         
               "-" , 8  ,          
               _   , 0            
          ]
        ];
      ];
    ];
    
    If[ Position[ lComm, 0 ] =!= {} , Return[False] ]; 
    iLevel = 0;               
    Do[        
      If[ Abs[lComm[[iI]]]===2 , iLevel += lComm[[iI]] ];
    ,{iI,Length[sCh]}];
    If[ iLevel =!= 0 , Return[ False ] ];  
   
    For[ iPos=iI=1, iI<=Length[sCh], iI++,  
      If[ lComm[[iI]] === 4 ,              
        sCh[[iI]] = "fBlock[" <> sCh[[iI]] <> "," <> ToString[iPos] <> "]"; 
        iPos++;                                  
        lComm[[iI]] = 1 ;       
      ];
    ];
   
    If[ (Print /. {opts})===On,
      Print["blocks: ",Thread[{sCh,lComm}]//MatrixForm];
    ];
     
    For[ iI=Length[sCh]-1, iI>0, iI--, 
      If[ lComm[[iI]]===8 ,
         If[ lComm[[iI+1]] == 1 ,
            (* single block *)
            sCh[[iI]]   = "fInv[" <> sCh[[iI+1]] <> "]";
            lComm[[iI]] = 1;
            sCh   = Drop[ sCh   , {iI+1} ];
            lComm = Drop[ lComm , {iI+1} ];
         , sCh[[iI]]   = "fInv[" <> sCh[[iI+1]];
            lComm[[iI]] = 2;
            sCh    = Drop[ sCh   , {iI+1} ];
            lComm  = Drop[ lComm , {iI+1} ];
            iLevel = 0;
            iJ     = iI+1;
            While[ iLevel>-2 ,   
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ++;
            ];
            sCh[[iJ-1]]   = sCh[[iJ-1]] <> "]";
         ];
      ];
    ];
  
    If[ (Print /. {opts})===On,
      Print["operation '-': ",Thread[{sCh,lComm}]//MatrixForm];
    ];
 
    For[ iI=2, iI<Length[sCh], iI++, 
      If[ lComm[[iI]]===6,
        iPrev = iI-1;              
        If[ lComm[[iI-1]]===-2 , 
            iLevel = 0;
            iJ     = iI-2;
            While[ iLevel<2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ--;
            ];
            iPrev = iJ+1 ;
        ];
        sCh   =  Insert[ sCh  , "fSpace[{", iPrev ];
        lComm =  Insert[ lComm, 2         , iPrev ]; 
        iNext = iI+1;  
       
        While[ iNext<=Length[sCh] && lComm[[iNext]]===6,
          sCh[[iNext]]   = ",";
          lComm[[iNext]] = 9  ;
          iNext++;              
          If[ lComm[[iNext]]===2 ,    
            iLevel = 0;
            iJ     = iNext+1;
            While[ iLevel>-2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ++;
            ];
            iNext = iJ-1;
          ];
          iNext++;   
        ];
        sCh   =  Insert[ sCh  , "}]" , iNext    ];
        lComm =  Insert[ lComm, -2   , iNext    ]; 
      ];
    ];
    
    If[ (Print /. {opts})===On,
      Print["operation ' ': ",Thread[{sCh,lComm}]//MatrixForm];
    ];
    
    For[ iI=2, iI<Length[sCh], iI++,  
      If[ lComm[[iI]]===7,
      
        iPrev = iI-1;              
        If[ lComm[[iI-1]]===-2 ,   
            iLevel = 0;
            iJ     = iI-2;
            While[ iLevel<2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ--;
            ];
            iPrev = iJ+1 ;
        ];
        sCh   =  Insert[ sCh  , "fComma[{", iPrev ];
        lComm =  Insert[ lComm, 2         , iPrev ]; 
         iNext = iI+1;  
       
        While[ iNext<=Length[sCh] && lComm[[iNext]]===7,
          sCh[[iNext]]   = ",";
          lComm[[iNext]] = 9  ;
          iNext++;              
          If[ lComm[[iNext]]===2 ,  
            iLevel = 0;
            iJ     = iNext+1;
            While[ iLevel>-2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ++;
            ];
            iNext = iJ-1;
          ];
          iNext++;   
        ];
        sCh   =  Insert[ sCh  , "}]" , iNext    ];
        lComm =  Insert[ lComm, -2   , iNext    ]; 
      ];
    ];
   
    If[ (Print /. {opts})===On,
      Print["operation ',': ",Thread[{sCh,lComm}]//MatrixForm];
    ];
    
    For[ iI=2, iI<=Length[sCh], iI++,
      If[ lComm[[iI]]===5,
         iPrev = iI-1;             
        If[ lComm[[iI-1]]===-2 ,   
            iLevel = 0;
            iJ     = iI-2;
            While[ iLevel<2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ--;
            ];
            iPrev = iJ+1 ;
        ];
        sCh   =  Insert[ sCh  , "fPlus[{", iPrev ];
        lComm =  Insert[ lComm, 2       , iPrev ];
        iNext = iI+1;   
        
        sCh[[iNext]]   = ",";
        lComm[[iNext]] = 9  ; 
        iNext++;            
        If[ lComm[[iNext]]===2 ,   
            iLevel = 0;
            iJ     = iNext+1;
            While[ iLevel>-2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ++;
            ];
            iNext = iJ-1;
        ];
        iNext++;    
        sCh   =  Insert[ sCh  , "}]" , iNext    ];
        lComm =  Insert[ lComm, -2   , iNext    ];
      ];       
    ];
   
    If[ (Print /. {opts})===On,
      Print["operation '+': ",Thread[{sCh,lComm}]//MatrixForm];
    ];

    sOut = StringJoin[ sCh ];
    If[ SyntaxQ[sOut]===False, Return[ False ]];    
    ToExpression[ sOut ]                            
     ]
     
 (*### Conversion from Conway Symbols for Polyhedral Knots and Links ####*)
     
fFormVList[Out_List]:= Module[{ll = {}, i},
    Do[Switch[Out[[i]],
               -1,ll=Append[ll,{4i-2,4i-1,4i,4i-3}],
                1,ll=Append[ll,{4i-3,4i,4i-1,4i-2}]]
         , {i, 1, Length[Out]}];
    ll
        ]
        
fMoja[{x_?NonNegative,y_?NonNegative}, br_]:={x+br,y+br}
fMoja[{x_?NonNegative,y_?Negative}, br_]:={x+br,y}
fMoja[{x_?Negative,y_?NonNegative}, br_]:={x,y+br}
fMoja[{x_,y_}, br_]:={x,y}

fConvertPoly[Conway_String]:=
  Module[{OutG=0,str,Con=Conway,Vert,br=0,SignList={},el,min,p,ll,pos,x6,x61,
      x8},pos=Flatten[StringPosition[Con,"*"]];
    If[pos=={},
      If[StringTake[Con,1]==".",OutG="x61";Con=StringDrop[Con,1],
        OutG="x6"],OutG="x"<>StringTake[Con,{1,pos[[1]]-1}];
      If[StringLength[Con]==pos[[1]],Con="",
        Con=StringTake[Con,{pos[[1]]+1,StringLength[Con]}]]];
    OutG=ToExpression[OutG];
    If[Head[OutG]==List,Vert=fFormVList[OutG[[3]]];
      If[Con!="",
        If[StringTake[Con,1]==".",Con=StringInsert[Con,"1",1]]];
      While[StringPosition[Con,".."]!={},
        Con=StringReplace[Con,".."->".1."]];
      If[Con=="",Do[Con=Con<>"1.",{i,1,Length[OutG[[3]]]}],
        pos=Union[Flatten[StringPosition[Con,"."]]];
        Do[Con=Con<>".1",{i,1,Length[OutG[[3]]]-1-Length[pos]}];
        Con=Con<>".";
        If[StringTake[Con,1]==".",Con=StringDrop[Con,1]];];
      pos=Union[Flatten[StringPosition[Con,"."]]];
      Do[str=StringTake[Con,{1,pos[[i]]-1}];
        Con=StringDrop[Con,{1,pos[[i]]}];
        pos=pos-pos[[i]];
        (*Print["str ",str];*)If[
          StringTake[str,-1]=="0"&&StringPosition[str," \
"]!={},
          el=fRef[fConvert[StringDrop[str,-2]]],el=fConvert[str]];
        el[[1]]+=br*4;el[[2]]+=br*4;
        br=br+Length[el[[3]]];
        SignList=Append[SignList,el[[3]]*OutG[[3,i]]];
        min=Min[Vert[[i]]]+4;
        p=4*Length[el[[3]]]-4;
        If[p!=0,ll=OutG[[2]];ll=Map[(#-min)&,ll];
          ll=Map[fMoja[#,p]&,ll];
          ll=Map[(#+min)&,ll];OutG[[2]]=ll];
        ll={Vert[[i]][[1]]->el[[2]][[1]][[1]],
            Vert[[i]][[2]]->el[[2]][[1]][[2]],
            Vert[[i]][[3]]->el[[2]][[2]][[2]],
            Vert[[i]][[4]]->el[[2]][[2]][[1]]};
        OutG[[2]]=ReplaceAll[OutG[[2]],ll];
        Do[Vert[[j]]+=p,{j,i+1,Length[Vert]}];
        OutG[[2]]=Join[OutG[[2]],el[[1]]],{i,1,Length[OutG[[3]]]}];
      OutG[[2]]=Sort[OutG[[2]]];
      OutG[[3]]=Flatten[SignList]];
    OutG]

(*########## #############*)
 
 (* ########## ############# *)
 
fFun[k_] := Module[{l},
    If[Mod[k, 4] == 1 \[Or] Mod[k, 4] == 2, l = IntegerPart[k/4] + 1, 
      l = -(IntegerPart[(k - 1)/4] + 1)];
    l]

fGraphHa[ConSym_, OutLinks_List, n_] :=
  	Module[{InLinks = {}, pom, p, i, ModSym, KLGraph, HL},
    ModSym = Map[Mod[#, 2] &, ConSym];
    Do[If[ModSym[[i]] == 1, 
        InLinks =  Union[Flatten[Map[Append[InLinks, #] &,
         {{4i - 3, 4i - 1}, {4i - 2, 4i}}],1]], 
        InLinks =  Union[Flatten[ Map[Append[InLinks, #] &,
         {{4i - 3, 4i},{4i - 2, 4i - 1}}],1]]]
     , {i, 1, Length[ModSym]}];
    KLGraph = FromUnorderedPairs[Join[InLinks, OutLinks]];
    HL = ConnectedComponents[KLGraph];
    If[Length[HL] == 1,(*KNOT*)HL = Map[fFun[#] &, HL[[1]]]; 
      HL = Take[HL, {1, Length[HL], 2}],
       (*LINK*)pom = {};
        Do[p = Map[fFun[#] &, HL[[i]]];
          	p = Take[p, {1, Length[p], 2}];pom = Append[pom, p]
       , {i, 1, Length[HL]}];
       HL = pom];
    {ModSym, HL} ]

fFormElem[p_, k_, q_] := Module[{i, LL = {}},
    If[p > 0, Do[LL = Append[LL, i], {i, k + 1, k + q}], 
      Do[LL = Append[LL, k + q + 1 - i], {i, 1, q}]];
    LL]
(* ##############fConSymPoly[Conway_String]###########*)
(*nalazi "pun" Konvejev simbol poliedra*)
(*ulaz Konvej*)
(*pomocna za DOWKERA *)

fConSymPoly[Conway_String]:=Module[{Conpom=Conway,pos,ConSym, i},
    If[Union[Flatten[StringPosition[Conpom,"*"]],
          Flatten[StringPosition[Conpom,"."]]]=={},
      (*nije poliedarski*)
      ConSym=ReadList[
          StringToStream[
            StringReplace[
              Conway,{"#"->" ",","->" ","+"->" ","-"->" ",
                "("->" ",")"->" "}]],Number]
      ,
      pos=Flatten[StringPosition[Conpom,"*"]];
      If[pos=={},
          If[StringTake[Conpom,1]==".",Conpom=StringDrop[Conpom,1]]; 
        If[Conpom!="",
                      
          If[StringTake[Conpom,1]==".",
            Conpom=StringInsert[Conpom,"1",1] ]] ;
               pol=6
           (*ako je . na prvom mestu ovo je ".1" koji je x6" *),  
            If[pos[[1]]==2,
                pol=ToExpression[StringTake[Conpom,1]],  
                 pol=ToExpression[StringTake[Conpom,{1,2}]] 
              ];
                (*sada nam je u pol oznaka poliedra tj. 
              broj njegovih presecnih tacaka*)
            
        If[StringLength[Conpom]==1,
                      Conpom="", 
            Conpom=StringTake[Conpom,{pos[[1]]+1,StringLength[Conpom]}]];
                   (*uzimamo ostatak posle * 
            ili prve . (ako ga ima) *)
            If[Conpom!="",
                      
          If[StringTake[Conpom,1]==".",
            Conpom=StringInsert[Conpom,"1",1] ]]
        ];
      
      (*ubacujemo tackice*)
               
      While[StringPosition[Conpom,".."]!={},
        Conpom=StringReplace[Conpom,".."->".1."]];
                    (*ubacujemo 1 u preskocenim temenima *)
               
      If[Conpom=="",
                      
        Do[Conpom=Conpom<>"1.",{i,1,Length[fConvertPoly[Conway][[3]]]}], 
        	     pos=Union[Flatten[StringPosition[Conpom,"."]]];
              Do[Conpom=Conpom<>".1",{i,1,pol-1-Length[pos]}];
                Conpom=Conpom<>".";
                
        If[StringTake[Conpom,1]==".",Conpom=StringDrop[Conpom,1]]];
                ConSym=ReadList[StringToStream[StringReplace[Conpom,
              {","->" ","+"->" ","-"->" ","("->" ",
                ")"->" ","."->" "," 0"->" "}]],Number]];
    ConSym
    ]
(* #############################################################*)
  (*############# DIRECT PRODUCT ###############*)
   (* 27.8.2003*)  
    (*############## DIRECT PRODUCT #################################*)

fRot[{U1_List,{{a_,b_},{a1_,b1_}},Z1_List}]:={U1,{{a,a1},{b,b1}},Z1}
(*NOVA 17.8.2003 *)

fDirectProduct[L_List]:=Module[{pomList,A1,B1,ind,pom,p,p1,d,d1},
    A1=L[[1]];B1=L[[2]]; 
    pomList=Rest[L];
    While[Not[SameQ[pomList,{}]],
                  ind=0;
                 If[Length[A1[[2]]]>2,ind=1];
                 If[Length[B1[[2]]]>2,If[ind==1,ind=3,ind=2]];
      If[ind==2, 
      (*zameni mesta prvom i drugom*)
       pom=A1;A1=B1;B1=pom];
       (*ako je samo 1 poly on je sada u A1*)
       (*ind oznacava da li ima poliedarskih:
       0-nema;1,2 ima jedan;3 oba poly*)
      p=A1[[2,1]];  d=A1[[2,2]];
      p1=B1[[2,1]];   d1=B1[[2,2]];
      If[ind==0, (*ne poly*)
        A1={Join[A1[[1]],B1[[1]],{{p[[2]],p1[[1]]},{d[[2]],d1[[1]]}}],
                                 {{p[[1]],d[[1]]},{p1[[2]],d1[[2]]}},
                                 Join[A1[[3]],B1[[3]]]},
                      (*jedan ili vise poliedarskih*)
                        \
 A1={Join[A1[[1]],B1[[1]]],
                                 Sort[Join[Rest[A1[[2]]],Rest[B1[[2]]],
                {Sort[{p[[1]],p1[[2]]}],Sort[{p[[2]],p1[[1]]}]}]],
                                  Join[A1[[3]],B1[[3]]]}];
                  (*u A1 je trenutni rezultat*)
                   
      pomList=Rest[pomList];
                   If[Not[SameQ[pomList,{}]],B1=First[pomList]]
      ] (*end of While *);
    A1
    ]
(*##########################*)
fConvertDirect[Conway_String]:=Module[{pos,pom,ind,llCon={},i,j=0},
    pos=Flatten[Map[Take[#,1]&,StringPosition[Conway,"#"]]];
    ind=SameQ[
        Union[StringPosition[Conway,"*"],StringPosition[Conway,"."]],{}];
    (* Dowker obezbedjuje da ovo nije prazna lista*)
    Do[If[i==1,
                pom=StringTake[Conway,{1,pos[[1]]-1}],
                 If[i==Length[pos]+1,  
                     
          pom=StringTake[Conway,{Last[pos]+1,StringLength[Conway]}],
                    pom=StringTake[Conway,{pos[[i-1]]+1,pos[[i]]-1}]
                 ]];
      If[Union[StringPosition[pom,"*"],StringPosition[pom,"."]]=={},
        If[ind,pom=fRot[fConvert[pom]],
          pom=fConvert[pom]],
        pom=fConvertPoly[pom]];
      pom={pom[[1]]+j,pom[[2]]+j,pom[[3]]};
         j=j+4*Length[pom[[3]]];
           llCon=Append[llCon,pom]
      ,{i,1,Length[pos]+1}];
    llCon=fDirectProduct[llCon];
    llCon
    ]
(*#################  fConSymDirect######################### *)
(*pravi "pun" Konvejev simbol za direktne proizvode *)
(*pomocna za Dowkera*)

fConSymDirect[Conway_String]:=Module[{pos,llCon={},i},
    pos=Flatten[Map[Take[#,1]&,StringPosition[Conway,"#"]]];
    Do[If[i==1,
                pom=StringTake[Conway,{1,pos[[1]]-1}],
                 If[i==Length[pos]+1,  
                pom=StringTake[Conway,{Last[pos]+1,StringLength[Conway]}],
                    pom=StringTake[Conway,{pos[[i-1]]+1,pos[[i]]-1}]
                 ]];
          pom=fConSymPoly[pom];
          llCon=Append[llCon,pom]
      ,{i,1,Length[pos]+1}];
    Flatten[llCon]
    ]
 (*####################### fGenSignDirProd #################*)
 
 fGenSignDirProd[Conway_String]:=Module[{pos,llCon={},i},
    pos=Flatten[Map[Take[#,1]&,StringPosition[Conway,"#"]]];
    Do[If[i==1,
                pom=StringTake[Conway,{1,pos[[1]]-1}],
                 If[i==Length[pos]+1,  
                   pom=StringTake[Conway,{Last[pos]+1,StringLength[Conway]}],
                    pom=StringTake[Conway,{pos[[i-1]]+1,pos[[i]]-1}]
                 ]];
       pom=fGenSign[pom];
          llCon=Append[llCon,pom]
      ,{i,1,Length[pos]+1}];
    Flatten[llCon]
    ]

(*################ fDToDDirect #####################*)
fDToDDirect[ConDir_String]:=Module[{ddow,l},
    ddow=Dowker[ConDir][[4]];
    l=Map[Length,ddow];
    If[SameQ[Union[l],{0}],l={Length[l]}];
    {l,Flatten[ddow]*fGenSignDirProd[ConDir]}]
 (* #############################################################*)
 
 (* #############################################################*)
 
 
    (*################## DOWKER ######################*)
   (* 16.8.2003*)
(*################# DOWKER ######################### *)

Dowker[Conway_String]:=
  Module[{ConList,Conpom,LinInd=0,L,h,n,ConSym=0,DL={},LDL={},LenList,LL={},
      PomList,lparova={},pol,LenDowList={},Ham,l,pos,k=0,pom=0,ll={},ll1={},
      ind,r,ListforSigns,q, i,j},
    Conpom=Conway;
    If[StringPosition[Conpom,"#"]!={},
            (*ne poledarski DIRECT PRODUCT*)
              
      ConList=fConvertDirect[Conpom];  ConSym=fConSymDirect[Conway]
      (*ConSym nije odgovarajuci za poliedre vidi crvemno dole *)
      ,
      If[Flatten[
            Union[Map[
                Flatten[StringPosition[Conpom,#]]&,{"*","."}]]]=={},
                (*ne poledarski Prime*)
                  
        ConList=fConvert[Conway];ConSym=fConSymPoly[Conway],
                 (*poliedri*)
                   ConList=fConvertPoly[Conway];
        ConSym=fConSymPoly[Conway]
          ]
      ];(*kraj prvog if-a od Direct*)
    (*Print["sym ",ConSym];
      Print["Convert ",ConList];*)
    (*ConSym={3,3};
      ConList={{{3,6},{4,5}},{{1,2},{7,8}},{1,1}};
      ConList=ReplacePart[ConList,{1,1,1,1,1,1,1,1,1},3];*)
    
    If[Head[ConSym]==List,n=Length[ConSym];
      L=fGraphHa[ConSym,Join[ConList[[1]],ConList[[2]]],n];
      DL={ConSym,L[[1]],L[[2]]};
      LenList=Map[Length[#]&,L[[2]]];
      Do[pom=pom+LenList[[i]];LL=Append[LL,pom],{i,1,Length[LenList]}];
      PomList=Flatten[Abs[L[[2]]]];k=0;
      Do[pom=0;
        Do[pom=pom+ConSym[[PomList[[k+j]]]],{j,1,LenList[[i]]}];
        LenDowList=Append[LenDowList,pom];
        k=k+LenList[[i]],{i,1,Length[L[[2]]]}];
      pom=0;
      Do[pom=pom+LenDowList[[i]];
        LDL=Append[LDL,pom],{i,1,Length[LenDowList]}];
      l=Length[LenList];
      If[Head[L[[2]][[1]]]==List,LinInd=1];
      Ham=Flatten[L[[2]]];
      k=0;
      Do[q=ConSym[[Abs[Ham[[i]]]]];
        ll=Append[ll,fFormElem[Ham[[i]],k,q]];
        k=k+q,{i,1,Length[Ham]}];
      pom=LenList[[1]];
      Ham=Abs[Ham];
      Do[ind=False;
        For[j=i+1,j<=Length[Ham]&&Not[ind],j++,
          ind=(Ham[[i]]==Ham[[j]]);
          If[ind,lparova=Append[lparova,{i,j}]];
          r=Position[LL,j-1];
          If[r!={},r=r[[1,1]]+1];
          
          If[LinInd==1&&IntegerQ[r]&&EvenQ[ll[[i,1]]-ll[[j,1]]]&&ind,
            Do[ll[[k]]=Map[(#-1)&,ll[[k]]],{k,LL[[r-1]]+1,LL[[r]]}];
            Part[ll,LL[[r-1]]+1,1]=LDL[[r]]]],{i,1,Length[Ham]-1}];
      Do[i=lparova[[h]][[1]];
        j=lparova[[h]][[2]];
        Do[
          Switch[Mod[ll[[i]][[t]],2],0,
            ll1=Union[ll1,{{ll[[j]][[t]],ll[[i]][[t]]}}],1,
            ll1=Union[ll1,{{ll[[i]][[t]],ll[[j]][[t]]}}]],{t,1,
            Length[ll[[i]]]}],{h,1,Length[lparova]}];
      ListforSigns=ll;
      ll1=Sort[Partition[Flatten[ll1,2],2]];
      ll1=Flatten[Map[Drop[#,1]&,Sort[ll1]]];PomList={};
      If[LinInd==1,LDL=Map[(#/2)&,LDL];
        Do[
          If[i==1,PomList=Append[PomList,Take[ll1,{1,LDL[[1]]}]],
            PomList=Append[PomList,Take[ll1,{LDL[[i-1]]+1,LDL[[i]]}]]],{i,1,
            Length[LDL]}];
        ll1=PomList];
      DL=Append[DL,ll1];
      DL=Append[DL,ListforSigns]];
    DL](*Vraca {} ako je greska*)

(*#########################################################*)
(* ######################################################*)

fGenSign[Conway_String] :=
   Module[{BlockSign={},i,DL,par,fConSignList,Ham,VertSign,Dow={}},
    DL=Dowker[Conway];
    If[DL!={},
    Dow=Flatten[DL[[4]]];
    VertSign=DL[[5]];
    Do[BlockSign = Append[BlockSign, 0], {i, 1, Length[DL[[1]]]}];
    If[Flatten[Union[Map[Flatten[StringPosition[Conway,#]]&,{"*","."}]]]=={},
      fConSignList = fConvert[Conway][[3]],
      fConSignList=fConvertPoly[Conway][[3]]];
    Ham=Flatten[DL[[3]]]; 
    Do[ par=Flatten[Position[Abs[Ham],i]];
        BlockSign[[i]]:=Sign[Ham[[par[[1]]]]Ham[[par[[2]]]]]*
        fConSignList[[i]];
        VertSign[[par[[1]]]]=BlockSign[[i]]*VertSign[[par[[1]]]];
        VertSign[[par[[2]]]]=BlockSign[[i]]*VertSign[[par[[2]] ]]
    ,{i,1,Length[BlockSign]}];VertSign=Flatten[VertSign];
  Dow=Map[(Sign[VertSign[[Flatten[Position[Abs[VertSign],#]][[1]]]]])&,Dow]
  ];
        Dow
    ](*Vraca {} ako je greska*)
 
fGen[Conway_String] := 
  Module[{Dow,OddDow={},proL = {}, ll={}, ulL = {}, izL = {}, pomL, 
  pomLL1,Gen={}, pom=0, i, k, p, n1, j},
    Dow=Dowker[Conway];
If[Dow!={},
    n1 = Length[Dow[[4]]];
    OddDow=Flatten[Dow[[4]]];
    OddDow=Sort[Map[(#-1)&,OddDow]];
    If[Head[Dow[[4]][[1]]] == List, 
     Do[If[i == 1, pom += Length[Dow[[4]][[i]]];
           ll = Append[ll, Take[OddDow,{1, pom}]],
           ll = Append[ll,Take[OddDow, 
           {pom+1,pom+Length[Dow[[4]][[i]]]}]]; 
        pom += Length[Dow[[4]][[i]]]], {i, 1, n1}];
     pomLL1 = {};
     Do[If[i == 1, 
        pomLL1 =Append[pomLL1,Flatten[{Dow[[4]][[n1]][[1]], 
                   Drop[Dow[[4]][[i]], 1]}]], 
        pomLL1 =Append[pomLL1,Flatten[{Dow[[4]][[i-1]][[1]],
                   Drop[Dow[[4]][[i]],1]}]]]
        ,{i,1,n1}];
       pomLL1 = Flatten[pomLL1]
    ];
    pomL = Flatten[Dow[[4]]];
    If[Head[Dow[[4]][[1]]] == Integer, pomLL1 = pomL];
    n = Length[pomL];
    Do[ proL = Append[proL, pomL[[i]]/2], {i,1,n}];
    Do[ For[j = 1, j <= n && (OddDow[[j]] != pomL[[i]] - 1), j++];
        ulL = Append[ulL, pomL[[j]]/2]  ;
        p = Mod[pomL[[i]] + 1, 2n];
        For[k = 1, k <= n && (OddDow[[k]] != p), k++];
        izL = Append[izL, pomLL1[[k]]/2]
      , {i, 1, n}];
    Gen = {};
    Do[  For[k = 1, k <= n && (proL[[k]] != i), k++];
         Gen = Append[Gen, ulL[[k]]];
         Gen = Append[Gen, izL[[k]]];
         Gen = Append[Gen, proL[[k]]]
      , {i, 1, n}]
      
      ];
       
    Gen
    ](*Vraca {} ako je greska*)
 
 Dow[Conway_String]:=Module[{dd,ll={},str="Dowker code: {"},
   dd=Dowker[Conway];
   If[dd!={},dd=dd[[4]];
   If[Head[dd[[1]]]==List,ll=Map[Length[#]&,dd]];
   If[ll=={}, str=str<>"{"<>ToString[Length[dd]]<>"},"
            <>""<>ToString[dd*fGenSign[Conway]]<>"}"
            ,str=str<>ToString[ll]<>","
            <>ToString[Flatten[dd]*fGenSign[Conway]]<>"}"
            ]
            
            ];
  If[dd=={},str="Input data is not correct. You have probably exceeded 
  the number of possible basic polyhedra for the given nuber of crossings."];
            
  str
   
    ]
(* ##########  Gauss Code - Extended Dowker Code ########### *)    
fGaussExt[Ul_] := 
  Module[{Dow, SL, LL, DowPair = {}, DowPair1, DowFl, i, k = 0, DL = {}, 
      pomL = {}},
    If[Head[Ul] == List,
            Dow = Abs[Ul[[2]]]; 
            If[Length[Ul[[1]]] == 1,
               LL = {},
               LL = 2*Ul[[1]]];
           SL = Flatten[Map[Sign[#] &, Ul[[2]]]]];
    If[Head[Ul] == String,
           Dow = Dowker[Ul];
           If[Dow != {},
               Dow = Dow[[4]];
               SL = fGenSign[Ul];
              If[SameQ[Flatten[Dow], Dow],
                  LL = {},
                  LL = Map[2*Length[#] &, Dow]]
        ]];
    If[Dow != {},
       DowFl = Flatten[Dow];
       Table[DowPair = Append[DowPair, {2i - 1, DowFl[[i]]}];
       DowPair = Append[DowPair, {DowFl[[i]], 2i - 1}], {i, Length[DowFl]}];
       DowPair = Union[Map[Sort[#] &, Sort[DowPair]]];
       DowPair1 =Map[Append[#, 
           SL[[Position[DowFl, Select[#, EvenQ][[1]]][[1, 1]]]]*
              Position[DowPair, #][[1, 1]]] &, DowPair];
      DL =Union[Map[Take[#, -2] &, DowPair1], 
          Map[{First[#], Last[#]} &, DowPair1]];
      DL = Flatten[Map[Take[#, -1] &, DL]];
      If[LL != {},
        Do[If[i == 1,
              k = LL[[1]]; pomL = Append[pomL, Take[DL, {1, LL[[1]]}]], 
              pomL = Append[pomL, Take[DL, {k + 1, k + LL[[i]]}]];
              k = k + LL[[i]]],
             {i, 1, Length[LL]}];  
        DL = pomL]];
    If[Flatten[DL] == DL, DL = {DL}];
    DL](*ako je greska vraca {}*)
    
(* #####  Linking No. #### *)   
(* April 3. 2005 *)   

LinkingNo[Ul_] := Module[{GE, LinNoList = {}, d, i, j, k, p,str},
     If[Not[SameQ[Head[Ul], String]],
     If[Length[Ul] != 1,
          If[SameQ[Select[Ul[[2]], OddQ[#1] &], {}],
             GE = fGaussExtSigns[Ul],
             GE = Ul],
         GE = Ul],
     GE = fGaussExtSigns[Ul]];
     d = Length[GE];
     If[SameQ[GE, Flatten[GE]],
        (*ako je u pitanju cvor GE nema podlista *)
       str="Knot",
       Do[Do[
          p = Map[Sign[#] &, Intersection[GE[[i]], GE[[j]]]];
         LinNoList = Append[LinNoList, Sum[p[[k]], {k, 1, Length[p]}]]
         , {j, i+1, d}], {i, 1, d - 1}];
        str="Linking number: "<>
            ToString[ 
        Sum[Abs[LinNoList[[k]]], {k, 1, Length[LinNoList]}]/2
            ]
      ];
      str
    ]
    
(* ########## UnKnot for Rational Knots and Links ############### *)

IsNotKnot[Conway_List]:=
  Module[{Con=Conway,Con1,Con2},
    If[SameQ[Union[Conway],{0}],(*Fractioni ne prolaze-same nule*)
      
      Con={},Con1=FromContinuedFraction[Con];
      Con2=Con;
      If[SameQ[Con1,ComplexInfinity],
      Con=FromContinuedFraction[Reverse[Con2]],
        Con=Con1];
      If[SameQ[Con,ComplexInfinity]&&
          SameQ[FromContinuedFraction[Reverse[Con2]],ComplexInfinity],Con=1,
        Con];
      Con=ContinuedFraction[Con];
      If[Length[Con]>1,
        If[Abs[First[Con]]==0,Con=Drop[Con,{1,2}],
          If[Abs[First[Con]]==1,Con[[2]]=Con[[1]]+Con[[2]];
            Con=Drop[Con,{1}]]]];
      Con=If[SameQ[Con,{}], {1},Con]
      ];
    (* Print["Is is not knot ",Con]; *) Con] (*25.9.2003*)
(*############## RatReduce-samo 1 konveja########*)

RatReduce[Conway_String]:=
  Module[{Con,Con1},
    Con=ReadList[
        StringToStream[
          StringReplace[
            Conway,{","->" ","+"->" ","("->" ",
              ")"->" "}]],Number];
    Con=IsNotKnot[Con];
    Con=Last[Sort[{Con,Reverse[Con]}]];
    Con=If[SameQ[Con,{0}],{1},Con];
    StringReplace[ToString[Con],{","->"","{"->"","}"->""}]
    ]

(*25.9.2003*)
(*#################*)


Redu[Conway_List]:=
    Module[{Con,i=1,j,l,ind=False,Con1,Res={}},
      Switch[Head[Conway[[1]]],Integer,l=1,List,l=Length[Conway]];
      While[i<=l&&Not[ind],
        Switch[Head[Conway[[1]]],Integer,Con=Conway,List,Con=Conway[[i]]];
          Con=IsNotKnot[Con];
          If[Con=={}||Con=={0}||Con=={1}||Con=={-1},
                 ind=True;Res={},  
                j=1;
             While[Not[MemberQ[Res,{0}]]&&j<=Length[Con],
                       
            Res=Union[
                Append[Res,ReplacePart[Con,Con[[j]]-2*Sign[Con[[j]]],j]]];
                (* Print[Res,MemberQ[Res,{}]]; *)
                       j++]
          
          ];
        i++];
      Res];
(* April 3. 2005 - puca zabog Japanaca za 1, 2, -2 *)
UnR[Conway_String]:=
  Module[{Con,br=-1},
    Con=ReadList[
        StringToStream[
          StringReplace[
            Conway,{","->" ","+"->" ","("->" ",
              ")"->" "}]],Number];
    Con=IsNotKnot[Con];
    If[Not[Con=={}||Con=={0}||Con=={1}||Con=={-1}], 
      	While[Con!={},br++; Con=Redu[Con]]
               ,br=0];
               br=If[fSignat[Conway]!=0,
                Range[IntegerPart[(fSignat[Conway]+1)/2],br],
  Range[1,br]];
    br]
    
 (* ######################## Millet ################################ *)  
  
NoToLett[Ken_, Br1_, Br2_, br_] := Module[{p, Ken11 = Ken},
      If[br == 3, p = "a"];
      If[br == 5, p = "b"]; If[br == 7, p = "c"]; If[br == 9, p = "d"];
      Ken11[[Br1]] = ReplacePart[Ken11[[Br1]], p, Br2];
      Ken11
      ];

fKenSym[Ken_, Br1_Integer, Br2_Integer] := 
  Module[{el = Ken[[Br1, Br2]], Ken11 = Ken, pos, pos1, ind1, ind2  },
    pos = Flatten[Position[Ken11[[Br1]], el]]; 
    If[Length[pos] == 1, 
           pos1 = Flatten[Position[Ken11[[el]], Br1]]; 
            Ken11 = NoToLett[Ken11, Br1, Br2 + 1, pos1[[1]]],
            pos1 = Flatten[Position[Ken11[[el]], Br1]];
            ind1 = pos[[2]] - pos[[1]] == 2;
            ind2 = pos1[[2]] - pos1[[1]] == 2;
            If[Or[ind1 && ind2, Not[ind1] && Not[ind2]],
                    Ken11 = NoToLett[Ken11, Br1, Br2 + 1, pos1[[2]]];     
                    Ken11 = NoToLett[Ken11, Br1, pos[[2]] + 1, pos1[[1]]],
                    Ken11 = NoToLett[Ken11, Br1, Br2 + 1, pos1[[1]]]; 
                    Ken11 = NoToLett[Ken11, Br1, pos[[2]] + 1, pos1[[2]]]
                ] ];
                If[Length[pos]==1,pos=Append[pos,-1]];
   { Ken11,pos[[2]]}
    ]

DowToKen[Conway_String] := 
  Module[{DL, Dow, VertSign, GenList, Ul, Iz,Pro, i, j,
   p, m, n,pom,Ind,ll={}, Ken},
    Switch[Conway, "1", Ken={},
                    "-1", Ken={},
                   "2",
Ken={{1,"-",2,d,2,c,2,b,2,a},{2,"-",1,d,1,c,1,b,1,a}},
                   "-2", 
Ken={{1,"+",2,b,2,a,2,d,2,c},{2,"+",1,b,1,a,1,d,1,c}},
                   "1 1",
Ken={{1,"-",2,d,2,c,2,b,2,a},{2,"-",1,d,1,c,1,b,1,a}},
                    "-1 -1", 
Ken={{1,"+",2,b,2,a,2,d,2,c},{2,"+",1,b,1,a,1,d,1,c}},
                   "1 -1", 
Ken={{1,"-",2,c,2,b,2,a,2,d},{2,"+",1,c,1,b,1,a,1,d}}, 
                    "-1 1", 
Ken={{1,"+",2,c,2,b,2,a,2,d},{2,"-",1,c,1,b,1,a,1,d}}];
  Con=StringReplace[Conway,"-"->""];
  If[Conway!="1"&&Conway!="-2"&&Conway!="2"&&Conway!="-1"
    &&Conway!="1 1"&&Conway!="1 -1"&&Conway!="-1 1"&&Conway!="-1 -1",
    DL = Dowker[Con];
    VertSign = fGenSign[Con];
    If[VertSign[[1]]==-1,VertSign=-VertSign];
    GenList = fGen[Con];
    Dow = Flatten[DL[[4]]];
    Do[ll=Append[ll,{Dow[[i]],VertSign[[i]]}], {i,1,Length[Dow]}];
    ll={};
    Do[ll=Append[ll,Dow[[i]]/2->i], {i,1,Length[Dow]}];
    GenList=ReplaceAll[GenList,ll];
    GenList=Flatten[Sort[Map[Reverse[#]&,Partition[GenList,3]]]];
    Pro= Take[GenList, {1, Length[GenList], 3}];
    Ul = Take[GenList, {2, Length[GenList], 3}];
    Iz = Take[GenList, {3, Length[GenList], 3}];
    Ken = Table[0, {i, Length[Dow]}, {j, 10}];
    
    Do[ Ken[[i, 1]] = Pro[[i]];
       Switch[VertSign[[i]]
       , 1, 
          Ken[[i, 2]] = "+", -1, Ken[[i, 2]] = "-"];
          Ken[[i, 3]] =Position[Iz,i][[1]][[1]];(*izlazni za i*)
          Switch[Ken[[i, 2]], "+", Ken[[i, 5]] = Iz[[i]], "-", 
          Ken[[i, 5]] = Ul[[i]]];
         Ken[[i, 7]] = Position[Ul,i][[1]][[1]];
         Switch[Ken[[i, 2]], "+", Ken[[i, 9]] = Ul[[i]], "-", 
                                  Ken[[i, 9]] = Iz[[i]]];
        , {i, 1, Length[Dow]}];   
    Ind=fGenSign[Conway]*
        fGenSign[StringReplace[Conway,"-"->""]]*
                      Flatten[DL[[4]]];
    Ind=Flatten[Map[Position[Ind,#]&,Select[Ind,#1<0&]]];
     Do[    m=Take[Ken[[Ind[[i]]]],2];
           Switch[m[[2]],"+",m[[2]]="-","-",m[[2]]="+"];
           n=Take[Ken[[Ind[[i]]]],{3,10}];
           pom=Ken[[Ind[[i]],5]];
           Ken[[Ind[[i]],5]]=Ken[[Ind[[i]],9]];
           Ken[[Ind[[i]],9]]=pom;
           If[First[fGenSign[StringReplace[Conway,"-"->""]]]>0,
              If[fGenSign[Conway][[Ind[[i]]]]>0,
                 n=RotateLeft[n,2],n=RotateRight[n,2]],
              If[fGenSign[Conway][[Ind[[i]]]]>0, 
				 n=RotateRight[n,2],n=RotateLeft[n,2]]] ;
           Ken[[Ind[[i]]]]=Join[m,n]
       ,{i,1,Length[Ind]}];
  Do[ m=-1;n=-1;
       For[j = 3, j <= 9, j = j + 2,
       If[j!=m&&j!=n,
              Ken=fKenSym[Ken, i, j][[1]];
              If[m==-1,n=fKenSym[Ken, i, j][[2]],n=fKenSym[Ken, i, j][[2]]]]]
      , {i, 1, Length[Dow]}]];
    Ken
    ]

fMillett[Conway_String] := 
  Module[{Ken, str, i, j, str1, f = "izlaz11.txt"},
       Ken = DowToKen[Conway];
       str = " " <> ToString[Length[Ken]] <> ".  00:\n";
       Do[str1 = "";
             Do[   str1 = str1 <> ToString[Ken[[i, j]]] ,
        {j, 1, Length[Ken[[i]]]}];
      str = str <> str1 <> "\n", {i, 1, Length[Ken]}];
    str = str <> "\n";
    OpenWrite[f];
     WriteString[f, str];
       Write[f];
    Close[f];
    Run["lmp1" <> " izlaz11.txt" <> 
        " < " <> "izlaz11.txt" <> " > " <> 
        "polinomi1.txt"];
    Import["polinomi1.txt"] 
    ]

    
(* ############ Connections with knotbycomp #################### *)

fCreatePData[Conway_String] := Module[{Dow, i, l, ind, pd = {}},
    Dow = Dowker[Conway];
 If[Dow!={},Dow=Dow[[4]];
    l = fGenSign[Conway];
    ind = l*fGenSign[StringReplace[Conway, "-" -> ""]];
    pd = Append[pd, If[IntegerQ[First[Dow]], 
          {Length[l]}, Map[Length, Dow]]];
    pd = Append[pd, Table[Sort[Table[
               If[Part[ind, i] == 1,{2i - 1, Flatten[Dow][[i]]*l[[i]]},
                                    {Flatten[Dow][[i]], (2i - 1)*l[[i]]}] 
                            ,{i, 1, Length[l]}]][[i]][[2]]
          , {i, 1, Length[l]}]];
    pd={pd[[1]],-pd[[2]]}
    
    ];
    
    If[Dow=={},pd={}];
    
    pd ]  
    
 fDowfromPD[pdata_List]:= 
  Module[{pd, pd1, pd2, ll = {}, i, j, z1, p = 0, q, z = {}},
    pd = Flatten[pdata[[2]]];
     pd1 = Complement[Table[i, {i, 1, 2Length[pd]}], Abs[pd]];
    pd1 = Table[{pd1[[i]], pd[[i]]}, {i, 1, Length[pd]}];
    pd1 = Join[pd1, Map[Reverse[#] &, pd1]];
    pd1 = Sort[Map[(#*Sign[#[[1]]]) &, pd1]];
    pd2 = Take[Flatten[pd1], {2, Length[Flatten[pd1]], 2}];
    z1 = Select[pd2, (#1 < 0) &];
    z1 = Map[(Abs[#] -> #) &, z1];
    pd1 = Abs[pd1];
    Do[p = p + 2pdata[[1, i]]; 
      ll = Append[ll, p], {i, 1, Length[pdata[[1]]]}];
    Do[z = {};
      If[OddQ[Flatten[pd1][[2ll[[i - 1]] + 2]]],
                 p = ll[[i - 1]] + 1; q = ll[[i]];
        	Table[ If[j == q, z = Append[z, q -> p],
                                                     
            z = Append[z, j -> j + 1]], {j, p, q}];
        pd1 = ReplaceAll[pd1, z]
         ]
      , {i, 2, Length[ll]}];
    pd1 = Sort[Map[If[EvenQ[#[[1]]], Reverse[#], #] &, pd1]];
    pd1 = ReplaceAll[pd1, z1];
    pd1 = Sort[Map[{Abs[#[[1]]], #[[2]]} &, Union[pd1]]];
    pd1 = Flatten[Map[Take[#, -1] &, pd1]];
    {pdata[[1]], pd1}
    ]
    
 (* ######### Unlinking or unknotting number ###############*)       
fGenerate[{LK_List, L_List}] := 
  Module[{pom, p, i, pl = L, l, rez = {}, pDa = {}, br = 1},
    Do[   
      If[MemberQ[Abs[pl], i] == False, 
        pDa = Append[pDa, {i*Sign[pl[[br]]], pl[[br]]}]; br++]
          , {i, 1, 2Length[pl]}];  
    Do[pom = pDa; l = {LK};
          pom[[i]] *= -1;
           pom[[i]] = Reverse[pom[[i]]];
          pom = Sort[pom, Abs[#1[[1]]] < Abs[#2[[1]]] &];
          (*p je oblika originalnog linka*)
            p = Flatten[Map[Take[#, -1] &, pom]];
              l = Append[l, p];
              rez = Append[rez, ReductionKnotLink[l]]
       , {i, 1, Length[pDa]}];
       rez;
     Union[rez]
    ]
    
(* calculates UnKLNo od Konveja, pdata i dowkera *)
(* April 3. 2005- puca za 1, 2, -2 zbog Japanaca*)


UnKnotLink[Ulaz_]:=
  Module[{pD=Ulaz,str="",R,NoGen=1,un0},
  
   If[SameQ[Head[Ulaz],String]&&
    SameQ[Union[StringPosition[Ulaz,","],StringPosition[Ulaz,"("],
        StringPosition[Ulaz,"*"],StringPosition[Ulaz,"."],
        StringPosition[Ulaz,"+"]],{}],NoGen=UnR[Ulaz],   
    If[SameQ[Head[Ulaz],String],pD=fCreatePData[Ulaz]];
    If[SameQ[Head[Ulaz],List]&&SameQ[Select[Ulaz[[2]],OddQ[#]&],{}],
       pD=fPDataFromDow[Ulaz]];
      un0=MemberQ[ReductionKnotLink[pD],{}];
    If[un0==True,
          str="Unknott",
                  R=fGenerate[pD];
     While[MemberQ[R,{{},{}}]==False&&
          MemberQ[R,{{0},{}}]==False,  
        	R=Flatten[Map[fGenerate[#]&,R],1];
        	NoGen++]];
     (*  NoGen=If[fSignat[pD]!=0,Range[IntegerPart[(fSignat[pD]+1)/2],NoGen],
  Range[1,NoGen]] *) ]; 
  NoGen    
    ](*31.8.2004 *)

    
 (* ######### Graphics for KnotPlot ###################*)       

fCreateGraphics[Conway_] := 
  Module[{pdata, k, i, j, pom, str = "", mm,
    filestr = "C:\\Program Files\\KnotPlot\\graphics.txt"},
    pdata = fCreatePData[Conway];
 If[pdata!={},
    Install["DrawKnot"];
    mm = AccountingForm[ N[GetDrawData[pdata[[2]], pdata[[1]]], 6]][[1]];
    Uninstall["DrawKnot"];
    Table[Do[pom = mm[[i, j]];
             Table[str = str <> StringReplace[
          ToString[AccountingForm[pom[[k]]]], {"(" -> "-", ")" -> ""}] 
<>
               " ", {k, 3}];
                        str = str <> "\n"  , {j, 1, Length[mm[[i]]]}];
                   str = str <> "\n" , {i, Length[mm]}];
    filestr = OpenWrite[filestr];
    WriteString[filestr, str];
    Write[filestr];
    Close[filestr],
    Print["Input data is incorrect"]
    ]
    ]
    
 (* #########  KnotLinkBase #########*)
 
 GetKnotLink[str_String,nn_Integer]:=Module[{s,rez},
   s=ToExpression[str];
   If[Head[s]==List,
    If[nn>Length[s],
     s="You have exceeded the number of knots and links from the data base "
     <>str<>".",rez=s[[nn]];
     s="Conway symbol: "<> ToString[rez];
      ]];
  If[Head[s]==Symbol,s="Data base "<>str<> " does not exist.";rez={}];
       Print[s];
  rez   
]

NumberOfKL[str_String]:=Module[{s="There are ",n},
    t=str;
   n=Read[StringToStream[StringDrop[t,1]],Number];
   If[Head[ToExpression[str]]==List,
   s=s<>ToString[Length[ToExpression[str]]];
   If[StringTake[str,1]=="a",s=s<>" alternating",s=s<>" non-alternating"];
   If[n==11,s=s<>" knots with ",s=s<>" knots and links with "];
   s=s<>ToString[n]<>" crossings."];
   If[Head[ToExpression[str]]==Symbol,s="Data base "<>str<> " does not \
exist."];
s]
   
 (* ########################################################*) 
 (*#########################################################*)
 (*############### K2KC1 ###################################*)
 
(*****************************************************
fProjections[string] (R.S. & M.S.)
30.05.2003. 
******************************************************
OPERATIONS:  (in the order of the execution)
-----------
o Operation "X n" --> "(X,1,1,...,1)" where "n" is a number 
  and there are "n" aces; "X" can be anything
  "(4,3) 2"  -->  "((4,3),1,1)"
  This also works with the negative number:
  "(4,3) -2"  -->  "((4,3),-1,-1)"
o Operation "X+n" --> "X,1,1,...,1" where "n" is a number and
  there are "n" aces; "X" can be anything
  "-(2,3)+3" --> "-(2,3),1,1,1"
  This also works with the negative number:  "3+-2"  -->  "3,-1,-1"
o Other combinations with " ", "+", "-" or "," are left unchanged
o Operation "X,1" (or "X,-1") --> makes multiples (different projections?)
  and "X" can be anything (separated by: a. "(", b. "," or c. beggining
  "2,1"      --> "2,1"  and  "1,2"
  "2,-1"     --> "2,-1" and  "-1,2"
  "2 (-2),1" --> "2 (-2),1"  and "1,2 (-2)"
o Operation "1,X" (or "-1,X") --> same thing ( "1,3" --> "1,3" and "3,1")

******************************************************)

(* ================================================ *)
(* = function fProjections[string]                = *)
(* = (previously fConvert2)                       = *)
(* = Output is the list of strings.               = *)
(* = The function is not checking all             = *)
(* = syntax errors that user may make!            = *)
(* = option: Print->On prints intermediate results= *)  
(* ================================================ *)
Options[ fProjections ] = { Print -> Off };
fProjKL[ sInput_String, opts___Rule ] :=
  Module[ 
    {
       sOut,sCh,lComm,iI,iJ,iK,iLevel,iPos,iPrev,iNext,iNumRead,
       sOne,tmp,ltmp,
       lAllCh,lAllComm        (* <--- we keep multiple outputs here *)
    } 
  ,

    (* we keep here splitted string, and later recognized commands *)
    sCh = Characters[ sInput ] ;

    (* ----------------------------------------------- *)
    (* - This part is recognizing the numbers and    - *)
    (* - letters, and we put the description in      - *)
    (* - the array lComm with symbols:               - *)
    (* - 0  --> not recognized (yet)                 - *)
    (* - 1  --> single block                         - *)
    (* - 2  --> open bracket    "("                  - *)
    (* - -2 --> closed bracket  ")"                  - *)
    (* - 3  --> operation "+"                        - *)
    (* - 4  --> number                               - *)
    (* - 5  --> letter "+"                           - *)
    (* - 6  --> letter " "                           - *)
    (* - 7  --> letter ","                           - *)
    (* - 8  --> letter "-"                           - *)
    (* - 9  --> separator after recognized operations- *)
    (* - -4 --> -number (e.g. "-2")                  - *)
    (* ----------------------------------------------- *)
    lComm = {};     (* array for code of symbols *)
    For[ iI=1, iI<=Length[sCh], iI++,  
    (* For loop since we change the length *)
      If[ sCh[[iI]]//DigitQ ,
        (* we have a digit *)
        While[ iI+1<=Length[sCh] && DigitQ[sCh[[iI+1]]] ,
           sCh[[iI]] = sCh[[iI]] <> sCh[[iI+1]];
           sCh       = Drop[ sCh , {iI+1} ];
        ];
        AppendTo[ lComm , 4 ];   (* we have a number *)
      ,
        (* we have something else *)
        AppendTo[ lComm ,
          Switch[ sCh[[iI]],
               "(" , 2  ,          (* open bracket   *)
               ")" , -2 ,          (* closed bracket *)
               "+" , 5  ,          (* letter "+"     *)
               " " , 6  ,          (* letter " "     *)
               "," , 7  ,          (* letter ","     *)
               "-" , 8  ,          (* letter "-"     *)
               _   , 0             (* something else *)
          ]
        ];
      ];
    ];
    (* ---------------------------------------------- *)

    (* ---------------------------------------------- *)
    (* - Check for unrecognized symbols and         - *)
    (* - missing brackets                           - *)
    (* ---------------------------------------------- *)
    If[ Position[ lComm, 0 ] =!= {} , Return[False] ];  (* other symbols? *)
    iLevel = 0;               
    Do[          (* missing brackets? *)
      If[ Abs[lComm[[iI]]]===2 , iLevel += lComm[[iI]] ];
    ,{iI,Length[sCh]}];
    If[ iLevel =!= 0 , Return[ False ] ];   (* error, return False *)
    (* ---------------------------------------------- *)

    (* ---------------------------------------------- *)
    (* - lets find "-", output: nothing             - *)
    (* - it just rearanges the arrays...            - *)
    (* ---------------------------------------------- *)
    For[ iI=Length[sCh]-1, iI>0, iI--,  
    (* For loop since we change the length *)
      If[ lComm[[iI]]===8 ,
         If[ lComm[[iI+1]] === 4 ,
            lComm[[iI]] = -4;   (* single block, i.e. -number *)
         ,
            lComm[[iI]] = 2;    (* bracket *)
         ];
         sCh[[iI]]   = "-" <> sCh[[iI+1]] ;
         sCh   = Drop[ sCh   , {iI+1} ];
         lComm = Drop[ lComm , {iI+1} ];
      ];
    ];
    (* ---------------------------------------------- *)
    If[ (Print /. {opts})===On,
       Print["blocks: ",Thread[{sCh,lComm}]//MatrixForm];
    ];
(* Other way of printing:
Print[ Thread[{Range[sCh//Length],sCh,lComm}]//Transpose//MatrixForm ];
*)

    (* ---------------------------------------------- *)
    (* - operation " n"  ("n" is number!)           - *)
    (* ---------------------------------------------- *)
    For[ iI=2, iI<Length[sCh], iI++,       (* For loop since we change the \
length *)
      If[ lComm[[iI]]===6 && Abs[ lComm[[iI+1]] ]===4 ,     (* " num" *)
        (* what is the prevoius element? *)
        iPrev = iI-1;               (* this assumes single block *)
        If[ lComm[[iI-1]]===-2 ,    (* or it is a bracket        *)
            iLevel = 0;
            iJ     = iI-2;
            While[ iLevel<2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ--;
            ];
            iPrev = iJ+1 ;
        ];
        sCh   =  Insert[ sCh  , "("       , iPrev ];
        lComm =  Insert[ lComm, 2         , iPrev ];
        iNext = iI+1;                                
        (* we inserted one element  *)
        (* now we search for other elements *)
        If[ iNext<=Length[sCh] && lComm[[iNext]]===6,     
        (* bilo While, sad je If *)
          sCh[[iNext]]   = ",";
          lComm[[iNext]] = 7  ; (* this is ordinary comma now *)
          iNext++;              (* next element          *)
          If[ lComm[[iNext]]===2 ,    (* is it a bracket? *)
            iLevel = 0;
            iJ     = iNext+1;
            While[ iLevel>-2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ++;
            ];
            iNext = iJ; (* go to the next element, if any *)
	  ,  (* no, it is a number *)
            iNumRead = ToExpression[ sCh[[iNext]] ];
	    If[ iNumRead < 0 ,
	       iNumRead = -iNumRead;
               sOne     = "-1";
            ,
               sOne     = "1";
            ];
	    sCh[[ iNext ]] = sOne; 
	    iNext++;
	    Do[
              sCh   =  Insert[ sCh  , "," , iNext ];
              lComm =  Insert[ lComm, 7   , iNext ]; (* comma *)
	      iNext++;
              sCh   =  Insert[ sCh  , sOne , iNext ];
              lComm =  Insert[ lComm, 4*ToExpression[sOne] , iNext ]; 
              (* number    *)
	      iNext++;
            ,{iK,iNumRead-1}];
          ];
        ];
        sCh   =  Insert[ sCh  , ")" , iNext    ];
        lComm =  Insert[ lComm, -2   , iNext    ]; 
        (* we consider this as a bracket *)
      ];
    ];
    (* ---------------------------------------------- *)
    If[ (Print /. {opts})===On,
      Print["operation ' ': ",Thread[{sCh,lComm}]//MatrixForm];
    ];

    (* ---------------------------------------------- *)
    (* - operation "+n"  ("n" is number!)           - *)
    (* ---------------------------------------------- *)
    For[ iI=2, iI<Length[sCh], iI++,       (* For loop since we change the \
length *)
      If[ lComm[[iI]]===5 && Abs[ lComm[[iI+1]] ]===4 ,     (* "+num" *)
        (* what is the prevoius element? *)
        iPrev = iI-1;               (* this assumes single block *)
        If[ lComm[[iI-1]]===-2 ,    (* or it is a bracket        *)
            iLevel = 0;
            iJ     = iI-2;
            While[ iLevel<2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ--;
            ];
            iPrev = iJ+1 ;
        ];
        iNext = iI;       
        (* now we search for other elements *)
        If[ iNext<=Length[sCh] && lComm[[iNext]]===5,     
        (* bilo While, sad je If *)
          sCh[[iNext]]   = ",";
          lComm[[iNext]] = 7  ; (* this is comma now *)
          iNext++;              (* next element          *)
          If[ lComm[[iNext]]===2 ,    (* is it a bracket? *)
            iLevel = 0;
            iJ     = iNext+1;
            While[ iLevel>-2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ++;
            ];
            iNext = iJ;         (* go to the next element, if any *)
	  ,  (* no, it is a number *)
            iNumRead = ToExpression[ sCh[[iNext]] ];
	    If[ iNumRead < 0 ,
	       iNumRead = -iNumRead;
               sOne     = "-1";
            ,
               sOne     = "1";
            ];
	    sCh[[ iNext ]] = sOne; 
	    iNext++;
	    Do[
              sCh   =  Insert[ sCh  , "," , iNext ];
              lComm =  Insert[ lComm, 7   , iNext ]; (* comma *)
	      iNext++;
              sCh   =  Insert[ sCh  , sOne , iNext ];
              lComm =  Insert[ lComm, 4*ToExpression[sOne] , iNext ]; (* \
number    *)
	      iNext++;
            ,{iK,iNumRead-1}];
          ];
        ];
      ];
    ];
    (* ---------------------------------------------- *)
    If[ (Print /. {opts})===On,
      Print["operation '+': ",Thread[{sCh,lComm}]//MatrixForm];
    ];
    
    lAllCh   = { sCh   } ;   (* mupltiple outputs *)
    lAllComm = { lComm } ;
    (* ---------------------------------------------- *)
    (* - operation ",1"    --> find multiples      -- *)
    (* ---------------------------------------------- *)
      For[ iK=1, iK<=Length[lAllCh], iK++ ,      
      (* loop over all multiples *)
      For[ iI=2, iI<Length[lAllCh[[iK]]] , iI++ ,    (* characters in one \
realization *)
        sCh   = lAllCh[[iK]];
	lComm = lAllComm[[iK]];      (* new realization to be generated *)
        If[ lComm[[iI]]===7 && (sCh[[iI+1]]==="1" || sCh[[iI+1]]==="-1"), 
        (* ",1" *)

          (* what is the prevoius element: a. begginig, b. "," or c. "(" *)
	  For[ iPrev=iI-1, iPrev>0 && lComm[[iPrev]]=!=7 && lComm[[iPrev]]=!=2, \
iPrev--,
             If[ lComm[[iPrev]]===-2 ,    (* or it is a bracket   ")"  *)
               iLevel = 0;
               While[ iLevel<2 ,
                 iPrev--;
                 If[ Abs[lComm[[iPrev]]]===2 , iLevel += lComm[[iPrev]] ];
               ];
             ];
          ];
          iPrev++;
      
          If[ !(iPrev===iI-1 && sCh[[iPrev]]===sCh[[iI+1]]) ,  
          (* no sence to swap "1,1" *)
            (* and swap them *)
	    tmp  = sCh[[iI+1]];
	    ltmp = lComm[[iI+1]];
	    Do[
                sCh[[ iJ ]]   = sCh[[ iJ-2 ]];
                lComm[[ iJ ]] = lComm[[ iJ-2 ]];
            ,{iJ,iI+1,iPrev+2,-1}];
            sCh[[iPrev]]     = tmp;
	    lComm[[iPrev]]   = ltmp;
            sCh[[iPrev+1]]   = ",";
	    lComm[[iPrev+1]] = 7;
 
	    (* do we have it alreadey in the list? *)
	    If[ !(MemberQ[lAllCh,sCh] && MemberQ[lAllComm,lComm]) ,	
	       AppendTo[ lAllCh, sCh ];
	       AppendTo[ lAllComm, lComm ];
            ];	  
          ];
        ];
      ];
    ];
    (* ---------------------------------------------- *)
    (* - operation "1,"    --> find multiples      -- *)
    (* ---------------------------------------------- *)
    For[ iK=1, iK<=Length[lAllCh], iK++ ,    (* loop over all multiples *)
      For[ iI=2, iI<Length[lAllCh[[iK]]], iI++ , (* characters in one \
realization *)
        sCh   = lAllCh[[iK]];
	lComm = lAllComm[[iK]];      (* new realization to be generated *)
        If[ lComm[[iI]]===7 && (sCh[[iI-1]]==="1" || 
        sCh[[iI-1]]==="-1"),     
 (* "1," *)

          (* what is the next element: a. end, b. "," or c. ")" *)
	  For[ iNext=iI+1, 
                 iNext<=Length[sCh] && lComm[[iNext]]=!=7 && \
lComm[[iNext]]=!=-2
          , iNext++,
             If[ lComm[[iNext]]===2 ,    (* or it is a bracket   ")"  *)
               iLevel = 0;
               While[ iLevel>-2 ,
                 iNext++;
                 If[ Abs[lComm[[iNext]]]===2 , iLevel += lComm[[iNext]] ];
               ];
             ];
          ];
          iNext--;
      
          If[ !(iNext===iI+1 && sCh[[iNext]]===sCh[[iI-1]]) ,  
          (* no sence to swap "1,1" *)
            (* and swap them *)
	    tmp  = sCh[[iI-1]];
	    ltmp = lComm[[iI-1]];
	    Do[
                sCh[[ iJ ]]   = sCh[[ iJ+2 ]];
                lComm[[ iJ ]] = lComm[[ iJ+2 ]];
            ,{iJ,iI-1,iNext-2}];
            sCh[[iNext]]     = tmp;
	    lComm[[iNext]]   = ltmp;
            sCh[[iNext-1]]   = ",";
	    lComm[[iNext-1]] = 7;

	    (* do we have it alreadey in the list? *)
	    If[ !(MemberQ[lAllCh,sCh] && MemberQ[lAllComm,lComm]) ,	
	       AppendTo[ lAllCh, sCh ];
	       AppendTo[ lAllComm, lComm ];
            ];	  
          ];
        ];
      ];
    ];
    (* ---------------------------------------------- *)

    sOut =  Map[ StringJoin , lAllCh ];   (* Join the strings! *)
    If[ (Print /. {opts})===On,
       Print["Input: \"" <> sInput <> "\"" ];
       Do[
             Print[ "  " <> ToString[iK] <> " --> \"" 
             <> sOut[[iK]] <> "\"" \

];
       ,{iK,lAllCh//Length}];
    ];
    If[ Head[ tmp=(Print /. {opts}) ] === OutputStream ,
       WriteString[tmp,"Input: \"" <> sInput <> "\"\n" ];
       Do[
         WriteString[tmp,"  " <> ToString[iK] <> " --> \"" 
         <> sOut[[iK]] <> "\
\"\n" ];
       ,{iK,lAllCh//Length}];
    ];

    sOut    (* Output *)
  ];
(* ================================================ *)


(**************************************************)
(*############# kraj fPROJECTIONS ################*)


(*############# fPROJECTIONS ADD ################*)

fPartL[LL_List, LL1_List] := Module[{ss = LL, pp = LL1, i},
    For[i = 1, i < Length[pp] + 1, ss = Part[ss, pp[[i]]]; i++];
    ss]

fPartList[LL_List, LL1_List, i_Integer] := Module[{ss = LL, pp = LL1},
    pp = Join[Drop[pp, -1], {i}];
    pp = fPartL[ss, pp];
    pp
    ]

fBlockOne[LL_List] := 
  Module[{dd2 = LL, mm, mm1, mm2, i}, 
    mm = Table[If[SameQ[dd2[[i]], 1], 1, 0], {i, Length[dd2]}];
    mm = Table[If[SameQ[dd2[[i]], 1], 1, 0], {i, Length[dd2]}];
    mm1 = Table[If[SameQ[Head[dd2[[i]]], List], 2, 0], {i, Length[dd2]}];
    mm2 = 
      Table[If[SameQ[NumberQ[dd2[[i]]], True] && dd2[[i]] >= 2, 3, 0], {i, 
          Length[dd2]}];
    mm = mm1 + mm2;
    mm = Split[mm];
    mm = Flatten[
        Table[If[MemberQ[mm[[i]], 2] && Length[mm[[i]]] > 1, 
            Partition[mm[[i]], 1], {mm[[i]]}], {i, Length[mm]}], 1];
    mm1 = 
      Flatten[Position[Table[If[MemberQ[mm[[i]], 2], 1, 0], 
      {i, Length[mm]}], 1]];
    mm2 = mm1 - 1;
    mm = Table[
        If[MemberQ[mm[[mm2[[i]]]], 0], Length[mm[[mm2[[i]]]]], 0], {i, 
          Length[mm2]}];
    mm
    ]


fChangePart[Ulaz_String] := Module[{ss, dep, dd, ff, dd1, dd2, dd3, k, i},
    ss = ToExpression[StringReplace[Ulaz, {"(" -> "{", ")" -> "}"}]];
    dep = Depth[ss];
    Do[
      dd = Depth[ss];
ff = dd;
      dd = Flatten[ss, dd - k];
      dd = Select[dd, SameQ[Head[#], List] &];
      dd1 = Flatten[Table[Position[ss, dd[[i]]], {i, Length[dd]}], 1];
       dd1 = Select[dd1, Length[#] == ff - k + 1 &];
      dd1 = Union[dd1]; 
      dd3 = fPartL[ss, Drop[dd1[[1]], -1]];
      (* ovo je deo na kome radimo *)
      ll = fBlockOne[dd3];
      dd2 = Table[fPartL[ss, dd1[[i]]], {i, Length[dd1]}];
      dd2 = 
        Table[If[OddQ[ll[[i]]], Reverse[dd2[[i]]], dd2[[i]]], {i, 
            Length[dd2]}];
      Do[ss = ReplacePart[ss, dd2[[i]], dd1[[i]]], {i, Length[dd1]}];
      ss, {k, 3, dep}];
    ss = ToString[ss];
    ss = StringReplace[ss, {"{" -> "(", "}" -> ")", " " -> ""}];
    ss
    ]
    
    fFixP[Ul_String] := Module[{pr = Ul, pp, i}, 
    pp = StringPosition[pr, " "];
    pr = StringJoin["{", 
        StringReplace[pr, {"(" -> "{", ")" -> "}", " " -> ","}], "}"];
    pr = ToExpression[pr];
    pr = Table[
        StringReplace[
          ToString[pr[[i]]], {" " -> "", "{" -> "(", "}" -> ")"}], {i, 
          Length[pr]}];
    pr = Map[fChangePart, pr];
    pr = Table[
        If[StringLength[pr[[i]]] == 3, StringDrop[StringDrop[pr[[i]], 1], \
-1],
           pr[[i]]], {i, Length[pr]}];
    pr = StringDrop[
        StringJoin[Table[StringJoin[pr[[i]], ","], {i, Length[pr]}]], -1];
    pr = StringReplacePart[pr, " ", pp];
    pr]
    

   (*############### fPPoly #######################*)

(*computes projections of polyhedral KL *)
(* 15. 05. 2005 *)
(*".2 2..2..2"="6*.2 2..2..2"and 
"2 2..2..2"="6*2 2..2..2"*)
(*computes projections of polyhedral KL *)

fPPoly[Conway_String]:=
  Module[{Con=Conway,pos1,pos2,res={},str,pom},
  If[SameQ[StringPosition[Con,"*"],{}]&&
  Not[SameQ[StringPosition[Con,"."],{}]]
&&Not[SameQ[StringPosition[Con,"."][[1]],{1,1}]],Con="6*"<>Con,Con];
    pos1=Flatten[StringPosition[Con,"*"]];
     If[pos1=={},Con=Con<>".",
      If[pos1[[1]]!=StringLength[Con],Con=Con<>"."]];
         pos2=Union[Flatten[StringPosition[Con,"."]]];
     If[pos1!={},str=StringTake[Con,pos1[[1]]];
      Con=StringDrop[Con,pos1[[1]]],str=StringTake[Con,pos2[[1]]];
      Con=StringDrop[Con,pos2[[1]]];
      pos2=Rest[pos2](*izbacili 1. tackicu*)];
        (*izdvojili smo pocetak i sad idemo od tackice do tackice tacnije 1. 
    Prvo je tackica:prepisemo je;
      2. prvo je broj:uzimamo deo do sledece tackice i radimo fProjections*)
      res={{str}};
    If[Con==".",Con=""];
    While[Con!="",
      If[StringTake[Con,1]==".",(*prva tackica*)Con=StringDrop[Con,1];
        res=Map[(#<>".")&,res],(*pom je parce do sledece tackice*)pom=
          StringTake[Con,StringPosition[Con,"."][[1,1]]-1];
        (*brisemo iz stringa parce*)Con=
          StringDrop[Con,StringPosition[Con,"."][[1,1]]];
        If[StringLength[pom]>2,
          If[SameQ[StringTake[pom,-2]," 0"],(*ako je nula na kraju*)
          pom=StringDrop[pom,-2];
            pom=fProjKL[pom];
            
            pom=Map[#<>" 0"&,pom],(*ako nema _0 i duzi od 3*)pom=
              fProjKL[pom]],(*ako nema _0 i kraci od 3*)pom=
            fProjKL[pom]];(*sve projekcije datog parceta*)
           pom = Flatten[Map[fFixP, pom]];
          Do[
          If[Con=="",str=Map[(res[[i]]<>#)&,pom],
            str=Map[(res[[i]]<>#<>".")&,pom]];
          res=ReplacePart[res,str,i],{i,1,Length[res]}];
        res=Flatten[res]];(*kraj ifa*)If[Con==".",Con=""]];
    res]
    
    
  
    (*############### fPROJECTIONS  ##################*)
    
    fProjections[Con_String] := 
  Module[{res, s, s1, pr},    
   If[Not[SameQ[StringTake[Con,-1],"*"]],   
    If[Union[Flatten[StringPosition[Con, "*"]], 
          Flatten[StringPosition[Con, "."]]] == {}, res = fProjKL[Con], 
      res = fPPoly[Con]];
      If[Union[Flatten[StringPosition[res[[1]], "*"]], 
          Flatten[StringPosition[res[[1]], "."]]] == {},
    If[Length[res] == 1,
      pr = res,
      pr = Flatten[Map[fFixP, res]]];
    pr = Flatten[pr];
    pr, pr=res], pr={Con}];
    pr]
    
    
(* fProjections[Con_String] := 
  Module[{res, s, s1, pr}, 
    If[Union[Flatten[StringPosition[Con, "*"]], 
          Flatten[StringPosition[Con, "."]]] == {}, res = fProjKL[Con], 
      res = fPPoly[Con]];
      Print[res];
    If[Length[res] == 1,
      pr = res,
      pr = Flatten[Map[fFixP, res]]];
    pr = Flatten[pr];
    pr]*)
     

(*################ F DIFFERENT PROJECTIONS ############*)
fDiffProjectionsAltKL[Con_String]:=Module[{ind,l,i,j,res={},pom},
    l=fProjections[Con];
    If[Con=="2",
        res={{"2",{{1, 1},{4, 2}}}},
         If[Con=="-2",
        res={"2",{{1, 1},{4, 2}}}, 
    If[Length[l]==1&&Head[l[[1]]]==List,l=l[[1]]];
    (*red samo za 6* i sl*)
    Do[ pom=Abs[MinDowProjAltKL[l[[i]]]];
              ind=False;j=1;
           While[  ind==False&&j<=Length[res],
                           ind=SameQ[pom,res[[j,2]]];j++];
       If[ind==False,res=Append[res,{l[[i]],pom}]]
      ,{i,1,Length[l]}]
    ]];res]
    
(*#############################################*)


(*############# fTorus KL ##########################*)

(*fTorusKL calculates for a torus KL[a,b] braid word,min.no.of
 crossings, unknotting no.,no.of components,brodge no.,
 Alexander polynomial and Murasugi signature*)
(**********)
fTorusKL[a_Integer,b_Integer]:=
  Module[{qq,rr,dd,www,bw,pd,pdr,cr,mcr,br,bwr,Alex},
    qq=Max[{a,b}];
    rr=Min[{a,b}];
    dd=GCD[qq,rr];
    www=FromCharacterCode[Range[rr-1]+96];
    bw=StringJoin[Table[www,{i,qq}]];
    Print["Braid word: ",bw];
    pd=KnotFromBraid[bw];
    pd=ReductionKnotLink[pd];
    Print["Pdata: ",pd];
    Print["Braid word: ",bw];
    If[dd>1,
        Print["Crossings: ",Length[pd[[2]]]];
        Print["Number of components: ",fComponentNo[pd]];
        (*Print["Unlinking number: ",UnKnotLink[pd]];*)
        Print[SkeinPolynomial[-1,pd]];
      	Print["Signature: ",fSignat[pd]],
      cr=Max[qq(rr-1),rr(qq-1)];
      mcr=Min[qq(rr-1),rr(qq-1)];
      Print["Crossings: ",mcr];
       Print["Number of components: ",fComponentNo[pd]];
      If[SameQ[dd,1],Print["Unknotting number: ",(qq-1)(rr-1)/2],
        Print[dd," -component link"]];
      br=Min[qq,rr];
      Print["Bridge number: ",br];
      Alex=PolynomialQuotient[(1-x)((1-x^(qq*rr)/dd)^dd),(1-x^qq)(1-x^rr),x];
      Print["Alexander polynomial: ",Alex];
      Print["Signature: ",fSignat[pd]];
      Print["Murasugi signature: ",
        If[SameQ[a,b],If[SameQ[Mod[a,2],1],(a^2-1)/2,(a^2-2)/2],
          fSigTor[a,b]]]
      ];
    ShowKnotfromPdata[pd];
    pd]
(*############### fAlexPoly ####################*)

(*calculates multivariable Alexander polynomial *)

(*############### fALEX ####################*)(*calculates multivariable \
Alexander polynomial*)
(*April 3, 2005 *)
fAlexPoly[Conway_String]:=
  Module[{lL0,res,lL,t,i,s,sl,j,p,MAT},
    If[Conway=="1",res=0,
      If[Conway=="2",res=1,lL0=fGenerators[Conway][[1]];
        lL=Flatten[lL0];
        t=ZeroMatrix[1,Length[lL0]][[1]];
        Do[
          If[i==1,t[[i]]=Length[lL0[[i]]]/3,
            t[[i]]=t[[i-1]]+Length[lL0[[i]]]/3],{i,1,Length[lL0]}];
        s=fGenerators[Conway][[2]];
        sl=Length[s];(*duzina*)MAT=ZeroMatrix[sl];
        j=1;p="a";
        For[i=1,i<=sl,i++,If[i>t[[j]],j++;
            p=FromCharacterCode[ToCharacterCode["a"]+j-1]];
            Switch[s[[i]],1,MAT[[i,Part[lL,3(i-1)+1]]]=1,-1,
            MAT[[i,Part[lL,3(i-1)+1]]]=-p];
          
          Switch[s[[i]],1,MAT[[i,Part[lL,3(i-1)+2]]]=-p,-1,
            MAT[[i,Part[lL,3(i-1)+2]]]=1];
          MAT[[i,Part[lL,3i]]]=p-1;];
          MAT=
          Apply[PolynomialGCD,
            Table[First[Part[Minors[MAT],i]],{i,1,Length[MAT]}]];
       If[SameQ[MAT,0],res=0,
        MAT=If[SameQ[MAT[[1]],-1],-MAT,MAT];
        res=Expand[Divide[MAT,PolynomialGCD[MAT[[1]],MAT]]]]] ];
        res]
(*########################################################*)
 
 (*## GAUSS CODE ####*)
(*##################################################*)
(*#### fGaussExtSigns ##############*)

(*# Input List-
      Dowker Code {{Lengths of Components},{code}} 
      ConwaySymbol_ String 
      PData
Output Extended Dowker=
    Gauss Code with signs #*)
    (*Neded for MinDowPr for Knots*)
    
  fGaussExtSigns[Ulaz_]:=
  Module[{DowL,SL,LL,DowPair={},DowPair1,DowFl,i,k=0,DL,
      pomL={},Ul=Ulaz},
       (*Input:Pdata*)
      If[SameQ[Head[Ul],List],
      If[MemberQ[Map[OddQ,Abs[Ul[[2]] ] ] ,True],
      Ul=fDowfromPD[Ul]]];
        (*sad nam je Ul iliKonvej ili dowker *)
      (*Input:Dowker*)
      If[SameQ[Head[Ul],List],DowL=Abs[Ul[[2]]];
      If[Length[Ul[[1]]]==1,LL={},LL=2*Ul[[1]]];
      SL=Flatten[Map[Sign[#]&,Ul[[2]]]]];
    (*Input:Conway*)
    If[Head[Ul]==String,DowL=Dowker[Ul];
      If[DowL!={},DowL=DowL[[4]];
        SL=fGenSign[Ul];
        If[SameQ[Flatten[DowL],DowL],LL={},LL=Map[2*Length[#]&,DowL]]]];
    (*Print[DowL,SL,LL];*)(*Form ext*)DowFl=Flatten[DowL];
    Table[DowPair=Append[DowPair,{2i-1,DowFl[[i]]}];
      DowPair=Append[DowPair,{DowFl[[i]],2i-1}],{i,Length[DowFl]}];
    DowPair=Union[Map[Sort[#]&,Sort[DowPair]]];
    (*Print["Ovde je",DowPair];*)DowPair1=
      Map[Append[#,
            SL[[Position[DowFl,Select[#,EvenQ][[1]]][[1,1]]]]*
              Position[DowPair,#][[1,1]]]&,DowPair];
    (*Print["Ili ovde ",DowPair1];*)DL=
      Union[Map[Take[#,-2]&,DowPair1],Map[{First[#],Last[#]}&,DowPair1]];
    DL=Flatten[
        Map[Take[#,-1]&,DL]];(*Print["Prvi prosireni ",
          DL];*)(*ako je link moramo da ga podelimo*)If[LL!={},
      Do[If[i==1,k=LL[[1]];pomL=Append[pomL,Take[DL,{1,LL[[1]]}]],
          pomL=Append[pomL,Take[DL,{k+1,k+LL[[i]]}]];
          k=k+LL[[i]]],{i,1,Length[LL]}];(*Print["DL ",DL];*)DL=
        pomL];(*Print["DP ",DowPair];*)DL]

(*#### fGaussExt ##############*)
(*# Input List-
      Dowker Code {{Lengths of Components},{code}} ConwaySymbol_ String \
Output Extended Dowker=Gauss Code without signs #*)
fGaussExt[Ul_]:=Module[{r},r=Abs[fGaussExtSigns[Ul]]]

(*############ Dowker from ExtGauss #################*)
fDowfromGaussExt[Ul_List]:=Module[{LL,ExtD={Ul},pom,k=0,i,Dow={}},
    If[SameQ[Length[Ul],Length[Flatten[Ul]]],
      LL={};
      ExtD=ExtD[[1]],
      ExtD=Flatten[ExtD,1];
      LL=Map[(Length[#]/2)&,ExtD];
      ExtD=Flatten[fFixLink[ExtD]]];
     pom=Union[Map[Append[Flatten[Position[ExtD,#]],Sign[#]]&,ExtD]];
    pom=Sort[Map[If[EvenQ[#[[1]]],{#[[2]],#[[1]],#[[3]]},#]&,pom]];
    pom=Map[(#[[2]]*#[[3]])&,pom];
    If[LL=={},LL={Length[pom]};Dow=pom,
      Do[If[i==1,Dow=Append[Dow,Take[pom,{1,LL[[1]]}]],
          Dow=Append[Dow,Take[pom,k+{1,LL[[i]]}]]];
        k=k+LL[[i]],{i,1,Length[LL]}]];
    {LL,Dow} ]
(*##################################################*)
(*##################################################*)

(*## MINIMIZATION ####*)
(*##################################################*)
(*# fMinComp ###########*)(*# Minimizes ExtDow of a knot #*)
(*# pomocna za fMinDowPr #*)
(*# Could be used for minimizing link components #*)

fMinComp[ExtDowComp_List]:=Module[{d1,d2,m,p,pL,i,pom, t},
    d1=Map[Flatten[Position[ExtDowComp,#]]&,ExtDowComp];
    d1=Union[
        Map[If[Length[#]==1,{#[[1]],
                Length[ExtDowComp]},{#[[1]],#[[2]]-#[[1]]}]&,d1],
        Map[If[Length[#]==1,{#[[1]],Length[ExtDowComp]},{#[[1]],
                Length[d1]-#[[1]]+#[[2]]}]&,Map[Reverse[#]&,d1]]];
    d1=Flatten[Map[Take[#,-1]&,d1]];
    d2=ReplaceAll[Length[d1]-d1,0->Length[d1]];
    m=Min[Join[d1,d2]];
    (*m je minimal distance in both lists*)
    p={Flatten[Position[d1,m]],Flatten[Position[d2,m]]};
   pL=
      Join[Map[{Take[RotateLeft[d1,#-1],{1,Length[d1],2}],#,0}&,p[[1]]],
        Map[{Take[Reverse[RotateRight[d2,Length[d2]-#]],{1,Length[d1],2}],#,
              1}&,p[[2]]]];
    (*Element pl je oblika:odgovarajuca lista razliak, 
      pos od koje smo krenuli i 0 ili 1 za smer*)
      pL=Sort[pL];
    pom=Flatten[Union[Map[Take[#,1]&,pL]],1];
    (*Biramo one koje imaju razlicite razlike bez 
    obzira na smer i poziciju*) 
    m={};
    Do[m=Append[m,Select[pL,SameQ[#[[1]],pom[[i]]]&][[1]]],{i,1,
        Length[pom]}];
    m;
    (* naredni deo dopisan -uzima samo minimalne *)
    t=Select[
        Sort[Table[{Take[Sort[m][[i]][[1]],2],Sort[m][[i,2]],
              Sort[m][[i,3]]},{i,Length[Sort[m]]}]],#[[1]]==
            Sort[Table[{Take[Sort[m][[i]][[1]],2],Sort[m][[i,2]],
                    Sort[m][[i,3]]},{i,Length[Sort[m]]}]][[1,1]]&];
       m=Table[
       If[SameQ[Take[m[[i]][[1]],2],t[[i]][[1]]],m[[i]],{}],{i,Length[t]}];
      m
    ]
    
(*#### fKnittDow ##############*)
(*# pomocna za fMinDowPr #*)
fKnittDow[ExtDow_List,UlList_List]:=Module[{pom,Ps},
    If[Head[UlList[[1]]]==Integer,
      If[UlList[[2]]==0,pom=RotateLeft[ExtDow,UlList[[1]]-1],
        pom=Reverse[RotateRight[ExtDow,Length[ExtDow]-UlList[[1]]]]];
      Ps=Map[Sign[#]&,Sort[Union[pom],Abs[#1]<Abs[#2]&]]];
    {pom,Ps}]

(*#### fFormDow ZAMENJENA ##############*)

fFormDow[KKK_List,Sig_List]:=Module[{LL=KKK,q, i,t,t1},
    LL=Union[Map[Flatten[Position[Abs[KKK],Abs[#]]]&,KKK]];
    t=Table[LL[[i]][[1]],{i,Length[KKK]/2}];
    t1=Table[
        Sign[KKK[[t[[i]]]]],{i,Length[Table[LL[[i]][[1]],{i,Length[LL]}]]}];
    q=Table[{LL[[i]],Sign[KKK[[t[[i]]]]]},{i,Length[t]}];
    q=Table[
        Sort[Table[
                If[EvenQ[q[[i,1,1]]],{q[[i,1,2]],
                    q[[i,1,1]]*q[[i]][[2]]},{q[[i,1,1]],
                    q[[i,1,2]]*q[[i]][[2]]}],{i,Length[LL]}]][[i]][[2]],{i,
          Length[LL]}];
    q
    ]
(*###################### ? ###################*)

(*## INCIDENCY GRAPH ####*)
(*##################################################*)
(*##########################*)
(*##### fFormGr#######*)
(*#radi za cvorove ili komponente linkova #*)
fFormGr[Ul_List]:=Module[{i,rez={}, j},
    If[SameQ[Length[Ul],Length[Flatten[Ul]]],Do[If[i==Length[Ul],
          rez=Append[rez,{Last[Ul],First[Ul]}],
          rez=Append[rez,{Ul[[i]],Ul[[i+1]]}];],{i,1,Length[Ul]}], 
      rez=Flatten[
          Table[Append[
              Table[{Ul[[j]][[i]],Ul[[j]][[i+1]]},{i,
                  Length[Ul[[j]]]-1}],{Last[Ul[[j]]],First[Ul[[j]]]}],{j,
              Length[Ul]}],1]];
    rez=Map[If[Abs[#[[1]]]>Abs[#[[2]]],{#[[2]],#[[1]]},#]&,rez]; 
    rez]
(*#####fMakeS #######*)
fMakeS[LL_List,SL_List]:=Module[{r},r=Map[(#*SL[[#]])&,LL];
    r]
(*##### fSignSum#######*)
fSignSum[LL_List]:=Module[{r,i},r=Sum[Sign[LL[[i]]],{i,1,Length[LL]}];
    r]
(*##### fGrInc #######*)
(*# radi graf incidencije od prosirenog dowkera #*)
(*# Input ExtDow Output List of UnOrderedPairs (Graph) #*)
(*# Needs:#*)
fGrInc[Ul_List]:=Module[{rez={},sl},
    If[SameQ[Length[Ul],Length[Flatten[Ul]]],
      rez=fFormGr[Ul],
      rez=Flatten[Map[fFormGr[#]&,Ul],1]];
      rez=Sort[rez,Abs[#1]<Abs[#2]&];
      sl=Sort[Union[Flatten[rez]],Abs[#1]<Abs[#2]&];
     sl=Map[Sign[#]&,sl];
    {Sort[Abs[rez]],sl}]
(*############ fGraphInc ###########################*)
(*# isto sto i fGrinc, samo direktno iz 
Conwaya ili Dowkera ili pdata #*)
fGraphInc[Ul_]:=fGrInc[fGaussExtSigns[Ul]]
(*##########################*)
(*## ORIENTED LINKS ####*)
(*###################################################*)

(*######### OldLinkingNo ##################*)
OldLinkingNo[Ul_]:=
  Module[{GE,LinNoList={},d,i,j,k,p,rez={}},
    If[Not[SameQ[Head[Ul],String]],
      If[Length[Ul]!=1,
        If[SameQ[Select[Ul[[2]],OddQ[#1]&],{}],
        GE=fGaussExtSigns[Ul],GE=Ul],GE=Ul],
      GE=fGaussExtSigns[Ul]];
    d=Length[GE];
    If[Length[GE]==1,Print["Knot"],
      Do[Do[p=Map[Sign[#]&,Intersection[GE[[i]],GE[[j]]]];
          
          LinNoList=Append[LinNoList,Sum[p[[k]],{k,1,Length[p]}]],{j,i+1,
            d}],{i,1,d-1}];
      rez=
        Append[rez,{Sum[LinNoList[[k]],{k,1,Length[LinNoList]}],
            Sum[Abs[LinNoList[[k]]],{k,1,Length[LinNoList]}]}]
      (*Print["Linking number of oriented link is:  ",
          Sum[LinNoList[[k]],{k,1,Length[LinNoList]}]];
        Print["Linking number of nonoriented link is:  ",
          Sum[Abs[LinNoList[[k]]],{k,1,Length[LinNoList]}]]*)];
    Flatten[rez,1]/2]
(*radi od prosirenig Dowkera*)

(*#######Pomocne za fOrientedLink ###############*)
(*#####fFixLink #######*)
(*# Rotation of Components due to oddity #*)
fFixLink[Ul_List]:=
  Module[{FlD=Flatten[Ul],p,rez={Ul[[1]]},i},
    Do[p=Flatten[Position[FlD,Ul[[i,1]]]];
      If[EvenQ[p[[2]]-p[[1]]],rez=Append[rez,RotateLeft[Ul[[i]]]],
        rez=Append[rez,Ul[[i]]]],{i,2,Length[Ul]}];
    rez]
(*##### fVarP#######*)
fVarP[n_Integer]:=Module[{r,i},r=Table[IntegerDigits[i,2,n],{i,0,2^n-1}];
    r]
    
     (*##### fVarNewP #######*)
    
    fVarNewP[n_Integer]:=Module[{r,i},
    r=Drop[Table[IntegerDigits[i,2,n],{i,0,2^n-1}],1];
    r=Sort[Table[{Count[r[[i]],1],r[[i]]},{i,Length[r]}]];
    r=Drop[
        Union[Table[
            If[r[[i,1]]>IntegerPart[(n+1)/2],{},r[[i,2]]],{i,Length[r]}]],1];
    r]
(*##### fChangeCommSign #######*)
(*# Menja znake onih elemenata ostalih 
komponenti koji pripadaju komponenti sa kojom radimo #*)
fChangeCommSign[LL_List,p_List]:=
  Module[{i,rez={},q},Do[q=Map[If[MemberQ[Abs[p],Abs[#]],-#,#]&,LL[[i]]];
      rez=Append[rez,q],{i,1,Length[LL]}];
    rez]
    
   
(*##### fChangeOne #######*)
(*Menja znake u komponenti koju trenutno radimo samo onima koje nisu \
samopreseci*)
fChangeOne[DL_List,CL_List,b_Integer]:=
  Module[{p,q},p=Flatten[Complement[DL,{DL[[b]]}]];
    p=Intersection[Abs[p],Abs[CL]];
    q=Map[If[MemberQ[Abs[p],Abs[#]],-#,#]&,CL];
    (*q je komponenta*)
    {Reverse[q],p}]
(*vraca datu komponentu kako treba i p njene zajednicke sa ostalima*)


(*##### fMakeOL #######*)
fMakeOL[DL_List,Ind_List]:=
  Module[{rez={},i,pom=DL},
    Do[If[Ind[[i]]==0,rez=Append[rez,First[pom]];
        pom=Drop[pom,1],p=fChangeOne[DL,First[pom],i];
        (*Print["el:",p];*)rez=fChangeCommSign[rez,p[[2]]];
        (*Print[rez];*)rez=Append[rez,p[[1]]];
        (*Print[rez];*)pom=fChangeCommSign[Drop[pom,1],p[[2]]];
        (*Print[pom]*)],{i,1,Length[DL]}];
    rez]

(*####### fOrientedLink ##################*)

fOrientedLink[Ul_]:=
  Module[{GE,VarL,EL,ExtList,rez={},i},GE=fGaussExtSigns[Ul];
    (*GE je prosireni Gauss dowker*)
    If[Not[SameQ[GE,Flatten[GE]]],
      VarL=fVarP[Length[GE]];
      ExtList=Map[fMakeOL[GE,#]&,VarL];
      ExtList=Map[fFixLink[#]&,ExtList];
      EL=Map[OldLinkingNo[#][[1]]&,ExtList];
      Do[rez=Append[rez,{ExtList[[i]],VarL[[i]],EL[[i]]}],    
           {i,1,Length[ExtList]}],
      Print["Error: Knot"]];
    rez]
(*##################   PLANARNI ##########################*)

(* fPlanarEmb[neu_List]: 
    ulaz-Conway ili Dow izlaz- ulazni graf, graf izomorfan sa ulaznim,  
   planar embedding izomorfnog grafa, izomorfizam i faces izomorfnog grafa *)
   
fFindOrientation[L_List,VerLabel_List]:=Module[{p,i},
    p=Drop[TranslateVertices[L,-L[[1]]],1];
    p=Map[(#[[1]]+I*#[[2]])&,p];
    p=Map[Arg[#]&,p]+Pi;
    Do[p=ReplacePart[p,{p[[i]],VerLabel[[i+1]]},i],{i,1,Length[p]}];
    p=Prepend[Map[#[[2]]&,Reverse[Sort[p]]],First[VerLabel]];
    p]
    
fFind4[Graf_List]:=Module[{l},
    l=Union[Flatten[Graf]];
     l=Union[Map[Count[Flatten[ Graf],#]&,l]];
    Select[l,#!=4&]=={}
    ]
(* Vraca False ako nisu sva temena cetvorovalentna*)
(*Znaci True ako je graf \
u redu*)


fOrientGr[GG_Graph,neUredjeni_List]:=
  Module[{l,g,rez={},CRD,sus,j,i,ind=True,GNew},
    l=Flatten[FromAdjacencyMatrix[GG[[1]]][[1]],1];
    Do[sus=Select[l,MemberQ[#,i]&];
          sus=Map[If[SameQ[#[[1]],i],#[[2]],#[[1]]]&,sus];
         CRD=Map[Extract[GG[[2]],#]&,sus];
          
      rez=Append[rez,
          fFindOrientation[Prepend[CRD,GG[[2,i]]],Prepend[sus,i]]],
      {i,1,Length[GG[[1]]]}];
    g=Flatten[
        Table[Table[{rez[[j]][[1]],rez[[j]][[i]]},{i,2,Length[rez[[j]]]}],{j,
            Length[rez]}],1];
    g=Union[Map[Sort,g]];
    g=Isomorphism[FromUnorderedPairs[Union[neUredjeni]],
    FromUnorderedPairs[g],All];
    j=1;
    While[j<=Length[g]&&ind,                         
      GNew=Sort[
          Map[Sort,ReplaceAll[neu,Table[i->g[[j,i]],{i,Length[g[[j]]]}]]]];
           If[fFind4[GNew],ind=False,j++]
      ];
    {GNew,g[[j]],rez}
    ]
    
    
(*GG je novi graf,
  g[[j]] odgovarajuci izomorfizam a rez je planar embeding dat ciklim u \
temenima*)



fOrientGr1[GG_Graph,neUredjeni_List]:=
  Module[{l,g,rez={},CRD,sus,j,i,ind=True,iso,GNew},
    l=Flatten[FromAdjacencyMatrix[GG[[1]]][[1]],1];
    Do[sus=Select[l,MemberQ[#,i]&];
          sus=Map[If[SameQ[#[[1]],i],#[[2]],#[[1]]]&,sus];
         CRD=Map[Extract[GG[[2]],#]&,sus];
          
      rez=Append[rez,
          fFindOrientation[Prepend[CRD,GG[[2,i]]],Prepend[sus,i]]],
      {i,1,Length[GG[[1]]]}];
    g=Flatten[
        Table[Table[{rez[[j]][[1]],rez[[j]][[i]]},{i,2,Length[rez[[j]]]}],{j,
            Length[rez]}],1];
    g=Union[Map[Sort,g]];(*novi graf- dobijen bez duplih ivica*)
    
    g=Isomorphism[FromUnorderedPairs[Union[neUredjeni]],
    FromUnorderedPairs[g], All];
    j=1;
    If[fFind4[neUredjeni],
      While[j<=Length[g]&&ind,
                   
        GNew=Sort[
            Map[Sort,
              ReplaceAll[neUredjeni,Table[i->g[[j,i]],{i,Length[g[[j]]]}]]]];
             If[fFind4[GNew],ind=False,j++];
        iso=g[[j]]
        ],
      (*ako graf nije 4-
          valentan onda daj trivijalan izomorfizam i sam neu kao nov*)
      
      GNew=neUredjeni;iso=Union[Flatten[neUredjeni]]
      ];
    {GNew,iso,rez}
    ]

(*GG je novi graf,
  g[[j]] odgovarajuci izomorfizam a rez je planar embeding dat ciklim u \
temenima*)
(*za PlanarEmbGraf mi treba samo GG i rez*)


    




  fNasCikl[LL_List]:=Module[{i,kk,l,p},kk=LL;
    (*Print["Nas cikl radi sa: ",kk];*)
    l=kk[[1]];
    kk=Rest[kk];
    Do[p=Position[kk,l[[i+1]]][[1]];
      l=Append[l,kk[[p[[1]],ReplaceAll[Mod[p[[2]]+1,2],0->2]]]];
      kk=Drop[kk,{p[[1]]}],{i,1,Length[kk]-1}];
    (*Print["Nas cikl vraca: ",l];*)
    l
    ]
    
(*pomocna za fPljosni *)
fDtEd[LL_List]:=Map[{LL[[1]],#}&,LL[[2]]]
fSve[LL_List,elL_List]:=Union[Map[MemberQ[LL,#]&,elL]]
fMakeVert[LL_List]:=Map[First[#]&,LL]
fMakeEd[FF_List]:=Module[{i,ed={}},
    Do[
      ed=Append[ed,Take[FF,{i,i+1}]],{i,1,Length[FF]-1}];
      ed=Append[ed,{Last[FF],First[FF]}];
    ed=Map[Sort[#]&,ed]
    ]
    
fMakeVert[LL_List]:=Map[First[#]&,LL]
(* radi od Konture spoljne ,
    adj mat- matrica incidencije
      grafa,UnPair -neuredjeni parovi sa dodatim digonima*)

fPljosni[Kontura_List,adjmat_List,graf_,gUnPair_List]:=
  Module[{DT,isogr,p,n,isoLvert={},i,pos,pom,res={},lepi,tro},
    isogr=Union[gUnPair];
    n=Length[Union[Flatten[isogr]]];
    isogr=Union[isogr, Map[Reverse[#]&,isogr]];
    Do[p=  Select[isogr,#[[1]]==i&];
       isoLvert=Append[isoLvert,
          Flatten[Map[ Take[#,{2}]&,p],1]]
       ,{i,1,n}];
    DT=DelaunayTriangulation[graf[[2]]];
    lepi=Flatten[Map[fDtEd[#]&,DT],1];
    pom=Map[fMakeEd[#[[2]]]&,DT];
      Do[ pom=ReplacePart[pom,Map[Prepend[#,i]&,pom[[i]]],i] ,{i,1,
        Length[pom]}];
    pom=Map[fMakeEd[#]&,Flatten[pom,1]];
    res=Select[pom,fSve[isogr,#]=={True}&];
    res=Map[Union[Flatten[#]]&,res];
    res=Union[Map[Sort[#]&,res]];
    pom=Complement[pom,res];
    pom=Union[Map[Sort[#]&,pom]];
    lepi=Union[Map[Sort[#]&,Complement[lepi,isogr]]];
    res=Map[fMakeEd[#]&,res];
    While[lepi!={},
                pos=Position[pom,lepi[[1]]];
        p=Select[Union[   pom[[pos[[1,1]]]],
                    pom[[pos[[2,1]]]]],
                   (#1!=lepi[[1]])&];
     If[fSve[isogr,p]=={True},
                             res=Append[res,p],
                             pom=ReplacePart[pom,p,pos[[1,1]]];
        pom=Drop[pom,{pos[[2,1]]}]           ];    
        lepi=Rest[lepi]
         ];
    res=Map[fNasCikl[#]&,res];
    res
   ]
     
     
fFacOrient[faces_List]:=Module[{pom,ind,t,newfac,ind1,f1,ind0,indold,i},
    pom={};
    ind=Join[{1},Table[0,{i,Length[faces]-1}]];
    While[MemberQ[ind,0],
      indold=ReplaceAll[ind,0->1];
      newfac=Table[
          If[ind[[i]]>=0,faces[[i]],Reverse[faces[[i]]]],{i,
            Length[ind]}];
         ind1=Abs[ind];
      pom=
        Union[Flatten[
            Table[If[SameQ[ind1[[i]],0],{},
                If[SameQ[ind1[[i]],1],
                  Join[Table[{newfac[[i,j]],newfac[[i,j+1]]},{j,
                        Length[newfac[[i]]]-1}],{{Last[newfac[[i]]],
                        First[newfac[[i]]]}}],
                  Map[Reverse,
                    Join[Table[{newfac[[i,j]],newfac[[i,j+1]]},{j,
                          Length[newfac[[i]]]-1}],{{Last[newfac[[i]]],
                          First[newfac[[i]]]}}]]]],{i,Length[ind]}],1]];
      pom=Map[Reverse,pom];
      f1=Table[
          Join[Table[{newfac[[j]][[i]],newfac[[j]][[i+1]]},{i,
                Length[newfac[[j]]]-1}],{{Last[newfac[[j]]],
                First[newfac[[j]]]}}],{j,Length[newfac]}];
      f1=Drop[f1,1];
      ind0=
        Table[If[
            SameQ[Intersection[f1[[i]],pom],{}]&&
              SameQ[Intersection[Map[Reverse,f1[[i]]],pom],{}],0,1],{i,
            Length[f1]}];
      ind1=
        Table[If[
            SameQ[If[
                  SameQ[Intersection[f1[[i]],pom],{}]&&
        SameQ[Intersection[Map[Reverse,f1[[i]]],pom],{}],0,2],2]&&
              Not[SameQ[Intersection[f1[[i]],pom],{}]],1,-1],{i,
            Length[faces]-1}];
      ind=Join[{1},ind0*ind1];
      ind=indold*ind;
      newfac=
        Table[If[ind[[i]]>=0,faces[[i]],Reverse[faces[[i]]]],{i,
            Length[ind]}];
      ];
    newfac
    ]  



(*#############################################*)
 (*############# fPlanar Emb ###########################*)



(* ##### CONVEX DRAW by M.Ochiai, N.Imafuji, N.Morimura ####### *)


makematrix[mtx_, len_,xy_]:=Block[{i,j,k,d,len1,amt,bmt,lumat},
	len1=Length[mtx];
	amt=Table[0,{i,len1-len},{j,len1-len}];
	bmt=Table[0,{i,len1-len},{j,len}];
	i=len+1;k=1;
	While[i<=len1,
		d=Sum[mtx[[i,j]],{j,len1}];
		For[j=1,j<=len1,j++,
			If[mtx[[i,j]]==1,
				If[j<=len,bmt[[k,j]]=1/d,amt[[k,j-len]]=1/d]]];
		k++;
		i++];
    	lumat=LUDecomposition[IdentityMatrix[len1-len]-amt];
    	Return[LUBackSubstitution[lumat,bmt . xy]]
    ]
rearrange[mtx_,cyl_]:=Block[{i,j,len1,len2,cset,tmm1,tmm2},
	len1=Length[mtx];len2=Length[cyl];
	cset=Complement[Table[i,{i,len1}],cyl];
	tmm1=mtx;
	tmm2=Table[0,{i,len1},{j,len1}];
	For[j=0,j<2,j++,
		For[i=1,i<=len2,i++,tmm2[[i]]=tmm1[[cyl[[i]]]]];
		For[i=len2+1,i<=len1,i++,tmm2[[i]]=tmm1[[cset[[i-len2]]]]];
		tmm1=Transpose[tmm2];
	];
	Return[tmm1]
]
circluarvertices[n_]:=Block[{i,s=Pi/2-2 Pi / n,x=N[2 Pi / n]},
	  Return[
      Chop[N[Map[({{Cos[s],-Sin[s]},{Sin[s],Cos[s]}} . #)&,
            Table[{(Cos[x  i]),(Sin[x  i])},{i,n}]]]]]
	]


Next[v_List]:=Part[v,2]


FindPlotRange[v_List] :=
  Module[{xmin=Min[Map[First,v]], xmax=Max[Map[First,v]],
			ipsmin=Min[Map[Next,v]], ipsmax=Max[Map[Next,v]]},
		{ {xmin - 0.05 Max[1,xmax-xmin], xmax + 0.05 Max[1,xmax-xmin]},
		  {ipsmin - 0.05 Max[1,ipsmax-ipsmin], 
        ipsmax + 0.05 Max[1,ipsmax-ipsmin]}}] 
showgraphics2D[vtx_,prs_]:=Show[Graphics[
			Join[{PointSize[0.02]},Map[Point,Chop[vtx]],
					Map[(Line[Chop[vtx[[#]]]])&,prs]]
		],{AspectRatio->1,PlotRange->FindPlotRange[vtx]}]
showgraphics[vtx_,prs_]:=
		Show[Graphics3D[
			Join[{PointSize[0.02]},Map[Point,Chop[vtx]],
					Map[(Line[Chop[vtx[[#]]]])&,prs]]],{AspectRatio->1,PlotRange->All}]
innerproduct[v1_,v2_]:=
  Sqrt[(First[v1]-First[v2])^2+(Next[v1]-Next[v2])^2+(Last[v1]-Last[v2])^2]
averagepos[vtx_,vw_]:={
	Apply[Plus,Map[(vtx[[#,1]])&,vw]],
  Apply[Plus,Map[(vtx[[#,2]])&,vw]],
  Apply[Plus,Map[(vtx[[#,3]])&,vw]]}/Length[vw];
centerpos[vtx_]:=Block[{sx=0,sy=0,sz=0,len,i},
    len = Length[vtx];
    For[i=1,i<=len,i++,sx += vtx[[i,1]];sy += vtx[[i,2]]; 
      sz += vtx[[i,3]]];
    Return[{sx/len,sy/len,sz/len}]]
drawgraph[mtx_, cyl_]:=
  Block[{len1,len2,xy0, xy,pairs,mx,xyz,minlen,cnt,ccentar},
	len1=Length[mtx];len2=Length[cyl];
         xy0=circluarvertices[len2];
    	If[len1==len2,xy=circluarvertices[len1];mx=mtx,
      			mx=rearrange[mtx,cyl];
      			xy=Join[xy0,makematrix[mx,len2,xy0]]
      	];
    	pairs=Select[Position[mx, _?(Function[n,n != 0])],OrderedQ];
    	xy=Join[xy0,Chop[normalizedmatrix[mx,xy,Length[cyl], xy0]]];
        showgraphics2D[xy,pairs];
        xyz=Map[(Append[#,1/2])&,xy];
        For[i=1,i<=len1,i++, 
            xyz[[i]]  *= 
            xyz[[i,3]]/(2 xyz[[i,1]]^2+2 xyz[[i,2]]^2+2 xyz[[i,3]]^2)];
             ccentar=centerpos[xyz];
        While[N[innerproduct[{0.0,0.0,0.0},ccentar],20]>0.0000001,
          For[i=1,i<=len1,i++,  
              xyz[[i]] -= ccentar;   
        xyz[[i]]/=innerproduct[{0.0,0.0,0.0},xyz[[i]]]];
        ccentar=centerpos[xyz]
      ];
          If[Length[len1]>30,cnt=10,cnt=5];
          For[j=1,j<cnt,j++,
      		For[i=1,i<=len1,i++,
        			
        xyz[[i]]=(xyz[[i]]+
                averagepos[xyz,
                  Flatten[Position[ mx[[i]], _?(Function[n,n != 0]) ]]])/2.0;
        			xyz[[i]]/=innerproduct[{0.0,0.0,0.0},xyz[[i]]]]];		
    	xyz =N[xyz,20];
    	showgraphics[xyz,pairs]; 
    Return[{mx,xyz}]
    ]	
Edgesnn[Graph[e_,_]]:=e
Vert[Graph[_,v_]]:=v
Ve[Graph[e_,_]]:=Length[e]
Me[Graph[g_,_],___]:=Apply[Plus,Map[(Apply[Plus,#])&,g]]/2
Me[Graph[g_,_],Directed]:=Apply[Plus,Map[(Apply[Plus,#])&,g]]
ChVertices[g_Graph,v_List]:=Graph[Edgesnn[g],v]
ChEdges[g_Graph,e_List]:=Graph[e,Vert[g]]
DFS[v_Integer]:=(dfi[[v]]=cnt++;
    AppendTo[visit,v];
    Scan[(If[dfi[[#]]==0,AppendTo[eedgs,{v,#}];DFS[#]])&,e[[v]]])
DepthFirstTrans[g_Graph,start_Integer,flag_:Vertex]:=
  Block[{visit={},e=ToAdjacencyL[g],eedgs={},dfi=Table[0,{Ve[g]}],cnt=1},
    DFS[start];
    If[flag===Edge,eedgs,visit]]
ToAdjacencyL[Graph[g_,_]]:=
  Map[(Flatten[Position[#,_?(Function[n,n!=0])]])&,g]
ComplGraph[0]:=Graph[{},{}]
ComplGraph[1]:=Graph[{{0}},{{0,0}}]
ComplGraph[n_Integer?Positive]:=CircGraph[n,Range[1,Floor[(n+1)/2]]]
CircGraph[n_Integer?Positive,l_List]:=
  Module[{i,r},r=Prepend[MapAt[1&,Table[0,{n-1}],Map[List,Join[l,n-l]]],0];
    Graph[Table[RotateRight[r,i],{i,0,n-1}],CircVertices[n]]]
EmptyGr[n_Integer?Positive]:=
  Module[{i},Graph[Table[0,{n},{n}],Table[{0,i},{i,(1-n)/2,(n-1)/2}]]]
ComplGraph[l__]:=
  Module[{ll=List[l],t,i,x,rroww,stages=Length[List[l]]},
      t=FoldList[Plus,0,ll];
      Graph[
        Apply[Join,
          Table[rroww=
              Join[Table[1,{t[[i-1]]}],Table[0,{t[[i]]-t[[i-1]]}],
                Table[1,{t[[stages+1]]-t[[i]]}]];
            Table[rroww,{ll[[i-1]]}],{i,2,stages+1}]],
        Apply[Join,
          Table[Table[{x,i-1+(1-ll[[x]])/2},{i,ll[[x]]}],{x,stages}]]]]/;
    TrueQ[Apply[And,Map[Positive,List[l]]]]&&(Length[List[l]]>1)
CircVertices[0]:={}
CircVertices[n_Integer]:=
  Module[{i,x=N[2 Pi/n]},Chop[Table[N[{(Cos[x i]),(Sin[x i])}],{i,n}]]]
CircVertices[Graph[g_,_]]:=Graph[g,CircVertices[Length[g]]]
showgraphics2D[vtx_,prs_]:=Show[Graphics[
			Join[{PointSize[0.02]},Map[Point,Chop[vtx]],
					Map[(Line[Chop[vtx[[#]]]])&,prs]]
		],{AspectRatio->1,PlotRange->FindPlotRange[vtx]}]
defese[v_Integer]:=(dfi[[v]]=cnt++;
    AppendTo[visit,v];
    Scan[(If[dfi[[#]]==0,AppendTo[eedgs,{v,#}];defese[#]])&,e[[v]]])
depthfirstsearch[mtx_List,start_Integer]:=
  Block[{visit={},e,eedgs={},dfi=Table[0,{Length[mtx]}],cnt=1},
    e=Map[(Flatten[Position[#,_?(Function[n,n!=0])]])&,mtx];
    defese[start];
    Return[eedgs]]
getalledges[adjmat_]:=
  Select[Position[ adjmat, _?(Function[n,n != 0]) ],OrderedQ]
gettreeedges[adjmat_,vtx_,stpt_]:=
  If[Length[vtx]==0,
    DepthFirstTrans[Graph[adjmat,Vert[ComplGraph[Length[adjmat]]]],stpt,
      Edge],
    DepthFirstTrans[Graph[adjmat,vtx],stpt,Edge]]
getleafedges[alledges_,treeedges_]:=
  
  Complement[alledges,
    Union[Select[treeedges, OrderedQ],
      Map[Reverse,Complement[treeedges,Select[treeedges, OrderedQ]]]]]
ShowTree[adjmat_,vtx_,stpt_]:=
  showgraphics2D[vtx,gettreeedges[adjmat,vtx,stpt]]
Dijks[g_Graph,start_Integer]:=First[Dijks[g,{start}]]
Dijks[g_Graph,l_List]:=
  Module[{x,start,e=ToAdjacencyL[g],i,p,pparentt,untraversed},
    p=Edgesnn[PathConditionGraph[g]];
    Table[start=l[[i]];
      pparentt=untraversed=Range[Ve[g]];
      distan=p[[start]]; distan[[start]]=0;
      Scan[(pparentt[[#]]=start)&,e[[start]]];
      While[untraversed!={},x=First[untraversed];
        Scan[(If[distan[[#]]<distan[[x]],x=#])&,untraversed];
        untraversed=Complement[untraversed,{x}];
        Scan[(If[distan[[#]]>distan[[x]]+p[[x,#]],
                distan[[#]]=distan[[x]]+p[[x,#]];
                pparentt[[#]]=x])&,e[[x]]];];
      {pparentt,distan},{i,Length[l]}]]
ShortPath[g_Graph,s_Integer,e_Integer]:=
  Module[{pparentt=First[Dijks[g,s]],i=e,lst={e}},
    While[(i!=s)&&(i!=pparentt[[i]]),
      PrependTo[lst,pparentt[[i]]];
      i=pparentt[[i]]];
    If[i==s,lst,{}]]
    
MakeUndir[Graph[g_,v_]]:=
  Module[{i,j,n=Length[g]},
    Graph[Table[
        If[g[[i,j]]!=0||g[[j,i]]!=0,1,0],{i,n},{j,n}],v]]
FromOrdP[{}]:=Graph[{},{}]
FromOrdP[l_List,v_List]:=
Graph[MapAt[1&,Table[0,{Length[v]},{Length[v]}],l],v]
FromUnordP[l_List]:=MakeUndir[FromOrdP[l]]
FromUnordP[l_List,v_List]:=MakeUndir[FromOrdP[l,v]]
ShortestPathSTree[g_Graph,s_Integer]:=
  Module[{pparentt=First[Dijks[g,s]],i},
    FromUnordP[Map[({#,pparentt[[#]]})&,Complement[Range[Ve[g]],{s}]],
      Vert[g]]]
AllPairsShortestP[g_Graph]:=
  Module[{p=Edgesnn[PathConditionGraph[g]],i,j,k,n=Ve[g]},
      Do[p=Table[Min[p[[i,k]]+p[[k,j]],p[[i,j]]],{i,n},{j,n}],{k,n}];
      p]/;Min[Edgesnn[g]]<0
AllPairsShortestP[g_Graph]:=Map[Last,Dijks[g,Range[Ve[g]]]]
RemoveSelfL[g_Graph]:=Module[{i,e=Edgesnn[g]},Do[e[[i,i]]=0,{i,Ve[g]}];
    Graph[e,Vert[g]]]
PathConditionGraph[Graph[e_,v_]]:=
  RemoveSelfL[Graph[ReplaceAll[e,0->Infinity],v]]
characteristiccyclematrix[adjmat_,vtx_,stpt_]:=
  Block[{i,len,len1,len2,at,al,alle,tttree,leaf,idm,bpm},
    len=Length[adjmat];
    alle=getalledges[adjmat];
    tttree=depthfirstsearch[adjmat,stpt];
    leaf=getleafedges[alle,tttree];
    len1=Length[tttree];
    len2=Length[leaf];
    at=Table[0,{i,1,len},{j,1,len1}];
    al=Table[0,{i,1,len},{j,1,len2}];
    For[i=1,i<=len1,i++,Map[(at[[#,i]]=1)&,tttree[[i]]]];
    For[i=1,i<=len2,i++,Map[(al[[#,i]]=1)&,leaf[[i]]]];
    idm=IdentityMatrix[len2];
    bpm=Mod[Transpose[Inverse[Delete[at,len],Modulus->2] . Delete[al,len]],
        2];
    Return[Map[(Flatten[Append[idm[[#]],bpm[[#]]]])&,Table[i,{i,1,len2}]]]
    ]
showgraphics2D[vtx_,prs_]:=Show[Graphics[
			Join[{PointSize[0.02]},Map[Point,Chop[vtx]],
					Map[(Line[Chop[vtx[[#]]]])&,prs]]
		],{AspectRatio->1,PlotRange->FindPlotRange[vtx]}]
ToOrdPairs[g_Graph]:=Position[Edgesnn[g],_?(Function[n,n!=0])]
TranslateVert[v_List,{x_,y_}]:=Map[(#+{x,y})&,v]
TranslateVert[Graph[g_,v_],{x_,y_}]:=Graph[g,TranslateVert[v,{x,y}]]
DilateVert[v_List,d_]:=(d*v)
DilateVert[Graph[e_,v_],d_]:=Graph[e,DilateVert[v,d]]
NormalizeVert[v_List]:=Module[{v1},v1=TranslateVert[v,{-Min[v],-Min[v]}];
    DilateVert[v1,1/Max[v1,0.01]]]
NormalizeVert[Graph[g_,v_]]:=Graph[g,NormalizeVert[v]]
PointsAndLines[Graph[e_List,v_List]]:=
  Module[{pairs=ToOrdPairs[Graph[e,v]]},
    Join[{PointSize[0.025]},Map[Point,Chop[v]],
      Map[(Line[Chop[v[[#]]]])&,pairs]]]
ShowLabelGraph[g_Graph]:=ShowLabelGraph[g,Range[Ve[g]]]
ShowLabelGraph[g1_Graph,labels_List]:=
  Module[{pairs=ToOrdPairs[g1],g=NormalizeVert[g1],v},v=Vert[g];
    Show[Graphics[
        Join[PointsAndLines[g],Map[(Line[Chop[v[[#]]]])&,pairs],
          GraphLabels[v,labels]]],{AspectRatio->1,
        PlotRange->FindPlotRange[v]}]]
GraphLabels[v_List,l_List]:=
  Module[{i},Table[Text[l[[i]],v[[i]]-{0.03,0.03},{0,1}],{i,Length[v]}]]


    
     (*############# fPlanar Emb ###########################*)

fFindPre[trazi_Integer, gde_Integer, CC_List] := 
  Module[{l = CC[[gde]], nadjen, pom},
    pom = Position[l, trazi][[1, 1]];
   If[pom == 1, nadjen = Last[l], nadjen = l[[pom - 1]]]; nadjen
    ]
    
fCloseFace[CC_List, Ed_List, PP_List] := 
  Module[{sledeci, izbaci = {}, poc = PP},
    sledeci = fFindPre[poc[[1]], poc[[2]], CC];
    While[sledeci != poc[[1]],
            izbaci = Append[izbaci, {Last[poc], sledeci}];
           poc = Append[poc, sledeci];
           sledeci = fFindPre[poc[[-2]], Last[poc], CC]
            ];
    izbaci = Append[izbaci, {Last[poc], sledeci}];
    {poc, Complement[Ed, izbaci]}
    ]
    
fFaces[CiklL_List] := Module[{Ul = CiklL, t, res, konacna = {}},
    (*Lista svih stranica*)
    t = Sort[
        Flatten[Table[
            Table[{i, Ul[[i, j]]}, {j, Length[Ul[[i]]]}], {i, Length[Ul]}], 
          1]];
        poc = First[t];  
        t = Rest[t];
      While[ t != {},
         res = fCloseFace[CiklL, t, poc];
        konacna = Append[konacna, res[[1]]];
        t = res[[2]];
      If[t != {}, poc = First[t]; t = Rest[t]]
      ];
    konacna
    ]  
    
fWrGraph[g_List, file_] := Module[{g1, edg, v, i, x, y, p},
    g1 = FromUnorderedPairs[g];
    p = Length[Union[Flatten[g]]];
    edg = ToAdjacencyMatrix[g1];
    edg = Table[Flatten[Position[edg[[i]], 1]], {i, Length[edg]}];
     v = Flatten[N[NormalizeVertices[g1[[2]]]], 1]; 
    OpenWrite[file];
    WriteString[file, "N=", ToString[p] <> "\n"];
    Do[WriteString[file, "", ToString[i - 1] <> ":"];
       {x, y} = Chop[v[[i]]];
      Scan[(WriteString[file, " ", ToString[# - 1]]) &, edg[[i]]] ;
      WriteString[file, " -1"];
      Write[file], {i, Length[v]}];
    Close[file];]
    
fPlanEmbedding[g_List] := Module[{vv, vv1, vv2, p,i},
    fWrGraph[g, "graph.txt"]; 
    Run["planarity" <> " graph.txt" <> " graph1.txt"];
    vv = Import["graph1.txt"];
    DeleteFile["graph.txt"];
    DeleteFile["graph1.txt"];
    vv = StringReplace[
        vv, {"\n" -> "", " -1" -> "}", ": " -> "{", " " -> ","}];
    vv1 = StringPosition[vv, "{"];
    vv2 = StringPosition[vv, "}"];
    p = Length[vv1];
    vv = StringDrop[vv, vv1[[1, 1]] - 1];
    For[i = 1, i < p, vv1 = Drop[StringPosition[vv, "{"], 1];
      vv2 = Drop[StringPosition[vv, "}"], -1]; 
      vv = StringReplacePart[vv, ",", {vv2[[i, 1]] + 1, vv1[[i, 1]] - 1}]; 
      i++];
    vv = ToExpression["{" <> vv <> "}"] + 1;
    vv
    ]
    
fPlanarEmb[Ul_] := Module[{vv, kk,i},
    kk = fGraphInc[Ul][[1]];
    vv = fPlanEmbedding[Union[kk]];
    vv = {kk, kk, Union[Flatten[kk]], 
        Table[Join[{i}, vv[[i]]], {i, Length[vv]}], 
        fFaces[vv]};
    vv
    ]


(*###########################################################*)
    
 fPlanarEmbKL[Ul_] := Module[{vv,kk,i},
    kk = fGraphInc[Ul][[1]];
    vv = fPlanEmbedding[Union[kk]];
    vv = {kk,Table[Join[{i}, vv[[i]]], {i, \
Length[vv]}],fFaces[vv]};
    vv
    ]
    
(*###########################################################*)    

(*########### fPlanarEmbGraph #############################*)

fPlanarEmbGraph[Ul_List] := Module[{vv,kk,i},
    kk = fPlanEmbedding[Union[Ul]];
    vv = Table[Join[{i}, kk[[i]]], {i, Length[kk]}];
    {Ul,vv,fFaces[kk]}
    ]
    
(*###########################################################*)
(*########### DrawPLanarEmbKL ######################*)



 
 Dijks[g_Graph,start_Integer]:=First[Dijks[g,{start}]]

Dijks[g_Graph,l_List]:=
  Module[{x,start,e=ToAdjacencyL[g],i,p,pparentt,untraversed},
    p=Edgesnn[PathConditionGraph[g]];
    Table[start=l[[i]];
      pparentt=untraversed=Range[Ve[g]];
      distan=p[[start]]; distan[[start]]=0;
      Scan[(pparentt[[#]]=start)&,e[[start]]];
      While[untraversed!={},x=First[untraversed];
        Scan[(If[distan[[#]]<distan[[x]],x=#])&,untraversed];
        untraversed=Complement[untraversed,{x}];
        Scan[(If[distan[[#]]>distan[[x]]+p[[x,#]],
                distan[[#]]=distan[[x]]+p[[x,#]];
                pparentt[[#]]=x])&,e[[x]]];];
      {pparentt,distan},{i,Length[l]}]]
    
findfacialcycle[adjmat_,vtx_,stpt_]:=
  Block[{i,j,len,len1,len2,sptreem,ssptr,al,alle,tttree,leaf,PPT,fpath,LPT,
      cccount,adjmat1,vtx0,vtx1},
    len=Length[adjmat];
    alle=getalledges[adjmat];
    tttree=depthfirstsearch[adjmat,stpt];
    leaf=getleafedges[alle,tttree];
    len1=Length[tttree];
    len2=Length[leaf];
    sptreem=Table[0,{i,1,len},{j,1,len}];
    For[i=1,i<=len1,i++,sptreem[[tttree[[i,1]],tttree[[i,2]]]]=1];
    sptreem=sptreem + Transpose[sptreem];
    If[Length[vtx]==0,vtx0=Vert[ComplGraph[len]],vtx0=vtx];
    ssptr=Graph[sptreem,vtx0];
    fpath={};
 For[i=1,i<=len2,i++,
      cccount=0;adjmat1=.;vtx1=.;
      PPT=ShortPath[ssptr,leaf[[i,1]],leaf[[i,2]]];
       LPT=Map[({#})&,PPT];
     For[j=1,j<=len2,j++,
        If[
          i!=j && MemberQ[PPT,leaf[[j,1]]] && 
            MemberQ[PPT,leaf[[j,2]]],cccount++];
        If[cccount==1,Break[]]];
      If[Length[adjmat]==Length[PPT],adjmat1=={};vtx1={},
        adjmat1=Delete[Transpose[Delete[adjmat,
                LPT]],LPT]];
      vtx1=Delete[vtx0,LPT];
     If[adjmat1=={} || vtx1=={} ,cccount++,
        If[Length[DepthFirstTrans[Graph[adjmat1,vtx1],1]]!=Length[adjmat1],
          cccount++]];
      If[cccount==0,If[Length[fpath]<Length[PPT],fpath=PPT]]
      ];
    Return[fpath]
    ]  
    

DrawPlanarEmbKL[Ul_]:=
  Module[{q=0,ind,Ul1,neu,neu1,adjmat,cikl,len1,len2,xy,xy0,mx,pairs},
    Ul1=Ul;
    If[SameQ[Head[Ul],String],
      ind=Length[
            ReadList[
              StringToStream[
                StringReplace[
                  Ul,{"*"->" 6 ","."->" 6 ","("->" 6 ",")"->" 6 ",
                    ","->" 6 "}]],Number]]==1
       ];
If[SameQ[Head[Ul],List],
      If[Length[Ul[[1]]]==1,q=Ul[[2]],q=iteratedTake[Ul[[2]],Ul[[1]]]];
          ind=SameQ[q,Dowker[ToString[Apply[Plus,Ul[[1]]]]][[4]]];
         If[ind,Ul1=ToString[Apply[Plus,Ul[[1]]] ]]
      ];
    If[ind==False,
           neu=fGaussExtSigns[Ul];
           neu=fGrInc[neu][[1]];
            neu1=neu;
           adjmat=ToAdjacencyMatrix[FromUnorderedPairs[Union[neu]]];
           cikl=findfacialcycle[adjmat,NULL,1];
           len1=Length[adjmat];len2=Length[cikl];
          If[len1==len2,xy=circluarvertices[len1];  mx=adjmat,
        			mx=rearrange[adjmat,cikl];
        			xy0=circluarvertices[len2];
        			xy=Join[xy0,makematrix[mx,len2,xy0]]
                         ];			
           pairs=Position[mx, _?(Function[n,n != 0])];
        ShowLabelGraph[Graph[mx,xy]],
      ShowLabeledGraph[FromUnorderedPairs[fMakeEd[Range[ToExpression[Ul1]]]]]
      ];
    ](*16.8.2003*)
 
 (*######### Draw Planar Embeding of the Graph ##############*)
 
DrawPlanarEmbGraph[Ul_]:=
  Module[{ind,neu,adjmat,cikl,len1,len2,xy,mx,xy0,pairs},
    neu=Union[Ul];
    ind=Union[Map[Count[Ul,#]&,Ul]];
    If[SameQ[ind,{4}]||SameQ[ind,{2}],
           (*slucaj 2 li -2 ili CIKL *)        
      ShowLabeledGraph[FromUnorderedPairs[neu]],
        adjmat=ToAdjacencyMatrix[FromUnorderedPairs[neu]];
      cikl=findfacialcycle[adjmat,NULL,1];
      len1=Length[adjmat];len2=Length[cikl];
      If[len1==len2,xy=circluarvertices[len1];
             mx=adjmat,
            mx=rearrange[adjmat,cikl];xy0=circluarvertices[len2];
            xy=Join[xy0,makematrix[mx,len2,xy0]]
        ];
      pairs=Position[mx,_?(Function[n,n!=0])];
      ShowLabelGraph[Graph[mx,xy]]
      ];
    ]
(* 16.08.2003 *)


(*########################################### *)  
(*########### fSignsKL ######################*)
(* fSignsKL racuna znake Alt KL. 
      Ulaz: Dow u smislu Knotscape. 
            Izlaz: Dow sa znacima. Dow za Alt KL unosimo bez znakova, 
  a za Nonalt samo znake tacaka koje su izmenjene u odnosu na Alt (sign \
change) *) 

UnsortedUnion[x_]:=Module[{f},f[y_]:=(f[y]=Sequence[];y);f/@x]

fSignsKL[Dow_List]:=
  Module[{DS,DowA,OldGaussExt,Emb,iso,i,j,CC,cik,q=1,QQQ={},uu,zakraj,eee},
    DS=Table[Sign[Dow[[2]][[i]]],{i,Length[Dow[[2]]]}];
    DowA=Abs[Dow];
    If[SameQ[Length[DowA[[1]]],1],OldGaussExt={fGaussExtSigns[DowA]},
      OldGaussExt=fGaussExtSigns[DowA]];
      eee=fPlanarEmb[DowA];
      Emb=eee[[4]];
     iso=eee[[3]];
    Emb=Sort[
        ReplaceAll[Emb,
          Table[iso[[i]]->i,{i,Length[iso]}]],#1[[1]]<#2[[1]]&];
    zakraj=Emb;
    Emb=Map[Rest,Emb];(*izomorfni PlanarEmb*)
    
    CC=Table[Position[OldGaussExt,i],{i,
          Length[Union[Flatten[OldGaussExt]]]}];
    CC=Map[If[OddQ[#[[1,2]]],#,Reverse[#]]&,CC];
    (*izgleda da pravimo novi Dow tj.ovih parova cemo uzimati druge \
koordinate za nove parove kod kojih je prvi neparan i po njima cemo \
sortirati*)

        CC=Table[{First[RotateLeft[OldGaussExt[[CC[[i,1,1]]]],
        CC[[i,1,2]]-2]],
          First[RotateLeft[OldGaussExt[[CC[[i,2,1]]]],CC[[i,2,2]]-2]],
          First[RotateLeft[OldGaussExt[[CC[[i,1,1]]]],CC[[i,1,2]]]],
          First[RotateLeft[OldGaussExt[[CC[[i,2,1]]]],CC[[i,2,2]]]]},{i,
          Length[CC]}];
    (*pravimo malu tabelicu ??? pomocu GaussExt i CC*)
    
    cik=Table[UnsortedUnion[CC[[i]]],{i,Length[CC]}];
    (*cik su orginalni CC bez dvojnih veza*)
    uu=Map[UnsortedUnion,Emb];
    (*uporedjujemo orijentaciju cik i uu.Ako neki 
    deo cik i Emb ima samo dva 
clana,pisemo 0*)
QQQ=Table[If[SameQ[Length[Union[cik[[j]]]],2],
          			0,
          			
          If[MemberQ[Table[RotateLeft[cik[[j]],i],{i,Length[cik[[j]]-1]}],
              uu[[j]]],
            		               1,
                                          -1]],
        		{j,Length[cik]}];
    QQQ=If[SameQ[Union[QQQ],{0}],Table[1,{i,Length[QQQ]}],QQQ];
    (*temena sa nulama cuvaju znak prethodnika*)
    (*Sad jos treba da \
odredimo u sta idu nule !!!!!*)
    (*u Emb trazimo kom lancu digona to teme \
sa 0 pripada i uzmemo taj znak*)
      pp=Flatten[Union[Position[QQQ,0]]];
    Do[ cik=Rest[zakraj[[pp[[i]]]]];
            cik=Select[cik,Not[SameQ[QQQ[[#]],0]]&][[1]];
        QQQ=ReplacePart[QQQ,QQQ[[cik]],pp[[i]]]
       ,{i,1,Length[pp]}];
    uu=Table[Position[OldGaussExt,i],{i,Length[QQQ]}];
    Do[OldGaussExt=
        ReplacePart[OldGaussExt,i*QQQ[[i]],{uu[[i,1]],uu[[i,2]]}],{i,
        Length[uu]}];
    OldGaussExt=fDowfromGaussExt[OldGaussExt];
    OldGaussExt={OldGaussExt[[1]],-Flatten[OldGaussExt[[2]]]*Flatten[DS]}]
(*3.12.2003*)


(*############################################################*)
(*************************)
(*pomocna za GraphKL*)
    
fNeighFac[LL_List,ff_List]:=Module[{l,p,dig,Lose,pom,i},
    dig=Union[Flatten[Select[LL,Length[#]==2&],1]];
    p=Position[LL,ff][[1,1]];
    l=Complement[Select[LL,Length[Intersection[#,ff]]>=1&],{ff}];
    Lose=Select[l,Length[#]!=2&];
    Lose=Map[{#,Complement[Intersection[#,ff],dig]}&,Lose];
    Lose=Select[Lose,#[[2]]=={}&];
    If[Lose!={},Lose=Map[Take[#,1]&,Lose]];
    Lose=Flatten[Lose,1];
    l=Complement[l,Lose];
    l=Map[Position[LL,#][[1,1]]&,l];
    Prepend[l,p]
    ]
fMakeCircEd[LL_List]:=Module[{i,res={}},
       Do[  res=Append[res,{LL[[1]],LL[[i]]}],{i,2,Length[LL]}];
    res]
  (*##############################################*)
(*############# fGraphKL ######################*)

fMakeEd[FF_List]:=Module[{i,ed={}},
    Do[
      ed=Append[ed,Take[FF,{i,i+1}]],{i,1,Length[FF]-1}];
      ed=Append[ed,{Last[FF],First[FF]}];
    ed=Map[Sort[#]&,ed]
    ]
    
(*fGraphKL:
    input Conway ili Dow,output:graf KL*)
   
fGraphKL[Ulaz_]:=
  Module[{ind,i,pe,faces,NewGr,pom,facesEd,pos,p},
    If[Head[Ulaz]==List,Ul=Abs[Ulaz];ind=Ul[[1]]];
    If[Head[Ulaz]==String,Ul=Ulaz;ind=fDToD[Ulaz][[1]]];
    If[SameQ[ind,{1,1}],
      (* slucaj "2","-2" *)
      pom={{1,2},{1,2}},
      (*ostali *)
      pe=fPlanarEmb[Ul];
      NewGr=pe[[2]];
      pom=Union[NewGr];
      faces=pe[[5]];(*Print[faces];*)      
      Do[If[Count[NewGr,pom[[i]]]==2,faces=Append[faces,pom[[i]]]],{i,1,\
          Length[pom]}];
      (*fac su faces+digons*)(*sad treba da nadjemo susedne*)
      facesEd=
        Map[fMakeEd,faces];
      nbs=Map[fNeighFac[facesEd,#]&,facesEd];
      (*cikli susednih oko faca*)pom=
        Map[{#,0,0}&,Range[Length[nbs]]];(*brojevi svih ivica*)pom=
        ReplacePart[pom,{1,0,1},1];(*bela je 0;crna je 1*)(*redni broj,
        nula za boju,nula za da li je do sada obojena*)While[
        Union[Flatten[Map[Take[#,-1]&,pom]]]!={1},
        Do[If[pom[[i,3]]==1,
            pom=Map[If[#[[3]]!=1&&MemberQ[nbs[[i]],#[[1]]],{#[[1]],
                      1-pom[[i,2]],1},#]&,pom]],{i,1,Length[nbs]}]];
      pom=
        Select[pom,#[[2]]==1&];(*pom je lista koj asdarzi crne \
pljosni*)
        pom=Flatten[Map[Take[#,1]&,pom]];
      (*redni brojevi crnih pljosni*)faces=
        Map[Prepend[#,Position[faces,#][[1,1]]]&,faces];
      faces=Select[faces,MemberQ[pom,#[[1]]]&];
      Do[faces=ReplacePart[faces,Flatten[{i,Rest[faces[[i]]]}],i],{i,1,
          Length[faces]}];
      faces=Flatten[Map[fMakeCircEd,faces],1];
      (*lista parova {stranica,teme} a mi cemo sad dobiti {str,str}*)n=
        Length[Union[Take[Flatten[faces],{2,Length[Flatten[faces]],2}]]];
      (*ovo su u sustini smao sva temena*)pom={};
      Do[p=Select[faces,#[[2]]==i&];
        pos=Map[Position[faces,#][[1,1]]&,p];
        pom=Append[pom,{p[[1,1]],p[[2,1]]}];
        faces=Drop[faces,{pos[[1]]}];
        faces=Drop[faces,{pos[[2]]-1}],{i,1,n}];
      ShowLabeledGraph[FromUnorderedPairs[pom]]
      ];
    pom]
    
(*#####################################################*)
(*############# fGenerators ########################*)
(*^^^^^^^^^^^^^^^^*)
(*checkArgs[s_,t_]:=
  If[Head[s]===List&&VectorQ[t,Head[#]===Integer&&#>=0&]&&
      Plus@@t<=Length[s],True,False]

iteratedTake[s_,t_]/;checkArgs[s,t]:=
  First/@Rest[FoldList[Through[{Take,Drop}[#1[[2]],#2]]&,{{},s},t]]*)

(*fGenerators[Ul_] od Ul_=
    Conway ili Ul_=
      Dow racuna generatore KL Izlaz:lista generatora UIP sa znacima*)
(*^^^^^^^^^^^^^^^^^^^^*)
fGenerators[Ul_]:=
  Module[{ges,t,ou,lpoz,pod,lneg,pro,qq,qq1,rr,ss,qq2,qq3,pro1,pro2,
  ul,iz,gen, qq4,sg, i, j},
    ges=fGaussExtSigns[Ul];
    (*ges nam odredjuje znake*)
    
    t=If[SameQ[ges,Flatten[ges]],Length[ges],Map[Length,ges]];
    ou=If[
        SameQ[Head[Ul],
          String],-Map[Sign,Flatten[ges]]*-Map[Sign,
              Flatten[fGaussExtSigns[StringReplace[Ul,"-"->""]]]]*
          Flatten[Abs[ges]]*Table[(-1)^i,{i,Length[Flatten[ges]]}],
        Map[Sign,Flatten[fGaussExtSigns[fDowfromPD[Ul]]]]*
          Map[Sign,Flatten[ges]]*Flatten[Abs[ges]]*
          Table[(-1)^i,{i,Length[Flatten[ges]]}]];
    t=If[SameQ[Head[t],Integer],{t},t];
    ou=iteratedTake[ou,t];
    ou=If[ou[[1,1]]>0,ou,-ou];
    (*ou=over-under iz koje izdvajamo generatore*)
    
    lpoz=Table[Select[ou[[i]],#>0&],{i,Length[ou]}];
    t=Map[Length,lpoz];
    lpoz=Flatten[lpoz];
    lpoz=Flatten[Table[Position[ou,lpoz[[j]]],{j,Length[lpoz]}],1];
    pod=iteratedTake[Table[ou[[lpoz[[i,1]],lpoz[[i,2]]]],{i,Length[lpoz]}],
        t];
    pod=Map[Length,pod];
    (*pod nam odredjuje raspodelu generatora po komponentama*)lneg=
      Table[Select[ou[[i]],#<0&],{i,Length[ou]}];
    pro=Flatten[
        Table[Append[
            Table[{lneg[[j,i]],lneg[[j,i+1]]},{i,Length[lneg[[j]]]-1}],{Last[
                lneg[[j]]],First[lneg[[j]]]}],{j,Length[lneg]}],1];
    (*pro je flatten lista negativnih uzastopnih u ou po komponentama*)qq=
      Table[Take[
          RotateLeft[ou[[Flatten[Position[ou,pro[[i,1]]]][[1]]]],
            Flatten[Position[ou,pro[[i,1]]]][[2]]],
          Flatten[Position[
                  RotateLeft[ou[[Flatten[Position[ou,pro[[i,1]]]][[1]]]],
                    Flatten[Position[ou,pro[[i,1]]]][[2]]],
                  pro[[i,2]]]][[1]]-1],{i,Length[pro]}];
    qq1=qq;
    rr=Complement[Union[Flatten[qq]],
        Flatten[If[
            SameQ[Union[
                  Table[If[Not[SameQ[qq[[i]],{}]],{qq[[i,1]]},{}],{i,
                      Length[qq]}]][[1]],{}],
            Drop[Union[
                Table[If[Not[SameQ[qq[[i]],{}]],{qq[[i,1]]},{}],{i,
                    Length[qq]}]],1],
            Union[Table[
                If[Not[SameQ[qq[[i]],{}]],{qq[[i,1]]},{}],{i,
                  Length[qq]}]]]]];
    ss=Flatten[Position[qq,{}]];
    Do[qq=ReplacePart[qq,{rr[[i]]},ss[[i]]],{i,Length[rr]}];
    qq2=Table[qq[[i,1]],{i,Length[qq]}];
    qq3=Select[qq1,Not[SameQ[#,{}]]&];
    pr=Flatten[
        Table[Table[qq3[[j,1]],{i,Length[qq3[[j]]]}],{j,Length[qq3]}]];
    qq3=-Flatten[qq3];
    pro1=Table[pro[[i,1]],{i,Length[pro]}];
    pro2=Table[pro[[i,2]],{i,Length[pro]}];
    ul=Flatten[
        Table[Flatten[qq2[[Flatten[Position[pro2,qq3[[i]]]]]]],{i,
            Length[qq]}]];
    iz=Flatten[
        Table[Flatten[qq2[[Flatten[Position[pro1,qq3[[i]]]]]]],{i,
            Length[qq]}]];
    gen=Flatten[Table[{ul[[i]],iz[[i]],pr[[i]]},{i,Length[ul]}]];
    gen=iteratedTake[gen,3*pod];
    ges=Union[Flatten[ges]];
    qq4=Abs[qq3];
    sg=Flatten[
        Table[Sign[ges[[Flatten[Position[Abs[ges],qq4[[i]]]]]]],{i,
            Length[ges]}]];
    {gen,sg}]
    (*23.11.2003*)
    
(*############# fColTest ############################*)
par[n_,1]:={{n}}
par[n_,r_]:=Module[{ans={},i},
    Do[
      AppendTo[ans,(Join[{i},#]&)/@par[n-i,r-1]],{i,0,n}];
    Flatten[ans,1]]
(* fColTab[stNo_Integer, cNo_Integer] 
formira sve izbore boja za KL sa stNo 
generatora i cNo boja *)
fColTab[stNo_Integer, cNo_Integer]:=Module[{tn,t, i, j},
    tn=Table[i,{i,cNo}];
    t=Union[Flatten[Table[Map[Sort,par[stNo,cNo]],{i,stNo}],1]];
    t=Union[
        Table[If[MemberQ[t[[i]],0]==False,t[[i]],{}],{i,Length[t]}]];
    t=If[SameQ[t[[1]],{}]==True,Drop[t,1],t];
    t=Table[
        Flatten[Table[Table[tn[[i]],{t[[j]][[i]]}],{i,Length[t[[j]]]}],1],{j,
          Length[t]}];
    t=Flatten[
        DistinctPermutations/@
          Union[If[SameQ[t[[1]],{}]==True,Drop[t,1],t]],1]
    ]
fChangeCol[L_List,CL_List]:=Module[{i},
    ReplaceAll[L,Table[i->CL[[i]],{i,Length[CL]}]]
    ]
    
    (*pomocna za fColTest- vraca 0,1,2,
  3 u zavisnosti da li je bojenje kan nemoguce,savrseno,
  standardno ili kaufmanovo*)
fColType[kan_List]:=Module[{pom,p,res},
    pom=Partition[kan,3];
    p=Map[Length[Union[#]]&,
        pom];(*broj razlicitih boje medju uip*)
    (*Print["broj boja",pom,
          p];*)
    If[SameQ[Union[p],{3}],
                 res=1(*Perfect*),
                  If[SameQ[Union[p],{1,3}],
                         res=2(*Standard*),
                         (*ili je Kauffman ili nije bojenje*)
                \
        pom=Map[MemberQ[{#[[1]],#[[2]]},#[[3]]]&,pom];
        	      If[Not[MemberQ[pom,True]],res=3(*Kauffman*),
                              res=0(*not coloring*)]]];
    (*Print[res];*)(*
      
      Switch[res,0,Print["nema"],1,Print["perfect"],2,Print["standard"],3,
          Print["Kauffman"]];*)
    res]
    
(* fColTest radi od Conwaya ili Dow i broja boja. 
Daje prvo nadjeno bojenje KL *)
fColTest[Ul_,cNo_Integer]:=
  Module[{g,s,stNo,ct,i,resP={},resS={},resK={},stop=True,maxCNo,ind},
    maxCNo=Length[fGenerators[Ul][[2]]];(*max no of colors*)
    If[cNo>maxCNo,
      (*previse boja*)
      s="Maximum number of colors is "<>ToString[l],
      (*trazimo bojenje sa cNo boja*)
      g=Flatten[Map[Reverse,
            
            Sort[Map[Reverse,
                Partition[Flatten[fGenerators[Ul][[1]]],3] ] ] ] ];
      Print["Generators: ",Partition[g,3]];
      (*g su  generatori u formi uip *)
      stNo=Length[g]/3;
      ct=fColTab[stNo,cNo];(*ct sadrzi sva moguca bojenja*)
      (*Print[" broj mogucih bojenja ",Length[ct]];*)
      i=1;ind=True;
      While[i<=Length[ct]&&ind&&stop,
             p=fChangeCol[g,ct[[i]]];(*kandidat*)
        (*Print[i,"Kandidat ",p];*)i=i+1;  
        ind=fColType[p];
             Switch[ind,
          	1,If[resP=={},resP=Append[resP,p]],
          	2,If[resS=={},resS=Append[resS,p]],
          	3,If[resK=={},resK=Append[resK,p]]]; 
        If[ind==1,stop=False];
             ind=Or[SameQ[resP,{}],SameQ[resS,{}],SameQ[resK,{}]]
          
          ]];
    If[Union[resP,resS,resK]=={},
          Print["No ",cNo,"-colorings"],
      If[Not[SameQ[resP,{}]],
        Print["First found perfect coloring: ",Partition[resP[[1]],3]] ;
        resK={};resS={},
        Print["No perfect colorings"];resP={0};
                   If[resS!={},
             Print["First found standard coloring: ",Partition[resS[[1]],3]];
          resK={},
             Print["No standard colorings"];resS={0};
            If[resK!={},
                  
            Print["First found Kauffman coloring: ",Partition[resK[[1]],3]],
                Print["No Kauffman colorings"]]]]];
    resP=Flatten[Join[resP,resS,resK]];
    If[Not[SameQ[resP,{}]],
    If[SameQ[resP[[1]],0],resP=Take[resP,-Length[resP]+Length[Position[resP,0]\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
]],resP];
    If[resP!={},resP=Take[resP,{3,Length[resP],3}];
      Do[resP=ReplacePart[resP,{i,resP[[i]]},i],{i,1,Length[resP]}]],resP];
     
    resP
    ]
(*##################################################*)

(*############# fPrimeKL ##########################*)

(*fPrimeKL checks that KL given by pData,
  Dowker code or Conway symbol is a prime or
  composite*)
(*vraca 1 ako je Prime 0 za direktan *)

fPrimeKL[Ul_]:=Module[{pd,Dow,gr, i},
       If[SameQ[Head[Ul],String],pd=fCreatePData[Ul],pd=Ul];
    If[SameQ[Union[Map[EvenQ,Abs[pd[[2]]]]],{True}],pd=fPDataFromDow[pd]];
    (*sada sigurno imamo pdata *)
    If[SameQ[pd[[1]],{1,1}],gr=1,
    pd=ReductionKnotLink[pd];
    Dow=fDowfromPD[pd];
    gr=fGaussExt[Dow];
    gr=fGrInc[gr][[1]];
    If[SameQ[Flatten[Union[Map[Bridges,
              Map[FromUnorderedPairs,
                Union[Table[Drop[gr,{i}],{i,Length[gr]}]]]]]],{}],
           gr=1, gr=0]];
    gr]
    
    (*############# fPrimeGraph ##########################*)
    (*fPrimeKLGraph checks that KL given by graph is a 
      prime or composite KL= direct product*)
	(*vraca 1 ako je Prime 0 za direktan *)

	fPrimeGraph[GG_List]:=Module[{gr, i},
   	 If[SameQ[Flatten[Union[
            Map[Bridges,
            Map[FromUnorderedPairs,
            Union[Table[Drop[GG,{i}],{i,Length[GG]}]]]]]],{}],
             gr=1,gr=0];
    gr]
(*####################################################*)
(* ############# Murasugi Sig ########################## *)

(* fSigTor calculates Murasugi 
signature for torus KL - Murasugi, pp.148 *)
fSigTor[m_Integer,n_Integer]:=Module[{a,b,r,p,s},
    p=1;
    s=0;
    r={p,{m,n},s};
    a=Max[m,n];
    b=Min[m,n];
    While[Not[SameQ[r[[2]],{0,0}]],
      p=r[[1]];
      a=Max[r[[2,1]],r[[2,2]]];
      b=Min[r[[2,1]],r[[2,2]]];
      s=r[[3]];
      r={p,{a,b},s};
      r=If[SameQ[Min[a,b],1],{p,{0,0},s},r];
      r=If[SameQ[Min[a,b],2],{p,{0,0},s+p*(Max[a,b]-1)},r];
      r=If[
          SameQ[Max[a,b],2*Min[a,b]]&&Not[SameQ[Min[a,b],1]]&&
            Not[SameQ[Min[a,b],2]],{{p},{0,0},s+p*(Min[a,b]^2-1)},r];
      r=If[
          Max[a,b]>2*Min[a,b]&&Not[SameQ[Max[a,b],2*Min[a,b]]]&&
            Not[SameQ[Min[a,b],1]]&&Not[SameQ[Min[a,b],2]]&&
            SameQ[Mod[Min[a,b],2],1],{p,Abs[{Max[a,b]-2*Min[a,b],Min[a,b]}],
            s+p*(Min[a,b]^2-1)},r];
      r=If[
          Max[a,b]>2*Min[a,b]&&Not[SameQ[Max[a,b],2*Min[a,b]]]&&
            Not[SameQ[Min[a,b],1]]&&Not[SameQ[Min[a,b],2]]&&
            SameQ[Mod[Min[a,b],2],0],{p,Abs[{Max[a,b]-2*Min[a,b],Min[a,b]}],
            s+p*Min[a,b]^2},r];
      r=If[
          Max[a,b]<2*Min[a,b]&&Not[SameQ[Max[a,b],2*Min[a,b]]]&&
            Not[SameQ[Min[a,b],1]]&&Not[SameQ[Min[a,b],2]]&&
            SameQ[Mod[Min[a,b],2],1],{-p,Abs[{Max[a,b]-2*Min[a,b],Min[a,b]}],
            s+p*(Min[a,b]^2-1)},r];
      r=If[
          Max[a,b]<2*Min[a,b]&&Not[SameQ[Max[a,b],2*Min[a,b]]]&&
            Not[SameQ[Min[a,b],1]]&&Not[SameQ[Min[a,b],2]]&&
            SameQ[Mod[Min[a,b],2],0],{-p,Abs[{Max[a,b]-2*Min[a,b],Min[a,b]}],
            s+p*(Min[a,b]^2-2)},r]];
     r[[3]]
    ]
(*#################################################################*)
(*######################*)

(*############# fPDataFromDow ##########################*)

(*Calculates pdata from Dow without signs for alternating KL, 
or from Dow with signs of changed crossings for nonalternating  KL *)

fPDataFromDow[Ul_List]:=Module[{dd,res,sc,s,tt, i},
    dd=fSignsKL[Ul];
    sc=Map[Sign,fSignsKL[Abs[Ul]]][[2]]*Map[Sign,dd][[2]];
    res=If[Not[MemberQ[sc,-1]],dd,
        s=Map[Sign,dd][[2]];
        sc=Map[Sign,Ul][[2]];
        tt=Table[{{2i-1,Abs[Ul[[2,i]]]},{s[[i]]}},{i,Length[Ul[[2]]]}];
        tt=
          Sort[Table[
              If[SameQ[sc[[i]],-1],{Reverse[{tt[[i,1,1]],tt[[i,1,2]]}],
                  tt[[i,2]]},tt[[i]]],{i,Length[tt]}]];
        tt={Ul[[1]],Flatten[Table[tt[[i,1,2]]*tt[[i,2]],{i,Length[tt]}]]}];
    res={res[[1]],-res[[2]]};
    res] (*NOVA VERZIJA 27.08.2004 *)
    
    
    (*############# fPDataFromDowker ##########################*)

(*Calculates pdata from Dowker code with signs*)

fPDataFromDowker[Dow_List] := Module[{dd1, dd2, cs, i},
    dd1 = Map[Sign, Dow];(*znaci cvora*)
    dd2 = Map[Sign, fSignsKL[Abs[Dow]]];(*znaci alternirajuceg*)
    cs = dd1*dd2;(*mesta gde se razlikuju*)
    dd2 = Table[{2i - 1, Abs[Dow[[2, i]]]}, {i, Length[Dow[[2]]]}];
    dd2 = 
      Table[If[cs[[2, i]] == 1, dd2[[i]], Reverse[dd2[[i]]]], {i, 
          Length[dd2]}];
    dd2 = 
      Sort[Table[{dd2[[i, 1]], dd1[[2, i]]*dd2[[i, 2]]}, {i, Length[dd2]}]];
    dd2 = {Dow[[1]], Table[dd2[[i, 2]], {i, Length[dd2]}]};
    dd2]


(*############# fComponentNo ##########################*)

(*fComponentNo for a knot given by Conway, Dow, 
  or pdata gives the number of components *)
  
fComponentNo[Ul_]:=Module[{pd,l},
    pd=If[SameQ[Head[Ul],String],Length[fCreatePData[Ul][[1]]],
        Length[Ul[[1]]]]
    ]

(*############# RATIONAL KL ##########################*)

(* compo pravi particije broja n *)
compo[0] := {{}}
compo[n_Integer] := compo[n] = Module[{i},
		Flatten[Table[(Join[{i}, #1] &) /@ compo[n - i], {i, n}], 1]]

(* vrsi racionalnu redukciju nealternirajuceg racionalnog *)

(* RationalKL computes number and list of Conway symbols
of KL with n crossings *)

RationalKL[n_Integer]:=Module[{pa,t,u,i},
    pa=compo[n];
    t=Union[
        Table[If[First[pa[[i]]]==1||Last[pa[[i]]]==1,{},
            pa[[i]]],{i,Length[pa]}]];
    t=If[t[[1]]=={},Drop[t,1],t];
    t=Union[
        Table[If[MemberQ[t,Reverse[t[[i]]]],
            Last[Sort[{t[[i]],Reverse[t[[i]]]}]],{}],{i,Length[t]}]];
    t=If[t[[1]]=={},Drop[t,1],t];
    t=StringReplace[
        Map[ToString,t],{"{"->"","}"->"",","->""}];
    Print[Length[t]];
    t]

(*############# RationalAmphiK ##########################*)

(*RationalAmphiK calculates number and list of 
Conway symbols of rational amphicheiral knots *)

RationalAmphiK[n_Integer]:=
  Module[{pa,t,u,i,pom},If[SameQ[EvenQ[n],True],pa=compo[n/2];
      t=Union[Table[If[First[pa[[i]]]==1,{},pa[[i]]],{i,Length[pa]}]];
      t=If[First[t]=={},Drop[t,1],t];
      t=Table[Join[t[[i]],Reverse[t[[i]]]],{i,Length[t]}];
      t=StringReplace[
          Map[ToString,t],{"{"->"","}"->"",","->""}];
      pom=Map[Dowker,t];
      pom=Table[Head[pom[[i]][[4]][[1]]],{i,Length[pom]}];
      pom=
        Sort[Union[
            Table[If[SameQ[pom[[i]],Integer]==True,t[[i]],{}],{i,
                Length[pom]}]]];
      t=Reverse[If[SameQ[pom[[-1]],{}]==True,Drop[pom,-1],pom]];
      Print[Length[t]];
      t,0]]
 (*###################################*)
 
 (*############# RationalAmphiK ##########################*)

(*RationalAmphiL calculates number and list of 
Conway symbols of rational amphicheiral links *)
 
 RationalAmphiL[n_Integer]:=
  Module[{pa,t,t1,u,i,pom},If[SameQ[EvenQ[n],True],pa=compo[n/2];
      t=Union[Table[If[First[pa[[i]]]\[Equal]1,{},pa[[i]]],{i,Length[pa]}]];
      t=If[First[t]\[Equal]{},Drop[t,1],t];
      t=Table[Join[t[[i]],Reverse[t[[i]]]],{i,Length[t]}];
      t=StringReplace[
          Map[ToString,t],{"{"\[Rule]"","}"\[Rule]"",","\[Rule]""}];
      pom=Map[Dowker,t];
      pom=Table[Head[pom[[i]][[4]][[1]]],{i,Length[pom]}];
      pom=
        Sort[Union[
            Table[If[SameQ[pom[[i]],Integer]\[Equal]True,t[[i]],{}],{i,
                Length[pom]}]]];
      t=Reverse[If[SameQ[pom[[-1]],{}]\[Equal]True,Drop[pom,-1],pom]];
      t1=Union[Table[If[First[pa[[i]]]\[Equal]1,{},pa[[i]]],{i,Length[pa]}]];
      t1=If[First[t1]\[Equal]{},Drop[t1,1],t1];
      t1=Table[Join[t1[[i]],Reverse[t1[[i]]]],{i,Length[t1]}];
      t1=StringReplace[
          Map[ToString,t1],{"{"\[Rule]"","}"\[Rule]"",","\[Rule]""}];
      t=Complement[t1,t];
      Print[Length[t]];
      t,0]]
 
 (*########### RacGenKL #############*)      

(* restrictedPartitions[m_Integer, r_Integer, 
      n_Integer] pravi ogranicene particije broja n, 
  na najvise r delova od kojih je svaki najvise m *)
  
partH[_, 0] := {{}}
partH[H_, n_] := {} /; First[H] > n
partH[H_, n_] := Module[{h, p, ans}, ans = Map[(h = #;
            p = partH[Select[H, # >= h &], n - h];
            Map[Join[{h}, #] &, p]) &, H];
    Flatten[ans, 1]]
partit[n_] := partH[Range[n], n]
partitionsOdd[n_] := partH[Range[1, n, 2], n]
partitionsDif[n_] := DeleteCases[partit[n], {___, x_, x_, ___}]

restrictedPartitions[m_Integer, r_Integer, n_Integer] := 
  Select[partH[Range[m], n], (Length[#] <= r) &]

(*########################*)

(* RacGenKL[n_Integer, m_Integer] su racionalni generating  za m = 
    3 i sourcelinks KL za m = 
      2 : racionalni generisuci KL ciji se zapis sastoji samo od 1, 2, 3, 
  koji direktno daju familije i racionalni source links m = 2 *) 

RatGenSourKL[n_Integer, m_Integer] := Module[{pa, t, u, i},
    pa = Flatten[Map[DistinctPermutations, 
    restrictedPartitions[m, n - 1, n]],1];
    t = Union[
        Table[If[First[pa[[i]]] == 1 || 
        Last[pa[[i]]] == 1, {}, pa[[i]]], {i, Length[pa]}]];
    t = Union[
        Table[If[MemberQ[t, Reverse[t[[i]]]], 
            Last[Sort[{t[[i]], Reverse[t[[i]]]}]], {}], {i, Length[t]}]];
    t = If[First[t] == {}, Drop[t, 1], t];
    Print[Length[t]];
    t = StringReplace[Map[ToString, t], {"{" -> "", "}" -> "", "," -> ""}] 
    ]
    
    (* RatSourceKLNo calculates the number of 
    rational source KL for k<=n crossings*)
RatSourceKLNo[n_Integer]:=Module[{b={1},i,p},
    If[n<4,1,
           If[n==4||n==5,
                  b={1},
                  b={1,1};
                 Do[
                       p=b[[i-1]]+b[[i-2]];
          	If[EvenQ[i],p=p-Fibonacci[(i-2)/2] ];
                      b=Append[b,p]
           ,{i,3,n-3}]
            ]
       ];
    Last[b]
    ]
 
(*########################*)

(* RatLinkU1 calculates number and Conway symbols of rational 
links with unlinking number 1 according to P.Kohn*)

RatLinkU1[n_Integer]:=Module[{pa,t,u,i},
    If[SameQ[OddQ[n],True],
      pa=compo[(n+1)/2];
      t=Union[
          Table[If[Last[pa[[i]]]!=1,{},pa[[i]]],{i,Length[pa]}]];
      t=If[First[t]=={},Drop[t,1],t];
      u=Table[Reverse[ReplacePart[t[[i]],t[[i]][[-2]]-1,-2]],{i,Length[t]}];
      t=Union[
          Table[If[MemberQ[Join[t[[i]],u[[i]]],0]==False,
              Join[t[[i]],u[[i]]],{}],{i,Length[t]}]];
      t=If[First[t]=={},Drop[t,1],t];
      t=Union[
          Table[If[First[t[[i]]]==1||Last[t[[i]]]==1,{},
              t[[i]]],{i,Length[t]}]];
      t=If[First[t]=={},Drop[t,1],t];
      t=StringReplace[
          Map[ToString,t],{"{"->"","}"->"",","->""}];
      Print[Length[t]];
      If[n==5,{"2 1 2"},t],0]]

(*########################*)

(*RatLinkU0[n_Integer]:n-odd =rational unlinks U=0*)

RatLinkU0[n_Integer]:=Module[{pa,t,u,i},
    If[SameQ[OddQ[n],True],
      pa=compo[(n+1)/2];
      t=Union[
          Table[If[Last[pa[[i]]]!=1,{},pa[[i]]],{i,Length[pa]}]];
      t=If[First[t]=={},Drop[t,1],t];
      u=Table[Reverse[ReplacePart[t[[i]],t[[i]][[-2]]-1,-2]],{i,Length[t]}];
      t=Union[
          Table[If[MemberQ[Join[t[[i]],u[[i]]],0]==False,
              Join[ReplacePart[t[[i]],-1,-1],u[[i]]],{}],{i,Length[t]}]];
      t=If[First[t]=={},Drop[t,1],t];
      t=Union[
          Table[If[First[t[[i]]]==1||Last[t[[i]]]==1,{},
              t[[i]]],{i,Length[t]}]];
      t=If[First[t]=={},Drop[t,1],t];
      t=StringReplace[
          Map[ToString,t],{"{"->"","}"->"",","->""}];
      Print[Length[t]];
      t,0]]

(*########################*)

(*MSigRat calculates Murasugi signature of a rational 
KL given by Conway symbol *)

MSigRat[Conway_String] := 
  Module[{Con, q, p, l, i}, 
    Con = ReadList[
        StringToStream[
          StringReplace[
            Conway, {"," -> " ", "+" -> " ", "(" -> " ", ")" -> " "}]], 
        Number];
    Con = FromContinuedFraction[Con];
    Con = ContinuedFraction[Con];
    (* Con = Last[Sort[{Con, Reverse[Con]}]]; *)
    Con1 = FromContinuedFraction[Con];
    Print[Con1];
    q = Denominator[Con1];
    p = Numerator[Con1];
    l = Table[(i - 1)*q, {i, p}];
    l = Table[
        If[Mod[l[[i]], 2p] > p, Mod[l[[i]], 2p] - 2p, Mod[l[[i]], 2p] ], {i, 
          Length[l]}];
    sig = Apply[Plus, Map[Sign, l]];
    sig
    ]
    
    (*########################*)
    
   (*RatKnotGenU1 generates rational knots with U=1 *)
   
fGenKnotU1[n_Integer]:=Module[{t,u,i},
    pa=compo[(n+1)/2];
    t=Union[Table[If[Last[pa[[i]]]!=1,{},pa[[i]]],{i,Length[pa]}]];
    t=If[First[t]=={},Drop[t,1],t];
    u=Table[Reverse[ReplacePart[t[[i]],t[[i]][[-2]]-1,-2]],{i,Length[t]}];
    t=Union[
        Table[If[MemberQ[Join[t[[i]],u[[i]]],0]==False,
            Join[t[[i]],u[[i]]],{}],{i,Length[t]}]];
    t=If[First[t]=={},Drop[t,1],t]
    ]


RatKnotU1[k_Integer,UL_List]:=Module[{pa,t,u,i},UL;
    t=Union[UL,Reverse/@UL];
    t=Union[
        Table[If[k==1,ReplacePart[t[[i]],t[[i]][[1]]+1,1],
            Join[{k},t[[i]]]],{i,Length[t]}]];
    t=Table[
        If[t[[i]][[-1]]==1,
          Drop[ReplacePart[t[[i]],t[[i]][[-2]]+1,-2],-1],t[[i]]],{i,
          Length[t]}];
    t=StringReplace[
        Map[ToString,t],{"{"->"","}"->"",","->""}]]

RatKnotGenU1[n_Integer]:=Module[{i=5,res={}},
    If[n<3,Print["To small No !"],
         If[n==3,Print["{3}"],
          While[i<=2IntegerPart[n/2],
                        
            res=Append[res,Map[RatKnotU1[n-i,{#}]&,fGenKnotU1[i]]];
                      i=i+2];
           res=Prepend[res,{{ToString[n-2]<>" 2"}}];
                    ]];
                    Flatten[res]
    ]
    
    (*########################*)
    
    (*RatKnotGenU0[n_Integer] generates rational 
    unknots with U=0. *)

fGenKnotU0[n_Integer]:=Module[{pa,t,u,i},pa=compo[(n+1)/2];
    t=Union[Table[If[Last[pa[[i]]]!=1,{},pa[[i]]],{i,Length[pa]}]];
    t=If[First[t]=={},Drop[t,1],t];
    u=Table[Reverse[ReplacePart[t[[i]],t[[i]][[-2]]-1,-2]],{i,Length[t]}];
    t=Union[
        Table[If[MemberQ[Join[t[[i]],u[[i]]],0]==False,
            Join[ReplacePart[t[[i]],-1,-1],u[[i]]],{}],{i,Length[t]}]];
    t=If[First[t]=={},Drop[t,1],t]
   ]


RatKnotU0[k_Integer,UL_List]:=Module[{pa,t,u,i},UL;
    t=Union[UL,Reverse/@UL];
    t=Union[
        Table[If[k==1,ReplacePart[t[[i]],t[[i]][[1]]+1,1],
            Join[{k},t[[i]]]],{i,Length[t]}]];
    t=Table[
        If[t[[i]][[-1]]==1,
          Drop[ReplacePart[t[[i]],t[[i]][[-2]]+1,-2],-1],t[[i]]],{i,
          Length[t]}];
    t=StringReplace[
        Map[ToString,t],{"{"->"","}"->"",","->""}]]

RatKnotGenU0[n_Integer]:=Module[{i=3,res={}},
    If[n<5,Print["To small n!"],
        While[i<=2IntegerPart[n/2],
                      res=Append[res,Map[RatKnotU0[n-i,{#}]&,fGenKnotU0[i]]];
                    i=i+2];
                ];
               Flatten[res]
    ]
(*#######################################################*)
(*#######################################################*)
(*########## ALL STATES RATIONAL ###############*)
(*Pomocna za fChange ConVar, AllStatesRat *)
(*Radi od Konvej Liste i \
odgovarajuce varijacije- 
    elemment varijacije je broj koliko puta treba izvrsiti promenu u toj \
tacki*)
fReduCon[Con_List,VVvar_List]:=
  Module[{p,i,l=Con,rr},Do[p=Con[[i]]-2*Sign[Con[[i]]]*VVvar[[i]];
      l=ReplacePart[l,p,i],{i,1,Length[VVvar]}];
    IsNotKnot[l]]
(*menja vrednosti u Konvej listi prema varijaciji*)

fChangeConVar[Con_List,vvSum_List]:=Module[{l=Con},
         l=fReduCon[l,vvSum];
         l=StringReplace[ToString[l],"{}"->"{1}"];
         l=StringReplace[l,{"{"->"",","->"","}"->"", "0"->"1"}];
        {Apply[Plus,vvSum],"Con: "<>l}
    ]
(*Ulaz je Konvej ili PData *)
(*Radi sve pojedinacne promene znaka na fix projekciji*)

AllStatesRational[Conway_String]:=
  Module[{Con1,Con,ConList,vv,i,res={},res1},
    ConList=ReadList[
        StringToStream[
          StringReplace[
            Conway,{","->" ","+"->" ","("->" ",
              ")"->" "}]],Number];
    Con1=FromContinuedFraction[ConList];
    If[SameQ[Con1,ComplexInfinity],
      Con=FromContinuedFraction[Reverse[ConList]],Con=Con1];
    Con=ContinuedFraction[Con];
    (*nasli smo ga ko je zaista*)(*sad ga razvezujemo*)
    
    vv=fVarP[Apply[Plus,Map[Abs[#]&,Con]]];
    vv=Map[iteratedTake[#,Abs[Con]]&,vv];
    Do[vv=ReplacePart[vv,Map[Apply[Plus,#]&,vv[[i]]],i],{i,1,Length[vv]}];
    vv=Union[vv];
    (*sad je svaka varijacija podeljena prema duzinama clanova konveja*)
    
    vv=Union[Map[fChangeConVar[Con,#]&,vv]];
    Con=Sort[Map[Reverse[#]&,vv]];
    vv=Union[Map[Rest[#][[1]]&,vv]];
    Do[res=Append[res,Reverse[First[Select[Con,SameQ[#[[1]],vv[[i]]]&]]]],{i,
        1,Length[vv]}];
    res=Sort[res];
   res1=Table[{res[[i,1]],StringReplace[res[[i,2]],"-"->""]},
   {i,Length[res]}];
    {Select[res1,#[[2]]=="Con: 1"&][[1,1]],res}]

(* ################# AllStatesRatFast ###################*)

AllStatesRatFast[Conway_String]:=
  Module[{Con1,Con,ConList,vv,i,res={},res1},
    ConList=ReadList[
        StringToStream[
          StringReplace[
            Conway,{","->" ","+"->" ","("->" ",
              ")"->" "}]],Number];
    Con1=FromContinuedFraction[ConList];
    If[SameQ[Con1,ComplexInfinity],
      Con=FromContinuedFraction[Reverse[ConList]],Con=Con1];
    Con=ContinuedFraction[Con];
    (*nasli smo ga ko je zaista*)(*sad ga razvezujemo*)
    
    vv=fVarNewP[Apply[Plus,Map[Abs[#]&,Con]]];
    vv=Map[iteratedTake[#,Abs[Con]]&,vv];
    Do[vv=ReplacePart[vv,Map[Apply[Plus,#]&,vv[[i]]],i],{i,1,Length[vv]}];
    vv=Union[vv];
    (*sad je svaka varijacija podeljena prema duzinama clanova konveja*)
    
    vv=Union[Map[fChangeConVar[Con,#]&,vv]];
    Con=Sort[Map[Reverse[#]&,vv]];
    vv=Union[Map[Rest[#][[1]]&,vv]];
    Do[res=Append[res,Reverse[First[Select[Con,SameQ[#[[1]],vv[[i]]]&]]]],{i,
        1,Length[vv]}];
    res=Sort[res];
     res1=Table[{res[[i,1]],StringReplace[res[[i,2]],"-"->""]},
     {i,Length[res]}
];
    {Select[res1,#[[2]]=="Con: 1"&][[1,1]],res}]


(*#################################### 25.09.2003*)



fVarNewPGap[n_Integer,u_Integer]:=
  Module[{r,i},r=Drop[Table[IntegerDigits[i,2,n],{i,0,2^n-1}],1];
    r=Sort[Table[{Count[r[[i]],1],r[[i]]},{i,Length[r]}]];
    r=Drop[
        Union[Table[
            If[r[[i,1]]>IntegerPart[(n+1)/2],{},r[[i,2]]],{i,Length[r]}]],1];
    r=Sort[Table[{Count[r[[i]],1],r[[i]]},{i,Length[r]}]];
    r=Select[r,#[[1]]==u&];
    r=Table[r[[i,2]],{i,Length[r]}];
    r]

(*####################################*)
(*############### ALLSTATES ########################*)
(*menja znake u pdata prema varijaciji*)
fCrChVarPD[PD_List,var_List]:=
  Module[{l,n,i,pos},n=Length[PD[[2]]];
    l=Complement[Range[2n],Abs[PD[[2]]]];
    (*svi brojevi do 2 duzine koji vec nisu u pData*)Do[
      l=ReplacePart[l,{l[[i]]*Sign[PD[[2,i]]],PD[[2,i]]},i],{i,1,n}];
    pos=Flatten[Position[var,1]];
    Do[l=ReplacePart[l,-1*Reverse[l[[pos[[i]]]]],pos[[i]]],{i,1,
        Length[pos]}];
    l=Sort[l,Abs[#1[[1]]]<Abs[#2[[1]]]&];
    l=Flatten[Map[Take[#,-1]&,l]];
    l={PD[[1]],l};
    l=ReductionKnotLink[l];
    (*{Count[var,1],var,l}*)
    
    If[MemberQ[{{{},{}},{{0},{}},{{},{0}}},l],l="pd: {{},{}}",
      l="pd: "<>ToString[l]];
    {Count[var,1],l}]


(*Ulaz je Konvej ili PData*)
(*Radi sve pojedinacne promene znaka na fix \
projekciji*)
fAllStatesProj[UL_]:=Module[{pd=UL,i,prvi,vv,res={}},
    If[SameQ[Head[UL],String],pd=fCreatePData[UL]];
    (*sad su pd pdata*)
    pd=ReductionKnotLink[pd];
    vv=Rest[fVarP[Length[pd[[2]]]]];
    pd=Union[Map[fCrChVarPD[pd,#]&,vv ]];
    prvi=Union[Map[Rest[#][[1]]&,pd]];
    pd=Sort[Map[Reverse[#]&,pd]];
    (*vv=Union[Map[Apply[Plus,#]&,vv]];Print[vv];*)
    Do[res=Append[res,
          Reverse[First[Select[pd,#[[1]]==prvi[[i]]&]]]
          ],
      
      {i,1,Length[prvi]}];
    res=Reverse[Sort[res]];
    {Select[res,#[[2]]=="pd: {{},{}}"&][[1,1]],Sort[res]}]
(*################################################*)


(*menja znake u pdata prema varijaciji*)
fSignVarPD[PD_List,var_List]:=
  Module[{l,n,i,pos},n=Length[PD[[2]]];
    l=Complement[Range[2n],Abs[PD[[2]]]];
    (*svi brojevi do 2 duzine koji vec nisu u pData*)Do[
      l=ReplacePart[l,{l[[i]]*Sign[PD[[2,i]]],PD[[2,i]]},i],{i,1,n}];
    pos=Flatten[Position[var,1]];
    Do[l=ReplacePart[l,-1*Reverse[l[[pos[[i]]]]],pos[[i]]],{i,1,
        Length[pos]}];
    l=Sort[l,Abs[#1[[1]]]<Abs[#2[[1]]]&];
    l=Flatten[Map[Take[#,-1]&,l]];
    l={PD[[1]],l};
    l=ReductionKnotLink[l];
    (*{Count[var,1],var,l}*)
    If[MemberQ[{{{},{}},{{0},{}},{{},{0}}},l],l=1,
      l=0]]

(*################# PERIOD KL ####################*)

PeriodProjAltKL[Ul_]:=Module[{G,sG,Aut,k,Aut1, i,j},If[SameQ[Ul,"1"],Aut={},
      G=fGaussExtSigns[Ul];
      sG=fGrInc[G][[2]];
      G=FromUnorderedPairs[fGrInc[G][[1]]];
      Aut=Drop[Automorphisms[G],1];
      Aut1=Aut;
      Aut=
        If[Aut!={},
          If[First[
                Union[Table[
                    If[SameQ[
                          Sign[Table[
                                Table[
                                  Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]],{j,
                                    Length[sG]}],{i,Length[Aut]}]][[k]],
                          sG]==True,Aut[[k]],{}],{k,
                      Length[Aut]}]]]=={},
            Drop[Union[
                Table[If[
                    SameQ[Sign[
                            Table[Table[
                                Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]],{j,
                                  Length[sG]}],{i,Length[Aut]}]][[k]],
                        sG]==True,Aut[[k]],{}],{k,Length[Aut]}]],1],
            Union[Table[
                If[SameQ[
                      Sign[Table[
                            Table[Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]],{j,
                                Length[sG]}],{i,Length[Aut]}]][[k]],
                      sG]==True,Aut[[k]],{}],{k,Length[Aut]}]]],Aut];
      (*trazimo cikle*)
      Aut=If[Aut!={},Map[ToCycles,Aut],Aut];
      Aut=Select[Aut,Count[Map[Length,#],1]<3&];
      Aut=Select[Aut,Length[Complement[Union[Map[Length,#]],{1}]]==1&];
      (*iz nekog razloga duzine cikla su periodi*)
      
      Aut=Union[Map[Length,Flatten[Aut,1]]];
      (*ako nije prazna izbaci prvu 1-identitet*)
      
      Aut=If[Aut!={}&&First[Aut]==1,Drop[Aut,1],Aut]];
    Aut=If[SameQ[Aut,{}],{0},Aut];
    Aut=If[Length[Aut1]==1,{2},Aut];
    Aut
    ]

(*###############################################################*)
(*Radi od Konveja, dowkera i pdata *)
PeriodAltKL[Conway_String]:=Module[{res,razne},If[SameQ[Conway,"1"],
       res={0},
      razne=fDiffProjectionsAltKL[Conway];
      razne=Map[#[[1]]&,razne];
      res=Union[Flatten[Map[PeriodProjAltKL[#]&,razne]]];
      If[SameQ[res,{}],res={0}]];
      If[SameQ[res[[1]],0],res=Drop[res,1],res];
    Print["Period: ",res]
    ]
    (*popravljano 16.8.2003 *)
(*###############################################################*)
 (*############### MINDOWALT ######################## *) 
 
 (*############### MINDOWALT ########################*)(*fOrderComp:
    Input Con_String priprema link za pletenje.
    Kao izlaz daje Gauss code kod 
koga su sve komponente sa samopresecima 
dovedene na optimalan pocetak.Drugi
clan izlaza je Ulc kod koji pokazuje 
koje su komponente neodlucive Ne radi za
"2"*)(*fOrderComp:
    Input Con_String priprema link za pletenje.
    Kao izlaz daje Gauss code kod
koga su sve komponente sa samopresecima 
dovedene na optimalan pocetak.Drugi
clan izlaza je Ulc kod koji pokazuje 
koje su komponente neodlucive Ne radi za
"2"*)

fOrderComp[Con_]:=Module[{c,p,EDL,fmc,tc, i, \
j},EDL=Sort[fGaussExtSigns[Con]];    
    c=Table[Length[EDL[[i]]],{i,Length[EDL]}];    
    p=Table[
        Table[Length[Intersection[Abs[EDL[[j]]],Abs[EDL[[i]]]]],{i,
            Length[EDL]}],{j,Length[EDL]}];
    p=Map[Sort,
        Table[Table[If[p[[j,i]]>0,c[[i]],c[[j]]],{i,Length[p[[j]]]}],{j,
            Length[p]}]];(*c[[j]] ili 2*c[[j]]*)
    fmc=Map[fMinComp,EDL];
   fmc=Sort[Table[{c[[i]],p[[i]],fmc[[i]][[1]],EDL[[i]]},{i,Length[EDL]}]];
    tc=fmc;
    fmc=Table[{fmc[[i]][[1]],fmc[[i]][[2]],Sort[fmc[[i]][[3,1]]]},{i,
          Length[fmc]}];
    fmc=Split[fmc];
    fmc=Table[Length[fmc[[i]]],{i,Length[fmc]}];
    Ulc=Range[fmc]+Prepend[Table[Sum[fmc[[i]],{i,i}],{i,Length[fmc]-1}],0];
    (*ne valja rotacija*)tc=
      Table[If[SameQ[tc[[i]][[3,3]],0],
          RotateLeft[tc[[i]][[4]],tc[[i]][[3,2]]-1],
          Reverse[RotateRight[tc[[i]][[4]],
              Length[tc[[i]][[4]]]-tc[[i]][[3,2]]]]],{i,Length[tc]}];
    {tc,Ulc}]


(*fEqChoices[
      Ulc_List] pomocna funkcija za permutovanje 
      ekvivalentnih komponenata.Iz
svakog podskupa ekvivalentnih 
elemenata bira reprezentanta i gradi sve
permutacije*)

fEqChoices[Ulc_List]:=Module[{dp,ec, i},dp=Map[DistinctPermutations,Ulc];
    ec=KSubsets[Flatten[dp,1],Length[dp]];
    ec=Union[
        Table[If[SameQ[Union[Flatten[ec[[i]]]],Flatten[Ulc]],ec[[i]],{}],{i,
            Length[ec]}]];
    If[SameQ[First[ec],{}],Drop[ec,1],ec]]
fPletiLink[sara_List,poc_List]:=Module[{res},res=Map[poc[[#]]&,sara]]

(*############ VARIJANTA A ####################*)

(*Sredjuje jednu komponentu u zavisnosti 
od toga da li ima dvojnih elemenata
ili nema:ind {0,1} iz fFixC*)

fForFixC[LL_List,OstL_List,ind_Integer]:=Module[{res={},pom,pom1, i},  
    If[ind==0,   
      (*Print["Indikator", ind];*)
            Do[If[OddQ[i],
                      If[SameQ[Position[OstL,LL[[i]]],{}],
                                            res=Append[res,{0,LL[[i]]}],     
                                             
            res=Append[res,{Position[OstL,LL[[i]]][[1,1]],LL[[i]]}]], 
                      
          res=Append[res,{Position[OstL,LL[[i]]][[1,1]],LL[[i]]}]]
                ,{i,1,Length[LL]}](*;
        Print["ind je 0 a ovo je rez ",res]*)
      ,
           Do[If[SameQ[Position[OstL,LL[[i]]],{}],
                         res=Append[res,{0,LL[[i]]}],
                          
          res=Append[res,{Position[OstL,LL[[i]]][[1,1]],LL[[i]]}]],
        {i,1,Length[LL]}]
      ];
    
    pom=Min[Map[Take[#,1]&,res]];
    poc=Select[res,SameQ[#[[1]],pom]&];
    poc=Union[Flatten[Map[Take[#,-1]&,poc],1]];
    (*Print["pocetni ",poc];*)
    pom=Map[Position[LL,#][[1,1]]&,poc];
    (*Print["pos, rotacije ", pom];*)
    If[ind==0,
      res=
        Map[{RotateLeft[res,#-1],Reverse[RotateRight[res,Length[res]-#]]}&,
          pom];
      res=Map[Flatten[#]&,Flatten[res,1]];
      (*Print["pre uzimanja :", Length[res],res];*)
      
      res=Map[Take[#,{2,Length[#],2}]&,res]
      (*Print["Uzeli smo sve ",res]*),
      res=fMin2Dvojne[LL]
      ];
    If[Length[res]==2,
          If[Count[pom[[1]],0]==2, 
                pom=Map[RotateLeft[#,Position[#,poc[[1]]][[2]]]&,res];
                pom=Map[Reverse[#]&,pom];
                res=Union[Flatten[Map[Append[res,#]&,pom],1]]
                 ],
        If[Count[pom[[1]],0]==2, 
              pom=RotateLeft[res[[1]],Position[res[[1]],poc[[1]]][[2]]];
              res=Append[res,Reverse[pom]]
           ]];
    res=Map[Prepend[OstL,#]&,res];
    
    res]

(*######### sredjuje jednu komponentu A ############*)

(*lista sa ExtDowkerima u 
fixiranim pozicijama tj.komponente su zauzele svoja
mesta*)
(*nn oznacava broj komponente sa kojom radimo*)

fMin2Dvojne[komp_List]:=Module[{poc,pom,pos,
posraz,minrast,i,res={},n,dd,dd1,
dd2,dd3},
    n=Length[komp];
    poc=Select[Union[komp],Count[komp,#]==2&];
   
    pos=Map[{#,Flatten[Position[komp,#]],
    Flatten[Position[Reverse[komp],#]]}&,
        poc];
    posraz=Map[{#[[1]],#[[2,1]],#[[2,2]]-#[[2,1]],
                    #[[3,1]],#[[3,2]]-#[[3,1]]}&,pos];
    posraz=Map[Insert[#,n-#[[3]],4]&,posraz];
    posraz=Map[Append[#,n-#[[6]]]&,posraz];
   (* Print["rastojanja:",posraz];*)
    
    minrast=Min[Flatten[Map[{#[[3]],#[[4]],#[[6]],#[[7]]}&,posraz]]];
    (*Print["Min ",minrast];*)
    Do[ dd=Position[komp,poc[[i]]][[2,1]];
    dd1=Position[Reverse[komp],poc[[i]]][[2,1]];
    dd2=Position[komp,poc[[i]]][[1,1]];
    dd3=Position[Reverse[komp],poc[[i]]][[1,1]];
     If[posraz[[i,3]]==minrast,
                    res=Append[res,RotateLeft[komp,dd2-1]]];
            If[posraz[[i,4]]==minrast,res=Append[res,RotateLeft[komp,dd-1]]];
            If[posraz[[i,6]]==minrast,
                  res=Append[res,RotateLeft[Reverse[komp],dd3-1]]];
             If[posraz[[i,7]]==minrast,
                   res=Append[res,RotateLeft[Reverse[komp],dd1-1]]]
      ,{i,1,Length[poc]}];
       (* Print[res];*)
    res
    ]

fFixC[tcL_List,nn_Integer]:=Module[{p,pre,ost,pom,res},
    p=tcL[[nn]];
    pom=Flatten[Map[Count[p,#]&,Union[p]]];
    pom=Count[pom,2];
    (*broj dvojnih *)
    (*Print["Komponenta ",nn," Broj Dvojnih ",pom];*)
  
      If[pom>=2,
      (*ako ima vise od dva dvojna onda je fiksirana *)
      
      pom=fMin2Dvojne[tcL[[nn]]];
      pre=Take[tcL,nn-1];
      ost=Drop[tcL,nn];
      If[pre=={},res=Map[{{#},ost}&,pom],
        If[ost=={},res=Map[{pre,{#}}&,pom],
          res=Map[{pre,{#},ost}&,pom]]];
      res=Map[Flatten[#,1]&,res],
      If[pom==1,
            res=fForFixC[p,Take[tcL,{nn+1,Length[tcL]}],1]];
      If[pom==0,
            res=fForFixC[p,Take[tcL,{nn+1,Length[tcL]}],0]];
      
         If[nn!=1,
               res=Map[Join[Take[tcL,{1,nn-1}],#]&,res]],
                res=fForFixC[p,Take[tcL,{nn+1,Length[tcL]}],0];
                If[nn!=1,res=Map[Join[Take[tcL,{1,nn-1}],#]&,res]]]
    ;
    res]
(*sredjuje sve komponente koje nemaju zajednickih elemenata sa svojim \
prethodnicima*)

fFixLAFirst[LinkL_List,indL_List]:=Module[{i,n,res,pre={}},n=Length[LinkL];
    Do[
         If[SameQ[i,1],(*1.komponenta*)
        res=fFixC[LinkL,1];
        pre=LinkL[[1]];
        ,(*ostale ako je 0 radimo inace ne*)
        If[SameQ[indL[[i]],0],
          res=Union[Flatten[Map[fFixC[#,i]&,res],1]];
          pre=Union[pre,Flatten[LinkL[[i]]]],
          pre=Union[pre,Flatten[LinkL[[i]]]]]],{i,1,n}];
    res]

(*############### VARIJANTA B ##############*)

(*ind iz FixLink AB 1*)
(*fSeePrevNext posmatra susede i uzima bolji*)
\
(*pomocna za fFixLinkB*)
(*nn=1 ima neparnih*)

fSeePrevNext[CompL_List,PosLL_List,nn_Integer]:=Module[{pom,res},
    Switch[nn,0,pom={Take[CompL,{2}][[1]],Last[CompL]},1,
      pom={Take[CompL,{1}][[1]],Take[CompL,{3}][[1]]}];
    pom=Flatten[{Union[Select[PosLL,SameQ[#[[2]],pom[[1]]]&]],
          Union[Select[PosLL,SameQ[#[[2]],pom[[2]]]&]]},1];
    pom={pom[[1,1]],pom[[2,1]]};
    If[nn==0,pom=Reverse[pom]];
    (*If[pom[[2]]>pom[[1]],
               res=CompL,
          Switch[nn,1,res=Reverse[RotateLeft[CompL,3]],0,
            res=Reverse[RotateLeft[CompL,1]]];
          If[SameQ[pom[[1]],pom[[2]]],res=Prepend[{res},CompL]]];*)
    
    res={CompL};
    Switch[nn,1,res=Append[res,Reverse[RotateLeft[CompL,3]]],
      0,res=Append[res,Reverse[RotateLeft[CompL,1]]]];
    res]


(*fFixCB sredjuje B komponentu*)
(*pomocna za fFixLinkB*)

fBFixC[CompL_List,PList_List,nn_Integer]:=Module[{pom,minel,minpos},
    If[SameQ[nn,1],pom=Select[PList,OddQ[#[[1]]]&],pom=PList];
    minel=First[Sort[pom]][[2]];
    minpos=Position[CompL,minel][[1]];
    Switch[nn,1,p=RotateLeft[CompL,minpos-2],0,p=RotateLeft[CompL,minpos-1]];
    p=fSeePrevNext[p,PList,nn];
    If[Head[p[[1]]]==Integer,p={p}];
    p]

fGlupa[LL_List,bb_Integer,n_Integer]:=Module[{pre,Li,PosL,PomLink},
    pre=Flatten[Take[LL,n-1]];
    PomLink=LL;
    (*ako je bb 1 onda je B i radimo je *)
    If[SameQ[bb,1],
      Li=PomLink[[n]];
      PosL=
        Map[If[SameQ[
                Position[pre,#],{}],{1000,#},{Position[pre,#][[1,1]],#}]&,
          Li];
      If[SameQ[Select[PosL,OddQ[#[[1]]]&],{}],
                         (*nema neparnih  0*)
                         
        pom=fBFixC[Li,PosL,0],
                         pom=fBFixC[Li,PosL,1]];
      PomLink=Map[Insert[Drop[PomLink,{n}],#,n]&,pom]];
    PomLink]

(*fFixLinkB sredjuje sve komponente koje su B u linku gde su A vec sredjene*)
\

fFixLBSecond[LLink_List,ind_List]:=Module[{i,n,res},
    res=LLink;n=Length[LLink];    
     Do[If[SameQ[Head[res[[1,1]]],List],
        res=Map[fGlupa[#,ind[[i]],i]&,res];
       res=Flatten[res,1];
        If[SameQ[ind[[i]],0],res=Partition[res,n]]];
      If[SameQ[Head[res[[1,1]]],Integer],res=fGlupa[res,ind[[i]],i]],{i,2,
        n}];
    res]

fLinkAB[LinkL_List]:=
  Module[{ind={},i,pre,res},
    Do[If[i==1,ind={0};pre=Union[LinkL[[1]]],
        If[SameQ[Intersection[pre,LinkL[[i]]],{}],(*A prethodne ne uticu-0*)
            ind=Append[ind,0],(*B prethodne uticu-1*)ind=Append[ind,1]];
        pre=Union[pre,LinkL[[i]]]],{i,1,Length[LinkL]}];
    res=fFixLAFirst[LinkL,ind];
   res=Flatten[Map[fFixLBSecond[#,ind]&,res],1]; 
   Print[res,"222"];
    res]

(*radi od Konveja,
  pdata ili dowkera sa znacima*)
(*!!!!Promena za cvorove u select*)

(*pomocna za MinDowProj*)

fPickEven[LL_List]:=Module[{res},res=Flatten[LL[[2]]];
    res=Union[Map[EvenQ[#]&,res]];
    SameQ[res,{True}]]
(****)

MinDowProjAltKL[Ul_]:=
  Module[{p,vrti,AllLinks,d,Con=Ul},
    If[SameQ[Ul,"2"]||SameQ[Ul,"-2"],p={{1,1},{4,2}},
      If[SameQ[Head[Ul],List]&&SameQ[Ul[[1]],{1,1}],p={{1,1},{4,2}},
        If[SameQ[Head[Con],List],
          If[SameQ[Union[Map[OddQ,Abs[Con[[2]]]]],{True}],(*input pdata*)Con=
              fDowfromPD[Con]]];
        (*sad imamo dowkera ili Konveja*)d=fGaussExtSigns[Con];
       If[SameQ[d,Flatten[d]],(*cvor*)AllLinks=fMinComp[d];
         AllLinks=Map[fKnittDow[d,Take[#,-2]]&,AllLinks];
         AllLinks=Map[fDowfromGaussExt[#[[1]]]&,AllLinks];
         p=Sort[Table[{Abs[AllLinks[[i]]],AllLinks[[i]]},
         {i,Length[AllLinks]}]][[1,2]],  
         (* KORIGOVANO 11.01.2004*)
          (*LINKOVI*)
          p=fOrderComp[Con];
          (*komponente sredjene po nekoliko kriterijuma,
            a drugi deo daje klase ekvivalencije*)
          
          vrti=Map[Flatten[#,1]&,fEqChoices[p[[2]]]];
          AllLinks=Map[fPletiLink[#,p[[1]]]&,vrti];
          
                 
          AllLinks=Flatten[Map[fLinkAB[#]&,AllLinks],1];
               
          Print[AllLinks,"Posle AB"];
          
          AllLinks=Map[fDowfromGaussExt[#]&,AllLinks];  
          
          Print[AllLinks,"Kraj"];
                               
          AllLinks=Union[Select[AllLinks,fPickEven[#]&]];
          
        
          
                  
          p=First[Sort[Abs[AllLinks]]];
          AllLinks=Union[Select[AllLinks,SameQ[Abs[#],p]&]]; 
              
          
          p=Last[Sort[AllLinks]]
          (*p=
              First[Sort[AllLinks,
                  Abs[#1[[2]]]<Abs[#2[[2]]]&]]*)]]];(*kraj prvog if-
        a*){p[[1]],Flatten[p[[2]]]}]


MinDowAltKL[Con_String]:=Module[{l,res},l=fProjections[Con];
    If[Length[l]==1&&Head[l[[1]]]==List,l=l[[1]]];
    res=Map[MinDowProjAltKL[#]&,l];
    
    l=First[Sort[Abs[res]]];
    res=First[Union[Select[res,SameQ[Abs[#],l]&]]];
    res={res[[1]],Flatten[res[[2]]]}]

(*#########################################################*)
\
(*#########################################################*)
(*#########################################################*)
(*28.8.2003  po abs- za buduce narastaje- srediti*)
(*#########################################################*)
(*#########################################################*)
(*#########################################################*)
 
 (*############### AMPHICHEIRALITY ######################## *)

(*AmphiProjAltKL[Ul_] testira amphicheiralnost projekcije date
Konvejevim simbolom dowkerom ili pdata  *)

      
  AmphiProjAltKL[Ul_]:=Module[{G,sG,sG1,Aut,wr,res, i,j,k},
    G=fGaussExt[Ul];
    sG=fGraphInc[Ul][[2]];
    wr=Apply[Plus,sG];
    If[wr==0,sG1=-sG;
      G=FromUnorderedPairs[fGraphInc[Ul][[1]]];
      Aut=Drop[Automorphisms[G],1];
      (*Print[Aut];*)(*ovde izdvajamo automorfizme koji cuvaju znak*)Aut=
        If[Aut!={},
          If[First[
                Union[Table[
                    If[SameQ[
                          Sign[Table[
                                Table[
                                  Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]],{j,
                                    Length[sG]}],{i,Length[Aut]}]][[k]],
                          sG1]==True,Aut[[k]],{}],{k,
                      Length[Aut]}]]]=={},
            Drop[Union[
                Table[If[
                    SameQ[Sign[
                            Table[Table[
                                Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]],{j,
                                  Length[sG]}],{i,Length[Aut]}]][[k]],
                        sG1]==True,Aut[[k]],{}],{k,Length[Aut]}]],1],
            Union[Table[
                If[SameQ[
                      Sign[Table[
                            Table[Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]],{j,
                                Length[sG]}],{i,Length[Aut]}]][[k]],
                      sG1]==True,Aut[[k]],{}],{k,Length[Aut]}]]],Aut];
      If[Aut!={},res=1,res=0];
      Aut=If[Aut!={},Map[ToCycles,Aut],Aut];
      Aut=Union[Map[Length,Flatten[Aut,1]]];
      If[Aut!={}&&First[Aut]==1,Drop[Aut,1],Aut];
      , res=0 ]; 
    res]
    

AmphiAltKL[Conway_String]:=Module[{p,i,res},
    p=Flatten[fProjections[Conway]];
    i=1;
    res=AmphiProjAltKL[p[[i]]];
    If[res==1,Print["Amphicheiral   ",p[[i]]],
               While[i<Length[p]&&SameQ[res,0],
                 i++; res=AmphiProjAltKL[p[[i]]]
                ];
  If[i<Length[p],Print["Amphicheiral   ", p[[i]]];
   Print["Amphicheiral projection: ",p[[i]]];
        ,Print["Non-amphicheiral"]]]
    ] 
(*#######################################################*)     
     
(*#######################################################*)
(*#######################################################*)
(*############### CUTTING NUMBER ######################## *)

(* ^^^^^^^^^^^^^^^
(*  koristi se za UInfty *)
fAddSign[LL_List,{el_Integer,elS_Integer}]:=Module[{z,i,p,l=LL},
    z=1-elS;
    p=Position[LL,{el,0}][[1,1]];
    Do[
       l=ReplacePart[l,{l[[i,1]],z},i] ;
      z=1-z
      ,{i,p,Length[LL]}];
    z=1-elS;
    Do[z=1-z;(*Print[l];*)
           l=ReplacePart[l,{l[[p-i,1]],z},p-i]
         
          ,{i,1,p-1}];
    l ]*)
    (*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*)
(*fGiveSign[LL_List,SL_List]:=Module[{p,i,l=LL,saznakom=SL,sz={},sa},
     If[SL=={},
            p=l[[1]]; 
            Do[   p=ReplacePart[p,{p[[i,1]],Mod[i,2]},i];
                  ,{i,1,Length[p]}];
            l=ReplacePart[l,Prepend[p,1],1];
             sz=p,
            (*ako vec imamo neke sa znacimo trazimo ih u drugim \
komponentama*)
\
                sa=Union[Flatten[Map[Take[#,1]&,saznakom]]];
         (*   Print["SSSSS",sa,saznakom];*)
            
      Do[  If[Length[l[[i,1]]]!=0,
                          (*znaci da ona nije oznacena*)
                     \
      presek=Select[l[[i]],MemberQ[sa,#[[1]]]&];
                          If[presek!={},
                                  
            presek=Select[saznakom,#[[1]]==presek[[1,1]]&][[1]];
                                     p=fAddSign[l[[i]],presek];
                                   l=ReplacePart[l,Prepend[p,1],i];
                                   p=Map[If[MemberQ[sa,#[[1]]],-1,#]&,p];
                               sz=Union[sz,Complement[p,{-1}]]
                             ]]
        ,{i,1,Length[l]}]
      ];
    {l,sz}
    ]*)
     (*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*)
(* Pravi alternirajuci kad se polazi od pData... *)

(*fMakeAlt[CoLe_List,LL_List]:=Module[{i,pom=0,
CoPom={},saznakom={},ind=True},

      l=iteratedTake[LL,CoLe];
      While[ind,
                 (*Print["SZ",saznakom];*)
                  
      saznakom=fGiveSign[l,saznakom];
                  (* 
        Print["posle fGiveSigna",saznakom[[1]],"\n \n",saznakom[[2]],
            "\n \n"];*)
                  l=saznakom[[1]];
      saznakom=Last[saznakom];
                 (* Print[Map[ Length[#[[1]]]&,l]];*)
                  
      ind=Union[Map[ Length[#[[1]]]&,l]]!={0}
            ];
     Map[Drop[#,1]&,l]
    ]*)
(*pomocna za fMakeGaussPD- menja odnose iznad ispod za neparne PData *)
(*^^^^^^^^^^^^^^^^^^^^^^^*)
fMenjaObePojavePD[LL_List,el_Integer]:=Module[{p,i,l=LL},
    p=Select[l,#[[1]]==el&];
    (*Print[p,p[[1,1]]];*)
    (*Print[Position[l,p[[1,1]]]];*)
    
    p=Map[Take[#,1][[1]]&,Position[l,p[[1,1]]]];
    Do[
      l=ReplacePart[l,{l[[p[[i]],1]],1-l[[p[[i]],2]]},p[[i]]],{i,1,2}];
    l
    ]
(*radi od Ext Gausa sa odnosima iznad ispod*)
(* 
  ako moze sredi parnost *)
(*NE RADI RETROAKTIVNO!!!! *)
(*^^^^^^^^^^^^^^*)
(*koristi je UInfty*)
(*fParno[LL_List]:=Module[{l,lGExt,i,p},
    l=LL;
    lGExt=Map[Flatten,l];
    lGExt=Map[Take[#,{1,Length[#],2}]&,lGExt];
    (*extGaus bez odnosa iznad ispod *)
    
    Do[p= Flatten[Position[Flatten[Take[lGExt,i-1]],lGExt[[i,1]] ]];
          (*pozicija prvog u sledecoj komponenti u 
          prethodnom delu gausa*)
   If[p!={},
               If[Not[EvenQ[p[[1]]]],
                     lGExt=ReplacePart[lGExt,RotateLeft[lGExt[[i]]],i];
                       l=ReplacePart[l,RotateLeft[l[[i]]],i]
                ]]  ,{i,2,Length[l]}];
    l]
(* ako se posle secenja za PD u istoj komponenti jave dva ista 
uzastopna broja (ciklicno) sa suprotnim drugim znakom 
npr. {-18,1},{-18,0} izbacujemo ih.*)
(*radi od gaus ext sa iznad ispod*)
(*koristi je UInfty*)
fIzbaciPetlju[GK_List]:=Module[{p,l,pom,res=GK, i},
    l=Map[Take[#,1]&,GK]; 
    el=Union[Select[l,Count[l,#]==2&]];
    p=el;
    If[p!={},
          p=Map[Union[Flatten[Position[l,#]]]&,p];  
            l=Map[If[#!={1,Length[GK]},RotateLeft[l,#[[2]]-1],l]&,p]; 
       p={};
         Do[ p=Append[p,Flatten[Position[l[[i]],el[[i]]]]]
          ,{i,1,Length[el]}];
                l=GK; 
         Do[If[p[[i]]=={1,Length[GK]},
                  l=Select[l,#[[1]]!=el[[i,1]]&]
               ]
        ,{i,1,Length[el]}];res=l
      ];
    res]*)
    (*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*)
    (*koristi je UInfty*)
fMakeGaussPD[PD_List]:=Module[{l,n,i,l1={},nepar},
      n=Length[PD[[2]]];
     l=Complement[Range[2n],Abs[PD[[2]]]];
    (*svi brojevi do 2 duzine koji vec nisu u pData *)
    
    Do[l=ReplacePart[l,{l[[i]]*Sign[PD[[2,i]]],PD[[2,i]]},i],{i,1,n}];
    (*Print["www",l]; *)(*uredjeni parovi sa sve znacima*)
      
    nepar=Select[l,Or[OddQ[#[[2]]],And[EvenQ[#[[1]]],EvenQ[#[[2]]]]]&];
       nepar=Flatten[Map[Take[#,-1]&,nepar]];
    (*Print[nepar];*)
    l1=Map[#-#&,Range[2n]];
    Do[ l1=ReplacePart[l1,{l[[i,2]],0},Abs[l[[i,1]]]];
            l1=ReplacePart[l1,{l[[i,2]],0},Abs[l[[i,2]]]],{i,1,n}];
    (*Print["???",2PD[[1]],
          l1];*)
    (*ubacili smo indikatore i sad treba da ga napravimo da \
je alternirajuci*)
    l1=fMakeAlt[2PD[[1]],l1];
    (*menjamo odnose onima koji su u pdata bili neparni a to su i ovde*)
    
    l1=Flatten[l1,1];
    (* Print[" pre kvarenja ",l1]; *)
    
    If[nepar!={},
      Do[l1=fMenjaObePojavePD[l1,nepar[[i]]]    ,{i,1,Length[nepar]}]];
    (* Print[" posle kvarenja ",l1]; *)
    l1
    ]
(*pomocna za fBrCo- izbacuje iz Gaus parova i
 zavisnosti od "ociscenog"Gausa*)
(*^^^^^^^^^^^^^^^^^^^^^^^^^^*)
fBezComp[GPairs_List, GG_List]:=Module[{ll,p},
    ll=GPairs;
    Do[ p=Select[ll[[i]],MemberQ[GG[[i]],#[[1]]]&];
               ll=ReplacePart[ll,p,i];
      ,{i,1,Length[GG]}];
    ll
    ]
fUredjenobezdvojnih[LL_List]:=Module[{i,redLL=Flatten[LL]},
    Do[  redLL=
        Join[Take[redLL,i],Select[Drop[redLL,i],#!=redLL[[i]]&]],{i,
        1,Length[redLL]/2}];
    redLL
    ]
f2Ista[l_List]:=Module[{i=1,res=0},
    If[Length[l]==2,
            If[l[[1,2]]==l[[2,2]],
                      res=2],
                   While[i<Length[l]-1&&l[[i,2]]!=l[[i+1,2]],i++];
                            If[i!=Length[l]-1,
                            res=i+1,
                              If[l[[i,2]]==l[[i+1,2]],res=i+1,res=0]]];
    res]

(* pomocna za fDataForPData - alternira komponentu *)

fMakeCompAlt[Komp_List,ind_]:=Module[{l=Komp},
    l=Map[{#[[1]],Mod[Position[l,#][[1,1]],2]}&,l];
     l={l,Map[{#[[1]],1-#[[2]]}&,l]};
    If[ind!=-1,l=Select[l,#[[1,2]]==ind&]];
    l
    ]
(*pomocna za fPDataCuttNo, radi pd ext Gausa koji je napravljen u fBreakCp*)
(* 
  rezultata je lista svih mogucih alternirajucih*)

fDataForPData[GE_List]:=Module[{res,p,i},
    res={fMakeCompAlt[GE[[1]],GE[[1,1,2]]]};
    Do[ p=fMakeCompAlt[GE[[i]],-1];
           res=Join[Map[Append[#,p[[1]]]&,res],Map[Append[#,p[[2]]]&,res]]
      ,{i,2,Length[GE]}];
    res
    ]
(*iz liste svih potencijalnih bira dobru *)
(*kako?*)
(* 
  Oduzme potencijalnu od pocetne*)
(* 
  dobije listu uredjenih parova koji oblika {broj,a}
      gde je a 0- ako su znaci isti
               1- ako su razliciti*)
(* ako ima 1 trazimo drugu
      ako nema 1- dobra je;0 su dobri a 2 pokvareni*)

fSelectForPData[LL_List,poc_List]:=Module[{i=1,j,ind=True,pom,pF},
              pF=Flatten[poc];
            While[ind,
                         pom=Flatten[LL[[i]]];j=2;
                 
      While[j<=Length[pom],
        pom=ReplacePart[pom,Abs[pom[[j]]-pF[[j]]],j];
                          j=j+2
            ];
      pom=Partition[pom,2];
      If[Length[pom]==2Length[Union[pom]],
           (*nasli smo pravi!!!*)
                ind=False;
            pom=Select[Union[pom],#[[2]]==1&];
             pom=Flatten[Map[Take[#,1]&,pom]]
        ];
      i++];
    pom
    ]
(*radi od gaus EXT sa "nekim" iznad ispod*)
(* 
  daje pData *)
(*pomocna za fBReak CO*)

fPDataCuttNo[GE_List]:=Module[{l=Flatten[GE],svi,pokvareni},
    (*pravimo se mogucnosti za alternirajuci*)
    svi=fDataForPData[GE];
    (*nalazimo prvu "dobru" alternirajucu i ona nam daje njene pokvarene*)
   
     pokvareni=fSelectForPData[svi,GE];
    l=Take[l,{1,Length[l],2}];
    (*za pdata bez iznad ispod *)
    {Map[Length[#]&,GE]/2,
      fPDizNiz[l,pokvareni]}]
(*Formira PData iz nase vrste gausa*)

fPDizNiz[niz_List,gr_List]:=Module[{n=Flatten[niz]},
    n=Union[Map[{#,Flatten[Position[n,#]]}&,n]];
    n=Map[If[EvenQ[#[[2,1]]],{#[[1]],{#[[2,2]],#[[2,1]]}},#]&,n];
    n=Map[If[MemberQ[gr,#[[1]]],
                          {#[[2,2]],Sign[#[[1]]]*#[[2,1]]},
                          { #[[2,1]],Sign[#[[1]]]*#[[2,2]]}   
                          ]&,n];
    n=Sort[n,#1[[1]]<#2[[1]]&];
    n=Take[Flatten[n],{2,2Length[n],2}];
    n]
fBrCo[Ul_,k_Integer]:=Module[{l,p,pd1,p1,pd,pd0,nepar,res,pom},
               pd=Ul;(* sad nam je pd PData *)
                
    pd=ReductionKnotLink[pd];
                (*sad moramo da napravimo odgovarajuceg Gausa i da ga iz \
alterniramo*)
                
    pd0=fMakeGaussPD[pd];(* 
      Gaus parovi sa iznad-ispod tj. 0-1 *)
                
    pd1=Flatten[
        Map[Take[#,1]&,pd0]]; (*Gaus bez iznad ispod*)
                
    pd0=iteratedTake[pd0,2pd[[1]]];        
    	(*gaus za seckanje podeljen na komponente*)
                
    pd1=iteratedTake[pd1,2pd[[1]]];
                
    p=pd0[[k]];(*k-ta tj. komponenta koju secemo*)
                 
    p1=Flatten[Map[Take[#,1]&,p]];(*p bez iznad ispod*)
                 
    pd1=Map[UnsortedComplement[#,p1]&,Drop[pd1,{k}]];
                
    pd0=Drop[pd0,{k}]; (*gaus parovi bez izbacene komponente*)
               \
 Do[ pd0=fBezComp[pd0,pd1] ,{i,1,Length[pd0]}];
                (*do petlja izbacuje one koji nisu preostali *)      
    	 pd0=UnsortedComplement[pd0,{{}}];
    	 (*sad nam treba njihov pravi poredak tj. 
          izbacujemo drugu pojavu*)
              
    pd1=fUredjenobezdvojnih[pd1];
              l=Map[Length,pd0];
              If[pd1=={},
                  (* Print["Direct product or 2-component link"]; *)
        res={{},{}},
      (*  Print["posle secenja za PD ",pd0];*)
      
      pd0=Map[fIzbaciPetlju[#]&,pd0];(*sredjujemo parnost- 
          ako ne moze da sredi dobicemo dva neparna i dva parna uparena*)
    
        pd0=fParno[pd0];
      (*Print["Posle parnosti ",pd0];*)
     If[SameQ[pd0,{{}}],
        res={{},{}},
      pd0=fPDataCuttNo[pd0];
      (*Print["Pre redukcije ",pd0];*)
      res=ReductionKnotLink[pd0]]];
    res
    ]
          
(*rezultat salje u fBreak Comp *)
fBreakComp[Ul_,k_Integer]:=Module[{pd={},i},
If[k>fComponentNo[Ul],
Print["Maximum number of components is: ",fComponentNo[Ul]],
    If[SameQ[fComponentNo[Ul],1],
              Print[ "Knot"],
                If[SameQ[Head[Ul],String],pd=fCreatePData[Ul],pd=Ul];
                (*ako je konvej pravimo PD ako ne sve je PD*)
                \
  pd=fBrCo[pd,k]]];
    pd]
(*radi od pdata i daje rezultat pri jednom secenu svake komponente ponaosob*)

BreakCoAll[Ul_]:=Module[{i=1,res={},p={3}},
    If[fComponentNo[Ul]==1,res={};Print["Knot"],
      While[p!={1}&&i<=fComponentNo[Ul],
                    (* Print["Ulaz: ",Ul];*)
                  
        p=fBreakComp[Ul,i];
                (*  
          Print["Isecena je ",i,"-ta komponenta: ",p];*)
                 
        If[Not[MemberQ[{{{},{}},{{0},{}},{{},{0}}},p]],
                       i++;res=Append[res,p],
                      p={1};res={0}]
             ];
        If[p!={1},
        res=Table[{res[[i,1]],Flatten[res[[i,2]]]},{i,Length[res]}]]];
    res
    ]
(*vraca {0} ako je razvezao a u suprotnom vraca listu ili {} ako je cvor *)

CuttNo[Ul_]:=Module[{pd,cNo,i,R,R1={},p={3},res="0"},
       If[SameQ[Head[Ul],String],pd=fCreatePData[Ul],pd=Ul];
    (*sad imamo pdata -odmah reduction*)
    
    cNo=MemberQ[ReductionKnotLink[pd],{}];
    If[fComponentNo[pd]==1,
         Print["Knot"];cNo=0,
          pd=ReductionKnotLink[pd];
        If[cNo==True,Print["Razvezan posle redukcije"];cNo=0,
        (*Ako je odmah redukovan- 
            cutting No je nula u suprotnom radimo dalje*)
        (*Print[cNo,
              res,pd];*)
        R=BreakCoAll[pd];
        (*Print["trenutna lista ",Length[R],R];*)
        
        cNo=1;(*postavljamo ga na 1 jer smo tacno jednu komponentu isekli*) 
        While[
          Not[MemberQ[R,0]],(*Print["RADI WHILE"];*)
                        
          R1={};i=0;
                         While[p!={0}&&i<Length[R],  i++;
                                     p=BreakCoAll[R[[i]]];
                                     
            R1=Append[R1,p](*Print[i,"tren ",
                R1]*)
                                   ];
                       If[p=={0},R={0},R=Flatten[R1,1]];
            (*  Print["TRENUTNA ",R];*)
             cNo++]
        ]];(*kraj prvog if-a*)
    cNo
    ]
(*######################################################*)
(*######################################################*)
(*############### REAL CUTTING  ######################## *)

(*radi od gausEXT- sva moguca "secenja sa fix krajevima"*)

fCuttComp[gg_List]:=Module[{res={},p=gg, i},
    Do[  p=RotateLeft[p];
            res=Join[res,{p},{Reverse[p]}]
      ,{i,1,Length[gg]}];
    Union[res]
    ]
(*vraca {a,b}
      a=1 ako jeste izomorfizam 0 inace
    i b=1 ako je sense preserving
        b=-1, ako je sense reversing*)

fCuteIso[L1_List,L2_List]:=Module[{l={},i},
    Do[l=Append[l,{L1[[i]],L2[[i]]}],{i,1,Length[L1]}];
    l=Sort[Union[l],Abs[#1[[1]]]<Abs[#2[[1]]]&];
    If[2Length[l]==Length[L1],i=1,i={0,0}];
    If[i==1,
             If[Union[Map[#[[1]]*#[[2]]>0&,l]][[1]],
                (*sense-preserving *)    i={i,1},
                  (*sense=reversing*)i={i,-1}]];
    i]
(*radi sa EXTGauss*)
(*iz njihove liste vraca predstavnike klasa izomorfizama \
i to ako je ind=1 onih koje cuvaju znake
 a -1 ako ne vodi racuna o tome*)

fCuteIsoClass[LL_List,ind_]:=
  Module[{l=LL,res,i,pom,SenseRev={},SensePres={},p,f},
    l=Drop[l,1];
    res={LL[[1]]};
    While[l!={},
                    p=Last[res];i=1;
                    While[i<=Length[l],
                      pom=fCuteIso[p,l[[i]]];
                    
                       If[ind==1,f={1,1},pom=pom[[1]];f=1];
                     (*ako je ind 1 onda moraju cuvati
                ako je -1 onda je smao prvi clan bitan da je 1*)
        	   
        If[pom==f,
                          l=Drop[l,{i}],i++]
                     ];
      (*  Print["Ovde mora biti prazna za 3 ",l];*)
           
      If[l!={}, res=Append[res,First[l]];
               l=Rest[l]];
            If[Length[l]==1,res=Append[res,l[[1]]];l={}]
      ];
    res
    ]
(*pravi linkove tako sto na el dodaje sve iz r*)
(* i sve to flatten*)

fJoinL[el_List,r_List]:=Module[{res={el}},
    res=Map[Flatten[Append[res,#],1]&,r];
    (*Print["duzina rezultata ",Length[res]];*)res
    ]
(*Radi od PData ili Konveja*)

fCuttRealKL[UL_]:=Module[{pd=UL,gaus,sPreserv,svi,brCo,len},
    If[SameQ[Head[UL],String],pd=fCreatePData[UL]];
    pd=ReductionKnotLink[pd];
    gaus=fGaussExtSigns[fDowfromPD[pd]];
    If[Head[gaus[[1]]]==Integer,gaus={gaus}];
    gaus=Map[fCuttComp[#]&,gaus];
    (*Print["gaus",
          Map[Length[#]&,gaus]];*)
          (*sad imamo sve moguce rasecene gause
        podliste su mogucnosti za svaku komponentu*)
    (*sredjujemo linkove- i oni postaju 1 komponenta*)
        brCo=Length[gaus];
    If[brCo!=1,
            len=Map[Length[First[#]]&,gaus];
   (*lista duzina svake komponente spremno za iterated take*)
    res=gaus[[1]];(*Print["res ",res];*)
    Do[ res=Flatten[Map[fJoinL[#,gaus[[i]]]&,res],1]
        ,{i,2,brCo}];
      gaus={res}
         ];
    (*Print["Gaus za POCETNIKE ",
          Length[gaus]];*)
    (*trazimo neizomorfne predstavnike *)
   (*cuvaju znake*)
    sPreserv=Flatten[Map[fCuteIsoClass[#,1]&,gaus],1];
    Print["No of classes under sign-preserving isomorphisms: \
",Length[sPreserv]];
    (*okrecu znake*)(*Print["******************"];*)
    
    svi=Flatten[Map[fCuteIsoClass[#,-1]&,gaus],1];
    Print["No of classes: ",Length[svi]];
    (*{sPreserv,svi}*)]
(*#############################################################*)
(*############### SPLITTING NUMBER  ######################## *)

(*radi od gauss Ext broj disjunktnih tj, razdvojenih komponenti *)

fDisjointComp[LL_List]:=Module[{i,br=0},
    If[Head[LL[[1]]]==Integer,
      br=1];(*ako je cvor onda je to smao 1 komponenta*)
    
    If[Head[LL[[1]]]==List,
      Do[
        If[Intersection[LL[[i]],Flatten[Drop[LL,{i}]]]=={},br++]
        ,{i,1,Length[LL]}]];
    br
    ]
(*Menja znak u pojedinacnoj tacki- radi od PData*)

fCrossChangePD[PD_List,pos_Integer]:=Module[{l,n,i},
      n=Length[PD[[2]]];
     l=Complement[Range[2n],Abs[PD[[2]]]];
    (*svi brojevi do 2 duzine koji vec nisu u pData *)
    
    Do[l=ReplacePart[l,{l[[i]]*Sign[PD[[2,i]]],PD[[2,i]]},i],{i,1,n}];
    (*Print["www",l];*)
    l=ReplacePart[l,-1*Reverse[l[[pos]]],pos];
    (*Print["www",l];*)
    l=Sort[l,Abs[#1[[1]]]<Abs[#2[[1]]]&];
    l=Flatten[Map[Take[#,-1]&,l]];
    l=ReductionKnotLink[{PD[[1]],l}];
    l={UnsortedComplement[l[[1]],{{}}],l[[2]]};
    l={UnsortedComplement[l[[1]],{0}],l[[2]]}
    ]
(*Radi od PD!-promeni znake u svim tackama pojedinacno*)

fCrChAllPD[PD_List]:=Module[{res,l},
    l=Range[Length[PD[[2]]]];
    res=Map[fCrossChangePD[PD,#]&,l]
    ]
fGaussPDForDisjoint[PD_List]:=Module[{l,i,res,n},
      n=Length[PD[[2]]];
     l=Complement[Range[2n],Abs[PD[[2]]]];
    (*svi brojevi do 2 duzine koji vec nisu u pData *)
    
    Do[l=ReplacePart[l,{l[[i]]*Sign[PD[[2,i]]],PD[[2,i]]},i],{i,1,n}];
    l=Abs[l];
    res=Table[0,{i,2n}];
    Do[res=ReplacePart[res,i,l[[i,1]]];
          res=ReplacePart[res,i,l[[i,2]]]
      ,{i,1,n}];
    iteratedTake[res,2PD[[1]]]
    ]
(* Radi od Konveja ili pdata*)

SplittNo[PData_]:=Module[{pd=PData,gaus,korak=0,CoNo,ind=True},
    If[SameQ[Head[pd],String],
              pd=fCreatePData[PData]];
    pd= ReductionKnotLink[pd];(*za svaki slucaj*)
    
    pd={UnsortedComplement[pd[[1]],{{}}],pd[[2]]};
    pd={UnsortedComplement[pd[[1]],{0}],pd[[2]]};
    CoNo=Length[pd[[1]]];(*broj komponenti i kriterijum za zaustavljanje*)
   
     If[CoNo==1,
             (*Knot *)korak=0,
              (*Link*)
                gaus=fGaussPDForDisjoint[pd];
                br=fDisjointComp[gaus];
                pd={pd};
                While[br<CoNo&&ind,
                          pd=Union[Flatten[Map[fCrChAllPD[#]&,pd],1]];
                           korak++;                        
        		
        ind=Not[MemberQ[
              Union[Map[MemberQ[pd,#]&,{{{},{}},{{0},{}},{{},{0}}}]],True]];
                            
        ind=ind&&Not[MemberQ[Map[Length[#[[1]]]==1&,pd],True]];
                            (*ind je True ako nema razvezanih  tj. 
              ako treba dalje splitati*)                
        		If[ind,
                              gaus=Map[fGaussPDForDisjoint[#]&,pd];
                               
          br=Max[Map[fDisjointComp[#]&,
                gaus]] ](*kraj if *)
          ](*kraj while*)
      ];(*kraj \
prvog if *)
    korak
    
    ]
(*##############################################################*)
(****************************************************************)
(*#############################################################*)

(*############### CONVERSION FUNCTIONS ########################*)
fDToD[Con_String]:=Module[{n,p,DD},
    DD=Dowker[Con][[4]];
    If[Head[DD[[1]]]==List,
      n=Map[Length,DD];
      p={n,Abs[Flatten[DD]]*fGenSign[Con]}];
    If[Head[DD[[1]]]==Integer,n=Length[DD];
      p={{n},Abs[DD]*fGenSign[Con]}];
    p]

fKnotscapeDow[Con_String]:=
  Module[{l},(*ukoliko nema-u Konveju Dowker je nas abs*)
    l=Abs[fDToD[Con]];
    If[StringPosition[Con,"-"]!={},(*ako ima-*)
         l=Abs[fDToD[StringReplace[Con,"-"->""]]];
      l={l[[1]],
          l[[2]]*fGenSign[Con]*fGenSign[StringReplace[Con,"-"->""]]}];
    l]
(*################################################*)
(*############### MID-EDGE GRAPH  ######################## *)

(* radimo sa poliedrima *)
fMidEdgeGraph[L_List]:=Module[{l,ll,res={},p,ind,i},
    ind=SameQ[Union[Map[Count[L,#]&,L]],{2}];
    l=Sort[Union[L]];
    ll=Map[{#,Position[l,#][[1]]}&,l];
    Do[ p=Select[ll,Length[Intersection[#[[1]],l[[i]]]]==1&];
           p=Map[Sort[{i,#[[2,1]]}]&,p];
       res=Append[res, p],{i,1,Length[l]}];
       res=Flatten[res,1];
    If[Not[ind],res={Union[res]}];
    Flatten[res,1]
    ]

(*############### fKLfromGraph  ######################## *)

(*pomocna za fFixInc-ubacuje simbol temena na narednu poziciju*)

fPosDouble[L_List,red_Integer,el_Integer]:=
  Module[{l=L,p},(*Print[l];*)p=
      Position[l[[red]],el][[1,1]];(*Print[red,el,"pos",p];*)
      If[p!=4,
      l=Insert[l,el,{red,p+1}],l=Insert[l,el,{red,4}]];
      (*Print[l];*)l]
(*Ubacuje u listu incidencije dvojne veze*)

fFixInc[IL_List,dvoL_List]:=
  Module[{res=IL, i},Do[res=fPosDouble[res,dvoL[[i,1]],dvoL[[i,2]]];
      res=fPosDouble[res,dvoL[[i,2]],dvoL[[i,1]]],{i,1,Length[dvoL]}];res]
(*koliko se puta u datom redu pojavljuje element sa pozicije el i u kom \
redosledu*)

fBrojPojave[L_List,red_Integer,el_Integer]:=Module[{p,res,r},p=L[[red,el]];
    r=Map[If[Length[#]==2,#[[1]],#]&,L[[red]]];
    If[Length[p]==2,p=p[[1]]];
    (*Print[p];*)
    If[Count[r,p]==1,res=0,
      If[el==5,res=2,
        If[el==2,res=1,Switch[p,r[[el-1]],res=2,r[[el+1]],res=1];]]];
    res]

fPrvi[inc_List]:=
  Module[{l=inc},
    l=Map[Rest[#]&,l];(*izbrisemo oznake reda*)l=
      Select[Flatten[l,1],Length[#]==0&];
    If[l=={},l=0,l=l[[1]]];
    l]
(*vraca 0 ako su svi iskoristeni a prvi sledeci ako takav postoji*)
UnsortedComplement[x_List,y__List]:=
  Replace[x,Dispatch[(#:>Sequence[])&/@Union[y]],1]
(*u redu RED trazi p-tu pojavu elementa el,pa element 2 pozicije udaljen*)

fTrazi[L_List, red_Integer, el_Integer, brojPojave_Integer] := 
  Module[{p, pom, inc = L},
    Switch[brojPojave, 0, p = 1, 1, p = 2, 2, p = 1];
    
    (*ako se vec pojavio (kao 1) onda uzimamo drugu pojavu, inace prvu*)
    
    pom = Position[Rest[L[[red]]], el];
    pom = Select[pom, Or[Length[#] != 2, Length[#] == 2 && #[[2]] != 2] &];
    
    pom = pom[[p, 1]];
    (*ne uzimamo u obzir prvi el. 
          jer je on oznaka reda i zato posle ide ++*)
    
    pom = pom + 1;
    
    inc = ReplacePart[L, {L[[red, pom]], 0}, {red, pom}];
    (* stavljamo mu 0 jer smo njega trazili*)
    If[pom <= 3, pom = pom + 2, pom = pom - 2];
    (* ovo je pozicija u red - u el. koji je sledeci clan komponente*)
    (*i sad jos da oznacimo da je dati upotrebljen*)
    
    If[Length[L[[red, pom]]] != 2,
      inc = ReplacePart[inc, {inc[[red, pom]], 1}, {red, pom}]];
    
    {L[[red, pom]], {red, pom - 1}, inc}]
(*vraca element, njegovu poziciju, nov inc*)


fSrediParnost[L_List]:=Module[{res={},p,i,j},res={L[[1]]};
    Do[j=1;p=L[[i]];
      While[Not[MemberQ[Flatten[res],p[[j]]]],j++];
      If[EvenQ[Position[res,p[[j]]][[1,2]]-j],
        res=Append[res,RotateLeft[L[[i]]]],res=Append[res,L[[i]]]],{i,2,
        Length[L]}];
    res]

fKLfromGraph[InG_List]:=
  Module[{p,GNew,incPoc,inc,dvojne,res={},komp={},
  el,red=1,brP,prvi,DUPLE={}},
    p=Map[Sort[#]&,InG];
    If[Map[Count[p,#]&,Union[p]]=={4},res={{1,1},{4,2}},
      If[Union[Map[Count[p,#]&,p]]=={2},(*ovo je cvor ili link "n"*)el=
          Length[p]/2;
        res=Dowker[ToString[el]][[4]];
        If[OddQ[el],res={{el},Flatten[res]},res={{el/2,el/2},Flatten[res]}],
        p=fPlanarEmbGraph[p];
        GNew=p[[1]];(*planar embeded graf-Unorderedpairs*)
        inc=p[[2]];(*lista incidencije*)dvojne=
          Union[Select[GNew,Count[GNew,#]==2&]];
         inc=fFixInc[inc,dvojne];
        el=inc[[1,2]];brP=fBrojPojave[inc,1,2];
        komp=Append[komp,el];
        inc=ReplacePart[inc,{inc[[1,2]],1},{1,2}];
        (*ovaj el smo vec iskoristili*)
        prvi=el;(*pocetni element u komponenti*)
      
        While[
          And[el!=0,
            Union[Map[Count[Flatten[res],#]&,Flatten[res]]]!={2}],
           
          While[And[Length[el]!=2,Count[komp,prvi]!=3],
            p=el;
            el=fTrazi[inc,el,red,brP];
           inc=el[[3]];
            If[Length[el[[1]]]!=2,red=p;
              brP=fBrojPojave[inc,el[[2,1]],el[[2,2]]+1];
              el=el[[1]];
              (*nov element i njegov broj pojave*)
              (*nalazi u redu el element 
jednak redu+2 ili-2 pozicije*)komp=
                Append[komp,
                  el],(*drugi-4 put na isti "nacin" prolazimo kroz jednu \
tacku-znaci kraj*)
          el=el[[1]]]];

         If[First[komp]==Last[komp],komp=Drop[komp,-1]];
          res=Append[res,komp];
          res=UnsortedComplement[res,{{}}];
          el=fPrvi[inc];
          If[el!=0,komp={el};
            (*posto pocinjemo novu komponentu vracamo inc na nekoristeno \
stanje*)red=Position[inc,el];
            red=Select[red,Length[#]==2&&#[[2]]!=1&][[1]];
            inc=ReplacePart[inc,{inc[[red[[1]],red[[2]]]],1},red];
            (*stavljamo daje iskoristen*)
           
            brP=fBrojPojave[inc,red[[1]],red[[2]]];
            red=red[[1]];]];(*kraj while*)
           res=fSrediParnost[res];
       res=fDowfromGaussExt[res];
        res={res[[1]],Flatten[res[[2]]]};
     If[fComponentNo[res]==1,res=fMinDowKnotPr[res]]]];(*kraj \
ifova*)
      res]
      
      
      
      
      (*############################## DOWKER CODES \
#########################*)


fDowkerCode[Ul_List] := 
  Module[{inc = Ul, res = {}, komp = {}, el, red = 1, brP, prvi, DUPLE = {}, \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
p},
    el = inc[[1, 2]]; brP = fBrojPojave[inc, 1, 2];
    komp = Append[komp, el];
    inc = ReplacePart[inc, {inc[[1, 2]], 1}, {1, 2}];
    (*ovaj el smo vec iskoristili*)
    prvi = el;(*pocetni element u komponenti*)
    While[
      And[el != 0, Union[Map[Count[Flatten[res], #] &, Flatten[res]]] != \
{2}],
       While[And[Length[el] != 2, Count[komp, prvi] != 3], p = el;
        el = fTrazi[inc, el, red, brP];
        inc = el[[3]];
        If[Length[el[[1]]] != 2, red = p;
          brP = fBrojPojave[inc, el[[2, 1]], el[[2, 2]] + 1];
          el = el[[1]];
          (*nov element i njegov broj pojave*)(*nalazi u redu el element \
jednak redu + 2 ili - 2 pozicije*)komp = 
            Append[komp, el],(*drugi - 
              4 put na isti "nacin" prolazimo kroz jednu tacku - znaci kraj*)
            el = el[[1]]]];
      If[First[komp] == Last[komp], komp = Drop[komp, -1]];
      res = Append[res, komp];
      res = UnsortedComplement[res, {{}}];
      el = fPrvi[inc];
      If[el != 0, komp = {el};
        (*posto pocinjemo novu komponentu vracamo inc na nekoristeno stanje*)
          red = Position[inc, el];
        red = Select[red, Length[#] == 2 && #[[2]] != 1 &][[1]];
        inc = ReplacePart[inc, {inc[[red[[1]], red[[2]]]], 1}, red];
        (*stavljamo daje iskoristen*)
        brP = fBrojPojave[inc, red[[1]], red[[2]]];
        red = red[[1]];]];(*kraj while*)res = fSrediParnost[res];
    res = fDowfromGaussExt[res];
    res = {res[[1]], Flatten[res[[2]]]};
    res=MinDowProjAltKL[res];
    res]
    
    fCorr[Ul_List] := Module[{pp, pp1, i},
    pp = Map[First, Ul];
    pp1 = Length[Union[pp]];
    pp1 = Flatten[Map[First, Table[Position[pp, i], {i, pp1}]]];
    pp = Table[Ul[[pp1[[i]]]], {i, Length[pp1]}];
    pp]
    
    fDowCodes[n_Integer] := Module[{vv, vv1, vvf, kk, kk1, mm, i, j, k},
    Run["plantri -apc2 " <> ToString[n] <> " " <> ToString[n] <> ".txt"];
    vv = Import[ToString[n] <> ".txt"];
    DeleteFile[ToString[n] <> ".txt"];
    vv = StringDrop[
          "{" <> StringReplace[
              vv, 
              {ToString[n] <> " " -> "", "\n" -> ",", "," -> ","}], -1] 
<>
         "}";
    vv = Partition[Map[ToCharacterCode, Map[ToString, ToExpression[vv]]] - \
96,
         n];
    vv = Union[
        Table[If[Max[Map[Length, vv[[i]]]] <= 4 , vv[[i]], {}], {i, 
            Length[vv]}]];
    vv = Partition[Flatten[If[SameQ[vv[[1]], {}], Drop[vv, 1], vv], 1], n];
    (* lista svih sa valencom <= 4 *)
    vv1 = 
      Table[Union[
          Map[Sort, 
            Union[Table[
                Table[Length[vv[[k, vv[[k, j, i]]]]], {i, 
                    Length[vv[[k, j]]]}], {j, Length[vv[[k]]]}]]]], {k, 
          Length[vv]}];
    vv = Union[
        Table[If[
            SameQ[Select[
                vv1[[i]], (Length[#] == 2 && MemberQ[#, 4]) && 
                    Union[Map[EvenQ, #]][[1]] &], {}], vv[[i]], {}], {i, 
            Length[vv1]}]];
    vv = Flatten[If[SameQ[vv[[1]], {}], Drop[vv, 1], vv], 1];
    vv = Partition[
        Table[If[
            Length[vv[[i]]] == 2, {vv[[i, 1]], vv[[i, 1]], vv[[i, 2]], 
              vv[[i, 2]]}, vv[[i]]], {i, Length[vv]}], n];
    (* ovde treba napuniti odgovarajuca trovalentna temena *)
    vv1 = Map[Split, Map[Sort, Map[Flatten, vv]]];
    vv1 = Map[Max, Table[Map[Length, vv1[[i]]], {i, Length[vv1]}]];
    vv = Union[Table[If[vv1[[i]] > 4, {}, vv[[i]]], {i, Length[vv1]}]];
    vv = If[SameQ[vv[[1]], {}], Drop[vv, 1], vv];
    (* svi obojivi sa dvojnim vezama samo u dvovalentnim temenima *)
    vvf = Select[vv, SameQ[Union[Map[Length, #]], {4}] &];
    (* vvf su vec obojene *)
    vv = Complement[vv, vvf];
    vv1 = 
      Table[Table[
          If[vv[[j, i, 1]] == 
              vv[[j, i, 2]], {{vv[[j, i, 1]], i}, {vv[[j, i, 3]], 
                i}}, {}], {i, Length[vv[[j]]]}], {j, Length[vv]}];
    vv1 = Map[Sort, Table[Flatten[vv1[[i]], 1], {i, Length[vv1]}]];
    kk = Table[Complement[Range[n], Map[First, vv1[[i]]]], {i, Length[vv1]}];
    kk = Table[
        Table[{kk[[j, i]], 0}, {i, Length[kk[[j]]]}], {j, Length[kk]}];
    vv1 = Table[Union[vv1[[i]], kk[[i]]], {i, Length[kk]}];
    vv1 = Map[fCorr, vv1];
    vv = Table[
        Table[If[vv1[[j, i, 2]] != 0 && vv[[j, i, 1]] != vv[[j, i, 2]], 
            Insert[vv[[j, i]], vv1[[j, i, 2]], 
              Position[vv[[j, i]], vv1[[j, i, 2]]]], vv[[j, i]]], {i, n}], \
{j,
           Length[vv]}];
    vv = Union[vv, vvf];
    vvf = Select[vv, SameQ[Union[Map[Length, #]], {4}] &];
    (* vec obojene koje dodajemo na kraju *)
    vv = Complement[vv, vvf];
    kk = Table[Apply[Plus, Map[Length, vv[[i]]]], {i, Length[vv]}];
    kk = Table[Divide[Times[4, n] - kk[[i]], 2], {i, Length[kk]}];
    vv1 = 
      Table[Union[
          Flatten[Table[
              Table[{{j, vv[[k, j, i]]}, 
              Length[vv[[k, vv[[k, j, i]]]]]}, {i, 
                  Length[vv[[k, j]]]}], {j, Length[vv[[k]]]}], 1]], {k, 
          Length[vv]}];
    vv1 = 
      Table[Union[Map[Sort, Map[First, Select[vv1[[i]], #[[2]] == 3 &]]]], \
{i,
           Length[vv]}];
    kk1 = Map[Split, Map[Sort, Map[Flatten, vv]]];
    kk1 = 
      Map[Union, 
        Map[Flatten, 
          Table[Select[kk1[[i]], Length[#] < 4 &], {i, Length[kk1]}]]];
    kk1 = 
      Table[Select[vv1[[i]], Length[Intersection[#, kk1[[i]]]] == 2 &], {i, 
          Length[vv]}];
    kk1 = Table[KSubsets[kk1[[i]], kk[[i]]], {i, Length[kk]}];
    kk1 = 
      Table[kk = 
          Union[Table[
              If[SameQ[Length[Flatten[kk1[[j, i]]]], 
                  Length[Union[Flatten[kk1[[j, i]]]]]], 
                  kk1[[j, i]], {}], {i, 
                Length[kk1[[j]]]}]];
        kk = 
          If[Not[SameQ[kk, {}]], If[SameQ[kk[[1]], {}], Drop[kk, 1], kk], 
            kk],
        {j, Length[kk1]}];
    vv1 = 
      Table[Table[
          Union[kk1[[j, i]], Map[Reverse, kk1[[j, i]]]], {i, 
            Length[kk1[[j]]]}], {j, Length[kk1]}];
    (* dobili smo sve kombinacije ivica koje dodajemo da bismo napravili 4 - 
        valentan graf *)
    mm = Table[Table[
          kk = vv1[[k, j]];
          
          kk1 = Table[
              Position[vv[[k, kk[[i, 1]]]], kk[[i, 2]]], {i, Length[kk]}];
          
          Table[{kk[[i, 1]], 
              Insert[vv[[k, kk[[i, 1]]]], kk[[i, 2]], kk1[[i]]]}, {i, 
              Length[kk]}],
          {j, Length[vv1[[k]]]}], {k, Length[vv1]}];
    vv1 = 
      Table[Table[{i, vv[[j, i]]}, {i, Length[vv[[j]]]}], {j, Length[vv]}];
    vv1 = Table[Select[vv1[[i]], Length[#[[2]]] == 4 &], {i, Length[vv1]}];
    vv1 = 
      Table[Table[Union[mm[[j, i]], vv1[[j]]], {i, Length[mm[[j]]]}], {j, 
          Length[mm]}];
    vv1 = 
      Table[Table[Map[Last, vv1[[j, i]]], {i, Length[vv1[[j]]]}], {j, 
          Length[vv1]}];
    vv1 = Flatten[Complement[vv1, {{}}], 1];
    vv1 = Union[vv1, vvf];
    vv1 = 
      Table[Table[Flatten[{i, vv1[[j, i]]}], {i, Length[vv1[[j]]]}], {j, 
          Length[vv1]}];
    vv1 = Union[Map[fDowkerCode, vv1]];
    vv1
    ]   
    
    
    
    
    
(* ################### fGenKL ###################### *)


fBasic[n_Integer] := Module[{del, del1, i, j, k},
    If[n < 4, del = {}, del = compo[n];
      del = Select[del, Not[MemberQ[#, 1]] &];
      del = Select[del, Length[#] >= 2 &];
      del1 = Sort[Table[{Sort[del[[i]]], del[[i]]}, {i, Length[del]}]];
      del = Union[Table[Sort[del[[i]]], {i, Length[del]}]];
      del = Table[Select[del1, #[[1]] == del[[i]] &], {i, Length[del]}];
      del = 
        Table[Table[del[[j, i, 2]], {i, Length[del[[j]]]}], {j, 
            Length[del]}];
      del = 
        Map[Last, 
          Flatten[Table[
              Union[Table[
                  Union[Table[
                      RotateLeft[del[[k, j]], i], {i, Length[del[[k, j]]]}], 
                    Map[Reverse, 
                      Table[RotateLeft[del[[k, j]], i], {i, 
                          Length[del[[k, j]]]}]]], {j, 
                    Length[del[[k]]]}]], {k, Length[del]}], 1]]];
    del = {del};
    del = 
      Flatten[Table[
          Map[StringReplace[ToString[#], 
          {"{" -> "", "}" -> "", " " -> ""}] 
&,
             del[[i]]], {i, Length[del]}]];
    del]
    
fRtan[n_Integer] := Module[{stel, pp, ll, ll1, i, j},
    stel = fBasic[n];
    stel = 
      Map[ReadList[StringToStream[StringReplace[#, "," -> " "]], Number] &, 
        stel];
    stel = Map[fMakePart[#] &, stel];(*svaki strem veci od 1 partitionira :*)

    
    stel = Flatten[Map[fPerm[#] &, stel], 1];
    stel = 
      Table[Table[
          ToExpression[StringReplace[stel[[j, i]], {" " -> ","}]], {i, 
            Length[stel[[j]]]}], {j, Length[stel]}];
    stel = 
      Map[Last, 
        Union[Table[
            Union[Table[RotateLeft[stel[[j]], i], {i, Length[stel[[j]]]}], 
              Map[Reverse, 
                Table[RotateLeft[stel[[j]], i], {i, Length[stel[[j]]]}]]], \
{j,
               Length[stel]}]]];
    stel = 
      Table[Map[
          StringReplace[ToString[#], {"{" -> "", "}" -> "", " " -> ""}] &, 
          stel[[i]]], {i, Length[stel]}];
    stel = Table[StringReplace[stel[[i]], "," -> " "], {i, Length[stel]}];
    stel = 
      Table[StringReplace[
          ToString[stel[[i]]], {"{" -> "", "}" -> "", ", " -> ","}], {i, 
          Length[stel]}];
    stel = Table[fStelString[stel[[i]]], {i, Length[stel]}];
    stel = Flatten[Table[Permutations[stel[[i]]], {i, Length[stel]}], 1];
    stel = Table[Map[ToString, stel[[i]]], {i, Length[stel]}];
    stel = 
      Table[StringReplace[stel[[i]], {"{" -> "", "}" -> "", ", " -> " "}], \
{i,
           Length[stel]}];
    stel = Map[Reverse, Union[Map[Sort, stel]]];
    pp = Map[Sort, Map[Split, Map[Sort, stel]]];
    ll = Map[Reverse, 
        Map[Sort, Table[Map[Length, pp[[i]]], {i, Length[pp]}]]];
    ll1 = Union[ll];
    pp = Map[Flatten, Table[Position[ll, ll1[[i]]], {i, Length[ll1]}]];
    stel = 
      Table[Table[stel[[pp[[j, i]]]], {i, Length[pp[[j]]]}], {j, 
          Length[pp]}];
    stel = Table[Map[Split, stel[[i]]], {i, Length[stel]}];
    stel = Table[Map[Sort, stel[[i]]], {i, Length[stel]}];
    ll = Table[Map[Length, stel[[i, 1]]], {i, Length[stel]}];
    stel = Table[Map[Flatten, stel[[i]]], {i, Length[stel]}];
    stel = {stel, ll};
    stel]

fRepres[n_Integer, pp_String] := 
  Module[{vv, ll, ll1, t, n1, stel, ss, i,j, ff, ff1, ppp},
    vv = fRtan[n];
    ll = vv[[2]];
    ll1 = 
      Map[Flatten, 
        Table[Table[
            Table[Length[ll[[j]]] + 1, {i, ll[[j, i]]}] - i + 1, {i, 
              Length[ll[[j]]]}], {j, Length[ll]}]];
    ll = Map[Length, ll1];
    t = iteratedTake[Map[ToString, Flatten[ll1]], ll];
    ppp = StringPosition[pp, "2"];
    If[ppp != {} && StringPosition[pp, "*"] != {},
      i = 1;
      While[ppp != {} && i <= Length[ppp],
        If [ppp[[i, 1]] < StringPosition[pp, "*"][[1, 1]], ppp = Rest[ppp], 
          i++]
        ]];
    n1 = Length[ppp];
    ff = Map[Length, t];
    t = Select[t, Length[#] == n1 &];
    stel = Table[Flatten[Map[Permutations, {t[[i]]}], 1], {i, Length[t]}];
    ss = Table[
        Flatten[Table[
            StringReplacePart[pp, Flatten[stel[[j, i]]], ppp], {i, 
              Length[stel[[j]]]}]], {j, Length[stel]}];
    ff = Table[Map[Abs, Map[MinDowAltKL, ss[[i]]]], {i, Length[ss]}];
    ff1 = Table[Union[ff[[i]]], {i, Length[ff]}];
    ff1 = 
      Map[Union, 
        Table[Table[Position[ff, ff[[j, i]]], {i, Length[ff[[j]]]}], {j, 
            Length[ff]}]];
    ff1 = Sort[Map[First, Flatten[ff1, 1]]];
    ff1 = Table[ss[[ff1[[i, 1]], ff1[[i, 2]]]], {i, Length[ff1]}];
    ff1
    ]


fDerive[b_Integer, pp_String] := 
  Module[{vv, ll, ll1, t, n1, stel, ss, ff, i, ff1, ppp, n,  j},
    n = fAdjustNo[b, pp];
    vv = fRtan[n];
    ll = vv[[2]];
    ll1 = 
      Map[Flatten, 
        Table[Table[
            Table[Length[ll[[j]]] + 1, {i, ll[[j, i]]}] - i + 1, {i, 
              Length[ll[[j]]]}], {j, Length[ll]}]];
    ll = Map[Length, ll1];
    t = iteratedTake[Map[ToString, Flatten[ll1]], ll];
    ppp = StringPosition[pp, "2"];
    If[ppp != {} && StringPosition[pp, "*"] != {},
      i = 1;
      While[ppp != {} && i <= Length[ppp],
        If [ppp[[i, 1]] < StringPosition[pp, "*"][[1, 1]], ppp = Rest[ppp], 
          i++]
        ]];
    n1 = Length[ppp];
    ff = Map[Length, t];
    t = Select[t, Length[#] == n1 &];
    stel = Table[Flatten[Map[Permutations, {t[[i]]}], 1], {i, Length[t]}];
    ss = Table[
        Flatten[Table[
            StringReplacePart[pp, Flatten[stel[[j, i]]], ppp], {i, 
              Length[stel[[j]]]}]], {j, Length[stel]}];
    ff = Table[Map[Abs, Map[MinDowAltKL, ss[[i]]]], {i, Length[ss]}];
    ff1 = Table[Union[ff[[i]]], {i, Length[ff]}];
    ff1 = 
      Map[Union, 
        Table[Table[Position[ff, ff[[j, i]]], {i, Length[ff[[j]]]}], {j, 
            Length[ff]}]];
    ff1 = Sort[Map[First, Flatten[ff1, 1]]];
    ff1 = Table[ss[[ff1[[i, 1]], ff1[[i, 2]]]], {i, Length[ff1]}];
    ll = ppp;
    ss = Table[
        Table[StringTake[ff1[[j]], ll[[i]]], {i, Length[ll]}], {j, 
          Length[ff1]}];
    ll = Map[Split, Map[Reverse, Map[Sort, ss]]];
    ll1 = Map[Flatten, ll];
    ll = Table[Map[Length, ll[[i]]], {i, Length[ll]}];
    ll = Flatten[Table[Position[vv[[2]], ll[[i]]], {i, Length[ll]}]];
    ll = Table[vv[[1, ll[[i]]]], {i, Length[ll]}];
    ll1 = 
      Map[Flatten, 
        Table[Map[First, 
            Table[Position[ll1[[j]], ss[[j, i]]], 
            {i, Length[ss[[j]]]}]], {j, 
            Length[ss]}], 1];
    ll1 = 
      Flatten[Table[
          Table[Table[ll[[k, j, ll1[[k, i]]]], {i, Length[ll1[[k]]]}], {j, 
              Length[ll[[k]]]}], {k, Length[ll1]}], 1];
    ll1 = Table[StringReplacePart[pp, ll1[[i]], ppp], {i, Length[ll1]}];
    ll1
    ]

fAdjustNo[n_Integer, pp_String] := Module[{p, i, ppp},
    p = Apply[Plus, fCreatePData[pp][[1]]];
    ppp = StringPosition[pp, "2"];
    If[ppp != {} && StringPosition[pp, "*"] != {},
      i = 1;
      While[ppp != {} && i <= Length[ppp],
        If [ppp[[i, 1]] < StringPosition[pp, "*"][[1, 1]], ppp = Rest[ppp], 
          i++]
        ]];
    n - p + 2Length[ppp]]

fGenKL[b_Integer, pp_String] := Module[{kk, n, res = {}, ppp, i},
    ppp = StringPosition[pp, "2"];
    If[ppp != {} && StringPosition[pp, "*"] != {},
      i = 1;
      While[ppp != {} && i <= Length[ppp],
        If [ppp[[i, 1]] < StringPosition[pp, "*"][[1, 1]], ppp = Rest[ppp], 
          i++]
        ]];
    If[Length[ppp] == 1,
            (*Only one digon *)
            n = fAdjustNo[b, pp]; ppp = Union[Flatten[ppp]];
           kk = Select[compo[n], Not[First[#] == 1] &];
            
      kk = Map[StringJoin, Table[Table[StringJoin[ToString[kk[[j, i]]], " "],
              {i, Length[kk[[j]]]}], {j, Length[kk]}]];
            kk = Table[StringDrop[kk[[i]], -1], {i, Length[kk]}];
        Do[ n = StringReplacePart[pp, kk[[i]], {ppp[[1]], ppp[[1]]}];
        res = Append[res, n],
        {i, Length[kk]}],
      (*if it has more than one digon*)
      res = fDerive[b, pp]];
    res]



(*##################### fMinDowKnot ################*)

(*#### fMinDowKnotPr##############*)
(*# Minimizes Dowker codes for KNOTs only and 
works with single knot projection #*)
(*# Input List- Dowker Code with signs
                {{Lengths of Components},{code}}
                ConwaySymbol_ String
    Output Minimized Dowker Code#*)
(*# Needs:fGaussExt,FMinComp,fKnittDow,fFormDow #*)
fMinDowKnotPr[Unos_]:=Module[{DL,DowM},DL=fGaussExtSigns[Unos];
    If[SameQ[DL,Flatten[DL]],
      (*Knot*)DowM=fMinComp[DL];(*Print["MC ",DowM];*)
      DowM=Map[fKnittDow[DL,Take[#,-2]]&,DowM];
      (*Print["MKnitt ",DowM];*)
      DowM=Map[fFormDow[#[[1]],#[[2]]]&,DowM];
      DowM=Sort[DowM,Abs[#1]<Abs[#2]&],
      Print["Link"]];
      (* Ako ne zelimo orijentaciju *)
      If[Negative[First[DowM[[1]]]],DowM=-DowM,DowM=DowM];
  {{Length[First[DowM]]}, First[DowM]}
    ]

(*#### fMinDowKnot ##############*)
(*# Minimizes Dowker codes for KNOTs #*)
(*# Input ConwaySymbol_ String Output Minimized Dowker Code #*)
(*# Needs:fMinDowPr,fGaussExt,fMinComp,fKnittDow,fFormDow,fProjections #*)
fMinDowKnot[ConUl_String]:=Module[{ConPr,pom,m={}, \
res},ConPr=fProjections[ConUl];
    (* ConPr=Map[{Abs[fMinDowKnotPr[#]],#}&,ConPr]; *)
    ConPr=Map[{fMinDowKnotPr[#],#}&,ConPr];
    pom=Union[Abs[Flatten[Union[Map[Take[#,1]&,ConPr]],1]]];
   (*Biramo one koje imaju razlicite Dow bez obzira na projekciju*)
     Do [m=Append[m,
              Select[ConPr,SameQ[Abs[#[[1]]],pom[[i]]]&][[1]]]
         ,{i,1,Length[pom]}];
        ConPr=m;
        {Table[ConPr[[i]][[1]][[2]],{i,Length[ConPr]}],
        Table[ConPr[[i]][[2]],{i,Length[ConPr]}],
        res=First[Table[ConPr[[i]][[1]],{i,Length[ConPr]}]]};
        res={res[[1]],Flatten[res[[2]]]}
    ]



(*############### fKLinGraph  ######################## *)

(* fKLinGraph: Ulaz- lista unorderud pairs; 
  Izlaz:  neizomorfni KL u datom grafu. 
 Ako je dat graf L kao lista unordered pairs:
  1) nalazimo sve njegove podskupove duzine
   4<=d<= Length[L]; 
  2) proveravamo cetvorovalentnost- svaki broj 
  se javlja po 4 puta u Flatten delu;
  3) proveravamo planarnost izdvojenog 4-
  valentnog podgrafa PlanarQ; 
  4) 4-valentne planarne selektujemo na osnovu
   izomorfizma- biramo neizomorfne;
  5) medju njima biramo prime. *)
(*valenca temena *)

fVerVal[LL_List,el_Integer]:=Module[{res},res=Count[Flatten[LL],el];res]
(* da li je graf cetvorovalentan*)

fFourValQ[GG_List]:=Module[{l=Union[Flatten[GG]]},
    Union[Map[fVerVal[GG,#]&,l]]=={4}]
(* prvi iz klase neizomorfnih u listi grafova istih duzina*)

fClassRep[LL_List]:=Module[{l,ll,i,lP},
    ll=Map[FromUnorderedPairs[#]&,LL];lP={First[LL]};
   Print[lP];
    If[ll!={},l={First[ll]}];ShowLabeledGraph[ll[[1]]];
    ll=Rest[ll];
    While[ll!={},
      Do[ 
        If[IsomorphicQ [Last[l],ll[[i]]],ll=ReplacePart[ll,0,i]],{i,1,
          Length[ll]}];
      ll=Complement[ll,{0}];
      If[ll!={},l=Append[l,ll[[1]]];
        lP=Append[lP,ToUnorderedPairs[ll[[1]]]];
        Print[ToUnorderedPairs[ll[[1]]]];
        ShowLabeledGraph[ll[[1]]]];
      If[ll!={},ll=Rest[ll]]
      ];
    lP
    ]
(*trazi cvorove i linkove u grafu*)

fKLinGraph[UOPair_List]:=Module[{l=Flatten[UOPair],val,res,i,svi={}},
    UnOrPair=UOPair;
    val=Union[Map[fVerVal[l,#]&,l]];
    If[Select[val,#>3&]=={},
          (*ne moze se naci u manje od 3-valentnom grafu*)
              
      Print["No 4-valent subgraphs"];res=0,
            Print["Input Graph"];
      	 ShowLabeledGraph[FromUnorderedPairs[UnOrPair]];
         Print["4-valent subgraphs"];
            Do[svi=Append[svi,KSubsets[UnOrPair,i]] ,{i,4,Length[UnOrPair]}];
              
      svi=Complement[Map[Select[#, fFourValQ[Flatten[#]]&]&,svi],{{}}];
              svi=Map[fClassRep[#]&,svi];
                res=Select[svi,PlanarQ[FromUnorderedPairs[#]]&]
      ];
    res
    ]

(*############### fAddDig  ######################## *)

fAddDig[Ul_List] := Module[{vv, vv1, vvf, kk, kk1, mm, i, j, k},
    vv = Union[Ul, Map[Reverse, Ul]];
    n = Length[Union[Flatten[Ul]]];
    vv = Table[Select[vv, #[[1]] == i &], {i, n}];
    vv = Table[Map[Last, vv[[i]]], {i, Length[vv]}];
    vv = {If[Last[Union[Map[Length, vv]]] <= 4, vv, {}]};
    vv1 = If[Not[SameQ[vv, {{}}]],
        vv = 
          Union[Table[
              If[Max[Map[Length, vv[[i]]]] <= 4 , vv[[i]], {}], {i, 
                Length[vv]}]];
        vv = 
          Partition[Flatten[If[SameQ[vv[[1]], {}], Drop[vv, 1], vv], 1], n];
        (* lista svih sa valencom <= 4 *)
        vv1 = 
          Table[Union[
              Map[Sort, 
                Union[Table[
                    Table[Length[vv[[k, vv[[k, j, i]]]]], {i, 
                        Length[vv[[k, j]]]}], {j, Length[vv[[k]]]}]]]], {k, 
              Length[vv]}];
        vv = 
          Union[Table[
              If[SameQ[
                  Select[vv1[[i]], (Length[#] == 2 && MemberQ[#, 4]) && 
                        Union[Map[EvenQ, #]][[1]] &], {}], vv[[i]], {}], {i, 
                Length[vv1]}]];
        vv = Flatten[If[SameQ[vv[[1]], {}], Drop[vv, 1], vv], 1];
        vv = 
          Partition[
            Table[If[
                Length[vv[[i]]] == 2, {vv[[i, 1]], vv[[i, 1]], vv[[i, 2]], 
                  vv[[i, 2]]}, vv[[i]]], {i, Length[vv]}], n];
        (* ovde treba napuniti odgovarajuca trovalentna temena *)
        vv1 = Map[Split, Map[Sort, Map[Flatten, vv]]];
        vv1 = Map[Max, Table[Map[Length, vv1[[i]]], {i, Length[vv1]}]];
        vv = Union[Table[If[vv1[[i]] > 4, {}, vv[[i]]], {i, Length[vv1]}]];
        vv = If[SameQ[vv[[1]], {}], Drop[vv, 1], vv];
        (* svi obojivi sa dvojnim vezama samo u dvovalentnim temenima *)
        vvf = Select[vv, SameQ[Union[Map[Length, #]], {4}] &];
        (* vvf su vec obojene *)
        vv = Complement[vv, vvf];
        vv1 = 
          Table[Table[
              If[vv[[j, i, 1]] == 
                  vv[[j, i, 2]], {{vv[[j, i, 1]], i}, {vv[[j, i, 3]], 
                    i}}, {}], {i, Length[vv[[j]]]}], {j, Length[vv]}];
        vv1 = Map[Sort, Table[Flatten[vv1[[i]], 1], {i, Length[vv1]}]];
        kk = 
          Table[Complement[Range[n], Map[First, vv1[[i]]]], {i, 
              Length[vv1]}];
        kk = 
          Table[Table[{kk[[j, i]], 0}, {i, Length[kk[[j]]]}], {j, 
              Length[kk]}];
        vv1 = Table[Union[vv1[[i]], kk[[i]]], {i, Length[kk]}];
        vv1 = Map[fCorr, vv1];
        vv = 
          Table[Table[
              If[vv1[[j, i, 2]] != 0 && vv[[j, i, 1]] != vv[[j, i, 2]], 
                Insert[vv[[j, i]], vv1[[j, i, 2]], 
                  Position[vv[[j, i]], vv1[[j, i, 2]]]], vv[[j, i]]], {i, 
                n}], {j, Length[vv]}];
        vv = Union[vv, vvf];
        vvf = Select[vv, SameQ[Union[Map[Length, #]], {4}] &];
        (* vec obojene koje dodajemo na kraju *)
        vv = Complement[vv, vvf];
        kk = Table[Apply[Plus, Map[Length, vv[[i]]]], {i, Length[vv]}];
        kk = Table[Divide[Times[4, n] - kk[[i]], 2], {i, Length[kk]}];
        vv1 = 
          Table[Union[
              Flatten[Table[
                  Table[{{j, vv[[k, j, i]]}, 
                      Length[vv[[k, vv[[k, j, i]]]]]}, {i, 
                      Length[vv[[k, j]]]}], {j, Length[vv[[k]]]}], 1]], {k, 
              Length[vv]}];
        vv1 = 
          Table[Union[
              Map[Sort, Map[First, Select[vv1[[i]], #[[2]] == 3 &]]]], {i, 
              Length[vv]}];
        kk1 = Map[Split, Map[Sort, Map[Flatten, vv]]];
        kk1 = 
          Map[Union, 
            Map[Flatten, 
              Table[Select[kk1[[i]], Length[#] < 4 &], {i, Length[kk1]}]]];
        kk1 = 
          Table[Select[vv1[[i]], 
              Length[Intersection[#, kk1[[i]]]] == 2 &], {i, Length[vv]}];
        kk1 = Table[KSubsets[kk1[[i]], kk[[i]]], {i, Length[kk]}];
        kk1 = 
          Table[kk = 
              Union[Table[
                  If[SameQ[Length[Flatten[kk1[[j, i]]]], 
                      Length[Union[Flatten[kk1[[j, i]]]]]], 
                    kk1[[j, i]], {}], {i, Length[kk1[[j]]]}]];
            
            kk = If[Not[SameQ[kk, {}]], 
                If[SameQ[kk[[1]], {}], Drop[kk, 1], kk], kk],
            {j, Length[kk1]}];
        vv1 = 
          Table[Table[
              Union[kk1[[j, i]], Map[Reverse, kk1[[j, i]]]], {i, 
                Length[kk1[[j]]]}], {j, Length[kk1]}];
             mm = Table[Table[
              kk = vv1[[k, j]];
              
              kk1 = Table[
                  Position[vv[[k, kk[[i, 1]]]], kk[[i, 2]]], {i, 
                    Length[kk]}];
              
              Table[{kk[[i, 1]], 
                  Insert[vv[[k, kk[[i, 1]]]], kk[[i, 2]], kk1[[i]]]}, {i, 
                  Length[kk]}],
              {j, Length[vv1[[k]]]}], {k, Length[vv1]}];
        vv1 = 
          Table[Table[{i, vv[[j, i]]}, {i, Length[vv[[j]]]}], {j, 
              Length[vv]}];
        vv1 = 
          Table[Select[vv1[[i]], Length[#[[2]]] == 4 &], {i, Length[vv1]}];
        vv1 = 
          Table[Table[Union[mm[[j, i]], vv1[[j]]], 
          {i, Length[mm[[j]]]}], {j, 
              Length[mm]}];
        vv1 = 
          Table[Table[Map[Last, vv1[[j, i]]], {i, Length[vv1[[j]]]}], {j, 
              Length[vv1]}];
        vv1 = Flatten[Complement[vv1, {{}}], 1];
        vv1 = Union[vv1, vvf];
        vv1 = 
          Table[Table[Flatten[{i, vv1[[j, i]]}], {i, Length[vv1[[j]]]}], {j, 
              Length[vv1]}], {}];
    vv1 = 
      Table[Flatten[
          Table[Table[{vv1[[k, j, 1]], vv1[[k, j, i]]}, {i, 2, 
                Length[vv1[[k, j]]]}], {j, Length[vv1[[k]]]}], 1], {k, 
          Length[vv1]}];
    vv1 = Table[Select[vv1[[i]], #[[1]] < #[[2]] &], {i, Length[vv1]}];
    If[Or[SameQ[vv1, {}], Length[vv1] == 1], vv1, vv1 = fClassRep[vv1]];
    vv1
    ]
 
(*############### RECOGNIZE CONWAY  ######################## *)

(*Radi od graph INc*)
fFindBP[gr_List]:=Module[{bp=gr,q,d=2},
    If [Select[gr,Count[gr,#]==2&]!={},
      While[d>=2,
        bp=Sort[Flatten[FromUnorderedPairs[bp][[1]],1]];
        q=Union[Select[bp,Count[bp,#]>=2&]][[1]];
        bp=Sort[Flatten[Contract[FromUnorderedPairs[bp],q][[1]],1]];
        d=
          If[SameQ[bp,{}],0,
            Max[Union[
                Table[Length[Flatten[Position[bp,bp[[i]]]]],{i,
                    Length[bp]}]]]]]];
    If[SameQ[bp,{}],1,bp]
    ]


fBasicPoly[Ul_List]:=Module[{pdata,bp,i=1,str="*",pom,duz},
    If[SameQ[Ul[[1]],{1,1}],bp="1*",
      If[MemberQ[Union[Map[OddQ,Ul[[2]]]],True],
      pdata=fPDataFromDow[Abs[Ul]],
        pdata=Ul];
     (*  pdata=ReductionKnotLink[pdata];*)
      bp=fDowfromPD[pdata];
      bp=fGraphInc[bp][[1]];
      bp=fFindBP[bp];
      If[Not[SameQ[Head[bp],List]],
        bp="1*",(*poliedarski*)
        pom=FromUnorderedPairs[bp];
        If[Max[Flatten[bp]]\[LessEqual]12,
          While[Not[IsomorphicQ[pom,FromUnorderedPairs[GraphBP[[i,2]]]]],
            i++];
          bp=GraphBP[[i,1]],str=ToString[Max[Flatten[bp]]]<>ToString[i]<>str;
          duz=StringLength[str]-2;          
          While[Not[IsomorphicQ[pom,
          FromUnorderedPairs[fGraphInc[str][[1]]]]],
            i++;str=StringTake[str,duz]<>ToString[i]<>"*"];
          bp=str] (*kraj If[Max....]*)]];
    (*kraj od If[Not[SameQ[Head[bp],List]]*)
    bp]

    
    
    fCompositePoly[Conway_String]:=Module[{p,q,q1, i},
    p=fGraphInc[Conway][[1]];
    q=KSubsets[p,4];
    q=Drop[
        Union[Table[
            If[Length[Union[Flatten[q[[i]]]]]>5,q[[i]],{}],{i,Length[q]}]],
        1];
    q=Table[Complement[p,q[[i]]],{i,Length[q]}];
    q1=Map[ConnectedComponents,Map[FromUnorderedPairs,q]];
    q1=Flatten[Position[Table[Length[q1[[i]]],{i,Length[q1]}],2]];
    q1=If[Not[SameQ[q1,{}]],1,0]]
    
    fPolyFlype[Conway_String]:=Module[{p},
    p=fGraphInc[Conway][[1]];
    p=First[
        Union[Table[
            EdgeConnectivity[DeleteVertex[FromUnorderedPairs[p],i]],{i,
              Length[p]/2}]]];
    p=Abs[3-p];
    p]
    
    
(*################################################### *)
(* ############ f FIND CON ##########################*)
(*radi od Liste Uredjenih parovakoj aima {Konvej, {lista njegivih \
neizomorfnih projekcija}} koji su grincovi
    i dovkera elementa kako se trazi*)

fPickCon[LL_List,DD_List]:=Module[{ind=False,p,i=1,j,res},
    While[ind!=True,
                  p=LL[[i,3]];j=1;
                  While[ind!=True&&j<=Length[p],
                                ind=IsomorphicQ[FromUnorderedPairs[p[[j]]],
                                                                       
            FromUnorderedPairs[fGraphInc[DD][[1]]]];
                               j++];
               If[ind==True,res=LL[[i,2]]];
             i++];
    res]
(*radi od Doukera ili PData *)

fFindCon[Ul_]:=Module[{DD,i=1,n},
    If[SameQ[Head[Ul],String],   DD=fDToD[Ul];
      ];(*ako je Konvej*)
     If[SameQ[Head[Ul],List] ,
           If[Select[Ul[[2]],OddQ[#]&]=={},
                 (*ako je ulaz Dowker *) DD=Ul,
                                DD= fDowfromPD[ReductionKnotLink[Ul]]]];
    (*sada nam je DD bas DOWKER*)
    n=Length[DD[[2]]];
    (*broj presecnih tacaka *)
    bp= fBasicPoly[DD];
    (*bazicni poliedar *)
    (*sada nam je DD graph spreman za izomorfizam- 
        gde da ga trazimo*)
    n=ToExpression["GraphA"<>ToString[n]];
    (*suzavamo izbor na one sa kojima iam isti bazicni*)
    
    n=Select[n, #[[1]]==bp&];
    fPickCon[n,DD]
    ]
(*#######################################################*)
(*#######################################################*)
(*############### SEIFERT MATRIX AND SIGNATURE  ######################## *)

(* The author of the functions SeifertMatrix and ssmW is S.Orevkov *)
fSeifert[Ul_]:=Module[{dd,sl,pd,sei},
      If[SameQ[Head[Ul],String],pd=fCreatePData[Ul],Ul];
      If[SameQ[Head[Ul],List],
      If[Union[Map[OddQ,Abs[Ul[[2]]]]]=={False},dd=fSignsKL[Abs[Ul]];
        sl=Map[Sign,Ul[[2]]]*Map[Sign,dd[[2]]];
        pd=
          If[SameQ[Ul[[1]],
                dd[[1]]]&&(SameQ[Ul[[2]],dd[[2]]]||SameQ[Ul[[2]],-dd[[2]]]),
            Ul,fPDataFromDow[{Abs[Ul[[1]]],Abs[Ul[[2]]]*sl}]],pd=Ul],Ul];
    bword=GetBraidRep[ReductionKnotLink[pd]];
    brd=fBrdFromBword[bword];
    sei=ssmW[brd[[1]],brd[[2]]];
    sei
    ]
    
fSignat[Ul_]:=Module[{dd,sl,pd,sig},
      If[SameQ[Head[Ul],String],pd=fCreatePData[Ul],Ul];
     If[SameQ[Head[Ul],List],
      If[Union[Map[OddQ,Abs[Ul[[2]]]]]=={False},dd=fSignsKL[Abs[Ul]];
        sl=Map[Sign,Ul[[2]]]*Map[Sign,dd[[2]]];
        pd=
          If[SameQ[Ul[[1]],
                dd[[1]]]&&(SameQ[Ul[[2]],dd[[2]]]||SameQ[Ul[[2]],-dd[[2]]]),
            Ul,fPDataFromDow[{Abs[Ul[[1]]],Abs[Ul[[2]]]*sl}]],pd=Ul],Ul];
    bword=GetBraidRep[ReductionKnotLink[pd]];
    brd=fBrdFromBword[bword];
    sig=ssmW[brd[[1]],brd[[2]]]; 
    sig=fSignatureKL[sig]; 
    sig
    ]
fBrdFromBword[bword_String]:=Module[{m,brd, i},
    m=Length[Union[ToCharacterCode[ToUpperCase[bword]]]]+1;
    brd=ToCharacterCode[bword]-64;
    brd=Table[
        If[brd[[i]]>26,brd[[i]]=-Mod[brd[[i]],32],brd[[i]]],{i,Length[brd]}];
    {m,brd}
    ]
SeifertMatrix[m_Integer,brd_List]:=
    Module[{n,e,V,X,q,c,i,j,h,a,b},
      a={{0,1,-1,0},{-1,0,1,0},{0,0,0,0},{1,-1,0,0}};
      b={{-1,1,0,0},{1,-1,0,0},{0,0,0,0},{0,0,0,0}};
      n=Length[brd];V=Table[0,{i,n},{j,n}];X=Table[{n},{h,m-1}];
      Do[h=Abs[brd[[q]]];e=Sign[brd[[q]]];
        c[1]=X[[h,1]];X[[h]]={c[2]=q};
        c[3]=If[h<m-1,X[[h+1,1]],n];c[4]=If[h>1,X[[h-1,1]],n];
        Do[Do[V[[c[i],c[j]]]+=a[[i,j]]+e*b[[i,j]],{i,4}],{j,4}],{q,n}];
      V=Delete[V,X];
      If[Length[V]>0,V=Transpose[Delete[Transpose[V],X]]/2];
      V];
      (*popravljena 15.08.2003 *)

ssmW[m_Integer,brd_List]:=
    Module[{bq,n,d,e,V,X,p=1,q,r=1,c,i,j,h,a,b},
      a={{0,0,-1,1,2},{0,0,1,-1,-2},{-1,1,0,0,0},{1,-1,0,0,0},{2,-2,0,0,0}};
      b={{-2,2,0,0},{2,-2,0,0},{0,0,0,0},{0,0,0,0}};
      d=n=Length[brd];Do[If[Not[IntegerQ[brd[[q]]]],d++;p++],{q,n}];
      V=Table[0,{i,d},{j,d}];X=Table[{d},{h,m-1}];
      Do[bq=brd[[q]];
        If[IntegerQ[bq],h=Abs[bq];e=Sign[bq],h=Abs[bq[[2]]];
          V[[r,r]]=2*bq[[1]]*Sign[bq[[2]]];c[5]=r++;];
        c[1]=X[[h,1]];X[[h]]={c[2]=p++};
        c[3]=If[h<m-1,X[[h+1,1]],d];c[4]=If[h>1,X[[h-1,1]],d];
        If[IntegerQ[bq],
          Do[Do[V[[c[i],c[j]]]+=a[[i,j]]+e*b[[i,j]],{i,4}],{j,4}],
          Do[Do[V[[c[i],c[j]]]+=a[[i,j]],{i,5}],{j,5}]],{q,n}];
       V=Delete[V,X];
      If[Length[V]>0,V=Transpose[Delete[Transpose[V],X]]/2];
      V];
      (*popravljena 15.08.2003 *)
      
(*fSignatureKL[S_List]calculates signature from a Seifert matrix S*)

fSignatureKL[S_List]:=Module[{s},s=S+Transpose[S];
    s=Abs[Apply[Plus,Map[Sign,Eigenvalues[N[S]]]]]]

(*##########################################*)
(*##########################################*)
(*############### SYMMETRY  ######################## *)
(*SymmKL[Ul] racuna automorfizme projekcije KL koji cuvaju znakove*)

Symm[Ul_]:= Module[{G,k,i,j, sG, Aut},
    G = fGaussExt[Ul];
    sG = fGrInc[G][[2]];
    (*Print[sG];*)G = FromUnorderedPairs[fGrInc[G][[1]]];
    Aut = Automorphisms[G];
    (*Print[Aut];*)(*ovde izdvajamo automorfizme koji cuvaju znak*)
    Aut = 
      If[Aut != {}, 
        If[First[
              Union[Table[
                  If[SameQ[
                        Sign[Table[
                              Table[Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]], {j, 
                                  Length[sG]}], {i, Length[Aut]}]][[k]], 
                        sG] == True, Aut[[k]], {}], {k, Length[Aut]}]]] == \
{},
           Drop[Union[
              Table[If[
                  SameQ[Sign[
                          Table[Table[
                              Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]], {j, 
                    Length[sG]}], {i, Length[Aut]}]][[k]], sG] == 
                    True, Aut[[k]], {}], {k, Length[Aut]}]], 1], 
          Union[Table[
              If[SameQ[
                    Sign[Table[
                          Table[Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]], {j, 
                              Length[sG]}], {i, Length[Aut]}]][[k]], sG] == 
                  True, Aut[[k]], {}], {k, Length[Aut]}]]], Aut];
       {Aut,If[Aut != {}, Map[ToCycles, Aut], Aut],Length[Aut]}
    ]

(*############### SYMMETRY  ######################## *)
(* nalazi max simetricnu projekciju alt KL od konveja*)

MaxSymmProjAltKL[Con_String]:=Module[{l},
    l=fDiffProjectionsAltKL[Con];
    l=Map[#[[1]]&,l];
    l=Map[{Length[Symm[#][[1]]],#}&,l];
    l=Reverse[Sort[l]][[1]];
    Print["Number of automorphisms: ",l[[1]]];
    Print["Max symmetric projection: ",l[[2]]];
    ShowKnotfromPdata[fCreatePData[l[[2]]]];
    ]

(*###################################################*)
(*############### UINFTYNO  ######################## *)

(*pomocna za fMirorDve *)
(*Revertuje znake jednostrukih elemenata*)
UnsortedComplement[x_List,y__List]:=
  Replace[x,Dispatch[(#:>Sequence[])&/@Union[y]],1]
fAddSign[LL_List,{el_Integer,elS_Integer}]:=Module[{z,i,p,l=LL},z=1-elS;
    p=Position[LL,{el,0}][[1,1]];
    Do[l=ReplacePart[l,{l[[i,1]],z},i];
      z=1-z,{i,p,Length[LL]}];
    z=1-elS;
    Do[z=1-z;(*Print[l];*)l=ReplacePart[l,{l[[p-i,1]],z},p-i],{i,1,p-1}];
    (*Print["rezultat iz fADDSIGN ",l];*)l]
fGiveSign[LL_List,SL_List]:=
  Module[{p,i,l=LL,saznakom=SL,sz={},sa},If[SL=={},p=l[[1]];
      Do[p=ReplacePart[p,{p[[i,1]],Mod[i,2]},i];,{i,1,Length[p]}];
      l=ReplacePart[l,Prepend[p,1],1];
      sz=p,(*ako vec imamo neke sa znacimo trazimo ih u drugim komponentama*)
        sa=Union[Flatten[Map[Take[#,1]&,saznakom]]];
      (*Print["SSSSS",sa,saznakom];*)Do[
        If[Length[l[[i,1]]]!=0,(*znaci da ona nije oznacena*)presek=
            Select[l[[i]],MemberQ[sa,#[[1]]]&];
          
          If[presek!={},
            presek=Select[saznakom,#[[1]]==presek[[1,1]]&][[1]];
            p=fAddSign[l[[i]],presek];
            l=ReplacePart[l,Prepend[p,1],i];
            p=Map[If[MemberQ[sa,#[[1]]],-1,#]&,p];
            sz=Union[sz,Complement[p,{-1}]]]],{i,1,Length[l]}]];
    {l,sz}]
(*Pravi alternirajuci kad se polazi od pData...*)

fMakeAlt[CoLe_List,LL_List]:=
  Module[{i,pom=0,CoPom={},saznakom={},ind=True},l=iteratedTake[LL,CoLe];
    While[ind,(*Print["SZ",saznakom];*)saznakom=fGiveSign[l,saznakom];
      (*Print["posle fGiveSigna",saznakom[[1]],"\n \n",saznakom[[2]],
            "\n \n"];*)l=saznakom[[1]];
      saznakom=Last[saznakom];
      (*Print[Map[Length[#[[1]]]&,l]];*)ind=
        Union[Map[Length[#[[1]]]&,l]]!={0}];
    Map[Drop[#,1]&,l]]
(*pomocna za fMakeGaussPD-menja odnose iznad ispod za neparne PData*)
checkArgs[s_,t_]:=
  If[Head[s]===List&&VectorQ[t,Head[#]===Integer&&#>=0&]&&
      Plus@@t<=Length[s],True,False]

iteratedTake[s_,t_]/;checkArgs[s,t]:=
  First/@Rest[FoldList[Through[{Take,Drop}[#1[[2]],#2]]&,{{},s},t]]
fParno[LL_List]:=Module[{l,lGExt,i,p},l=LL;
    lGExt=Map[Flatten,l];
    lGExt=Map[Take[#,{1,Length[#],2}]&,lGExt];
    (*extGaus bez odnosa iznad ispod*)Do[
      p=Flatten[Position[Flatten[Take[lGExt,i-1]],lGExt[[i,1]]]];
      (*pozicija prvog u sledecoj komponenti u prethodnom delu gausa*)If[
        p!={},
        If[Not[EvenQ[p[[1]]]],
          lGExt=ReplacePart[lGExt,RotateLeft[lGExt[[i]]],i];
          l=ReplacePart[l,RotateLeft[l[[i]]],i]]],{i,2,Length[l]}];
    l]
(*ako se posle secenja za PD u istoj komponenti 
jave dva ista uzastopna broja 
(ciklicno) sa suprotnim drugim znakom npr.{-18,1},{-18,
      0} izbacujemo ih.*)
(*radi od gaus ext sa iznad ispod*)
fIzbaciPetlju[GK_List]:=Module[{p,l,pom,res=GK},l=Map[Take[#,1]&,GK];
    el=Union[Select[l,Count[l,#]==2&]];
    p=el;
    If[p!={},p=Map[Union[Flatten[Position[l,#]]]&,p];
      l=Map[If[#!={1,Length[GK]},RotateLeft[l,#[[2]]-1],l]&,p];
      p={};
      Do[p=Append[p,Flatten[Position[l[[i]],el[[i]]]]],{i,1,Length[el]}];
      l=GK;
      Do[If[p[[i]]=={1,Length[GK]},
          l=Select[l,#[[1]]!=el[[i,1]]&]],{i,1,Length[el]}];res=l];
    res]
(*############### UINFTYNO ########################*)(*pomocna za \
fMirorDve*)(*Revertuje znake jednostrukih elemenata*)
  fRevertSigns[ll_List,rev_List]:=Module[{p=rev,pom=ll},p=Map[Take[#,1]&,p];
    p=Flatten[Select[p,Count[p,#]==1&]];
    pom=Map[If[MemberQ[p,#[[1]]],{-#[[1]],#[[2]]},#]&,pom]]
(*pomocna za fMirorComp,fMiror1*)
(*izbaci jednu dvojnu tacku i da gaus EXt \
sa iznad ispod bez nje*)
(*radi od gausEXt i rednog broja komponent koju \
menjamoo i el koji izbacijemo*)

fMirorDve[ExtGaus_List,redni_Integer,el_Integer]:=
  Module[{pom=ExtGaus,len,komp,pos,p},komp=pom[[redni]];
    pos=Position[Map[Take[#,1]&,komp],el];
    pos=Map[#[[1]]&,pos];
    (*pozicije na kojima se dvojni nalazi*)p=
      Reverse[Take[komp,{pos[[1]]+1,pos[[2]]-1}]];
    (*p sadrzi deo komponente koji se revertovao*)komp=
      Join[Take[komp,pos[[1]]-1],p,Drop[komp,pos[[2]]]];
    pom=ReplacePart[pom,komp,redni];
    len=Map[Length[#]&,pom];
    pom=Select[Flatten[pom,1],#[[1]]!=el&];
    (*moze i bez prenumeracije*)pom=fRevertSigns[pom,p];
    pom=iteratedTake[pom,len];
    (*Print["PMiror DVE Pre Izbaci petlju  ",pom];
      pom=Map[fIzbaciPetlju[#]&,pom];Print["posle petlji ",pom];*)pom]
(*pomocna za fDvojneQ*)
(*vraca True ako ima samo jednostruke elemente*)

fDvojneKomp[Ulaz_List]:=Module[{res},
    res=Ulaz;
    res=Union[Map[Count[res,#]&,res]];
    If[res=={1},res=True,res=False];
    res](*pomocna za fMakeGaussPD-menja odnose iznad ispod za neparne PData*)
fMenjaObePojavePD[LL_List,el_Integer]:=
  Module[{p,i,l=LL},p=Select[l,#[[1]]==el&];
    (*Print[p,p[[1,1]]];*)(*Print[Position[l,p[[1,1]]]];*)p=
      Map[Take[#,1][[1]]&,Position[l,p[[1,1]]]];
    Do[l=ReplacePart[l,{l[[p[[i]],1]],1-l[[p[[i]],2]]},p[[i]]],{i,1,2}];
    l]
(*radi od Ext Gausa sa odnosima iznad ispod*)
(*ako moze sredi parnost*)
(*NE \
RADI RETROAKTIVNO!!!!*)

fMakeGaussPD[PD_List]:=Module[{l,n,i,l1={},nepar},n=Length[PD[[2]]];
    l=Complement[Range[2n],Abs[PD[[2]]]];
    (*svi brojevi do 2 duzine koji vec nisu u pData*)Do[
      l=ReplacePart[l,{l[[i]]*Sign[PD[[2,i]]],PD[[2,i]]},i],{i,1,n}];
    (*Print["www",l];*)(*uredjeni parovi sa sve znacima*)nepar=
      Select[l,Or[OddQ[#[[2]]],And[EvenQ[#[[1]]],EvenQ[#[[2]]]]]&];
    nepar=Flatten[Map[Take[#,-1]&,nepar]];
    (*Print[nepar];*)l1=Map[#-#&,Range[2n]];
    Do[l1=ReplacePart[l1,{l[[i,2]],0},Abs[l[[i,1]]]];
      l1=ReplacePart[l1,{l[[i,2]],0},Abs[l[[i,2]]]],{i,1,n}];
    (*Print["???",2PD[[1]],
          l1];*)(*ubacili smo indikatore i sad treba da ga napravimo da je \
alternirajuci*)l1=fMakeAlt[2PD[[1]],l1];
    (*menjamo odnose onima koji su u pdata bili neparni a to su i ovde*)l1=
      Flatten[l1,1];
    (*Print[" pre kvarenja ",l1];*)If[nepar!={},
      Do[l1=fMenjaObePojavePD[l1,nepar[[i]]],{i,1,Length[nepar]}]];
    (*Print[" posle kvarenja ",l1];*)l1]
(*pomocna za fMirore*)
(*radi od podeljenog Gauss Ext sa odnosima iznad-ispod*)
fMakeGaussPD[PD_List]:=Module[{l,n,i,l1={},nepar},n=Length[PD[[2]]];
    l=Complement[Range[2n],Abs[PD[[2]]]];
    (*svi brojevi do 2 duzine koji vec nisu u pData*)Do[
      l=ReplacePart[l,{l[[i]]*Sign[PD[[2,i]]],PD[[2,i]]},i],{i,1,n}];
    (*Print["www",l];*)(*uredjeni parovi sa sve znacima*)nepar=
      Select[l,Or[OddQ[#[[2]]],And[EvenQ[#[[1]]],EvenQ[#[[2]]]]]&];
    nepar=Flatten[Map[Take[#,-1]&,nepar]];
    (*Print[nepar];*)l1=Map[#-#&,Range[2n]];
    Do[l1=ReplacePart[l1,{l[[i,2]],0},Abs[l[[i,1]]]];
      l1=ReplacePart[l1,{l[[i,2]],0},Abs[l[[i,2]]]],{i,1,n}];
    (*Print["???",2PD[[1]],
          l1];*)(*ubacili smo indikatore i sad treba da ga napravimo da je \
alternirajuci*)l1=fMakeAlt[2PD[[1]],l1];
    (*menjamo odnose onima koji su u pdata bili neparni a to su i ovde*)l1=
      Flatten[l1,1];
    (*Print[" pre kvarenja ",l1];*)If[nepar!={},
      Do[l1=fMenjaObePojavePD[l1,nepar[[i]]],{i,1,Length[nepar]}]];
    (*Print[" posle kvarenja ",l1];*)l1]
(*pomocna za fBrCo-izbacuje iz Gaus parova i zavisnosti od "ociscenog"Gausa*)
fDvojneQ[GausExt_List]:=Module[{res=GausExt,len},len=Map[Length[#]&,res];
    res=Flatten[Take[Flatten[res],{1,Length[Flatten[res]],2}]];
    res=iteratedTake[res,len];
    res=Map[fDvojneKomp[#]&,res];
    If[Union[res]=={True},res=True,res=False];
    res]
(*vraca True ako nijedna komponenta nema dvojnih*)
(*pomocna za fMirorComp*)
\
(*daje nove PData kad je izbacena 1 dvostruka*)
fMiror1[Ul_,redni_Integer,el2_Integer]:=
  Module[{pom,ind},pom=fMirorDve[Ul,redni,el2];
    (*sredjujemo parnost-
        ako ne moze da sredi dobicemo dva neparna i dva parna uparena*)If[
      pom=={{}},(*posle svega i izbacivanja ptelji doboili smo prazno-
          kraj*)pom={-1},pom=fParno[pom];
      pom=fPDataCuttNo[pom];
      pom=ReductionKnotLink[pom];
      If[MemberQ[{{{},{}},{{0},{}},{{},{0}}},pom]||
          fDvojneQ[iteratedTake[fMakeGaussPD[pom],2pom[[1]]]],pom={-1}]];
    (*ako u ovom razvezivanju nema dvojnih ili je unknott 
    onda je to kraj i 
vraca {-1}*)(*inace vraca neke pdata*)pom]
(*jedna komponenta->
    lista rezultata*)
(*Radi od pdata ili Konveja i broja komponente cije \
cemo dvostruke da izbacujemo*)
(*radi samo za komponente koje imaju dvojnih*)

(*vraca listu rezultata*)

fMirorKomp[Ul_,k_Integer]:=
  Module[{pd,pd0,pd1,p,i,dvojne,res={},brzi=True},
    If[SameQ[Head[Ul],String],pd=fCreatePData[Ul],pd=Ul];
    pd=ReductionKnotLink[pd];
    (*sad moramo da napravimo odgovarajuceg Gausa i da ga iz \
alterniramo*)pd0=
      fMakeGaussPD[pd];
    (*Gaus parovi sa iznad-ispod tj.0-1*)pd1=
      Flatten[Map[Take[#,1]&,pd0]];(*Gaus bez iznad ispod*)pd0=
      iteratedTake[pd0,2pd[[1]]];
    (*gaus za seckanje podeljen na komponente*)pd1=
      iteratedTake[pd1,2pd[[1]]];
    p=pd0[[k]];(*k-ta tj.komponenta koju secemo*)p1=
      Flatten[Map[Take[#,1]&,p]];(*p bez iznad ispod*)(*Print["Gaus",p1];*)
      dvojne=Union[Select[p1,Count[p1,#]==2&]];
    If[fDvojneQ[pd0],res={},Do[If[brzi,brzi=fMiror1[pd0,k,dvojne[[i]]];
          res=Append[res,brzi];
          If[brzi=={{1}},brzi=False,brzi=True]],{i,1,Length[dvojne]}];
      res];
    If[MemberQ[res,{{-1}}],res={-1}];
    res]
(*ako je res={-1} vise ne treba da radimo u suprotnom imamo pdata nekih \
clanova nove generacije*)

SamePdata[LL_List]:=
  Module[{l,res={LL[[1]]}},If[Length[LL]!=1,l=Rest[LL];
      While[l!={},l=Select[l,Not[SamePDQ[Last[res],#]]&];
        If[l!={},res=Append[res,First[l]];
          If[Length[l]==1,l={},l=Rest[l]]]]];
    res]
(*kad primrnimo fMirorKomp na sve komponente i to skupimo i jednu listu*)
\
(*to je lista pdata-pravimo gause-
    ako ima neki bez dvojnih to je kraj*)
(*radi od pdata*)

MirorAllKomp[Ul_]:=
  Module[{i=1,res={},p={2}},
    If[fComponentNo[Ul]==1,res=fMirorKomp[Ul,1],
      While[p!={-1}&&i<=fComponentNo[Ul],
        p=fMirorKomp[Ul,i];
        (*Print["Sredjena je ",i,"-ta komponenta: ",p];*)If[
          p!={-1}&&
            p!={{-1}},(*ako je lista pdata dodamo je u nove rez*)
            i++;
          If[p!={},res=Append[res,p]]
          (*ako je {} onda ne doprinosi rezultatu*),(*ako je {-1} ne treba \
vise da radimo*)p={-1};res={-1}]]];
    If[SameQ[res,{{}}]||SameQ[res,{}],res={-1}];
    (*Print[res];*)If[Not[SameQ[res,{-1}]],
      If[Length[res]==1,res=res[[1]]];
      res=SamePdata[res];
      res={res}];
    res]
(*vraca {-1} ako je razvezao a u suprotnom vraca listu pdata*)
fNSc[pd_List]:=Module[{NscNo=0,ind=True,gen={}},
    If[fDvojneQ[iteratedTake[fMakeGaussPD[pd],2pd[[1]]]],
      (*nema dvojnih*)(*Print["odma"];*)NscNo=0,
      (*ima dvojnih*)gen={pd};
      While[ind==True,NscNo++;
        gen=Flatten[Map[MirorAllKomp[#]&,gen],1];
        If[gen!={},gen=Flatten[gen,1]];
        ind=
          Not[MemberQ[gen,{-1}]||MemberQ[gen,-1]||SameQ[gen,{-1}]||
              SameQ[gen,{}]||SameQ[gen,{{}}]]]];
    
    NscNo]

(*###########################################*)
NoSelfCrossNo[Ul_]:=Module[{pd},
    If[SameQ[Head[Ul],String],pd=fCreatePData[Ul],pd=Ul];
    (*sad imamo pdata-odmah reduction*)pd=ReductionKnotLink[pd];
    Min[fNSc[pd],fNSc[GetMirrorImageKnot[pd]]]
    ]
(*###########################################*)
(*###########################################*)
(* ###################  SAME KNOT LINK #########################*)
(*uporedjuje dva KL data Konvejima *)

SameAltConKL[Con1_String,Con2_String]:=
  Module[{d1,d2,res},
    d1=MinDowAltKL[Con1];
      d2=MinDowAltKL[Con2];
    If[SameQ[Abs[d1],Abs[d2]],
          res=1 (*same*),
        res=0 (*not same*)];
    res]

(*uporedjuje dva KL data Konvejima *)

SameAltProjKL[Con1_,Con2_]:=
  Module[{d1,d2,res},
    d1=MinDowProjAltKL[Con1];
    d2=MinDowProjAltKL[Con2];
    If[SameQ[Abs[d1],Abs[d2]],
          res=1 (*same*),
        res=0 (*not same*)];
    res]
    
(*################################################### *)

(*################## GETBRAIDREPRESENT ################### *)



fGaussExtSignsBraid[Ulaz_] := Module[{Ul = Ulaz, n, pp, pp1, 
pp2, ppr, ppi, i},
    (*Input : Pdata*)
    If[SameQ[Head[Ul], List], If[MemberQ[Map[OddQ, Abs[Ul[[2]]]], True],
        Ul], 
      If[SameQ[Head[Ul], String], Ul = fCreatePData[Ul], 
        Ul = fPDataFromDow[Ul]]];
    (*sad nam je Ul iliKonvej ili dowker*)
    n = Apply[Plus, Ul[[1]]];
    pp = Complement[Range[Apply[Plus, 2*Ul[[1]]]], Map[Abs, Ul[[2]]]];
    pp1 = Map[Sign, Ul[[2]]];
    pp2 = Table[If[EvenQ[Ul[[2, i]]], 1, -1], {i, n}];
    pp = Table[{pp[[i]], Abs[Ul[[2, i]]]}, {i, n}];
    ppr = Map[Reverse, pp];
    pp = Union[Table[{pp[[i]], {pp1[[i]], pp2[[i]]}}, {i, n}], 
        Table[{ppr[[i]], {pp1[[i]], pp2[[i]]}}, {i, n}]];
    pp1 = Map[First, pp];
    pp = Table[
        If[Map[First, Position[pp1, pp1[[i, 2]]]][[1]] == pp[[i, 1, 1]], 
          pp[[i]] , {Reverse[pp[[i]][[1]]], pp[[i]][[2]]}], {i, 2n}];
    pp = Table[{pp[[i, 1, 1]], pp[[i, 2]]}, {i, 2n}];
    ppr = Union[Map[First, pp]];
    pp = Table[{Flatten[Position[ppr, pp[[i, 1]]]][[1]], pp[[i, 2]]}, {i, 
          2n}];
    ppr = Map[Last, Union[Table[{pp[[i, 1]], pp[[i, 2, 1]]}, {i, 2n}]]];
    ppi = Table[-pp[[i, 2, 2]]*(-1)^i, {i, 2n}];
    pp = Map[First, pp];
    pp = Table[pp[[i]]*ppi[[i]], {i, 2n}];
    pp = {iteratedTake[pp, 2Ul[[1]]], ppr};
    pp
    ]
    
 fGetBraidRepresent[Ul_] := Module[{f = "br.txt", pp = Ul, pp1, pp2},
    pp = fGaussExtSignsBraid[pp];
    pp1 = ToString[pp[[1]]];
    pp2 = ToString[pp[[2]]];
    pp1 = 
      StringReplace[
        StringReplace[
          pp1, {"}, {" -> ";", "," -> "", "{" -> "", "}" -> ""} ], {";" -> 
            ", "}];
    pp2 = 
      StringReplace[
        pp2, {"1" -> "+", "-1" -> "-", "{" -> "", "}" -> "", "," -> " "}];
    pp = pp1 <> " / " <> pp2;
    OpenWrite[f];
    WriteString[f, pp];
    Write[f];
    Close[f];
    Run["fromdos" <> " br.txt"];
    Run["braid -v" <> " br.txt" <> " brout.txt"];
    pp = Import["brout.txt"];
    pp1 = StringPosition[pp, "="][[1, 1]];
    pp2 = StringLength[pp];
    pp = StringTake[pp, {pp1 + 1, pp2 - 2}]
    ]  
    
    

(*################## BraidWReduce###################*)
(*redukuje braid words pomocu Herves Iberovog programa
hr koji ucitavamo*)

BraidWReduced[PD_List]:=Module[{pp}, 
    pp=ReductionKnotLink[PD];
    pp=GetBraidRep[pp];
    pp
        ]

(*###################################################*)
(*##################fClassicToCon ###################*)
fClassicToCon[Ulaz_String]:=Module[{brPres,pos,pom},
    pos=StringPosition[Ulaz,"_"][[1,1]];
    brPres=ToExpression[StringTake[Ulaz,pos-1]];
    (*Print[Head[brPres],Select[CCtoC,SameQ[#[[1]],10]&]];*)
    
    pom=Flatten[Select[CCtoC,SameQ[#[[1]],brPres]&],1];
    pom=Rest[pom];
    (*Print[pom];*)
    Select[pom, SameQ[#[[1]],Ulaz]&][[1,2]]
    
    ]
(*################################################### *)
(*################### EDMONDS############*)



fGrCompl[Ul_List] := Module[{res,i},
    res = Union[Ul, Map[Reverse, Ul]];
    res = Table[Select[res, SameQ[#[[1]], i] &], {i, Max[Flatten[Ul]]}];
    res
    ]

fSelParts[LL_List] := Module[{res, i, j, L},
    res = 
      Flatten[Table[
          Join[{LL[[1, i]]}, {LL[[2, j]]}], {i, Length[LL[[1]]]}, {j, 
            Length[LL[[2]]]}], 1];
    Do[L = {res, LL[[i]]}; 
      res = Flatten[
          Table[Join[L[[1, i]], {L[[2, j]]}], {i, Length[L[[1]]]}, {j, 
              Length[L[[2]]]}], 1]; i = i + 1, {i, 3, Length[LL]}];
    res
    ]

fAllRot[Ul_List] := Module[{gg,i,j},
    gg = fGrCompl[Ul];
    gg = Table[
        Table[Join[{gg[[j, 1]]}, 
            DistinctPermutations[
                Table[gg[[j, i]], {i, 2, Length[gg[[j]]]}]][[i]]], {i, 
            Length[DistinctPermutations[
                Table[gg[[j, i]], {i, 2, Length[gg[[j]]]}]]]}], {j, 
          Length[gg]}];
    gg = fSelParts[gg];
    gg
    ]

fFindCyc[UlL_List, poc_List] := Module[{pp, pom, uu, uu1, vv1,vv, i},
    pom = {poc};
    uu = pom[[1]];
    While[Not[SameQ[uu1, pom[[1]]]], vv = Reverse[uu];
      vv1 = Position[UlL, Reverse[uu]];
      uu = 
        If[vv1[[1, 2]] + 1 <= Length[UlL[[1]]], 
          UlL[[vv1[[1, 1]], vv1[[1, 2]] + 1]], UlL[[vv1[[1, 1]], 1]]];
      uu1 = uu;
      pom = Join[pom, {uu1}]]; pp = Complement[Flatten[UlL, 1], pom];
    pom = Flatten[Map[#[[1]] &, Drop[pom, -1]]];
    If[pp != {}, 
      pp = Table[Select[pp, SameQ[#[[1]], i] &], {i, Max[Flatten[pp]]}]];
    {pom, pp}
    ]

fFindCycAll[UlL_List, poc_List] := 
  Module[{res = {}, ind = True, komplementi = {}, novpoc = poc, pom},
    While[Or[ind, komplementi != {}], ind = False;
      pom = fFindCyc[UlL, novpoc];
      If[poc == novpoc, komplementi = Flatten[pom[[2]], 1], 
        komplementi = Intersection[komplementi, Flatten[pom[[2]], 1]]];
      If[komplementi != {}, novpoc = komplementi[[1]];
        komplementi = Drop[komplementi, 1]];
      res = Append[res, pom[[1]]]];
     res]

fEdmonds[UlL_List] := Module[{gg, res = {}, i,j,rr1,rr2},
    gg = fAllRot[UlL];
    Do[res = Append[res, fFindCycAll[gg[[i]], First[gg[[i, 1]]]]]  , {i, 1, 
        Length[gg]}];
    res = Table[{gg[[i]], res[[i]]}, {i, Length[gg]}];
    rr1 = 
      Table[{Sort[Map[Length, res[[i, 2]]]], 
          Sort[Map[Length, Map[Split, Map[Sort, res[[i, 2]]]]]]}, {i, 
          Length[res]}];
     rr2 = Union[rr1];
     rr2 = 
      Flatten[Map[First, Table[Position[rr1, rr2[[i]]], {i, Length[rr2]}]]];
    res = Table[res[[rr2[[i]]]], {i, Length[rr2]}];
    res = 
      Table[{Table[Map[Last, res[[j, 1, i]]], {i, Length[res[[j, 1]]]}], 
          res[[j, 2]], 
          Max[Flatten[UlL]] + Length[res[[j, 2]]] - Length[UlL], 
          Divide[2 - Max[Flatten[UlL]] - Length[res[[j, 2]]] + Length[UlL], 
            2]}, {j, Length[res]}];
    res]


    
(*###################################################*)
(*####### STELAR #######################*)


    fStellarBasic[n_Integer] := Module[{del, del1, i, j, k},
    If[n < 6, del = {}, del = compo[n];
      del = Select[del, Not[MemberQ[#, 1]] &];
      del = Select[del, Length[#] >= 3 &];
      del1 = Sort[Table[{Sort[del[[i]]], del[[i]]}, {i, Length[del]}]];
      del = Union[Table[Sort[del[[i]]], {i, Length[del]}]];
      del = Table[Select[del1, #[[1]] == del[[i]] &], {i, Length[del]}];
      del = 
        Table[Table[del[[j, i, 2]], {i, Length[del[[j]]]}], {j, 
            Length[del]}];
      del = 
        Map[Last, 
          Flatten[Table[
              Union[Table[
                  Union[Table[
                      RotateLeft[del[[k, j]], i], {i, Length[del[[k, j]]]}], 
                    Map[Reverse, 
                      Table[RotateLeft[del[[k, j]], i], {i, 
                          Length[del[[k, j]]]}]]], {j, 
                    Length[del[[k]]]}]], {k, Length[del]}], 1]]];
    del = {del};
    del = 
      Flatten[Table[
          Map[StringReplace[ToString[#], 
          {"{" -> "", "}" -> "", " " -> ""}]&,
             del[[i]]], {i, Length[del]}]];
    del]
(*Changed 9.03.2005*)


   
fMakePart[LL_List]:=Module[{res,pom,i=1},
    res=Map[compo[#]&,LL];
    While[i<=Length[res],
                pom=Select[res[[i]],First[#]!=1&];
      pom=Flatten[Map[DistinctPermutations,pom],1];
      pom=Select[pom, First[#]!=1&];
                
      If[pom=={},res=Drop[res,{i}],res=ReplacePart[res,pom,i]];
      i++
      ];
    Do[  pom=Map[StringReplace[ToString[#],", "->" "]&,res[[i]]];
             res=ReplacePart[res,pom,i]
      ,{i,1,Length[res]}];
    res
    ]
    
fCombine[res_List,NewL_List]:=Module[{i,res1,j,pom,combRes={}},
    res1=res;
    Do[
         Do[   
               pom=Join[res1[[j]],{NewL[[i]]}];
                combRes=Append[combRes,pom] ,
             {j,1,Length[res]}],
      {i,1,Length[NewL]}];
    combRes
    ]
    
fPerm[LL_List]:=Module[{i,perm={},res={}},
    perm=Map[{#}&,LL[[1]]];
    Do[perm=fCombine[perm,LL[[i]]],{i,2,Length[LL]}];
    perm
    ]

fStellar[n_Integer] := Module[{stel, i, j},
    stel = fStellarBasic[n];
    stel = 
      Map[ReadList[StringToStream[StringReplace[#, "," -> " "]], Number] &, 
        stel];
    stel = Map[fMakePart[#] &, stel];(*svaki strem veci od 1 partitionira :*)

    
    stel = Flatten[Map[fPerm[#] &, stel], 1];
    stel = 
      Table[Table[
          ToExpression[StringReplace[stel[[j, i]], {" " -> ","}]], {i, 
            Length[stel[[j]]]}], {j, Length[stel]}];
    stel = 
      Map[Last, 
        Union[Table[
            Union[Table[RotateLeft[stel[[j]], i], {i, Length[stel[[j]]]}], 
              Map[Reverse, 
                Table[RotateLeft[stel[[j]], i], 
                {i, Length[stel[[j]]]}]]], {j,Length[stel]}]]];
    stel = 
      Table[Map[
          StringReplace[ToString[#], {"{" -> "", "}" -> "", " " -> ""}] &, 
          stel[[i]]], {i, Length[stel]}];
    stel = Table[StringReplace[stel[[i]], "," -> " "], {i, Length[stel]}];
    stel = 
      Table[StringReplace[
          ToString[stel[[i]]], {"{" -> "", "}" -> "", ", " -> ","}], {i, 
          Length[stel]}];
    stel
    ]
    (*Changed 9.03.2005*)
    
    
    fStellarPlus[n_Integer]:=Module[{st, i, j},
    st=If[n<7,Print["n>6"],
        st=Table[fStellar[i],{i,6,n-1}];
        st=Flatten[
            Table[Table[
                StringJoin[st[[i,j]],"+",ToString[n-5-i]],{j,
                  Length[st[[i]]]}],{i,Length[st]}]]];
    st
    ]
    
    (*####### STELARNALT #######################*)
    
    fStelString[Ul_String] := Module[{ss1, ss2, sp, ss = Ul, i},
    ss1 = Length[Union[Flatten[StringPosition[ss, ","]]]];
    ss2 = Length[Union[Flatten[StringPosition[ss, " "]]]];
    sp = Divide[Union[Flatten[StringPosition[ss, ","]]], 2];
    ss = Flatten[
        Map[ReadList[StringToStream[StringReplace[#, "," -> " "]], 
              Number] &, {ss}]];
    ss = iteratedTake[ss, 
        Flatten[Prepend[{Length[ss] - sp[[-1]]}, 
            Flatten[Append[{sp[[1]]}, 
                Table[sp[[i]] - sp[[i - 1]], {i, 2, Length[sp]}]]]]]];
    ss]
    
    fCompl[Ul_List] := Module[{ss = Ul},
    ss = ReplaceAll[
        Table[If[ss[[i, -1]] == 1, 
            Drop[ReplacePart[ss[[i]], (ss[[i, -2]] + 1), -2], -1], 
            Append[ReplacePart[ss[[i]], ss[[i, -1]] - 1, -1], 1]], {i, 
            Length[ss]}], {1, 1} -> {2}];
    {Ul, ss}]
    
    fStelRot[Ul_List] := Module[{ppr, ppr1, ppr2, ppr3},
    ppr = Union[Table[RotateLeft[Ul[[1]], i], {i, Length[Ul[[1]]]}]];
    ppr1 = Map[Reverse, ppr];
    ppr2 = Union[Table[RotateLeft[Ul[[2]], i], {i, Length[Ul[[2]]]}]];
    ppr3 = Map[Reverse, ppr2];
    ppr = Union[ppr, ppr1, ppr2, ppr3];
    ppr
    ]
    
    fMakeStel[Ul_List] := Module[{stel = Ul, i},
    stel = 
      Table[Map[
          StringReplace[ToString[#], {"{" -> "", "}" -> "", " " -> ""}] &, 
          stel[[i]]], {i, Length[stel]}];
    stel = Table[StringReplace[stel[[i]], "," -> " "], {i, Length[stel]}];
    stel = 
      Table[StringReplace[
          ToString[stel[[i]]], {"{" -> "", "}" -> "", ", " -> ","}], {i, 
          Length[stel]}];
    stel]
    
    fStellarNalt[n_Integer] := Module[{stel, d, m, st1, stelpar, stelnep, i, \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
j, k},
    stel = fStellar[n];
    d = IntegerPart[n/2];
    m = IntegerPart[d/2];
    st1 = 
      Table[StringJoin[stel[[i]], "+-1"], {i, 
          Length[stel]}];   (*svi sa jednim minusom*)
    stel = 
      Table[Select[stel, Length[StringPosition[#, ","]] == i &], {i, 3, 
          d - 1}];
    stelpar = Table[stel[[2i - 1]], {i, m - 1}];
    stelnep = Reverse[Complement[stel, stelpar]];
    stelnep = 
      Flatten[Table[
          Table[StringJoin[stelnep[[j, i]], "+-", ToString[k]], {k, 2, 
              j + 1}, {i, Length[stelnep[[j]]]}], {j, 
            Length[stelnep]}]];     (*dodati svi minusi neparnim *)
    stelpar1 = 
      Flatten[Table[
          Table[StringJoin[stelpar[[j, i]], "+-", ToString[k]], {k, 2, j}, \
{i,
               Length[stelpar[[j]]]}], {j, 
            Length[stelpar]}]];  
            (*dodati svi minusi parnim izuzev zabranjenog broja *)
    stelpar = Table[Map[fStelString, stelpar[[i]]], {i, Length[stelpar]}];
    stelpar = Table[Map[fCompl, stelpar[[i]]], {i, Length[stelpar]}];
    stelpar = 
      Map[Union, Table[Map[fStelRot, stelpar[[i]]], {i, Length[stelpar]}]];
    stelpar = Table[Map[First, stelpar[[i]]], {i, Length[stelpar]}];
    stelpar = Table[Map[Reverse, stelpar[[i]]], {i, Length[stelpar]}];
    stelpar = Table[fMakeStel[stelpar[[i]]], {i, Length[stelpar]}];
    stelpar = 
      Flatten[Table[
          Table[StringJoin[stelpar[[j, i]], "+-", ToString[k]], {k, 2, 
              j + 1}, {i, Length[stelpar[[j]]]}], {j, Length[stelpar]}]];
    stel = Join[Reverse[Sort[st1]],Reverse[Sort[stelnep]], \
Reverse[Sort[stelpar]]];
    stel
    ]



(*################### NMoveRat ############################*)

ReduMove[n_Integer,Conway_List]:=
    Module[{Con,i=1,j,l,ind=False,Con1,Res={}},
      Switch[Head[Conway[[1]]],Integer,l=1,List,l=Length[Conway]];
      While[i<=l&&Not[ind],
        Switch[Head[Conway[[1]]],Integer,Con=Conway,List,Con=Conway[[i]]];
        Con=IsNotKnot[Con];
        If[
          Con=={}||Con=={0}||Con=={1}||Con=={-1}||
            Con=={2}||Con=={-2},ind=True;Res={},j=1;
          
          While[Not[MemberQ[Res,{0}]]&&j<=Length[Con],
            Res=Union[
                Append[Res,ReplacePart[Con,Con[[j]]-n*Sign[Con[[j]]],j]]];
            (*Print[Res,MemberQ[Res,{}]];*)j++]];
        i++];
      Res];

NMoveRat[n_Integer,Conway_String]:=
  Module[{Con,br=-1},
    Con=ReadList[
        StringToStream[
          StringReplace[
            Conway,{","->" ","+"->" ","("->" ",
              ")"->" "}]],Number];
    Con=IsNotKnot[Con];
    If[Not[
        Con=={}||Con=={0}||Con=={1}||Con=={-1}||
          Con=={2}||Con=={-2}],
      While[Con!={},br++;Con=ReduMove[n,Con]],br=0];
    br]
    
(*###########################################################*)

(*############### fJablanPoly ####################*)
(*calculates multivariable Jablan polynomial*)

fJablanPoly[Conway_String]:=
  Module[{lL0,res,lL,t,i,s,sl,j,p,MAT},
    If[Conway=="1",res=0,
      If[Conway=="2",res=1,lL0=fGenerators[Conway][[1]];
        lL=Flatten[lL0];
        t=ZeroMatrix[1,Length[lL0]][[1]];
        Do[
          If[i==1,t[[i]]=Length[lL0[[i]]]/3,
            t[[i]]=t[[i-1]]+Length[lL0[[i]]]/3],{i,1,Length[lL0]}];
        s=fGenerators[Conway][[2]];
        sl=Length[s];(*duzina*)
        MAT=ZeroMatrix[sl];
        j=1;p="a";
        For[i=1,i<=sl,i++,If[i>t[[j]],j++;
            p=FromCharacterCode[ToCharacterCode["a"]+j-1]];
          
          Switch[s[[i]],1,MAT[[i,Part[lL,3(i-1)+1]]]=1,-1,
            MAT[[i,Part[lL,3(i-1)+1]]]=-1];
          
          Switch[s[[i]],1,MAT[[i,Part[lL,3(i-1)+2]]]=-1,-1,
            MAT[[i,Part[lL,3(i-1)+2]]]=1];
          MAT[[i,Part[lL,3i]]]=p;];
        res=Det[MAT]
        ]];
    res]
(*########################################################*)

(*######### Unlinking or unknotting number of fixed projection \
###############*)

fGener[{LK_List,L_List}]:=
  Module[{pom,p,i,pl=L,l,rez={},pDa={},br=1},
    Do[If[MemberQ[Abs[pl],i]==False,
        pDa=Append[pDa,{i*Sign[pl[[br]]],pl[[br]]}];br++],{i,1,2Length[pl]}];
    Do[pom=pDa;
      l={LK};
      pom[[i]]*=-1;
      pom[[i]]=Reverse[pom[[i]]];
      pom=Sort[pom,Abs[#1[[1]]]<Abs[#2[[1]]]&];
      (*p je oblika originalnog linka*)
      
      p=Flatten[Map[Take[#,-1]&,pom]];
      l=Append[l,p];
      rez=Append[rez,l],{i,1,Length[pDa]}];
    rez;
    Union[rez]]

 (*################################################################*)
 
    
    (*Ulaz je Konvej ili PData*)
(*Radi sve pojedinacne promene znaka na fix \
projekciji*)
fUnKLFixed[UL_,n_Integer]:=Module[{pd=UL,i,prvi,vv,res={}},
    If[SameQ[Head[UL],String],
    pd=fCreatePData[UL]];
    (*sad su pd pdata*)
    pd=ReductionKnotLink[pd];
    vv=fVarNewPGap[Length[pd[[2]]],n];
    pd=Union[Map[fCrChVarPD[pd,#]&,vv ]];
    prvi=Union[Map[Rest[#][[1]]&,pd]];
    pd=Sort[Map[Reverse[#]&,pd]];
    (*vv=Union[Map[Apply[Plus,#]&,vv]];Print[vv];*)
    Do[res=Append[res,
          Reverse[First[Select[pd,#[[1]]==prvi[[i]]&]]]
          ],      
      {i,1,Length[prvi]}];
    res=Reverse[Sort[res]];
    res=Flatten[Select[res,#[[2]]=="pd: {{},{}}"&]];
    res=If[SameQ[res,{}],{},res[[1]]]]
    
    
    (*###########################################################*)
    
    (* calculates UnKLNo od Konveja, pdata i dowkera *)
UnKLNo[Ulaz_]:=
  Module[{pD=Ulaz,str="",R,NoGen=1,un0},
    If[SameQ[Head[Ulaz],String],pD=fCreatePData[Ulaz]];
    If[SameQ[Head[Ulaz],List]&&SameQ[Select[Ulaz[[2]],OddQ[#]&],{}],
       pD=fPDataFromDow[Ulaz]];
      un0=MemberQ[ReductionKnotLink[pD],{}];
    If[un0==True,
          str="Unknott",
                  R=fGenerate[pD];
           
      While[MemberQ[R,{{},{}}]==False&&
          MemberQ[R,{{0},{}}]==False,  
        	R=Flatten[Map[fGenerate[#]&,R],1];
        	NoGen++];];
   NoGen](*16.8.2003 *)
   
   (*######################fGap##################*)
   
   fGap[Conway_String]:=Module[{k,v,l},
    k=UnKLNo[Conway];
    v=k;
    l=fUnKLFixed[Conway,k];
    While[SameQ[fUnKLFixed[Conway,k],{}],k++];
    {k-v,{v,k},Conway}] (*27.08.2004*)
    
  (* fGapKnotsc[Ul_List]:=Module[{k,dd,l,v},
    dd=fPDataFromDow[Ul];
    k=UnKLNo[dd];
    v=k;
    l=fUnKLFixed[dd,k];
    While[SameQ[fUnKLFixed[dd,k],{}],k++];
    {k-v,{v,k},Ul}] *)
   
       
    (*################fUnRFixProj###########################*)
    
    (* tests that a fixed projection of rational 
    KL given by Conway can be unlinked by k crossing changes. Makes
    crossing changes directly in Conway symbol *)
    
    fUnRFixProj[Conway_String,k_Integer]:=Module[{Con,n,r,p,i,res},
    Con=ReadList[
        StringToStream[
          StringReplace[
            Conway,{","->" ","+"->" ","("->" ",
              ")"->" "}]],Number];
    n=Apply[Plus,Con];
    r=Length[Con];
    If[SameQ[r,1]&&EvenQ[n],res=If[SameQ[k,n/2],1,0],
    p=par[k,r];
    p=Select[p,Max[#]<=Max[Con]&];
    p=Flatten[Union[Table[DistinctPermutations[p[[i]]],{i,Length[p]}]],1];
    p=Table[Con-p[[i]],{i,Length[p]}];
    p=Select[p,Min[#]>=0&];
    p=Table[Con-p[[i]],{i,Length[p]}];
    p=Table[Con-2p[[i]],{i,Length[p]}];
    p=Map[IsNotKnot,p];
    res=If[MemberQ[p,{1}]||MemberQ[p,{-1}]||MemberQ[p,{0}],1,0]];
    res]
    (*24.6.2004*)
    
    (*######################fGapRat#############################*)
    
    (* calculates unlinking gap, unlinking no. and
    fixed projection unlinking no. for a rational KL *)
    
    fGapRat[Conway_String]:=Module[{k,l,m,v},
    k=Max[UnR[Conway]];
    v=k;
    l=fUnRFixProj[Conway,k];
    While[SameQ[fUnRFixProj[Conway,k],0],k++];
    {k-v,{v,k},Conway}]
    
    (*###################### AUTOMATS #############################*)
    (*################### fAutoSigInp ############################*)
    
    (*Gr_List={outgiong edges of graph, signlist,list of inputs}, NAND=1, \
AND=-1*)
    
fAutoSigInp[Gr_List]:=Module[{gr,s,e,n,ul,iz,i,j,k,v,res,st,inp,in},
    gr=Gr[[1]];
    s=Gr[[2]];
    inp=Gr[[3]];
    gr=Flatten[Table[Table[{i,gr[[i,j]]},{j,Length[gr[[i]]]}],
    {i,Length[gr]}],1];
    gr=Union[gr];
    e=Length[gr];
    v=fVarP[e];
    n=Length[Union[Flatten[gr]]];
    s=If[SameQ[s,{}],Table[1,{i,n}],s];
    in=Table[1,{i,n}];
    inp=If[SameQ[inp,{}],inp,Table[{inp[[i]]},{i,Length[inp]}]];
    inp=If[SameQ[inp,{}],in,ReplacePart[in,0,inp]];
    e=Table[{gr[[i]],i},{i,Length[gr]}];
    ul=Table[Select[e,#[[1,2]]==i&],{i,n}];
    ul=Table[Table[ul[[j,i,2]],{i,Length[ul[[j]]]}],{j,Length[ul]}];
    iz=Table[Select[e,#[[1,1]]==i&],{i,n}];
    iz=Table[Table[iz[[j,i,2]],{i,Length[iz[[j]]]}],{j,Length[iz]}];
    st=iz;
    ul=Flatten[
        Table[Table[
            Table[v[[i,ul[[j,k]]]],{k,Length[ul[[j]]]}],{j,Length[ul]}],{i,
            Length[v]}],1];
    iz=Flatten[
        Table[Table[
            Table[v[[i,iz[[j,k]]]],{k,Length[iz[[j]]]}],{j,Length[iz]}],{i,
            Length[v]}],1];
    e=Map[Max,Partition[Map[Length,Map[Union,iz]],n]];
    ul=Partition[Map[Min,ul],n];
    ul=Table[inp*ul[[i]],{i,Length[ul]}];
    iz=Partition[Map[Min,iz],n];
    iz=Table[(s+1)/2-s*iz[[i]],{i,Length[iz]}];
    res=Union[
        Table[If[SameQ[ul[[i]],iz[[i]]]&&SameQ[e[[i]],1],v[[i]],{}],{i,
            Length[v]}]];
    res=If[SameQ[res[[1]],{}],Drop[res,1],res];
    e=Map[Length,st];
    e=Flatten[
        Map[Union,
          Flatten[Union[Table[iteratedTake[res[[i]],e],{i,Length[res]}]],
            1]]];
    e=Partition[e,Length[st]];
    Print[res];
    e]
    
   
    (*################### fAutoKL ############################*)
    
    (* fAutoKL[Ul_List]: Ul={Conway,signs,inputs}; 
    If[signs={},signs=signs(Con) *)
    
    
fAutoKL[Ul_List]:=Module[{Con,gg,gg1,i,j,s,ss,inp},
    Con=Ul[[1]];
    inp=Ul[[3]];
    ss=Ul[[2]];
    gg=fGaussExtSigns[Con];
    Print["Gauss code:", gg];
    s=If[SameQ[ss,{}],Map[Sign,Reverse[Union[Flatten[gg]]]],ss];
    Print["Signs:", s];
    gg1=Union[Flatten[gg]];
    gg1=Sort[Table[{Abs[gg1[[i]]],Sign[gg1[[i]]]},{i,Length[gg1]}]];
    gg1=Table[gg1[[i,2]],{i,Length[gg1]}];
    gg=If[SameQ[gg,Flatten[gg]],
        Sort[Join[Abs[{{gg[[-1]],gg[[1]]}}],
            Table[Abs[{gg[[i]],gg[[i+1]]}],{i,Length[gg]-1}]]],
        Sort[Flatten[
            Table[Join[Abs[{{gg[[j,-1]],gg[[j,1]]}}],
                Table[Abs[{gg[[j,i]],gg[[j,i+1]]}],
                {i,Length[gg[[j]]]-1}]],{j,Length[gg]}],1]]];
                Print["Graph:", gg];
    gg=Table[Select[gg,#[[1]]==i&],{i,Length[Union[Flatten[gg]]]}];
    gg=Table[Table[gg[[i,j,2]],{j,Length[gg[[i]]]}],{i,Length[gg]}];
    gg=fAutoSigInp[{gg,s,inp}];
    gg]
    
(*################# LiangPoly #######################*)
    

LiangPoly[Conway_String]:=Module[{MM,gg,ss,n,i,j},
    gg=fGraphInc[Conway];
    ss=gg[[2]];gg=gg[[1]];
    n=Length[ss];
    MM=DiagonalMatrix[Map[Power[t,#]&,ss]];(*Print[n,MM];*)
    Do[Do[
        p=Count[gg,{i,j}];
        MM=ReplacePart[MM,p,{i,j}];
        MM=ReplacePart[MM,p,{j,i}];
        ,{j,i+1,n}],{i,1,n-1}];
    (*Print[MM];*)res=Det[MM];
    Print["Amphicheiral ",
      SameQ[res,Det[ReplaceAll[MM,{t\[Rule]1/t,1/t\[Rule]t}]]]];
    res
    ]
(*################## fBoolean ########################*)

(*pomocna za fBrisi stepen: u pol sve stepene xx mkonvertuje u xx *)

fImaX[Pol_,xx_]:=Module[{pom=Pol,ind=True},
    If[MemberQ[Variables[pom],xx],
      While[ind,
        ind=MemberQ[Variables[pom],xx]&&Not[MemberQ[Variables[pom/xx],xx]];
        If[ind,ind=False,pom=pom/xx;ind=True]]];
    pom]
    
(*stepene svih varijabli konvertuje u same varijable*)

fBrisiStepen[Pol_]:=Module[{pom=Pol,i,var},
    var=Variables[Pol];
    Do[pom=Map[ fImaX[#,var[[i]]]&,pom]
      ,{i,1,Length[var]}];pom
    ]
    
(*racuna polinome za logicke funkcije*)

fBoolean[LL_]:=Module[{p},
    Unprotect[And];
    Unprotect[Or];
    Unprotect[Not];
    Unprotect[Nand];
    Unprotect[Nor];
    Unprotect[Implies];
    Unprotect[Xor];
    Unprotect[Equal];
    a_&&b_:=a b;
    a_||b_:=a+b-a b;
    Not[a_]:=1-a;
    Nor[a_,b_]:=1-a-b+a b;
    Nand[a_,b_]:=1-a b;
    Implies[a_,b_]:=1-a+a b;
    Xor[a_,b_]:=a+b-2a b;
    Equal[a_,b_]:=1-a-b+2a b;
    p=Expand[LL];
    Clear[And];
    Clear[Or];
    Clear[Not];
    Clear[Nand];
    Clear[Nor];
    Clear[Implies];
    Clear[Xor];
    Clear[Equal];
    p=fBrisiStepen[p];
    p
    ]
    
    

(*################## fDiffSeq ########################*)

  
(* fPerBasic produces material for periodic sequences 
without dual material
or constant parts *)
fPerBasic[n_Integer,p_Integer]:=Module[{v,u,vr,vp,vd, i, j, k},
    v=Sort[fVarP[n]];
    v=KSubsets[v,p];
    u=Table[Table[Table[v[[k,i,j]],{i,p}],{j,n}],{k,Length[v]}];
    v=Union[
        Table[If[Length[Union[Map[Length,Map[Union,u[[i]]]]]]\[Equal]1,
            v[[i]],{}],{i,Length[v]}]];
    v=Map[Sort,If[v[[1]]\[Equal]{},Drop[v,1],v]];
    vd=Flatten[
        Table[{i,
            Position[v,Sort[ReplaceAll[v[[i]],{0\[Rule]1,1\[Rule]0}]]]},{i,
            Length[v]}]];
    vd=Flatten[Union[Map[Sort,Partition[vd,2]]]];
    vd=Table[vd[[2i-1]],{i,Length[vd]/2}];
    v=Map[Sort,Table[v[[vd[[i]]]],{i,Length[vd]}]];
    vp=Union[
        Map[Sort,Table[Table[Table[v[[k,i,j]],{i,p}],{j,n}],{k,Length[v]}]]];
    v=Table[
        Table[Table[vp[[k,i,j]],{i,Length[vp[[k]]]}],{j,p}],{k,Length[vp]}]
    ]
(* fPerm produces permutations different with regard to automorphisms of \
dihedral group - rotations and reverse *)
fPerm[n_Integer]:=Module[{p,q,l, i, j},
    p=Permutations[Range[n]];
    p=Map[First,
        Union[Table[Sort[Table[RotateLeft[p[[i]],j],{j,n}]],{i,Length[p]}]]]
    ]
fPerSeq[n_Integer,p_Integer]:=Module[{v,s, i, j, k},
    If[2^n\[GreaterEqual]p>2^(n-2),
      v=fPerBasic[n,p];
      s=fPerm[p];
      s=Flatten[
          Table[Table[Table[v[[k,s[[j,i]]]],{i,p}],{j,Length[s]}],{k,
              Length[v]}],1];
      s=Union[
          Table[Sort[Table[Table[s[[k,i,j]],{i,p}],{j,n}]],{k,Length[s]}]];
      s=Union[Table[Table[Table[s[[k,i,j]],{i,n}],{j,p}],{k,Length[s]}]],
      Print["2^(n-2)<p<=2^n"]]
    ]
fKauffAlg[n_Integer,p_Integer]:=Module[{l,t, i, j, k},
    l=fPerSeq[n,p];
    t=Table[
        Join[{Table[
              If[l[[k,1,j]]\[Equal]0,l[[k,-1]],{}],{j,Length[l[[k,1]]]}]},
          
          Table[Table[
              If[l[[k,i,j]]\[Equal]0,l[[k,i-1]],{}],{j,Length[l[[k,i]]]}],{i,
              2,Length[l[[k]]]}]],{k,Length[l]}];
    t=Table[Map[Union,Table[Table[t[[k,i,j]],{i,p}],{j,n}]],{k,Length[l]}];
    t=Table[
        Table[If[t[[j,i,1]]\[Equal]{},Drop[t[[j,i]],1],t[[i,1]]],{i,
            Length[t[[j]]]}],{j,Length[t]}]
    ]
(*pomocna za fBrisi stepen: u pol sve stepene xx mkonvertuje u xx *)

fImaX[Pol_,xx_]:=Module[{pom=Pol,ind=True},
    If[MemberQ[Variables[pom],xx],
      While[ind,
        ind=MemberQ[Variables[pom],xx]&&Not[MemberQ[Variables[pom/xx],xx]];
        If[ind,ind=False,pom=pom/xx;ind=True]]];
    pom]
(*stepene svih varijabli konvertuje u same varijable*)

fBrisiStepen[Pol_]:=Module[{pom=Pol,i,var},
    var=Variables[Pol];
    Do[pom=Map[ fImaX[#,var[[i]]]&,pom]
      ,{i,1,Length[var]}];pom
    ]

fOp[L_List]:=Module[{l,k=1,i,j,n,pom,n1=Length[L]},
    l=L;
    n=Length[l[[1,1]]];
    (*n je broj varijabli koje koristimo*)
    (*Print[l,n];*)
    Do[
      pom=Symbol[FromCharacterCode[96+k]];
      Do[Do[
            
          If[l[[i,j,k]]\[Equal]1,l=ReplacePart[l,pom,{i,j,k}],
            l=ReplacePart[l,1-pom,{i,j,k}]]
          ,{j,1,Length[L[[i]]]}],{i,1,n1}]
      ,{k,1,n}];
    (*Print["LLL ",l];*)
    
    l=Table[Expand[
          Times@@Table[1-Expand[Times@@l[[j,i]]],{i,Length[l[[j]]]}]],{j,
          Length[l]}];
    (*Print["Polinomi",l];*)
    l=Map[fBrisiStepen,l]]
(*pomocna za fDiffSeq*)
(*pomocna za fPraviZamene i fBalanced*)

fPZ[PP_List,CP_List,OP_Integer]:=Module[{res={},i},
    Do[  If[SameQ[OP,0],res=Append[res,{PP[[i]]\[Rule]CP[[i]]}],
                                                    
        res=Append[res,{PP[[i]]==CP[[i]]}]],{i,1,Length[CP]}];
    Flatten[res,1]]
(*pomocna za fDiffSeq*)
(*za n simbole pravi listu svih mogucih zamena na \
osnovu Distinct permutations*)
fPraviZamene[n_Integer]:=Module[{l,p},
    l=DistinctPermutations[
        Map[Symbol[FromCharacterCode[#]]&,Range[97,96+n]]];
    l=Map[fPZ[First[l],#,0]&,l]
    ]
(*pomocna za fDiffSeq*)
(*svakom elementu liste na odgovarajuce mesto \
dopisuje uredjen par sa odgovarajucom promenljivom- simbolom*)

fDodajSlovo[n_Integer, LL_List]:=Module[{l,res=LL},
    l=Map[Symbol[FromCharacterCode[#]]&,Range[97,96+n]];
    Do[  res=Map[ReplacePart[#,{l[[i]],#[[i]]},i]&,res],{i,1,n}];res]

fDiffSeq[n_Integer,p_Integer]:=Module[{l,lll={},ll,pp,qq,zam, i},
    l=fPerSeq[n,p];
    pp=fKauffAlg[n,p];
    qq=pp;
    pp=Map[fOp,pp];
    pp=fDodajSlovo[n,pp];
    (*svakom elementu dodaje odgovarajucu promenljivu*)
    
    zam=fPraviZamene[n];
    (*zamene na osnovu svih permutacija od n-simbola*)
    
    Do[pom=Sort[Map[Sort[ReplaceAll[pp[[i]],#]]&,zam]];
      (*Print["PPPPP",pom];*)
      lll=Append[lll,pom],{i,1,Length[pp]}];
    (*Print[Length[pp],"LLL ",Length[lll]];*)
    ll=Union[lll];
    (*Print["KKK",Length[ll],Length[lll]];*)
    
    pp=Sort[Table[Position[lll,ll[[i]]],{i,Length[ll]}]];;
    pp=Table[pp[[i,1,1]],{i,Length[pp]}];(*Print[pp,Length[lll]];*)
    
    ll=Table[lll[[pp[[i]]]],{i,Length[pp]}];
    (*Print[Length[l],Max[pp]];*)
    l=Table[l[[pp[[i]]]],{i,Length[pp]}];
    qq=Table[qq[[pp[[i]]]],{i,Length[pp]}];
    ll=Table[{ll[[i,1]][[1,2]],ll[[i,1]][[2,2]],ll[[i,1]][[3,2]]},{i,
          Length[ll]}];
    ll=Table[{l[[i]],qq[[i]],ll[[i]]},{i,Length[l]}];
    Print["There are ", Length[ll]," different periodic sequences for ",n, 
      " & ",p];
    ll
    ]
    
fBalanced[Diff_List,redBr_Integer]:=Module[{n,poly,l1,l2},
    n=Length[Diff[[1,1,1]]];
    If[redBr>Length[Diff],
         Print["There are only ", Length[Diff],
        " different periodic sequences. Please pick lesser number"]
      ,
      poly=Diff[[redBr,3]];
      l2=Map[Symbol[FromCharacterCode[#]]&,Range[97,96+n]];
      l1=fPZ[poly,l2,1]
      ];
    Print["Periodic sequence: ", Diff[[redBr,1]]];
    Print["Kauffman's code: ", Diff[[redBr,2]]];
    Print["Polynomials: ",poly];
    Print["Balanced states: ",Solve[l1,l2]];    
    ]

(*################## fKnotFind ########################*)


fForknotFind[Conway_String]:=Module[{l,DL,str,f="ff.txt"},
    DL=Dowker[Conway][[4]];
    OpenWrite[f];
     l=fGenSign[Conway]*fGenSign[StringReplace[Conway,{"-"\[Rule]""}]];
     DL*=l;
    str=ToString[Length[DL]]<>" 1 "<>" "<>
        StringReplace[ToString[DL],{"{"->"","}"->"",","->""}];
     WriteString[f,str];
       Write[f];
    Close[f]
    ]
    
fKnotFind[Conway_String]:=Module[{res,k,s},
    If[fComponentNo[Conway]\[Equal]1,
      res=Dowker[Conway][[4]];
      k=fForknotFind[Conway];
      s=fGenSign[Conway];
      s=s*fGenSign[StringReplace[Conway,{"-"\[Rule]""}]];
      Run["fromdos.exe","ff.txt"];
      Run["knotfind.exe","ff.txt","rez.txt"];
      res=Import["rez.txt"];
      res= If[Not[SameQ[StringTake[res,1], "u"]]==True,
      DeleteFile["ff.txt"];
      DeleteFile["rez.txt"] ;
      res=StringReplace[res,"  "->" "];
      res=StringReplace[res,"  "->" "];
      res=StringDrop[StringReplace[res,"  "->" "],3];
      res=Drop[ToCharacterCode[res],-1];
      res=If[First[res]\[Equal]32,Rest[res],res];
      res=ReadList[StringToStream[FromCharacterCode[res]],Number];
      res=Prepend[{res},{Length[res]}],1];
      res,
      Print["Link"]]]
      

RedKauffmanPolynomial[Ul_List] := Module[{kk, ff},
    kk = KauffmanPolynomial[Ul];
    If[SameQ[Position[kk[[1]], a], {}], kk,
      kk = ReplaceAll[kk, xa -> x*a]];
    ff = FactorList[kk][[2]];
    If[SameQ[Position[ff, x], {}],
      ff = Expand[Divide[kk, Power[ff[[1]], ff[[2]]]]], ff = kk];
    If[ff == 1, ff = kk];
     ff]
     
     (* ############# NONALGEBRAIC TANGLES ############ *)
     
PDataSumProdNonA[L_List, tangle_String] := 
  Module[{BrPreseka, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, 
      x21, x22, res, str, f = "NonABase.txt"},
    BrPreseka = Divide[Length[L[[2]]], 2];
    Switch[BrPreseka, 10, x10 = L, 11, x11 = L, 12, 
    x12 = L, 13, x13 = L, 14, 
      x14 = L, 15, x15 = L, 16, x16 = L, 17, x17 = L, 18, x18 = L, 19, 
      x19 = L, 20, x20 = L, 21, x21 = L, 22, x22 = L];
     str = StringJoin["x", ToString[BrPreseka], "=", ToString[L], ";"];
     OpenWrite[f];
    WriteString[f, str];
    Write[f];
    Close[f];
    Get["NonABase.txt"];
    res = fCreatePData[StringJoin[ToString[BrPreseka], "*", tangle]]
    ]
    
fProdTangles[s1_String, s2_String, tt_String] := Module[{pol, pol1, L, LL},
    L = ToExpression[StringDrop[s1, -1]];
    LL = ToExpression[StringDrop[s2, -1]];
    pol = {{}, Union[L[[2]], LL[[2]] + 2Length[L[[2]]]], 
        Join[L[[3]], -LL[[3]]], L[[4]], LL[[4]] + 2Length[L[[2]]]};
    pol1 = 
      Map[Sort, {{pol[[4, 1, 1]], pol[[5, 1, 2]]}, {pol[[4, 2, 1]], 
            pol[[5, 1, 1]]}, {pol[[4, 1, 2]], 
            pol[[5, 2, 2]]}, {pol[[4, 2, 2]], pol[[5, 2, 1]]}}];
    pol = {{}, 
        Union[Complement[pol[[2]], Map[Sort, Union[pol[[4]], pol[[5]]]]], 
          pol1], pol[[3]]};
    PDataSumProdNonA[pol, tt]]
    
fSumTangles[s1_String, s2_String, tt_String] := Module[{pol, pol1, L, LL},
    L = ToExpression[StringDrop[s1, -1]];
    LL = ToExpression[StringDrop[s2, -1]];
    pol = {{}, Union[L[[2]], LL[[2]] + 2Length[L[[2]]]], 
        Join[-L[[3]], -LL[[3]]], L[[4]], LL[[4]] + 2Length[L[[2]]]};
    pol1 = 
      Map[Sort, {{pol[[4, 1, 1]], pol[[5, 2, 1]]}, {pol[[4, 2, 1]], 
            pol[[5, 1, 1]]}, {pol[[4, 1, 2]], 
            pol[[5, 2, 2]]}, {pol[[4, 2, 2]], pol[[5, 1, 2]]}}];
    pol = {{}, 
        Union[Complement[pol[[2]], Map[Sort, Union[pol[[4]], pol[[5]]]]], 
          pol1], pol[[3]]};
    PDataSumProdNonA[pol, tt]]
    
     (* ############# VIAE ############ *)
     
     (*HAMCUBGR.M-Generate Hamilton cubic graphs*)
     (*13-DEC-1995 TAMARA BERTOK*)
     
     AdmissibleEdge[x_,y_,n_,c_]:=
  If[Abs[x-y]=!=1&&Abs[x-y]=!=n-1&&y-x+1\[GreaterEqual]c,True,False]
  
ListOfOneFactors[n_] := 
  Module[{koncni, seznam, s, x, k, l, tmp, tmp1, tmp2, tmp3, c}, 
    s = Range[n];
    koncni = {};
    seznam = {};
    k = 2;
    While[k <= Length[s], If[s =!= {}, x = s[[1]]];
      If[x == 1 && k > n/2 + 1, k = Length[s] + 1, 
        Do[If[koncni =!= {}, c = koncni[[1, 2]], c = 0];
          
          If[AdmissibleEdge[x, s[[k]], n, c], 
            Do[koncni = Append[koncni, {x, s[[k]]}];
              s = Complement[s, {x, s[[k]]}];
              k = 2], 
            Do[If[k < Length[s], k = k + 1, 
                If[koncni =!= {}, Do[tmp = Last[koncni];
                    s = Union[s, tmp];
                    koncni = Complement[koncni, {tmp}];
                    l = Position[s, tmp[[2]]][[1, 1]];
                    k = l + 1], k = Length[s] + 1]]]];
          If[s == {}, seznam = Append[seznam, koncni];
            tmp1 = Last[koncni];
            tmp2 = koncni[[Length[koncni] - 1]];
            koncni = Complement[koncni, {tmp1, tmp2}];
            s = Union[s, tmp1, tmp2];
            l = Position[s, tmp2[[2]]][[1, 1]];
            k = l + 1;];
          If[k > Length[s] && Length[s] < n, tmp3 = Last[koncni];
            koncni = Complement[koncni, {tmp3}];
            s = Union[s, tmp3];
            l = Position[s, tmp3[[2]]][[1, 1]];
            k = l + 1;
            While[tmp3[[1]] =!= 1 && l == Length[s], Do[tmp3 = Last[koncni];
                s = Union[s, tmp3];
                koncni = Complement[koncni, {tmp3}];
                l = Position[s, tmp3[[2]]][[1, 1]];
                k = l + 1]]]]] (*If*)];(*While*)
          Return[seznam]
    ]
    
    
 fDiffViae[n_Integer] := 
  Module[{ g = {}, i, j, g1 = {}, m = 0, gnova = {}, pot, d, d1, dd},
    gnova = ListOfOneFactors[2n];
    g = Table[
        Union[gnova[[i]], Map[Reverse, gnova[[i]]]], {i, Length[gnova]}];
    g = Flatten[
        Table[Table[g[[j, i, 2]] - g[[j, i, 1]], {i, 2n}], {j, Length[g]}]];
    g = Partition[
        Table[If[Abs[g[[i]]] > n, -Sign[g[[i]]] 2n + g[[i]], g[[i]]], {i, 
            Length[g]}], 2n];
    g = ReplaceAll[g, -n -> n];
    g = Table[
        First[Union[
            ReplaceAll[Table[RotateLeft[g[[i]], j], {j, 2n}], -n -> n], 
            ReplaceAll[
              Map[Reverse, Table[RotateLeft[g[[i]], j], {j, 2n}]], -n -> n], 
            ReplaceAll[-Table[RotateLeft[g[[i]], j], {j, 2n}], -n -> n], 
            ReplaceAll[-Map[Reverse, 
                  Table[RotateLeft[g[[i]], j], {j, 2n}]], -n -> n]]], {i, 
          Length[gnova]}];
    g1 = Union[g];
    g = Table[{g[[i]], gnova[[i]]}, {i, Length[gnova]}];
    g = Table[Select[g, #[[1]] == g1[[i]] &], {i, Length[g1]}];
    g = Sort[Table[First[g[[i]]][[2]], {i, Length[g]}]];
    d = Union[Table[{i, i + 1}, {i, 2n - 1}], {{1, 2n}}];
    g = Union[
        Table[If[PlanarQ[FromUnorderedPairs[Union[g[[i]], d]]] == True, 
            g[[i]], 0], {i, Length[g]}]];
    g = If[SameQ[g[[1]], 0], Drop[g, 1], g];
    g = Union[
        Table[If[VertexConnectivity[FromUnorderedPairs[Union[g[[i]], d]]] > \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
2,
             g[[i]], {}], {i, Length[g]}]];
    g = If[SameQ[g[[1]], {}], Drop[g, 1], g];
    g
    ] 
    
    
    fViaToKL[LL_List]:=Module[{d,gg},
    d=Length[LL];
    gg=Join[Table[{i,i+1},{i,2d-1}],{{1,2d}}];
    gg=Sort[
        Map[Sort,
          Table[{Position[LL,gg[[i,1]]][[1,1]],
              Position[LL,gg[[i,2]]][[1,1]]},{i,Length[gg]}]]];
    gg=fKLfromGraph[gg];
    gg
    ]
    
     (* ############# COMBINATORICA 42 ############ *)
     
     Unprotect[ConnectedComponents]
      
ConnectedComponents[g_Graph] :=
        Block[{untraversed=Range[V[g]], visit, comps={}, \
e=ToAdjacencyLists[g],
              dfi=Table[0,{V[g]}],cnt=1, $RecursionLimit = Infinity},
              While[untraversed != {},
                    visit = {}; edges = {};
                    DFS[First[untraversed]];
                    AppendTo[comps,visit];
                    untraversed = Complement[untraversed,visit]        
              ];
              comps
        ] /; UndirectedQ[g] 
        
      
        
   
       Unprotect[DFS]
        
DFS[v_Integer] :=
	( dfi[[v]] = cnt++;
	  AppendTo[visit,v];
	  Scan[ (If[dfi[[#]]==0,AppendTo[edges,{v,#}];DFS[#] ])&, e[[v]] ] )
     
           
End[]      
EndPackage[]






















































































































































































































































































































































































