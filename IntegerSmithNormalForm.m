
(* :Title: IntegerSmithNormalForm, Version 1.0, 1994*)

(* :Author: David Jabon
	    Eastern Washington University
	    jabon@ewu.edu

  

*)

(* :Summary:
   A package that defines functions that return the
   Smith normal form of an integral matrix, the
   invariant factors of such a matrix and the 
   elementary divisors of  such a matrix.
   An extended version is included also which 
   returns unimodular matrices which transform
   the given matrix into its Smith normal form.
*)

(* :Context: IntegerSmithNormalForm` *)

(* :Package Version: 1.0 *)

(* :Copyright: Copyright 1994 *)

(* :History:
  Created by David Jabon at Eastern Washington
  University, 1993-4
*)

(* :Keywords: 

  Smith normal form, invariant factors, elementary divisors
*)

(* :Mathematica Version: 2.2 *)

(* :Limitations: 

This package can only be used with integral matrices.

*)

BeginPackage["IntegerSmithNormalForm`"]

SmithForm::usage =
     "SmithForm[A] gives,for an integral matrix A, 
     the Smith normal form of A. "
     
ExtendedSmithForm::usage =    
     "ExtendedSmithForm[A] gives, for an integral
      matrix A, {D, {P,Q}} where 
      D is in Smith normal form and P and Q are 
      matrices such that P A Q = D. "      
 
InvariantFactors::usage=     
     "InvariantFactors[A] gives, for an integral 
      matrix A, the invariant factors of A."
      
ElementaryDivisors::usage=
      "ElementaryDivisors[A] gives, for an integral 
      matrix A, the elementary divisors of A."    
     
Begin["`Private`"]    
             
NormalizedExtendedGCD[a_,b_]:=
          If[ Mod[b,a]==0, Return[{Abs[a],{Sign[a],0}}],Return[
          ExtendedGCD[a,b]] ]     

R1[i_,u_,A_]:= If[Depth[A]==4, Return[Map[ R1[i,u,#] &, A]],Return[
      Array[Which[ #==i, u A[[i]] , True, A[[#]] ] &, Length[A] ] ]]    ;
C1[i_,u_,A_]:=If[Depth[A]==4, Return[Map[C1[i,u,#] &, A]], 
            Return[Transpose[R1[i,u, Transpose[A]]]]];
R2[c_,i_,j_,A_]:=If[Depth[A]==4, 
               Return[Map[R2[c,i,j,#] &, A]],
              Return[Array[Which[ # == j, A[[j]]+c A[[i]], 
             True, A[[#  ]]] &, Length[A]]] ];
C2[c_, i_,j_,A_]:=If[Depth[A]==4, 
                Return[Map[C2[c,i,j,#] &, A]],
                Return[ Transpose[R2[c,i,j,Transpose[A] ]]]] ;
R3[i_,j_, A_]:= If[Depth[A]==4, Return[Map[R3[i,j, #]&, A]],
              Return[
               Array[Which[ # == i, A[[j]], # ==j, 
                A[[i]], True, A[[#  ]]] &, Length[A] ]] ];              
C3[i_,j_, A_]:= If[Depth[A]==4, Return[Map[C3[i,j, #]&, A]],
                Transpose[R3[i,j,Transpose[A]]]];
R4[i_,j_,d_,a_,b_,x_,y_,A_] := If[Depth[A] ==4, 
                     Return[Map[R4[i,j,d,a,b,x,y,#] &,A]],
                     Return[Array[Which[ # == i, 
                     x A[[i]]+ y A[[j]],
                     # == j, -(b/d) A[[i]] +
                     (a/d) A[[j]], True, A[[#]] ] &,
                     Length[A]] ] ];
                            
C4[i_,j_,d_,a_,b_,x_,y_,A_] := If[Depth[A] ==4, 
                     Return[Map[C4[i,j,d,a,b,x,y, #] &,A]],
                     Transpose[R4[i,j,d,a,b,x,y, Transpose[A]]]];
             
 EliminateCol[A_]:=
            Module[{nA=A,ExtGCD},
            
            Do[If[nA[[m,j]] != 0,
            (ExtGCD=NormalizedExtendedGCD[nA[[m,m]], nA[[m,j]] ];
            nA=C4[m,j,ExtGCD[[1]],nA[[m,m]], nA[[m,j]],
            ExtGCD[[2,1]], ExtGCD[[2,2]], nA]
            ), (Continue[ ])],
             {j,m+1, cols}];
            Return[nA]];   

            
 EliminateRow[A_]:=
            Module[{nA=A, ExtGCD},
            
            Do[If[nA[[j,m]] != 0,
            (ExtGCD=NormalizedExtendedGCD[nA[[m,m]], nA[[j,m]] ];
            nA=R4[m,j,ExtGCD[[1]],nA[[m,m]], nA[[j,m]],
            ExtGCD[[2,1]], ExtGCD[[2,2]], nA]
            ), (Continue[ ])],
             {j,m+1, rows}];
            Return[nA]];            
                       
            
            
ThirdStep[A_]:=Module[{nA=A,InterchangeQ=False},
                (Do[ ( Do[If[ Mod[ A[[i, j]] , A[[m,m]] ] == 0, 
                (Continue[ ]),
                     (nA=R2[1,i, m,nA];
                     
                     nA=EliminateCol[nA];
                     
                     InterchangeQ=True;
                      Break[ ]) ],
                     {j, m+1, Dimensions[A][[2]]}]),
                 
                  {i,m+1,Dimensions[A][[1]]}];
                  
                  If[InterchangeQ , 
                  nA= EliminateRow[nA]
                  ];
                  Return[nA])];             

SmithForm[A_] := Block[ {nA=A,
                        rows, cols,MinEntry=0, Loc,
                        ZeroMatrixQ,NonZeroRowQ=True,NonZeroColQ=
                        True},
          (
                
          (* Some Initialization*)
          {rows, cols}= Dimensions[A];
           
          
                
           (* Now we start the big loop on m*)
           Do[(
           
          (*We first find the min  non-zero value in the remaining
          principal submatrix*)     
          Loc={m,m};ZeroMatrixQ=True;
          Do[ Which[nA[[i,j]] != 0 && ZeroMatrixQ, 
                       (MinEntry=Abs[ nA[[i,j]] ];Loc={i,j};
                        ZeroMatrixQ=False), 
                       nA[[i,j]] != 0 && Abs[ nA[[i,j]] ]< MinEntry,
                       (MinEntry=Abs[ nA[[i,j]] ];
                        Loc={i,j} ), 
                       True, Continue[ ] ], 
                      {i,m,rows},{j, m,cols} 
                      ]; 
           If[ZeroMatrixQ,  Break[ ]];  (*If the remaining
           principal submatrix has all zeros we are done*)        
           
           
           (*First Step: Put min value in (m,m) position *)
           
           If[Loc[[1]] != m, nA=R3[m,Loc[[1]],nA] ];
           If[Loc[[2]] != m, nA=C3[m, Loc[[2]],nA]];
           
           If[nA[[m,m]] < 0, nA=R1[m,-1,nA] ];
           
           
           (*Second Step:  Eliminate row and column *)
           
           If[ m< Min[rows, cols], (
           
           
           
           While[(Do[(
           If[ nA[[i,m]] !=0, (NonZeroRowQ=True;Break[ ]), (Continue[ ])]),
           {i,m+1, rows}];
           If[!NonZeroRowQ,
           Do[(If[ nA[[m,j]] !=0, (NonZeroColQ=True;
           Break[ ]), 
           (Continue[ ])]),
           {j,m+1,cols}]];
           NonZeroRowQ  || NonZeroColQ),
           
           
           If[NonZeroRowQ,nA = EliminateRow[nA]];
           
           
           nA = EliminateCol[ nA];
           
           NonZeroRowQ=False;NonZeroColQ=False];
           
           (*Third Step: Make sure m,m entry divides the rest*)
           nA = ThirdStep[nA];
           
       ),
            (
            Which[rows < cols, 
            
           nA = EliminateCol[nA], 
           rows> cols, 
           nA = EliminateRow[nA]
           , rows==cols, (
           Break [ ] )])]
                
                
                ),{m,1, Min[rows, cols]}];
                
                Return[nA]
                ) 
                ];               

 ExtendedEliminateCol[P_, A_, Q_]:=
            Module[{nA=A,nP=P,tA=A, nQ=Q,ExtGCD},
            
            Do[If[nA[[m,j]] != 0,
            (ExtGCD=NormalizedExtendedGCD[nA[[m,m]], nA[[m,j]] ];
            {nA, nQ}=C4[m,j,ExtGCD[[1]],nA[[m,m]], nA[[m,j]],
            ExtGCD[[2,1]], ExtGCD[[2,2]], {nA, nQ}]
            ), (Continue[ ])],
             {j,m+1, cols}];
            Return[{nP, nA, nQ}]];   

            
 ExtendedEliminateRow[P_, A_, Q_]:=
            Module[{nA=A,nP=P,nQ=Q, ExtGCD},
            
            Do[If[nA[[j,m]] != 0,
            (ExtGCD=NormalizedExtendedGCD[nA[[m,m]], nA[[j,m]] ];
            {nA, nP}=R4[m,j,ExtGCD[[1]],nA[[m,m]], nA[[j,m]],
            ExtGCD[[2,1]], ExtGCD[[2,2]], {nA, nP}]
            ), (Continue[ ])],
             {j,m+1, rows}];
            Return[{nP, nA, nQ}]];            
                       
            
            
ExtendedThirdStep[P_, A_, Q_]:=
                Module[{nA=A,nP=P, nQ=Q,InterchangeQ=False},
                (Do[ ( Do[If[ Mod[ A[[i, j]] , A[[m,m]] ] == 0, 
                (Continue[ ]),
                     ({nA,nP}=R2[1,i, m,{nA,nP}];
                     
                     {nP, nA, nQ}=ExtendedEliminateCol[nP,nA, nQ];
                     
                     InterchangeQ=True;
                      Break[ ]) ],
                     {j, m+1, Dimensions[A][[2]]}]),
                 
                  {i,m+1,Dimensions[A][[1]]}];
                  
                  If[InterchangeQ , 
                  {nP,nA,nQ}= ExtendedEliminateRow[nP,nA,nQ]
                  ];
                  Return[{nP,nA,nQ}])];             

ExtendedSmithForm[A_] := Block[ {nP, nA=A, nQ, 
                        rows, cols,MinEntry=0, Loc,
                        ZeroMatrixQ,NonZeroRowQ=True,NonZeroColQ=
                        True},
          (
                
          (* Some Initialization*)
          {rows, cols}= Dimensions[A];
           nP= IdentityMatrix[rows];
           nQ= IdentityMatrix[cols];
           
          
                
           (* Now we start the big loop on m*)
           Do[(
           
          (*We first find the min  non-zero value in the remaining
          principal submatrix*)     
          Loc={m,m};ZeroMatrixQ=True;
          Do[ Which[nA[[i,j]] != 0 && ZeroMatrixQ, 
                       (MinEntry=Abs[ nA[[i,j]] ];Loc={i,j};
                        ZeroMatrixQ=False), 
                       nA[[i,j]] != 0 && Abs[ nA[[i,j]] ]< MinEntry,
                       (MinEntry=Abs[ nA[[i,j]] ];
                        Loc={i,j} ), 
                       True, Continue[ ] ], 
                      {i,m,rows},{j, m,cols} 
                      ]; 
           If[ZeroMatrixQ,  Break[ ]];  (*If the remaining
           principal submatrix has all zeros we are done*)        
           
           
           (*First Step: Put min value in (m,m) position *)
           
           If[Loc[[1]] != m,{nA,nP}=R3[m,Loc[[1]],{nA, nP}] ];
           If[Loc[[2]] != m, {nA,nQ}=C3[m, Loc[[2]],{nA,nQ}]];
           
           If[nA[[m,m]] < 0, {nA,nP}=R1[m,-1,{nA,nP}] ];
           
           (*Second Step:  Eliminate row and column *)
           
           If[ m< Min[rows, cols], (
           
           
           
           While[(Do[(
           If[ nA[[i,m]] !=0, (NonZeroRowQ=True;Break[ ]), (Continue[ ])]),
           {i,m+1, rows}];
           If[!NonZeroRowQ,
           Do[(If[ nA[[m,j]] !=0, (NonZeroColQ=True;
           Break[ ]), 
           (Continue[ ])]),
           {j,m+1,cols}]];
           NonZeroRowQ  || NonZeroColQ),
           
           
           If[NonZeroRowQ,{nP, nA, nQ} = ExtendedEliminateRow[nP, nA, nQ]];
           
           
           {nP, nA, nQ} = ExtendedEliminateCol[nP, nA, nQ];
           
           NonZeroRowQ=False;NonZeroColQ=False];
           
           (*Third Step: Make sure m,m entry divides the rest*)
           {nP, nA, nQ} = ExtendedThirdStep[nP, nA, nQ]
       ),
            (
            Which[rows < cols, 
            
           {nP, nA, nQ} = ExtendedEliminateCol[nP, nA, nQ], 
           rows> cols, 
           {nP, nA, nQ} = ExtendedEliminateRow[nP, nA, nQ]
           , rows==cols, (
           Break [ ] )])]
                
                
                ),{m,1, Min[rows, cols]}];
                
                Return[{nA, {nP, nQ}}]
                ) 
                ];
               
  InvariantFactors[A_]:=Module[{B,del={1},t,dims=Dimensions[A]},
  If[Max[dims]< 7,

(Do[(t=Apply[GCD, Union[Flatten[Minors[A, k]]]];
AppendTo[del,t]; j=k-1; If[t==0, (Break[ ])] ), {k,1,
Min[dims]}];
Return[Table[ If[k<j+3,del[[k]]/del[[k-1]],0],
{k,2,Min[dims]+1}]]),  
  
(B=SmithForm[A];Return[Table[B[[i,i]],{i,1,Min[dims] }]])] ] ;

ElementaryDivisors[A_]:=Module[{nA=A,T},(T=Sort[Flatten[
 FactorInteger[InvariantFactors[nA]],1]];
  Table [ T[[i,1]]^T[[i,2]], {i, 1, Length[T]}]) ]             

End[ ]

EndPackage[ ] 

