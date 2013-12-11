(*
fEdmondsNew[UlL_List] := Module[{gg, edg, ver}, gg = fAllCycles[UlL];
    gg = If[Length[gg] == 1, {gg[[1]], gg[[1]]}, gg];
    gg = fSelectCycles[UlL, gg];
    gg = Table[
        Table[If[gg[[2, j, i]] > 2, 
            Select[KSubsets[gg[[1, i]], gg[[2, j, i]]], 
              Max[Union[Apply[Plus, #]]] <= 2 &], 
            KSubsets[gg[[1, i]], gg[[2, j, i]]]], {i, 
            Length[gg[[2, j]]]}], {j, Length[gg[[2]]]}];
    gg = Map[fSelParts, gg];
    gg = Complement[gg, {{}}];
     gg = Complement[Union[Table[fSelAll[gg, i], {i, Length[gg]}]], {{}}];
    
    gg = Table[
        Map[Rest, 
          Map[FindCycle, 
            Map[FromUnorderedPairs, 
              Table[Complement[
                  Union[Table[
                      If[SameQ[gg[[k, j, i]], 1], UlL[[i]], {}], {i, 
                        Length[gg[[k, j]]]}]], {{}}], {j, 
                  Length[gg[[k]]]}]]]], {k, Length[gg]}];
    edg = Length[UlL];
    ver = Last[Union[Flatten[UlL]]];
    gg = Table[{Sort[Table[gg[[j, i]], {i, Length[gg[[j]]]}]], 
          ver + Length[gg[[j]]] - edg, 
          Divide[2 - ver - Length[gg[[j]]] + edg, 2]}, {j, Length[gg]}];
    gg]
*)
