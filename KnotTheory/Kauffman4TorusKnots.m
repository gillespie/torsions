BeginPackage["KnotTheory`Kauffman4TorusKnots`", {"KnotTheory`"}]
Message[KnotTheory::loading, "Kauffman4TorusKnots`"]

Begin["`Private`"]

KOs = {
  TorusKnot[3, 2], TorusKnot[5, 2], TorusKnot[7, 2], TorusKnot[4, 3],
  TorusKnot[9, 2], TorusKnot[5, 3], TorusKnot[11, 2], TorusKnot[13, 2],
  TorusKnot[7, 3], TorusKnot[5, 4], TorusKnot[15, 2], TorusKnot[8, 3]
}

Is = {-a^(-4) - 2/a^2 + z/a^5 + z/a^3 + z^2/a^4 + z^2/a^2, 
 2/a^6 + 3/a^4 + z/a^9 - z/a^7 - (2*z)/a^5 + z^2/a^8 - (3*z^2)/a^6 - 
  (4*z^2)/a^4 + z^3/a^7 + z^3/a^5 + z^4/a^6 + z^4/a^4, 
 -3/a^8 - 4/a^6 + z/a^13 - z/a^11 + z/a^9 + (3*z)/a^7 + z^2/a^12 - 
  (2*z^2)/a^10 + (7*z^2)/a^8 + (10*z^2)/a^6 + z^3/a^11 - (3*z^3)/a^9 - 
  (4*z^3)/a^7 + z^4/a^10 - (5*z^4)/a^8 - (6*z^4)/a^6 + z^5/a^9 + z^5/a^7 + 
  z^6/a^8 + z^6/a^6, -a^(-10) - 5/a^8 - 5/a^6 + (5*z)/a^9 + (5*z)/a^7 + 
  (10*z^2)/a^8 + (10*z^2)/a^6 - (5*z^3)/a^9 - (5*z^3)/a^7 - (6*z^4)/a^8 - 
  (6*z^4)/a^6 + z^5/a^9 + z^5/a^7 + z^6/a^8 + z^6/a^6, 
 4/a^10 + 5/a^8 + z/a^17 - z/a^15 + z/a^13 - z/a^11 - (4*z)/a^9 + z^2/a^16 - 
  (2*z^2)/a^14 + (3*z^2)/a^12 - (14*z^2)/a^10 - (20*z^2)/a^8 + z^3/a^15 - 
  (3*z^3)/a^13 + (6*z^3)/a^11 + (10*z^3)/a^9 + z^4/a^14 - (4*z^4)/a^12 + 
  (16*z^4)/a^10 + (21*z^4)/a^8 + z^5/a^13 - (5*z^5)/a^11 - (6*z^5)/a^9 + 
  z^6/a^12 - (7*z^6)/a^10 - (8*z^6)/a^8 + z^7/a^11 + z^7/a^9 + z^8/a^10 + 
  z^8/a^8, 2/a^12 + 8/a^10 + 7/a^8 - (8*z)/a^11 - (8*z)/a^9 - z^2/a^12 - 
  (22*z^2)/a^10 - (21*z^2)/a^8 + (14*z^3)/a^11 + (14*z^3)/a^9 + 
  (21*z^4)/a^10 + (21*z^4)/a^8 - (7*z^5)/a^11 - (7*z^5)/a^9 - (8*z^6)/a^10 - 
  (8*z^6)/a^8 + z^7/a^11 + z^7/a^9 + z^8/a^10 + z^8/a^8, 
 -5/a^12 - 6/a^10 + z/a^21 - z/a^19 + z/a^17 - z/a^15 + z/a^13 + (5*z)/a^11 + 
  z^2/a^20 - (2*z^2)/a^18 + (3*z^2)/a^16 - (4*z^2)/a^14 + (25*z^2)/a^12 + 
  (35*z^2)/a^10 + z^3/a^19 - (3*z^3)/a^17 + (6*z^3)/a^15 - (10*z^3)/a^13 - 
  (20*z^3)/a^11 + z^4/a^18 - (4*z^4)/a^16 + (10*z^4)/a^14 - (41*z^4)/a^12 - 
  (56*z^4)/a^10 + z^5/a^17 - (5*z^5)/a^15 + (15*z^5)/a^13 + (21*z^5)/a^11 + 
  z^6/a^16 - (6*z^6)/a^14 + (29*z^6)/a^12 + (36*z^6)/a^10 + z^7/a^15 - 
  (7*z^7)/a^13 - (8*z^7)/a^11 + z^8/a^14 - (9*z^8)/a^12 - (10*z^8)/a^10 + 
  z^9/a^13 + z^9/a^11 + z^10/a^12 + z^10/a^10, 
 6/a^14 + 7/a^12 + z/a^25 - z/a^23 + z/a^21 - z/a^19 + z/a^17 - z/a^15 - 
  (6*z)/a^13 + z^2/a^24 - (2*z^2)/a^22 + (3*z^2)/a^20 - (4*z^2)/a^18 + 
  (5*z^2)/a^16 - (41*z^2)/a^14 - (56*z^2)/a^12 + z^3/a^23 - (3*z^3)/a^21 + 
  (6*z^3)/a^19 - (10*z^3)/a^17 + (15*z^3)/a^15 + (35*z^3)/a^13 + z^4/a^22 - 
  (4*z^4)/a^20 + (10*z^4)/a^18 - (20*z^4)/a^16 + (91*z^4)/a^14 + 
  (126*z^4)/a^12 + z^5/a^21 - (5*z^5)/a^19 + (15*z^5)/a^17 - (35*z^5)/a^15 - 
  (56*z^5)/a^13 + z^6/a^20 - (6*z^6)/a^18 + (21*z^6)/a^16 - (92*z^6)/a^14 - 
  (120*z^6)/a^12 + z^7/a^19 - (7*z^7)/a^17 + (28*z^7)/a^15 + (36*z^7)/a^13 + 
  z^8/a^18 - (8*z^8)/a^16 + (46*z^8)/a^14 + (55*z^8)/a^12 + z^9/a^17 - 
  (9*z^9)/a^15 - (10*z^9)/a^13 + z^10/a^16 - (11*z^10)/a^14 - 
  (12*z^10)/a^12 + z^11/a^15 + z^11/a^13 + z^12/a^14 + z^12/a^12, 
 5/a^16 + 16/a^14 + 12/a^12 - (16*z)/a^15 - (16*z)/a^13 - (10*z^2)/a^16 - 
  (76*z^2)/a^14 - (66*z^2)/a^12 + (60*z^3)/a^15 + (60*z^3)/a^13 + 
  (6*z^4)/a^16 + (138*z^4)/a^14 + (132*z^4)/a^12 - (78*z^5)/a^15 - 
  (78*z^5)/a^13 - z^6/a^16 - (122*z^6)/a^14 - (121*z^6)/a^12 + 
  (44*z^7)/a^15 + (44*z^7)/a^13 + (55*z^8)/a^14 + (55*z^8)/a^12 - 
  (11*z^9)/a^15 - (11*z^9)/a^13 - (12*z^10)/a^14 - (12*z^10)/a^12 + 
  z^11/a^15 + z^11/a^13 + z^12/a^14 + z^12/a^12, 
 a^(-18) + 9/a^16 + 21/a^14 + 14/a^12 - z/a^19 - (8*z)/a^17 - (28*z)/a^15 - 
  (21*z)/a^13 - z^2/a^18 - (22*z^2)/a^16 - (91*z^2)/a^14 - (70*z^2)/a^12 + 
  (14*z^3)/a^17 + (84*z^3)/a^15 + (70*z^3)/a^13 + (21*z^4)/a^16 + 
  (154*z^4)/a^14 + (133*z^4)/a^12 - (7*z^5)/a^17 - (91*z^5)/a^15 - 
  (84*z^5)/a^13 - (8*z^6)/a^16 - (129*z^6)/a^14 - (121*z^6)/a^12 + z^7/a^17 + 
  (46*z^7)/a^15 + (45*z^7)/a^13 + z^8/a^16 + (56*z^8)/a^14 + (55*z^8)/a^12 - 
  (11*z^9)/a^15 - (11*z^9)/a^13 - (12*z^10)/a^14 - (12*z^10)/a^12 + 
  z^11/a^15 + z^11/a^13 + z^12/a^14 + z^12/a^12, 
 -7/a^16 - 8/a^14 + z/a^29 - z/a^27 + z/a^25 - z/a^23 + z/a^21 - z/a^19 + 
  z/a^17 + (7*z)/a^15 + z^2/a^28 - (2*z^2)/a^26 + (3*z^2)/a^24 - 
  (4*z^2)/a^22 + (5*z^2)/a^20 - (6*z^2)/a^18 + (63*z^2)/a^16 + 
  (84*z^2)/a^14 + z^3/a^27 - (3*z^3)/a^25 + (6*z^3)/a^23 - (10*z^3)/a^21 + 
  (15*z^3)/a^19 - (21*z^3)/a^17 - (56*z^3)/a^15 + z^4/a^26 - (4*z^4)/a^24 + 
  (10*z^4)/a^22 - (20*z^4)/a^20 + (35*z^4)/a^18 - (182*z^4)/a^16 - 
  (252*z^4)/a^14 + z^5/a^25 - (5*z^5)/a^23 + (15*z^5)/a^21 - (35*z^5)/a^19 + 
  (70*z^5)/a^17 + (126*z^5)/a^15 + z^6/a^24 - (6*z^6)/a^22 + (21*z^6)/a^20 - 
  (56*z^6)/a^18 + (246*z^6)/a^16 + (330*z^6)/a^14 + z^7/a^23 - (7*z^7)/a^21 + 
  (28*z^7)/a^19 - (84*z^7)/a^17 - (120*z^7)/a^15 + z^8/a^22 - (8*z^8)/a^20 + 
  (36*z^8)/a^18 - (175*z^8)/a^16 - (220*z^8)/a^14 + z^9/a^21 - (9*z^9)/a^19 + 
  (45*z^9)/a^17 + (55*z^9)/a^15 + z^10/a^20 - (10*z^10)/a^18 + 
  (67*z^10)/a^16 + (78*z^10)/a^14 + z^11/a^19 - (11*z^11)/a^17 - 
  (12*z^11)/a^15 + z^12/a^18 - (13*z^12)/a^16 - (14*z^12)/a^14 + z^13/a^17 + 
  z^13/a^15 + z^14/a^16 + z^14/a^14, -7/a^18 - 21/a^16 - 15/a^14 + 
  (21*z)/a^17 + (21*z)/a^15 + (21*z^2)/a^18 + (126*z^2)/a^16 + 
  (105*z^2)/a^14 - (105*z^3)/a^17 - (105*z^3)/a^15 - (21*z^4)/a^18 - 
  (294*z^4)/a^16 - (273*z^4)/a^14 + (189*z^5)/a^17 + (189*z^5)/a^15 + 
  (8*z^6)/a^18 + (346*z^6)/a^16 + (338*z^6)/a^14 - (157*z^7)/a^17 - 
  (157*z^7)/a^15 - z^8/a^18 - (222*z^8)/a^16 - (221*z^8)/a^14 + 
  (65*z^9)/a^17 + (65*z^9)/a^15 + (78*z^10)/a^16 + (78*z^10)/a^14 - 
  (13*z^11)/a^17 - (13*z^11)/a^15 - (14*z^12)/a^16 - (14*z^12)/a^14 + 
  z^13/a^17 + z^13/a^15 + z^14/a^16 + z^14/a^14}

Is = Function /@ (Is /. {a -> #1, z-> #2})

MapThread[(Kauffman[#1] = #2)&, {KOs, Is}]

Clear[KOs, Is]

End[]; EndPackage[]

