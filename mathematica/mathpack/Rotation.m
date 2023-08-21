(* ::Package:: *)

(*************************************************************************)
(*                               Rotation.m                              *)
(*************************************************************************)
(*                                                                       *)
(* Author: Renato Maia Matarazzo Orsino                                  *)
(* 2020-04-16                                                            *)
(*                                                                       *)
(*************************************************************************)



Rotation = Function @ Module[{x},
	x["TransformList"] = List[##] /. {
		"x" -> (RotationTransform[#, {1,0,0}]&),
		"y" -> (RotationTransform[#, {0,1,0}]&),
		"z" -> (RotationTransform[#, {0,0,1}]&),
		"X" -> (RotationTransform[#, {1,0,0}]&),
		"Y" -> (RotationTransform[#, {0,1,0}]&),
		"Z" -> (RotationTransform[#, {0,0,1}]&)
		};
	Function[(TransformationMatrix @ (Simplify @ Inner[(#1 @ #2)&, x["TransformList"], List[##], Dot]))[[1;;3,1;;3]]]
	];

SkewToVec = If[ And @@ (Flatten @ PossibleZeroQ[# + Transpose[#]]), {#[[3,2]], #[[1,3]], #[[2,1]]}]&;

OrthoToVec  = {#[[3, 2]], #[[1, 3]], #[[2, 1]]}&;

VecToSkew = {{0, -#[[3]], #[[2]]}, {#[[3]], 0, -#[[1]]}, {-#[[2]], #[[1]], 0}}&;	

QuatToRot = {
	{#1^2-#2^2-#3^2+#4^2,2 #1 #2-2 #3 #4,2 #1 #3+2 #2 #4},
	{2 #1 #2+2 #3 #4,-#1^2+#2^2-#3^2+#4^2,2 #2 #3-2 #1 #4},
	{2 #1 #3-2 #2 #4,2 #2 #3+2 #1 #4,-#1^2-#2^2+#3^2+#4^2}
	}&;	