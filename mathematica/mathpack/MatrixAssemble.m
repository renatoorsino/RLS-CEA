(* ::Package:: *)

(*************************************************************************)
(*                            MatrixAssemble.m                           *)
(*************************************************************************)
(*                                                                       *)
(* Author: Renato Maia Matarazzo Orsino                                  *)
(* 2020-04-14                                                            *)
(*                                                                       *)
(*************************************************************************)



MatrixPartition[xM_List, xp_Integer: 2]:= AssociationThread[ 
	Flatten[Outer[List, #, #]& @ Range[xp], 1],
	Flatten[Partition[xM, {#,#}& @ (Mean @ Dimensions[xM])/xp], 1]
	];

MatrixDuplicateDiagonal[xM_List, xp_Integer: 2]:= ArrayFlatten[
	(Outer[List, #, #] & @ Range[xp]) /. (MatrixPartition[xM] /. (xL_List :> ArrayFlatten[{{xL, 0}, {0, xL}}]))
	];

MatrixDiagonal2[xL_List]:= ArrayFlatten[{{xL, 0}, {0, xL}}]

MatrixDiagonalAssemble[xM1_List, xM2_List, xp_: 1/2]:= Module[{zdim, zM1, zM2},
	zdim = (Mean @ Dimensions[xM2])*(1-xp);
	zM1 = ArrayFlatten[{{xM1, 0}, {0, ConstantArray[0, {zdim, zdim}]}}];
	zdim = (Mean @ Dimensions[xM1]) - zdim;
	zM2 = ArrayFlatten[{{ConstantArray[0, {zdim, zdim}], 0}, {0, xM2}}];
	zM1+zM2
	];

FEMMatrixAssemble[nM_, ne_]:= Fold[
	MatrixDiagonalAssemble, 
	MatrixDuplicateDiagonal /@ (nM[#] & /@ Range[ne])
	];

MatrixDiagonal[xM1_List, xM2_List, xp_: 1/2]:= Module[{zdim, zM1, zM2},
	zdim = (Last @ Dimensions[xM1]) + (Last @ Dimensions[xM2])*(1-xp);
	zM1 = PadRight[#, zdim]& /@ xM1;
	zM2 = PadLeft[#, zdim]& /@ xM2;
	Join[zM1, zM2]
	];

JacobiMatrixAssemble[nM_, ne_]:= Fold[
	MatrixDiagonal, 
	(nM[#] & /@ Range[ne])
	];	

VectorPartition[xM_List, xp_Integer: 2]:= AssociationThread[ 
	Range[xp], 
	Partition[xM, (Length[xM])/xp]
	];	

VectorInter[xM1_List, xM2_List, xp_Integer: 2]:= 
	Flatten @ Values @ MapThread[List, {VectorPartition[xM1,xp], VectorPartition[xM2,xp]}]

VectorAssemble[xv1_List, xv2_List, xp_: 1/2]:= Module[{zdim, zv1, zv2},
	zdim = Length[xv1] + Length[xv2]*(1-xp);
	zv1 = PadRight[xv1, zdim];
	zv2 = PadLeft[xv2, zdim];
	zv1+zv2
	];	

FEMVectorAssemble[nM_, ne_]:= Fold[
	VectorAssemble, 
	(nM[#] & /@ Range[ne])
	];