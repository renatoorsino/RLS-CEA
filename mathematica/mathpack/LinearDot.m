(* ::Package:: *)

(*************************************************************************)
(*                              LinearDot.m                              *)
(*************************************************************************)
(*                                                                       *)
(* Author: Renato Maia Matarazzo Orsino                                  *)
(* 2020-04-16                                                            *)
(*                                                                       *)
(*************************************************************************)



LinearToCoefficientArrays[xX_List, xRules_List: {}] :=
	If[ArrayDepth[xX] === 1, ColumnCoefficientArrays[xX, xRules], MatrixCoefficientArrays[xX, xRules]]

CoefficientArraysToLinear[xX_List] :=
	((First @ xX)[1] + Inner[Times, (Last @ xX), (First @ xX) /@ (Last @ xX)])

LinearDot[xX_List, xY_List, xToLinear_String: "Y"] :=
	Module[{xA, xB, xVariables},
		xA = LinearToCoefficientArrays[xX];
		xB = LinearToCoefficientArrays[xY];
		xVariables = Union[Last @ xA, Last @ xB];
		If[xToLinear === "Y", Expand @ CoefficientArraysToLinear[#], #]& @ {
			Association[ Union@@{
				{1-> (* Chop @ *) Expand[(First @ xA)[1].(First @ xB)[1]]},
				(#-> (* Chop @ *) Expand[((First @ xA)[#] //. Missing[xSo__] -> Array[0 &, Dimensions @ ((First @ xA)[1])]).(First @ xB)[1] 
					+ (First @ xA)[1].((First @ xB)[#] //. Missing[xSo__] -> Array[0 &, Dimensions @ ((First @ xB)[1])])])	& /@ xVariables
				}],
			xVariables
			}
		]
