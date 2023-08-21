(* ::Package:: *)

(*************************************************************************)
(*                           GetAllVariables.m                           *)
(*************************************************************************)
(*                                                                       *)
(* Author: Renato Maia Matarazzo Orsino                                  *)
(* 2020-04-14                                                            *)
(*                                                                       *)
(*************************************************************************)



HeadList = {
	Or, And, 
	Equal, Unequal, Inequality
	Less, LessEqual, 
	Greater, GreaterEqual
	};

GetAllVariables[xNumber_?NumericQ] := 
	Sequence[]

GetAllVariables[{}] := 
	Sequence[]

GetAllVariables[xRelationalOperator_] /; MemberQ[HeadList, xRelationalOperator] := 
	Sequence[]

GetAllVariables[x_List] := 
	DeleteDuplicates @ (Flatten @ (Union @ (GetAllVariables[#] & /@ x))) // Quiet

GetAllVariables[Derivative[xNumber_Integer][xFunction_][xArgument__]] :=
	Module[{xVariable},
   		If[MemberQ[Attributes[xFunction], NumericFunction] || MemberQ[HeadList, xFunction],
			(*-TRUE-*)
			xVariable = GetAllVariables[{xArgument}],
			(*-FALSE-*)
    		xVariable = Derivative[xNumber][xFunction][xArgument]
    		];
   		xVariable
		] // Quiet

GetAllVariables[xFunction_Symbol[xArgument__]] :=
	Module[{xVariable},
		If[MemberQ[Attributes[xFunction], NumericFunction] || MemberQ[HeadList, xFunction],
			(*-TRUE-*)
			xVariable = GetAllVariables[{xArgument}],
			(*-FALSE-*)
			xVariable = xFunction[xArgument]
			];
	xVariable
	]  // Quiet

GetAllVariables[xOther_] := 
	xOther