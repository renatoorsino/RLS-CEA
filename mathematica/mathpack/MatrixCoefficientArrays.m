(* ::Package:: *)

(*************************************************************************)
(*                       MatrixCoefficientArrays.m                       *)
(*************************************************************************)
(*                                                                       *)
(* Author: Renato Maia Matarazzo Orsino                                  *)
(* 2020-04-16                                                            *)
(*                                                                       *)
(*************************************************************************)



MatrixCoefficientArrays[xA_List, xRules_List: {}] :=
	Module[{xxMatrix, xxVariables, xxCoefficientMatrices},
		xxMatrix = xA //. xRules;
		xxVariables = Union @ GetVariables[xxMatrix];
		xxCoefficientMatrices = If[
			xxVariables == {}, 
			CoefficientArrays[xxMatrix, EMPTY],
			CoefficientArrays[xxMatrix, xxVariables]
			];

		{
		Association[ Union@@{
			{1-> Normal@Part[xxCoefficientMatrices,1]},
			MapThread[ 
				(#1-> Normal@Part[xxCoefficientMatrices,2,All,All,#2])&, 
				{#, Range@Length@#}, 
				1
				]& @ xxVariables
			}],
		xxVariables
		}
		]

ColumnCoefficientArrays[xA_List, xRules_List: {}] :=
	Module[{xxColumn, xxVariables, xxCoefficientMatrices},
		xxColumn = xA //. xRules;
		xxVariables = Union @ GetVariables[xxColumn];
		xxCoefficientMatrices = If[
			xxVariables == {}, 
			CoefficientArrays[xxColumn, EMPTY],
			CoefficientArrays[xxColumn, xxVariables]
			];
			
		{
		Association[ Union@@{
			{1-> Normal@Part[xxCoefficientMatrices,1]},
			MapThread[ 
				(#1-> Normal@Part[xxCoefficientMatrices,2,All,#2])&, 
				{#, Range@Length@#}, 
				1
				]& @ xxVariables
			}],
		xxVariables
		}
		]	