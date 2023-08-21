(* ::Package:: *)

(*************************************************************************)
(*                       LinearizedPseudoInverse.m                       *)
(*************************************************************************)
(*                                                                       *)
(* Author: Renato Maia Matarazzo Orsino                                  *)
(* 2020-04-16                                                            *)
(*                                                                       *)
(*************************************************************************)



LinearizedPseudoInverse[xLinMatrix_, xE_Real:2. 10^-4, xNZero_Rational:1 10^-8] :=
	Module[{x, xLinMatrixCoefficients, xCoordinates, xNA1, xNC1, xNCq, xc3},
		
		{xLinMatrixCoefficients,xCoordinates} = 
			MatrixCoefficientArrays[xLinMatrix];	

		xNA1 = xLinMatrixCoefficients[1];
		xNC1 = Re @ ((* Round[#, xNZero]& @ *) PseudoInverse[xNA1]);

		xc3 = {3./4., -3./20., 1./60.} / xE;
		
		xNCq = Association[
			(#-> Re @ ((* Round[#, xNZero]& @ *)(
				+ xc3[[1]] * PseudoInverse[xNA1 + xE * xLinMatrixCoefficients[#]]
				- xc3[[1]] * PseudoInverse[xNA1 - xE * xLinMatrixCoefficients[#]]	
				+ xc3[[2]] * PseudoInverse[xNA1 + 2 * xE * xLinMatrixCoefficients[#]]
				- xc3[[2]] * PseudoInverse[xNA1 - 2 * xE * xLinMatrixCoefficients[#]]	
				+ xc3[[3]] * PseudoInverse[xNA1 + 3 * xE * xLinMatrixCoefficients[#]]
				- xc3[[3]] * PseudoInverse[xNA1 - 3 * xE * xLinMatrixCoefficients[#]]
				))
				)& /@ xCoordinates
			];
		
		(xNC1 + Inner[Times, xCoordinates, xNCq /@ xCoordinates])
		]



LinearizedLinearSolve[xM_, xB_] :=
	Module[{xxCoefficientMatrices, xCoordinates, xCM, xCB, xCY1, xCY},
		
		xCoordinates = Union[GetVariables[xM], GetVariables[xB]];


		xxCoefficientMatrices = CoefficientArrays[xM, xCoordinates];

		xCM = Association[ Union@@{
			{1-> Normal@Part[xxCoefficientMatrices,1]},
			MapThread[ 
				(#1-> Normal@Part[xxCoefficientMatrices,2,All,All,#2])&, 
				{#, Range@Length@#}, 
				1
				]& @ xCoordinates
			}];

		xxCoefficientMatrices = CoefficientArrays[xB, xCoordinates];

		xCB = Association[ Union@@{
			{1-> Normal@Part[xxCoefficientMatrices,1]},
			MapThread[ 
				(#1-> Normal@Part[xxCoefficientMatrices,2,All,#2])&, 
				{#, Range@Length@#}, 
				1
				]& @ xCoordinates
			}];	

		xCY1 = LinearSolve[xCM[1], xCB[1]];

		xCY = Association[ Union@@{
			{1-> xCY1},
			MapThread[ 
				(#1-> LinearSolve[xCM[1], xCB[#1]-xCM[#].xCY1])&, 
				{#, Range@Length@#}, 
				1
				]& @ xCoordinates
			}];	
		
		(xCY1 + Inner[Times, xCoordinates, xCY /@ xCoordinates])
		]