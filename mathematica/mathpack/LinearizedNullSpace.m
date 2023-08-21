(* ::Package:: *)

(*************************************************************************)
(*                         LinearizedNullSpace.m                         *)
(*************************************************************************)
(*                                                                       *)
(* Author: Renato Maia Matarazzo Orsino                                  *)
(* 2020-04-16                                                            *)
(*                                                                       *)
(*************************************************************************)



LinearizedNullSpace[xLinMatrix_, (* xE_Real:1. 10^-3, *) xNZero_Rational:1 10^-30] :=
	Module[{x, xLinMatrixCoefficients, xCoordinates, xNA1, xNC1, xNCq, xc3, xnr},
		
		{xLinMatrixCoefficients,xCoordinates} = 
			MatrixCoefficientArrays[xLinMatrix];	

		xNA1 = xLinMatrixCoefficients[1];
		xNC1 = Re @ ((* Round[#, xNZero]& @ *) NullSpace[xNA1, Tolerance -> 10^-30]);

		xc3 = {3./4., -3./20., 1./60.} / xE;
		
		xNCq = Association[
			(#-> 
				Re @ ((* Round[#, xNZero]& @ *)(
				+ xc3[[1]] * NullSpace[xNA1 + xE * xLinMatrixCoefficients[#], Tolerance -> 10^-30]
				- xc3[[1]] * NullSpace[xNA1 - xE * xLinMatrixCoefficients[#], Tolerance -> 10^-30]	
				+ xc3[[2]] * NullSpace[xNA1 + 2 * xE * xLinMatrixCoefficients[#], Tolerance -> 10^-30]
				- xc3[[2]] * NullSpace[xNA1 - 2 * xE * xLinMatrixCoefficients[#], Tolerance -> 10^-30]	
				+ xc3[[3]] * NullSpace[xNA1 + 3 * xE * xLinMatrixCoefficients[#], Tolerance -> 10^-30]
				- xc3[[3]] * NullSpace[xNA1 - 3 * xE * xLinMatrixCoefficients[#], Tolerance -> 10^-30]
				))
				)& /@ xCoordinates
			];
		
		xnr = MaximalBy[
			Join[{Dimensions @ xNC1}, Dimensions /@ (xNCq // Values)], 
     		First
     		] // DeleteDuplicates // Flatten;
		
		(xNC1 + Inner[Times, xCoordinates, (xNCq /@ xCoordinates)])
		]


LinearizedOrthogonalComplement[xLinMatrix_] :=
	Module[{xLinMatrixCoefficients, xCoordinates, xNA1, xNC1, xNCq},
		
		{xLinMatrixCoefficients,xCoordinates} = 
			MatrixCoefficientArrays[xLinMatrix];	

		xNA1 = xLinMatrixCoefficients[1];
		xNC1 = Transpose @ Re @ ((* Round[#, xNZero]& @ *) NullSpace[xNA1, Tolerance -> 10^-30]);
		
		xNCq = Association[
			(#-> - Re @ PseudoInverse[xNA1] . xLinMatrixCoefficients[#] . xNC1
				)& /@ xCoordinates
			];
		
		(xNC1 + Inner[Times, xCoordinates, (xNCq /@ xCoordinates)])
		]		
