(* ::Package:: *)

(*************************************************************************)
(*                          LinearizedInverse.m                          *)
(*************************************************************************)
(*                                                                       *)
(* Author: Renato Maia Matarazzo Orsino                                  *)
(* 2020-04-16                                                            *)
(*                                                                       *)
(*************************************************************************)



LinearizedInverse[xLinMatrix_, xE_Real:1. 10^-3, xNZero_Rational:1 10^-8] :=
	Module[{x, xLinMatrixCoefficients, xCoordinates, xNA1, xNC1, xNCq, xc3},
		
		{xLinMatrixCoefficients,xCoordinates} = 
			MatrixCoefficientArrays[xLinMatrix];	

		xNA1 = xLinMatrixCoefficients[1];
		xNC1 = Re @ ((* Round[#, xNZero]& @ *) Inverse[xNA1]);

		xc3 = {3./4., -3./20., 1./60.} / xE;
		
		xNCq = Association[
			(#-> Re @ ((* Round[#, xNZero]& @ *)(
				+ xc3[[1]] * Inverse[xNA1 + xE * xLinMatrixCoefficients[#]]
				- xc3[[1]] * Inverse[xNA1 - xE * xLinMatrixCoefficients[#]]	
				+ xc3[[2]] * Inverse[xNA1 + 2 * xE * xLinMatrixCoefficients[#]]
				- xc3[[2]] * Inverse[xNA1 - 2 * xE * xLinMatrixCoefficients[#]]	
				+ xc3[[3]] * Inverse[xNA1 + 3 * xE * xLinMatrixCoefficients[#]]
				- xc3[[3]] * Inverse[xNA1 - 3 * xE * xLinMatrixCoefficients[#]]
				))
				)& /@ xCoordinates
			];
		
		N[xNC1 + Inner[Times, xCoordinates, xNCq /@ xCoordinates]]
		]