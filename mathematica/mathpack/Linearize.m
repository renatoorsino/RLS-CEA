(* ::Package:: *)

(*************************************************************************)
(*                              Linearize.m                              *)
(*************************************************************************)
(*                                                                       *)
(* Author: Renato Maia Matarazzo Orsino                                  *)
(* 2020-10-04                                                            *)
(*                                                                       *)
(*************************************************************************)


LinearExpansion[xE_] = {
	Derivative[2][xX_][t] :> Derivative[2][xX]["R"] + xE Derivative[2][xX][t], 
	Derivative[1][xX_][t] :> Derivative[1][xX]["R"] + xE Derivative[1][xX][t], 
	xX_[t] :> xX["R"] + xE xX[t]
	};

Linearize[xL_, xReferenceMotion_: {}] :=
	(Series[(((xL/.LinearExpansion[xE]) /. (xReferenceMotion /. {t -> "R"})) /. 
		{xX_["R"]-> 0}), 
		{xE,0,1}] // Normal) /. {xE-> 1};	



(* -- VERSION 2020-04-16 -- *)

(* 
LinearExpansion[xE_] = {
	Derivative[2][Subscript[Subscript[xX_, xId__], xId2__]][t] ->
		Superscript[Subscript[Subscript[Overscript[xX, ".."], xId], xId2], \[EmptySmallCircle]]
		+ xE Derivative[2][Subscript[Subscript[xX, xId], xId2]][t], 
	Derivative[1][Subscript[Subscript[xX_, xId__], xId2__]][t] ->
		Superscript[Subscript[Subscript[Overscript[xX, "."], xId], xId2], \[EmptySmallCircle]]
		+ xE Derivative[1][Subscript[Subscript[xX, xId], xId2]][t], 
	Derivative[xD_][Subscript[Subscript[xX_, xId__], xId2__]][t] /; (xD > 2) ->
		Superscript[Subscript[Subscript[Superscript[xX, "(" <> (ToString @ xD) <> ")"], xId], xId2], \[EmptySmallCircle]]
		+ xE Derivative[xD][Subscript[Subscript[xX, xId], xId2]][t], 
	Subscript[Subscript[xX_, xId__], xId2__][t] ->
		Superscript[Subscript[Subscript[xX, xId], xId2], \[EmptySmallCircle]]
		+ xE Subscript[Subscript[xX, xId], xId2][t], 
	Derivative[2][Subscript[xX_, xId__]][t] ->
		Superscript[Subscript[Overscript[xX, ".."], xId], \[EmptySmallCircle]]
		+ xE Derivative[2][Subscript[xX, xId]][t], 
	Derivative[1][Subscript[xX_, xId__]][t]->
		Superscript[Subscript[Overscript[xX, "."], xId], \[EmptySmallCircle]]
		+ xE Derivative[1][Subscript[xX, xId]][t], 
	Derivative[xD_][Subscript[xX_, xId__]][t] /; (xD > 2) ->
		Superscript[Subscript[Overscript[xX, "(" <> (ToString @ xD) <> ")"], xId], \[EmptySmallCircle]]
		+ xE Derivative[xD][Subscript[xX, xId]][t], 
	Subscript[xX_, xId__][t] ->
		Superscript[Subscript[xX, xId], \[EmptySmallCircle]]
		+ xE Subscript[xX, xId][t], 
	Derivative[2][Subscript[xX_, xId__]][t] ->
		Superscript[Subscript[Overscript[xX, ".."], xId], \[EmptySmallCircle]]
		+ xE Derivative[2][Subscript[xX, xId]][t], 
	Derivative[1][Subscript[xX_, xId__]][t] ->
		Superscript[Subscript[Overscript[xX, "."], xId], \[EmptySmallCircle]]
		+ xE Derivative[1][Subscript[xX, xId]][t], 
	Derivative[xD_][Subscript[xX_, xId__]][t] /; (xD > 2) ->
		Superscript[Subscript[Overscript[xX, "(" <> (ToString @ xD) <> ")"], xId], \[EmptySmallCircle]]
		+ xE Derivative[xD][Subscript[xX, xId]][t], 
	Subscript[xX_, xId__][t] ->
		Superscript[Subscript[xX, xId], \[EmptySmallCircle]]
		+ xE Subscript[xX, xId][t], 
	Derivative[2][xX_][t] ->
		Superscript[Overscript[xX, ".."], \[EmptySmallCircle]]
		+ xE Derivative[2][xX][t], 
	Derivative[1][xX_][t] ->
		Superscript[Overscript[xX, "."], \[EmptySmallCircle]]
		+ xE Derivative[1][xX][t], 
	Derivative[xD_][xX_][t] /; (xD > 2) ->
		Superscript[Overscript[xX, "(" <> (ToString @ xD) <> ")"], \[EmptySmallCircle]]
		+ xE Derivative[xD][xX][t], 
	xX_[t] ->
		Superscript[xX, \[EmptySmallCircle]]
		+ xE xX[t]
	};

Linearize[xL_, xReferenceMotion_: {}] :=
	(Series[(((xL/.LinearExpansion[xE]) /. xReferenceMotion) /. 
		{Superscript[xX_,\[EmptySmallCircle]]-> 0}), 
		{xE,0,1}] // Normal) /. {xE-> 1}; 
*)	