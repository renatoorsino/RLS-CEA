(* ::Package:: *)

(*************************************************************************)
(*                             GetVariables.m                            *)
(*************************************************************************)
(*                                                                       *)
(* Author: Renato Maia Matarazzo Orsino                                  *)
(* 2020-04-14                                                            *)
(*                                                                       *)
(*************************************************************************)



GetVariables[xX_List, xExcept_List: {}] := 
	Complement[ DeleteDuplicates @ Cases[xX, xVariable_[t], Infinity], xExcept]

GetVariables[xX_Association, xExcept_List: {}] := 
	Complement[ DeleteDuplicates @ Cases[xX["Matrix"], xVariable_[t], Infinity], xExcept]
