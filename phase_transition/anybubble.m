(* ::Package:: *)

(* Find bubble nucleation instantons *)

(* Copyright 2016, 2017, 2018, 2019 Ali Masoumi, Ken D. Olum, and Benjamin Shlaer

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*)

Needs["PowellHybrid`", "powell.m"]
BeginPackage["AnyBubble`"]

AnyBubbleVersion = "2.0.4"

FindBubble::usage = "FindBubble[potential, fieldname, truevacuum, falsevacuum] computes instanton profiles for multidimensional fields"

Options[FindBubble] = {Verbose -> False,
		       SpaceTimeDimension -> 4,
		       StartAnalyticFraction -> 1/100,
		       EndAnalyticFraction -> 1/100,
		       MaxIntervalGrowth -> 30,
		       MaxReadjustments -> 20,
		       InitialProfilePoints-> {},
		       PowellVerbosity -> 1, (* Here and below passed to PowellHybrid *)
		       MaxIterations -> 500,
		       AccuracyGoal -> Automatic,
		       WorkingPrecision -> MachinePrecision}

SpaceTimeDimension::usage = "SpaceTimeDimension is an option for FindBubble, giving the dimension of spacetime"
StartAnalyticFraction::usage = "StartAnalyticFraction is an option for FindBubble, giving the fraction of the way from the true to the false vacuum in field space that can be safely handled analytically"
EndAnalyticFraction::usage = "EndAnalyticFraction is an option for FindBubble, giving the fraction of the way from the false to the true vacuum in field space that can be safely handled analytically"
MaxIntervalGrowth::usage = "MaxIntervalGrowth is an option for FindBubble, giving the maximum exponential growth allowed in a single shooting interval"
MaxReadjustments::usage = "MaxReadjustments is an option for FindBubble, giving the maximum number of times that they position of the shooting intervals may be readjusted"
InitialProfilePoints::usage = "InitialProfilePoints is an option for FindBubble, giving a list of points in field space through which the initial guess should go from the true vacuum to the false"
PowellVerbosity::usage = "powellVerbosity is an option for FindBubble.  It is passed as the Verbosity option to PowellHybrid"

(* Errors *)
FindBubble::npos = "The true vacuum must be lower than the false vacuum"
FindBubble::negativeEigenvalue = "An eigenvalue of the Hessian at the starting position was negative"
FindBubble::tooManyReadjustments = "Failed to find solution after `` readjustments"
FindBubble::readjustmentFailed = "Attempting to readjust parameters left us with too few intermediate points"


Begin["`Private`"]
$ContextPath = Prepend[$ContextPath, "PowellHybrid`"]

(* We are solving the equation f''(r) + ((d-1)/r) f'(r) + dU/df = 0, where r is the "time" coordinate,
   representing the radius of the bubble, f is the field, and U is the potential *)

(* potential is a function of variables named field[1], field[2], etc. *)
FindBubble[potential_, fieldname_, inputTrueVacuum_, inputFalseVacuum_, OptionsPattern[]] :=
 (* Install many dynamic variables *)
  Block[{verbose = OptionValue[Verbose],
	 rDim = OptionValue[SpaceTimeDimension],
	 startAnalyticFraction = OptionValue[StartAnalyticFraction],
	 endAnalyticFraction = OptionValue[EndAnalyticFraction],
	 maxIntervalGrowth = OptionValue[MaxIntervalGrowth],
	 initialProfilePoints = OptionValue[InitialProfilePoints],
	 maxReadjustments = OptionValue[MaxReadjustments],
	 powellVerbosity = OptionValue[PowellVerbosity],
	 accuracyGoal = OptionValue[AccuracyGoal],
	 maxIterations = OptionValue[MaxIterations],
	 workingPrecision = OptionValue[WorkingPrecision],
	 fieldDim = Length[inputFalseVacuum],
	 inputU, (* functional form of given potential, with floating-point numbers *)
	 Uf,	 (* internal potential at false vacuum *)
	 fieldScale, (* scale we multiply by to externalize field *)
	 uScale, (* scale we multiply by to externalize potential *)
	 middleShooting, initialMotion,	 (* See note below on dynamic caching *)
	 (* Variables installed by setupTransform *)
	 falseVacuum, trueVacuum, U, internalizeField, externalizeField, internalizeU, externalizeU, initialPath,
	 (* Variables and functions installed by setupDerivatives *)
	 gradient, hessian, falseVacuumEigenvalues, falseVacuumEigenvectors, drMax,
	 bareRi,				 (* Divisions between shooting segments before offsetting *)
	 deltaScale,				 (* scale factor for offsetting r_i *)
	 guess1, answer, profile}, 
    If[accuracyGoal == Automatic, accuracyGoal = N[workingPrecision]/2 - 1];
    initialMotion[r1_, field1_] := initialMotion1[r1, field1];
    middleShooting[r1_, r2_, field1_, field1Prime_] := middleShooting1[r1, r2, field1, field1Prime];
    setupInputU[potential, fieldname];
    setupTransform[inputTrueVacuum, inputFalseVacuum];
    Uf = rePrecision[U[falseVacuum]];
    If[Uf <= 0, Message[FindBubble::npos]];
    setupDerivatives[];
    {bareRi, guess1} = rePrecision[initialGuess[]];
    deltaScale = Total[bareRi]/Length[bareRi];
    guess1 = Join[guess1,{0}];	(* Add delta parameter *)
    (* For debugging:
    Begin["AnyBubble`Private`"];
    Print["Debugging dialog: you are in the AnyBubble`Private` context"];
    Dialog[];
    End[];*)
    answer = findSolution[guess1];
    If[answer === "Failed", "Failed",
       profile = profileFunction[answer];
       {calculateAction[profile, answer], inputProfileFunction[profile]} (* Return action and profile function *)
    ]]

(* Note on dynamic caching for functions initialMotion and middleShooting:
   These functions cache their results, so that when they are called as part of Jacobian calculations
   they can just returned a cached result, rather than needing to do their computations again.
   I think it's better programming practice to store cached data (such as U) in dynamically scoped variables
   rather than globals, and for the same reason, I think cached function calls should be stored dynamically.
   The Block above makes these functions dynamically bound, so cached values are discarded when the block
   exits.  But it also means that global definitions are not accessible inside the Block.  So the actual
   code for some function F is defined below as function F1.  Inside the block, F with arbitrary
   arguments is defined to call F1.  Then F1 is defined to memoize the call in F, so that if F is
   called again with the same arguments, the result will be returned immediately rather than calling F1 *)

(* Set prototype definition for defineFieldFunction once and for all now.  This makes sure that field can never
   be localized into field$ or anything like that, which might happen if the prototype inside
   defineFieldFunction *)
defineFieldFunctionPrototype = Hold[defineFieldFunctionName[field_List] := defineFieldFunctionRhs]

(* This defines a function with the given name and one argument called field.  It operates only when field is a list,
   which could be either a numerical field point or a list of symbols such as {f[1], f[2], ...}.  The body of the
   function is the form  with entries such as field[1] replaced by field[[1]].  In the case of a list of symbols,
   the result is returned symbolically and can be differentiated.  Numerical quantities in the form are
   turned into floating point numbers with precision workingPrecision *)
defineFieldFunction[name_, form_] :=
  (* form is an expression involving things such as field[1].  The declaration below keeps this from turning into
     field[1.0].  We put the name of the function we're trying to define and the form into the prototype.
     Then we change field[1] to field[[1]], etc., in the held definition, and perform it.
     We can't just use "=" in the function definition because that would try to evaluate the right hand
     side first including taking parts of the symbol "field".  Thus := and Hold *)
  ReleaseHold[defineFieldFunctionPrototype /. {defineFieldFunctionName -> name,
					       defineFieldFunctionRhs -> N[form, workingPrecision]}
	      /. field[n_] :> field[[n]]]

SetAttributes[field, NHoldAll]	(* See defineFieldFunction *)

(* Define the function inputU.  It accepts a list of numbers giving a field point and returns the potential there.
   Alternatively, it can accept the list of symbols such as {f[1], f[2], ...} and return the expression in
   terms of those, which can be symbolically differentiated.  It does not, however, do anything unless the argument
   is a list.  Thus a function call which has not yet been evaluated such as path[lambda] will not be taken
   apart and used in the potential. *)
setupInputU[potential_, fieldname_] :=
  defineFieldFunction[inputU, potential /. fieldname -> field] (* Change name of field to our symbol named "field" *)

fieldList := Table[field[i], {i, 1, fieldDim}] (* Get list of field variables *)
fieldPrimeList := Table[field[i]', {i, 1, fieldDim}] (* And their derivatives *)
fieldListr := Map[#[r]&, fieldList] (* Get list of field variables as functions of r *)

(* Put numbers in our precision *)
rePrecision[x_] := SetPrecision[x, workingPrecision]

(* Instead of the input potential in terms of the input field space, we will construct a potential
   with the following properties:
   The true vacuum is at 0 and the false vacuum is at (1,0,0...).
   The potential at the true vacuum is 0 and the maximum value reached by the potential on the straight line
   or the user-provided path from true to false vacuum is 1.

   These choices mean that the typical field value and field derivative are of order 1
   
   Sets up dynamic variables:
   U -- new potential
   fieldScale -- what to multiply by to externalize field: the distance between the input vacua
   externalizeField -- TransformationFunction that takes a field vector in the internal coordinates and
     returns vector in the input coordinates
   internalizeField -- TransformationFunction that takes a field vector in input coordinates and gives internal version
   uScale -- what to multiply by to externalize potential: the range of the input U on the initial path
   externalizeU -- function that takes potential value in internal format and gives external version
   internalizeU -- function that takes potential value has given by user and returns internal version
   trueVacuum -- location of true vacuum = (0,0,0,...)
   falseVacuum -- location of false vacuum = (1,0,0,...)
   initialPath -- function on 0..1: straight line from true to false vacuum, or spline through user's points
   *)
setupTransform[inputTrueVacuum_, inputFalseVacuum_] :=
  Block[{$MinPrecision = workingPrecision}, (* All intermediate calculations at this precision *)
    Module[{offset = inputFalseVacuum - inputTrueVacuum, (* Vector from true to false *)
	    e1 = UnitVector[fieldDim, 1],
	    inputUt = inputU[inputTrueVacuum], (* Value of potential at true vacuum *)
	    initialPoints, 
	    path,
	    uMax},
      (*  Create transform that takes internal field coordinates as decribed above to external.  It maps
	  e1 into inputFalseVacuum *)
      fieldScale = Norm[offset];
      externalizeField = AffineTransform[{If[fieldDim==1, {offset}, (* 1x1 matrix to transform 1 to offset *)
					RotationMatrix[{e1, offset}] (* Matrix that rotates e1 to dir of offset *)
					* fieldScale], (* Scale to length of offset *)
				     inputTrueVacuum}]; (* Shift origin to inputTrueVacuum  *)
      internalizeField = InverseFunction[externalizeField]; (* Transform that takes external coordinates to internal *)
      trueVacuum = Table[0, {fieldDim}]; (* True vacuum at origin in simplified coordinate system *)
      falseVacuum = e1; (* Location of false vacuum in simplified coordinate system *)
      (* Now construct potential transform *)
      initialPoints = (* Points in initial profile in the rotated coordinates.
			 By default it is a straight line from true to false vacuum *)
	   Join[{trueVacuum},
		Map[internalizeField, initialProfilePoints], (* Convert user's points to internal coordinate *)
		{falseVacuum}];
      (* Interpolated path, a pure function with 0 giving the true vacuum and 1 the false vacuum *)
      initialPath = Interpolation[MapIndexed[{(#2-1)/(Length[initialPoints]-1), (* x = 0, 1/(N-1), 2/(N-1), ... 1 *)
				     #1}&, (* y = vector value at this point *)
				    initialPoints],
			   InterpolationOrder -> Min[3, Length[initialPoints]-1]]; (* Cubic spline unless too few pts *)
      (* Range of input potential along initial path.  *)
      (* It is not important to compute this accurately, since the point is only to rescale things in a rough manner.
	 If you set $MinPrecision you have to set workingPrecision in FindMaximum to at least that value, and there's
	 no need for that.  So we use low precision here and increase the precision of the resulting number so that
	 transformations are done precisely. *)
      Block[{$MinPrecision = MachinePrecision},
         uScale = rePrecision[FindMaximum[{inputU[externalizeField[initialPath[lambda]]], 0 < lambda < 1}, lambda,
					  PrecisionGoal -> 3][[1]]
			      - inputUt]];
      (* Convert user's potential to internal form with values 0<U<1 on the initial path *)
      internalizeU = Function[(#-inputUt)/uScale];
      externalizeU = Function[# uScale + inputUt]; (* Convert internal potential to user's form *)
      (* U is a function giving the scaled and shifted potential as a function of internal coordinates.  *)
      U = Function[internalizeU[inputU[externalizeField[#]]]];
    ]]

(* Distance from true to false vacuum is just 1 due to rescaling *)
falseVacuumNorm = 1

nr := Length[bareRi]		(* Number of intervals *)

(* Give ri from bareRi by offsetting by delta, the last parameter in the list *)
offsetRi[parameters_] := rePrecision[bareRi + deltaScale parameters[[-1]]]

(* Set up derivatives, the eigensystem of the Hessian at the true and false vacua, and a guess at the
   limit on interval size *)
setupDerivatives[] :=
  (* We reprecision the function itself, which means that any number stored there are turned into high precision.
     This avoids trouble from NDSolve. *)
  Module[{trueVacuumEigenvalues},
    defineFieldFunction[gradient, D[U[fieldList], {fieldList}]]; (* Gradient of the potential *)
    defineFieldFunction[hessian, D[gradient[fieldList], {fieldList}]]; (* Hessian matrix *)
    {falseVacuumEigenvalues, falseVacuumEigenvectors} = Eigensystem[hessian[falseVacuum]];
    trueVacuumEigenvalues = Eigenvalues[hessian[trueVacuum]];
    (* Maximum r-interval width to prevent solutions that grow as exp(sqrt{ev} r) from growing by more
       then given limit in one interval.  Of course we don't know whether we will encounter a larger
       eigenvalue somewhere else in the potential landscape *)
    drMax = Log[maxIntervalGrowth]/Sqrt[Max[Join[falseVacuumEigenvalues,trueVacuumEigenvalues]]]]

(* Figure out a guess at the initial profile by solving a modified one-dimensional problem on the initial path.
   Returns a list of the ri, the initial parameters, and the value for scaling delta (which is just the thin wall
   radius, and thus a typical value for the r_i *)
initialGuess[]:=
  Module[{uAlongThePath, (* Potential along the path as a function of the parameter NOT the arclength *)
	  oneDFalseVacuum, (* Path length distance of the false vacuum from the true vacuum *)
	  thinWallR, guessR, guessRPrime, r1, rN, n,
	  ri,			(* In this function, ri is a list of points *)
	  lambda0 = 0, 		(* lambda value corresponding to r=0 *)
	  centerLambda,
	  lambda1, lambdaN, n2,
	  pathLength
	 },
    (* lambda = 0..1 along the path. This is not the pathlength *)
    oneDFalseVacuum = NIntegrate[Norm[initialPath'[lambda]], {lambda, 0, 1}]; (* Path length to false vacuum *)
    (* Set centerLambda to the lambda position in the center of the path by length *)
    NDSolve[{pathLength'[lambda]==Norm[initialPath'[lambda]], pathLength[0]==0, 
	     WhenEvent[pathLength[lambda] == oneDFalseVacuum / 2, {centerLambda=lambda ,"StopIntegration"}]},
	    pathLength, {lambda, 0, 1}];
    (* Adjust potential so it has equal minima *)
    uAlongThePath[lambda_?NumericQ]:= U[initialPath[lambda]] + 12 Uf (lambda^4/4 - lambda^3/3);
    (* In thin-wall approximation, surface tension can be computed by integration.  Here we integrate with respect
      to path length. *)
    thinWallR = Module[{surfaceTension = NIntegrate[Sqrt[2 Max[0, uAlongThePath[lambda]]] Norm[initialPath'[lambda]],
						   {lambda, 0, 1}, AccuracyGoal-> 3]},
		 (* Determine R by setting surface energy and volume energy equal *)
		 (rDim-1) surfaceTension / Uf];
    lambda1 = startAnalyticFraction/2 falseVacuumNorm/Norm[initialPath'[0]]; (* Only go down to this value *)
    lambdaN = 1 - endAnalyticFraction/2 falseVacuumNorm/Norm[initialPath'[1]]; (* and only up to this value *)
    (* Find r[phi] for offset problem by integration.  We integrate in both directions from centerLambda *)
    guessR = NDSolveValue[{r'[lambda] == Norm[initialPath'[lambda]]/Sqrt[2 uAlongThePath[lambda]],
			   r[centerLambda] == thinWallR, (* Choose center of profile at thin wall radius *)
			   (* Stop if r reaches 0 and set field0 to field value there *)
			   WhenEvent[r[lambda] == 0, {lambda0=lambda,"StopIntegration"}]},
			  r,	(* Return function giving r[phi] *)
			  {lambda, lambda1, lambdaN}];
    (* Possibly reset field1 if we reached r=0 above *)
    lambda1 = lambda0 + startAnalyticFraction/2 falseVacuumNorm/Norm[initialPath'[0]];
    r1 = guessR[lambda1];
    rN = guessR[lambdaN];
    n2 = Ceiling[(rN-r1) / drMax]; (* Number of division points exclusively between 1 and N. *)
    ri = Table[r1 + (rN-r1) (i-1)/(n2+1), {i, n2+2}]; (* Evenly spaced r points *)
    {ri,
     Flatten[{initialPath[lambda1],	(* phi_1 *)
	      (* phi_i, phi'_i at n-3 intermediate points *)
              Table[Module[{lambda = lambda /. FindRoot[ri[[i]] == guessR[lambda],
							{lambda,
							 (* Start with linear guess *)
							 lambda1 + (ri[[i]]-r1)/(rN-r1) (lambdaN-lambda1),
							 lambda1, lambdaN}]}, (* Stay in range *)
			   {initialPath[lambda], (* phi at r_i *)
			    initialPath'[lambda]/guessR'[lambda]}], (* dphi/dr = (dphi/dlambda)/(dr/dlambda) *)
		    {i, 2, n2}],			      (* Internal parameters at N-3 = n2-1 points *)
	      initialPath[lambdaN]}]}]				      (* phi_N *)



(* Near r = 0, we use an analytic approximation rather than soling the differential equations
   numerically, for two reasons.  First, the appearance of 1/r in the differential equation leads to
   numerical problems for r near 0.  Second, since the solution leaves exponentially from the true
   vacuum, if we want to stay near the true vacuum for significant time, we have to start exponentially
   close. *)

(* The common part of the various functions below.  If derivativeQ is True, then the order of the Bessel
   function is increased by 1 for the derivative *)
initialBessel[r_, derivativeQ_, eigenvalue_] :=
  r^(1-rDim/2) If[eigenvalue<0,
		    If[derivativeQ, -1, 1] (* The derivative formula has a different sign for J *)
		    * BesselJ[rDim/2-If[derivativeQ, 0, 1], Sqrt[-eigenvalue] r],
		  BesselI[rDim/2-If[derivativeQ, 0, 1], Sqrt[eigenvalue] r]]

(* Special case for r = 0 and no derivative to give initial value without errors.  The case with B < 0
   and so J instead of K is asymptotically the same except for the sign in the square root *)
initialBessel0[eigenvalue_] := (Sqrt[Abs[eigenvalue]]/2)^(rDim/2-1) / Gamma[rDim/2]

(* Automatically thread these functions over lists of eigenvalues *)
SetAttributes[initialBessel, Listable]
SetAttributes[initialBessel0, Listable]

(* Given r_1 and phi_1, we return the value of phi_0 (but note that this may be have numerical error),
   a function giving the field from 0 to r1, the derivative at r1, and
   a Jacobian matrix giving the dependence of the derivative on the function value.  For the Jacobian,
   we imagine that the potential is expanded to quadratic order around the initial field point, once and for all,
   and so ignore changes to the eigenvectors and eigenvalues due to a different initial point.
   We cache the result to use for finding the Jacobian later.  See note above about dynamic caching.  *)
initialMotion1[r1_, field1_] :=
  initialMotion[r1, field1] =
  Module[{grad, eigenvectors, eigenvalues, solution, derivativepart,
	  r},			(* Avoid possible premature evaluation of r in Evaluate below *)
    {eigenvalues, eigenvectors} = Eigensystem[hessian[field1]];
    grad = eigenvectors . gradient[field1]; (* Gradient of potential at r1 in eigenvector coordinates *)
    (* Function giving field(r) - field(r1) in eigenvector coordinates *)
    solution = Function[r,
		Evaluate[grad/eigenvalues
			 *(initialBessel[r, False, eigenvalues]/initialBessel[r1, False, eigenvalues] - 1)]];
    (* Get ratio of derivative to field value entering into derivative and Jacobian *)
    derivativepart = initialBessel[r1, True, eigenvalues]/initialBessel[r1, False, eigenvalues];
    {Transpose[eigenvectors] . (grad/eigenvalues
				(initialBessel0[eigenvalues]/initialBessel[r1, False, eigenvalues] - 1))
       + field1, (* phi_0 *)
     Function[r, Evaluate[Transpose[eigenvectors] . solution[r] + field1]],
     (* The argument to the Bessel function is always sqrt{|B|}, and the prefactor is A/B.  Differentiation
	brings out the sqrt{|B|}, and by keeping this in the form sqrt{|B|}/B, we get the right sign *)
     Transpose[eigenvectors] . (Sqrt[Abs[eigenvalues]] grad/eigenvalues derivativepart), (* derivative at r_1 *)
     Transpose[eigenvectors] . (Sqrt[Abs[eigenvalues]] derivativepart eigenvectors)} (* Jacobian *)
  ]

(* Near the false vacuum, we use analytic solutions based on the approximation
   U(phi) = Uf - H(phi-falsevacuum), where H is the Hessian at the false vacuum.   Our goal is to arrive at
   the false vacuum at r = infinity.  In one dimension, these solutions would be the false vacuum plus
   a multiple of f(r) = r^(1-d/2) K_{d/2-1}(sqrt{H} r), where K is the Bessel function.  In coordinates
   where the Hessian is diagonal, we simply have such solutions in each variable.

   What we need to know here is the ratio f'(r)/f(r), so that if we have chosen the value of f, we
   can determine the corresponding value of f'.  That is given by
   -sqrt{H} K_{d/2}(sqrt{H} r)/K_{d/2-1}(sqrt{H} r).
   
   With more dimensions in the basis that diagonalizes the Hessian, each coordinate is simply multiplied by
   the above function based on the corresponding eigenvalue.  In a general basis, we first take then
   inner product of f with the eigenvector matrix (whose rows are eigenvectors) to transform to the
   diagonal basis, multiply by the appropriate numbers and transform back with the transpose matrix.
   *)

(* Return the asymptotic solution for the given field at the boundary as a function of r *)
(* This does not check whether fieldBoundary is actually near the false vacuum, but if it is too far, the ri will
   be readjusted later *)
asymptoticSolution[rBoundary_, fieldBoundary_] :=
  (* Convert desired field at the boundary into coefficients for the Bessel functions in the
     eigenvector basis *)
  Module[{r},			(* Make sure r is not evaluated inside function body *)
   Function[r, Evaluate[Transpose[falseVacuumEigenvectors]
			. (MapThread[If[#2 == 0, 0,
					#2 / (rBoundary^(1-rDim/2) BesselK[rDim/2-1, Sqrt[#1] rBoundary])
					r^(1-rDim/2) BesselK[rDim/2-1, Sqrt[#1] r]]&,
				     {falseVacuumEigenvalues, falseVacuumEigenvectors . (fieldBoundary - falseVacuum)}])
			+ falseVacuum]]]

(* Create matrix that converts boundary fields, with the false vacuum subtracted, into boundary field
   derivatives *)
asymptoticSolutionMatrix[r_] :=
  Dot[Transpose[falseVacuumEigenvectors],
      Map[- Sqrt[#] BesselK[rDim/2, Sqrt[#] r]/BesselK[rDim/2-1, Sqrt[#] r] &, falseVacuumEigenvalues]
      * falseVacuumEigenvectors]

(* Shooting procedure. *)

(* Give the differential equations that we are trying to solve in terms of functions named field[i]
   and independent variable r.  Arguments are the starting value of r and the initial fields and
   their derivatives *)
fieldEquations[r1_, field1_, field1Prime_] :=
  Block[{functionList = Table[field[i][r], {i, 1, fieldDim}]}, (* vector of functions of r *)
   Table[{field[i]''[r] + (rDim-1)/r field[i]'[r] == gradient[functionList][[i]],
	  field[i][r1] == rePrecision[field1[[i]]],
	  (D[field[i][r],r]/.r->r1) == rePrecision[field1Prime[[i]]]},
	 {i,1,fieldDim}]]

(* Basic shooting function.  Shoot from one r value to another.  We accept the values of the field
   and its derivatives at the initial point, and return the field at r2, its derivative there,
   and a function giving the fields as functions of r.
   We can shoot either toward increasing or decreasing r.
   We cache our results to use for finding the Jacobian later.  See note above about dynamic caching.  *)
middleShooting1[r1_, r2_, field1_, field1Prime_] :=
  middleShooting[r1, r2, field1, field1Prime] =
  Module[{solution = Null, field, derivative},
    If[verbose, Print["shooting from r =", r1, " to r = ", r2, ", starting field ", field1,
		    ", starting derivative ", field1Prime]];
    Catch[
      With[{fieldFunctions = fieldListr}, (* Don't figure this out each time *)
        (* Solve equations, return solution for fields and solution for their derivatives *)
	{field, derivative, solution} =
	   Check[NDSolveValue[fieldEquations[r1, field1, field1Prime],
			      (* Return field at r2, derivative at r2, field functions *)
			      {Map[#[r2]&,fieldList], Map[#[r2]&,fieldPrimeList], fieldList},
			      {r,r1,r2},
			      StepMonitor :> checkTooFar[fieldFunctions, field1],
			      WorkingPrecision -> workingPrecision],
		 {Null,Null,Null}, (* Set everything to Null if solution fails *)
		 {NDSolveValue::ndsz, NDSolveValue::precw}]],
      shootingTooFar];		(* Throw here if we get too far away *)
    If[solution === Null,
       If[verbose, Print["Shooting went too far away"]];Throw[Null, shootingFailed],
       {field, derivative,
	Function[r, Evaluate[Map[#[r]&, solution]]]}]] (* Convert list of solutions into single function *)

doNDSolveValue[args__] := (Print["NDSolveValue[",Map[InputForm,{args}],"]"]; NDSolveValue[args])

(* If this step has taken us more than twice the distance from the true to the false vacuum, abort *)
checkTooFar[field_, field1_] :=
  If[Norm[field - field1] > falseVacuumNorm, Throw[0, shootingTooFar]]

(* As above, but shoot outward from the center.  Accepts the field value at the first internal point r1,
   computes the derivative there, and then shoots outward from this point to r2.  Returns field at r2,
   derivative at r2, function giving fields from r1 to r2, function giving fields from 0 to r1. *)
(* This does not check whether field1 is actually near field0, but if it is too far, the ri will
   be readjusted later *)
firstShooting[r1_, r2_, field1_] :=
  Module[{initialSolution, initialDerivative, initialJacobian, field0},
    If[verbose, Print["shooting from center out to r =", r2, ", starting field ", field1]];
    {field0, initialSolution, initialDerivative, initialJacobian} =
	 rePrecision[initialMotion[r1, field1]];
    Append[middleShooting[r1, r2, field1, initialDerivative], (* Shoot the rest of the way *)
	   initialSolution]]

(* As above, but shoot inward from the boundary.  We accept the value of the field at the boundary position r,
   and compute the value of the field derivative there, then shoot back to r2.  *)
lastShooting[rBoundary_, r2_, fieldBoundary_] :=
  (If[verbose, Print["shooting from boundary at r = ", rBoundary, " in to ", r2,
		   ", starting field ", fieldBoundary]];
   middleShooting[rBoundary, r2, fieldBoundary,
		  rePrecision[asymptoticSolutionMatrix[rBoundary] . (fieldBoundary - falseVacuum)]])

(* Find the Jacobian giving the effect of changing the starting field and derivatives in shooting on the
   ending field and derivatives, by solving a differential equation.  Accepts the same arguments as
   middleShooting and returns the Jacobian as a 2D*2D matrix. *)
middleJacobian[r1_, r2_, field1_, field1Prime_] :=
  (* First get solution, usually cached. *)
  Module[{solution = middleShooting[r1, r2, field1,field1Prime][[3]], 
	  aux},
    Block[{QMatrix = ArrayFlatten[{{ConstantArray[0, {fieldDim,fieldDim}], IdentityMatrix[fieldDim]},
				   {hessian[solution[r]], -(rDim-1)/r IdentityMatrix[fieldDim]}}],
	   AuxFields = Array[aux[#1,#2][r]&, {2 fieldDim, 2 fieldDim}], (* Dummy variables for NDSolve *)
	   RHandSideMatrix, JacobianEqns, Jsol},
      RHandSideMatrix = QMatrix . AuxFields;
      (* Equations for the evolution of the Jacobian *)
      JacobianEqns = Flatten[Array[{D[aux[#1,#2][r], r] == RHandSideMatrix[[#1,#2]], (* evolution *)
				    aux[#1,#2][r1] == KroneckerDelta[#1,#2]}&, (* initial condition *)
				   {2 fieldDim, 2 fieldDim}]];
      (* Solve equations and return Jacobian evaluated at final position *)
      NDSolveValue[JacobianEqns, AuxFields /. r->r2, {r,r1,r2}, Method->{"EquationSimplification"->"Solve"},
		   WorkingPrecision -> workingPrecision]]]

(* Return a matrix of 2D rows and D columns giving the dependence of phi_2 and phi'_2 on phi_1. *)
firstJacobian[r1_, r2_, field1_] :=
  Module[{field0, solution, derivative, jacobian, j},
    (* Get cached initial solution *)
    {field0, solution, derivative, jacobian} = initialMotion[r1, field1];
    (* Get dependence of phi_2 and phi'_2 on phi_1 and phi_1' *)
    j = middleJacobian[r1, r2, field1, derivative];
    (* Now include dependence of phi_1' on phi_1 *)
    j . Join[IdentityMatrix[fieldDim], jacobian]]

(* Return a matrix of 2D rows and D columns giving the dependence of phi_{N-1} and phi'_{N-1} on
   the boundary value phi_N. *)
lastJacobian[rBoundary_, r2_, fieldBoundary_] :=
  Module[{j},
    (* Get dependence of phi_{N-1} and phi'_{N-1} on phi_N and phi_N' *)
    j = middleJacobian[rBoundary, r2, fieldBoundary,
		       rePrecision[asymptoticSolutionMatrix[rBoundary] . (fieldBoundary - falseVacuum)]];
    (* Add direct dependences on phi_N with dependence on phi_N' times dependence of that on phi_N *)
    j . Join[IdentityMatrix[fieldDim], asymptoticSolutionMatrix[ri[[nr]]]]]

(* Calculate result vector from parameters.
   Let N (nr in the code) be the number of division points r_i.  The minimum value is 3.
   The first and last points, r_1 and r_N must be in the intervals near the true and false vacuum
   where analytic approximations can be used.  But r_2 and r_{N-1} should not be in these regions.
   If r_1 is too small, they will not be possible to represent phi_1 accurately enough to avoid
   trouble with exponential growth on the way from r_1, and similarly at the other end.  Thus
   the final solution will be given in terms of N+1 intervals, the first and last of which are
   analytic approximations, while the N-1 intervals are computed numerically.
   We solve 2N-4 equations in 2N-4 unknown parameters. Each parameter and equation is a d-dimensional vector.
   
   The parameters are:
   The field value phi_1.  We don't need phi_0, and we compute phi'_1 from phi_1 analytically.
   The starting point and starting derivative phi_i, phi'_i at the first N-3 intermediate points (if any)
   The ending (boundary) point phi_N, where we transition to an analytic solution.  This
   analytic solution gives phi'_N in terms of phi_N.
   
   There are 2 equations associated with each of the N-2 intermediate points:
   At the first N-3 intermediate points, phi_i and phi'_i computed from the left should agree with
   the parameters phi_i and phi'_i.
   At the last intermediate point, phi_(N-1) and phi'_(N-1) computed by shooting from the right should
   agree with the same values computed by shooting from the left.
   *)

(* Accepts a list of 2N-4 vector-valued parameters as above, returns 2N-4 vector-valued results that
   we're trying to make zero, and caches data in the function shootingSaves to compute the Jacobian. *)
shooting1[parameters_] :=
  (* Make a list of fields and derivatives coming from the left at points 2...N-1 *)
  (Flatten[Table[If[i==1, firstShooting[ri[[1]], ri[[2]], parameters[[1]]], (* phi_1 gives phi_2, phi'_2 *)
		   middleShooting[ri[[i]], ri[[i+1]], (* From r_i to r_{i+1} *)
				    parameters[[2i-2]], (* field at r_i *)
				    parameters[[2i-1]]]] (* derivative at r_i *)
		 [[1;;2]],				 (* Keep only field and derivative *)
		 {i, nr-2}],
	  1]			(* Flatten first level to get one list with both fields and derivatives *)
   (* Now subtract solutions and derivatives at points 2...N-2 given by parameters, and those at N-1
      coming from the right. *)
   - Join[If[nr>3, parameters[[2;;2nr-5]], {}],
	  lastShooting[ri[[nr]], ri[[nr-1]], parameters[[2nr-4]]][[1;;2]]])

(* This is the function actually given to Powell's method, so it takes 2DN scalar parameters and returns
   2DN scalar values.  The condition keeps it from trying to run with symbolic input when called from
   Plot *)

shooting[flatparameters_ /; Apply[And, Map[NumberQ, flatparameters]]] :=
  Block[{ri = offsetRi[flatparameters]}, (* Offset by delta parameter *)
    Catch[Flatten[shooting1[Partition[flatparameters[[1;;-2]], fieldDim]]],
	  shootingFailed]]

(* Return a symbol defined to be the function giving the profile as a function r in our internal coordinates.
   The function triggers only on real number arguments, so it can be used for integration.  We can't use a pure
   function, because those don't have conditions. *)
profileFunction[flatparameters_] :=
  Module[{ri = offsetRi[flatparameters], (* Store values of ri in unique symbol *)
	  parameters = Partition[flatparameters[[1;;-2]], fieldDim],
	  functions, i,
	  nr = nr,		(* Preserve outside binding of bareRi *)
	  profile},		(* Unique name for function *)
      With[{field0 = initialMotion[ri[[1]], parameters[[1]]][[1]]}, (* Cache r=0 answer *)
	(* Make a list of functions for the N+1 intervals delimited by the N values of ri *)
        functions = Join[firstShooting[ri[[1]], ri[[2]], parameters[[1]]][[{4,3}]], (* 0 to 2 *)
			Table[middleShooting[ri[[i]], ri[[i+1]], parameters[[2i-2]], parameters[[2i-1]]][[3]],
			      {i, 2, nr-2}],
			{lastShooting[ri[[nr]], ri[[nr-1]], parameters[[2nr-4]]][[3]], (* n-1 to n *)
			 asymptoticSolution[ri[[nr]], parameters[[2nr-4]]]}]; (* n to infinity *)
	profile[r_/;Element[r,Reals]] := (* Don't try If unless a real number *)
	  If[r == 0, field0,	     (* Don't get error from initial analytic function *)	
     Catch[Do[If[r < ri[[i]], Throw[functions[[i]][r]]], {i, nr}];
		   functions[[nr+1]][r]]]];
     profile]			(* Return symbol *)
	 
(* Return a symbol defined to give the profile function in the input coordinates *)
inputProfileFunction[internalProfile_] :=
  Module[{profile,		(* Name of function *)
	  transform = externalizeField,
	  rescaleR = Sqrt[uScale] / fieldScale}, (* Must stretch profile to account for mods to potential *)
     profile[r_/;Element[r,Reals]] := transform[internalProfile[rescaleR r]];
     profile]
 
(* Return the Jacobian matrix giving the dependence of the results on the parameters.  This depends on
   the solutions in the shooting intervals, and those should be cached.
   We construct a matrix made of DxD blocks.  The columns are the
   parameters phi_1, phi_2, phi'_2, ... phi_N and the rows are the functions, which we can write
   phi_2, phi'_2, ... phi_N-1, phi'_N-1.  The matrix for N=6 thus looks like
   Am000000 where A and B give the dependence of phi_2, phi'_2 on phi_1,
   B0m00000 	  C and D that of phi_3 on phi_2 and phi'_2 while E and F give the dependence of phi'_3 on the same
   0CDm0000	  G, H, I, J are the same with dependence of quantities at 4 on those at 2
   0EF0m000       K, L, M, N are the same with dependence of quantities at 5 on those at 3
   000GHm00       P is minus the dependence of phi_5 coming from the right on phi_6, and 
   000IJ0m0       Q is minus the dependence of phi'_5 coming from the right on phi_6.    
   00000KLP       and m is minus the identity matrix.
   00000MNQ *)
shootingJacobian1[parameters_] :=
  Module[{z22 = ConstantArray[0, {2 fieldDim, 2 fieldDim}], (* 2Dx2D array of zeroes *)
	  z21 = ConstantArray[0, {2 fieldDim, fieldDim}],   (* 2Dx1D *)
	  z12 = ConstantArray[0, {fieldDim, 2 fieldDim}],   (* 1Dx2D *)
	  nid = -IdentityMatrix[2 fieldDim]},		    (* 2Dx2D negative identity matrix *)
    ArrayFlatten[		(* Assemble blocks that give the structure above, then flatten *)
     Table[Join[If[i == 1, {firstJacobian[ri[[1]], ri[[2]], parameters [[1]]]}, (* Left part up through AB *)
		   Join[{z21}, Table[z22, {i-2}],	      (* or DF, HJ, etc. *)
			{middleJacobian[ri[[i]], ri[[i+1]], parameters[[2i-2]], parameters[[2i-1]]]}]],
		If[i < nr-2, Join[{nid}, Table[z22, {nr-i-3}], {z21}], (* Right part *)
		   {-lastJacobian[ri[[nr]], ri[[nr-1]], parameters[[2nr-4]]]}]], (* Right part, bottom rows *)
	   {i, 1, nr-2}]]]

(* Assemble parameters into vectors *)
shootingJacobian[flatparameters_] :=
  Block[{ri = offsetRi[flatparameters]}, (* Offset by delta parameter *)
    Module[{parameters = Partition[flatparameters[[1;;-2]], fieldDim],
	   (* {phi_1, phi_2, phi'_2, phi_(nr-2), phi'_(nr-2), phi'_N} *)
	   jac, endPointData, endPointPhiDDPrime, f1, f2, firstTerm, lastTerm, middleTerms, jac1,
	   allTerms, initPointPhiDDPrime},
      (* phi_r1'(phi_r1). The derivative at phi1 is determined as a function of phi_1 *)
      f1 = rePrecision[initialMotion[ri[[1]],parameters[[1]]][[2]]'[ri[[1]]]]; 
      (* phi_nr'(phi_nr). The derivative at rn is determined as a function of phi_nr *)
      f2 = rePrecision[asymptoticSolution[ri[[nr]],parameters[[-1]]]'[ri[[nr]]]];
      (* This is a set of parameters obtained from shooting in the form
	 {{phi2, phi'2},{phi3,phi'3}, ...,{phi_(N-1)(Left), phi'_(N-1)(Left)},{phi_(N-1)(Right), phi'_(N-1)(Right)}} *)
      endPointData = Table[If[i==1, firstShooting[ri[[1]],ri[[2]],parameters[[1]]],
			     If[i==nr-1, lastShooting[ri[[nr]], ri[[nr-1]],parameters[[-1]]],
				middleShooting[ri[[i]],ri[[i+1]], parameters[[2i-2]], parameters[[2i-1]]]]] [[1;;2]],
			   {i,1,nr-1}];
      initPointPhiDDPrime = (* These are {phi''_r2, phi''_r3, phi''_r(N-2)} *)
	Table[gradient[parameters[[2i]]]- (rDim-1)/ri[[i+1]] parameters[[2i+1]], {i,nr-3}];
      endPointPhiDDPrime = (* {phi''_2, phi''_3, phi''_(N-1) Left, phi''_(N-1)(Right)} *)
        Table[gradient[endPointData[[i,1]]] - (rDim-1)/If[i!= nr-1, ri[[i+1]], ri[[i]]] endPointData[[i,2]], {i,nr-1}];
      jac1 = shootingJacobian1[parameters]; (* Jacobian without extra column *)
      jac = Partition[jac1, {fieldDim,fieldDim}]; (* This is the N*N Jacobian *)
      firstTerm = {endPointData[[1,2]] - jac[[1,1]] . f1, endPointPhiDDPrime[[1]] - jac[[2,1]] . f1};
      lastTerm = {endPointData[[-1,2]] + jac[[2nr-5,2nr-4]] . f2, endPointPhiDDPrime[[-1]] + jac[[2nr-4,2nr-4]] . f2};
      middleTerms = rePrecision[Table[{endPointData[[i,2]] - jac[[2i-1,2i-2]] . parameters[[2i-1]]
				       - jac[[2i-1,2i-1]] . initPointPhiDDPrime[[i-1]], 
				       endPointPhiDDPrime[[i]]- jac[[2i,2i-2]] . parameters[[2i-1]]
				       - jac[[2i,2i-1]] . initPointPhiDDPrime[[i-1]]}
				      - If[i!=nr-2, 0, lastTerm]
				      ,{i,2,nr-2}]];
      allTerms = Flatten[Join[firstTerm,middleTerms]];
      Join[jac1, deltaScale Transpose[{allTerms}], 2] (* Append additional column *)
    ]]

(* Find the solution using Powell's method *)
findSolution[firstGuess_] :=
  Module[{guess = firstGuess,
	  answer = Null},
    Catch[
      Do[guess = Catch[		(* If parameters need adjusting, exit this catch with new parameters *)
		     answer = PowellHybrid[shooting, guess, shootingJacobian, Verbosity -> powellVerbosity,
					   MaxIterations -> maxIterations, AccuracyGoal -> accuracyGoal,
					   WorkingPrecision -> workingPrecision,
					   StepMonitorFunction -> checkRi];
		     Break[],	(* Exit Do on success *)
		   readjusted],	(* Throw to this tag after readjusting to reenter PowellHybrid *)
	 {maxReadjustments+1}];	(* Readjust given number of time, plus once more to find solution *)
      (* Here on success or if too many readjustments *)
      If[answer === Null, Message[FindBubble::tooManyReadjustments, maxReadjustments]];
      answer,
     findSolutionFailed]]

(* If Ri need to be readjusted, throw out of Powell loop to do that *)
checkRi[parameters_] :=
  With[{newparameters = maybeReadjust[parameters]}, (* Return new parameters or Null if no changes needed *)
    If[newparameters =!= Null, Throw[newparameters, readjusted]]] (* If something to do, exit Powell loop *)

(* See if readjustment is needed.  If it is, change Ri and return new parameters.  If not, don't change global state
   and return Null *)
maybeReadjust[flatparameters_]:=
  Block[{ri = offsetRi[flatparameters]}, (* Offset dynamic variable by delta parameter *)
    maybeReadjust1[Partition[flatparameters[[;;-2]], fieldDim]]]

maybeReadjust1[parameters_] :=	   
  Module[{field1 = parameters[[1]], (* Current phi1 *)
	  fieldN = parameters[[-1]], (* Current phiN *)
	  startNorm = startAnalyticFraction falseVacuumNorm, (* Distance in field space for start analytic region *)
	  endNorm = endAnalyticFraction falseVacuumNorm, (* Same for ending region *)
	  field0, function1, derivative1,
	  readjustStart, readjustEnd},
    {field0, function1, derivative1} = initialMotion[ri[[1]], field1][[1;;3]]; (* Get info about current profile *)
    (* See what needs readjusting *)
    readjustStart = !(startNorm > Norm[field1 - field0] > startNorm/10); (* Start too close or too far *)
    readjustEnd = !(endNorm > Norm[fieldN - falseVacuum] > endNorm/10); (* Too close or too far from false vacuum *)
    If[!(readjustStart || readjustEnd), (* Nothing to do? *)
       Null,				(* Just give Null *)
      (* We're going to readjust and then CheckRi will throw out of Powell loop.  The variable ri is
	 dynamically bound to the offset ri values.  We will modify it and store it in bareRi *)
      Module[{functionN, r1, rN, newParameters = parameters,
	      data, nr = nr, (* Use variable that we can modify *)
	      rPosition}, (* Position of inserted r in the old ri *)
	If[verbose, Print["Readjusting.  Old ri = ", ri, ", old parameters =  ", parameters ]];
	If[readjustStart,
	   r1 = Abs[r] /. (* The function is even in r, so we allow it to get on the wrong branch and fix it later *)
		 FindRoot[Norm[function1[r] - field0] - startNorm/2, {r,ri[[1]]}]; (* Find new r1 value *)
	   (* Calculation does not need to be done with high accuracy, but now we need a precise number to use *)
	   r1 = rePrecision[r1];
	   If[r1<ri[[1]],  		(* New r value goes before old r1? *)
	      ri = Join[{r1},ri];	(* Yes.  Prepend new r1 *)
	      newParameters = Join[{function1[r1]}, (* New phi1 *)
				   {field1,	   (* Old phi1, now phi2 *)
				    derivative1},	   (* Old phi1', now phi2' *)
				   newParameters[[2;;]]], (* Rest of parameters *)
	     (* New r value comes after old r1.  We will discard all the old r values less that the new r1 *)
	     rPosition = FirstPosition[ri, r_/;r>r1][[1]]; (* Position of first to keep *)
	     (* If only one r is left, we would only have 2 points after adjusting.  We do not attempt to handle the
		case where we would add a new rN below and so create a third point. *)
	     If[rPosition >= nr, Message[FindBubble::readjustmentFailed]; Throw["Failed", findSolutionFailed]];
	     ri = Join[{r1}, ri[[rPosition;;-1]]]; (* r1, then other kept ri *)
	     newParameters = Join[{function1[r1]}, (* New phi1 *)
				  newParameters[[(2 rPosition - 2);;]]];	(* Parameters from remaining ri *)
	     nr = Length[ri] (* reset variable to account adjustment *)
	   ]]; 
	(* Now consider issues near false vacuum *)
	If[readjustEnd,
	   functionN = asymptoticSolution[ri[[nr]], newParameters[[-1]]];
	   rN = r /. FindRoot[Norm[functionN[r] - falseVacuum] - endNorm/2, {r,ri[[nr]]}]; (* choose new rN *)
	   rN = rePrecision[rN];
	   If[rN>ri[[nr]],	  (* New rN goes after old rN? *)
	      ri = Join[ri,{rN}]; (* Yes, append it *)
	      newParameters =
		Join[newParameters[[;;-2]], (* Parameters remaining after previous step, except last *)
		     lastShooting[ri[[nr]], ri[[nr-1]], fieldN] [[1;;2]], (* Add field and derivative at old r_(N-1) *)
		     {functionN[rN]}], (* Finally, new phiN *)
	     (* New rN goes before old rN *)
	     rPosition = FirstPosition[ri, r_/;r>rN][[1]]; (* New rN goes in place of this element and rest of list *)
	      (* We must have 3 points to proceed *)
	     If[rPosition <= 2, Message[FindBubble::readjustmentFailed]; Throw["Failed", findSolutionFailed]];
	     ri = Join[ri[[;;rPosition-1]],{rN}];
	     newParameters = Join[newParameters[[;; (2 rPosition-5)]],
				  {functionN[rN]}]]];
	If[verbose, Print["New ri = ", ri, ", new parameters =  ", newParameters]];
	bareRi = ri;			 (* Store new ri as bare values *)
	deltaScale = Total[bareRi]/Length[bareRi]; (* Readjust scaling of offset parameter *)
	Join[Flatten[newParameters],{0}] (* Convert to single list, add delta = 0 *)
	]]]

calculateAction[profile_, answer_] :=
  Block[{ri = offsetRi[answer]}, (* Offset dynamic variable by delta parameter *)
     Module[{f},
       f[r_?NumberQ] := - r^(rDim-1) profile[r] . (gradient[profile[r]]-gradient[falseVacuum]);
       Pi^(rDim/2)/Gamma[1+rDim/2] NIntegrate[f[r], Join[{r, 0}, ri, {Infinity}],
					      WorkingPrecision -> workingPrecision,
					      PrecisionGoal -> N[workingPrecision]/2,
					      MaxRecursion -> N[workingPrecision]/2]
       * fieldScale^rDim / uScale^(rDim/2-1)  (* Undo rescalings to give action of input potential *)
     ]]
(* Check that all symbols are bound in Block so they don't have values after exit *)
(* Names returns a set of strings.  But if we pass the string to Symbol, the result will be evaluated.
   Thus the need for string manipulation *)
checkBlock[] :=
  (Begin["AnyBubble`Private`"]; (* In the package at runtime, so we get short names *)
   Scan[checkName, Names["AnyBubble`Private`*"]] 
  End[];)

(* Name is a string in the AnyBubble`Private` package.  See if the symbol has a leftover value. *)
checkName[name_] :=
  If[ToExpression[StringJoin["ValueQ[", name, "]"]],
     (* Check list of and variables generated by module in profile functions *)
     If[!StringMatchQ[name,
		      RegularExpression["nr|fieldList|fieldListr|fieldPrimeList|falseVacuumNorm|defineFieldFunctionPrototype|ri\$[0-9]*|nr\$[0-9]*|functions\$[0-9]*|transform\$[0-9]*"]],
	Print[name, " has a value"]]]

End[]
EndPackage[]
