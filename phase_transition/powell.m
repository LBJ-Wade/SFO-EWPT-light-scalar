(* Solve multiple non-linear equations using Powell's hybrid method.
   See M. J. D. Powell, "A Hybrid Method for Nonlinear Equations", in
   Numerical Methods for Nonlinear Algebraic Equations by Philip Rabinowitz *)

(* Copyright 2014-2017 Ken D. Olum

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

(* initialParameters is a list of initial values for the parameters.
   function is a function that accepts a list of parameters and returns a list of values.
     There can be fewer values than parameters.  In this case we try to find some solution
     nearby the initial guess.
   Jacobian is a function that accepts a list of parameters and returns the Jacobian matrix there.
   The code for computing the Jacobian numerically has not been written *)

BeginPackage["PowellHybrid`"]

PowellHybridVersion = "1.1.1"

(* 1.0.1: Add more info about combining steps
   1.1: working precision
   1.1.1: return "Failed" on initial failure *)


Options[PowellHybrid] = {AccuracyGoal -> Automatic, InitialStep -> PureGradient, MaxStep -> Automatic,
			 Verbosity -> 0, MaxIterations->1000, MinFailStep -> 10^-7, MinStep -> 10^-7,
			 StepMonitorFunction -> Null, WorkingPrecision -> MachinePrecision}

PowellHybrid::minimum = "PowellHybrid found a local minimum with squared error ``"
PowellHybrid::tooManyEquations = "More equations than unknowns in PowellHybrid"
PowellHybrid::tooManyIterations = "PowellHybrid failed to complete after `` iterations"
PowellHybrid::functionFails = "PowellHybrid failed because a step of length `` led to invalid parameters"
PowellHybrid::initialFailure = "PowellHybrid function failed on initial parameters"
PowellHybrid::tooSmall = "PowellHybrid reduced stepsize to too-small value ``"
PowellHybrid::verbosityInvalid = "Verbosity option should be a number or Infinity"

MaxStep::usage = "MaxStep is an option to PowellHybrid, giving the threshold for the suggested gradient stepped to indicate that a minimum has been found"
InitialStep::usage = "InitialStep is an option to PowellHybrid, giving the length of the initial step or one of the special symbols PureNewton or PureGradient"
MinFailStep::usage = "MinFailStep is an option to PowellHybrid, giving the minimum distance to a failing value which will lead to giving up."
MinStep::usage = "MinStep is an option to PowellHybrid.  When repeated backing off produces this stepsize, PowellHybrid will give up."
PureNewton::usage = "PureNewton is a value for the InitialStep option to PowellHybrid, saying to take a Newton's method step in the first iteration"
PureGradient::usage = "PureGradient is a value for the InitialStep option to PowellHybrid, saying to take an unscaled gradient step in the first iteration"
StepMonitorFunction::usage = "StepMonitorFunction is option to PowellHybrid.  If given, it is a function to be called with a list of the parameters after a successful step."

(* The reason this is not called Verbose is that that is a symbol in the System package *)
Verbosity::usage = "Verbosity is an option to PowellHybrid giving the the desired verbosity level. 0: no messages, 1: print a dot for each step taken, 2: one line for each step taken, 3: everything but Jacobian, 4 or Infinity: everything"

Begin["`Private`"]

PowellHybrid::usage = "PowellHybrid[function, initialParameters, jacobianFunction] uses Powell's hybrid method to find a simultaneous zero of N functions depending on (at least) N parameters.  The function will be called with a list of parameters, starting with the initial guess in initialParameters.  The jacobianFunction will be called with the same parameter list and should return a Jacobian matrix.  Each column gives the dependence of the functions on one of the parameters."

PowellHybrid[function_, initialParameters_, jacobianFunction_, OptionsPattern[]] :=
  Module[{parameters = initialParameters, (* Current parameters *)
	  values = Null,		  (* Function values.  Not known yet. *)
	  verbosity = OptionValue[Verbosity],
	  stepsize = OptionValue[InitialStep], (* Step size, or PureNewton or PureGradient *)
	  maxStep = OptionValue[MaxStep],
	  accuracyGoal = OptionValue[AccuracyGoal],
	  maxIterations = OptionValue[MaxIterations],
	  minFailStep = OptionValue[MinFailStep],
	  minStep = OptionValue[MinStep],
	  stepMonitorFunction = OptionValue[StepMonitorFunction],
	  workingPrecision = OptionValue[WorkingPrecision],
	  tolerance,
	  previousscale = 1,	(* Scaling recommended on previous iteration *)
	  iterations = 0,
	  sqerror,		(* Squared error *)
	  jacobian,
	  gradient, gnorm,
	  step, nstep, gstep, nstepsize, gstepsize,
	  newparameters, newvalues, newsqerror,
	  stepped = True,	(* Flag for whether we took a step *)
	  dmult, ss, sp, lambda, mu	    (* Local variables for adjusting stepsize *)
	 },
    (* verbosity must be something that we can compare with numbers and get a definite answer *)
    If[verbosity>0, Null, Null, Message[PowellHybrid::verbosityInvalid]; verbosity = 0]
    If[accuracyGoal == Automatic, accuracyGoal = workingPrecision/2];
    If[maxStep == Automatic, maxStep = 10^(workingPrecision/2)];
    tolerance = 10^(-2 accuracyGoal); (* When squared error less than this, declare success *)
    parameters = SetPrecision[parameters, workingPrecision]; (* Reprecision starting point *)
    If[verbosity>=2, Print["Powell's hybrid method"]];
    If[verbosity>=3, Print["Starting from ", parameters]];
    values = function[parameters];
    If[Head[values] =!= List, (* Failed on initial parameters? *)
       Message[PowellHybrid::initialFailure];
       Return["Failed"]];	(* Return without doing anything *)
    If[Length[values] > Length[parameters], Message[PowellHybrid::tooManyEquations]];
    sqerror = values . values;
    If[verbosity>=3, Print["Values ", values]];
    If[verbosity>=2, Print["Initial squared error ", sqerror]];
    While[True,
      If[verbosity>=3, Print["----------------------------------------------------------------------"]];
      If[stepped,		(* If we took a step, check for completion and redo Jacobian *)
	 If[sqerror < tolerance, If[verbosity>=2, Print["Success"]];Break[]]; (* Success *)
	 If[iterations >= maxIterations, Message[PowellHybrid::tooManyIterations, iterations];Break[]];
	 jacobian = jacobianFunction[parameters]; (* Find Jacobian *)
	 If[verbosity>=4, Print["Jacobian ", jacobian]];
	 stepped = False];
      (* One trial step of Powell's hybrid method *)
      If[verbosity>=3, Print["Desired step size ", stepsize]];
      nstep = NewtonStep[values, jacobian]; (* Step recommended by Newton's method *)
      nstepsize = Norm[nstep];
      If[verbosity>=3, Print["Newton step ", nstep]];
      gradient = - 2 values . jacobian; (* Negative gradient of squared error *)
      gnorm = Norm[gradient];
      If[sqerror > maxStep gnorm, Message[PowellHybrid::minimum, sqerror];Break[]];
      gstep = gradientStep[jacobian, gradient, gnorm]; (* Step recommended by gradient method *)
      gstepsize = Norm[gstep];
      If[verbosity>=3, Print["Gradient step ", gstep]];
      (* We never take a step larger than stepsize.  If Newton recommends a shorter step, we take it.
	 If Newton recommends a longer step, but gradient recommends shorter step, then we mix them together.
	 If both recommend longer steps, we go distance stepsize in the gradient direction *)
      Which[stepsize === PureNewton, stepsize = nstepsize, (* Deal with special requests *)
	    stepsize === PureGradient, stepsize = gstepsize];
      Which[nstepsize <= stepsize, step = nstep; stepsize = nstepsize; If[verbosity>=3, Print["Using Newton"]],
	    stepsize <= gstepsize, Module[{scale = stepsize/Norm[gstep]},
					  step = scale gstep;
					  If[verbosity>=3, Print["Scaling gradient by ", scale, " gives ", step]]],
	    True, step = combineSteps[gstep, gstepsize, nstep, nstepsize, stepsize, verbosity]];
      step = SetPrecision[step, workingPrecision];
      newparameters = SetPrecision[parameters + step, workingPrecision];	(* Try step.  Don't lose precision. *)
      If[verbosity>=3, Print["Trying ", newparameters]];
      newvalues = function[newparameters];
      If[verbosity>=3,Print["Values ", newvalues]];
      If[Head[newvalues] =!= List, (* function failed *)
	 If[verbosity>=3, Print["Failed"]];
	 If[stepsize < minFailStep, Message[PowellHybrid::functionFails, stepsize];Break[],
	    stepsize = stepsize/2; previousscale = 1],
	 newsqerror = newvalues . newvalues; (* function succeeded *)
	 If[newsqerror > sqerror,
	    (* Worse: decrease stepsize by 2, don't take step *)
	    If[verbosity>=3, Print["Worse"]]; stepsize = stepsize/2; previousscale = 1;
	    If[stepsize < minStep, Message[PowellHybrid::tooSmall, stepsize];Break[]],
	    (* Determine the "quality" of the step just taken, i.e., how much did it reduce the error
	       compared to what we expected *)
	    Module[{predictedvalues = values + jacobian . step, (* First-order model of what values should be now *)
		    predictederror, quality},
	      (* Ratio of actual improvement in square error to expected improvement *) 
	      predictederror = predictedvalues . predictedvalues;
	      quality = (sqerror - newsqerror) / (sqerror - predictederror);
		   If[quality < 0, Print["Negative quality: ", sqerror//FullForm, " ", newsqerror//FullForm, " ", predictederror//FullForm]]
	      If[verbosity>=3, Print["Improvement ", quality, " of predicted"]];
	      If[quality < 0.1,
		 stepsize = stepsize/2; previousscale = 1; (* Not even 10% of what we hoped for: decrease by 2 *)
		 If[stepsize < minStep, Message[PowellHybrid::tooSmall, stepsize];Break[]],
		 dmult = sqerror - newsqerror - (sqerror - predictederror)/10;
		 sp = Apply[Plus, MapThread[#1(#1-#2)&, {newvalues, predictedvalues}]];
		 ss = Norm[newvalues - predictedvalues]^2;
		 lambda = If[ss == 0, 2, (* If prediction exact, double and avoid dividing by zero *)
			     Sqrt[1 + dmult/Sqrt[sp^2 + dmult ss]]];
		 (* Use only the lesser of the increases recommended in the previous 2 iterations *)
		 mu = Min[2, lambda, previousscale]; (* Never more than double *)
		 stepsize = Min[mu stepsize, maxStep];
		 previousscale = lambda / mu (* Amount of increase left for next time *)
	      ]];
	    parameters = newparameters; values = newvalues; sqerror = newsqerror; (* Take step *)
	    iterations++;
	    If[verbosity>=2, Print["Took step #", iterations, ".  Squared error now ", sqerror],
	       If[verbosity>=1, WriteString[$Output,"."]]];
	    If[stepMonitorFunction =!= Null,
	       stepMonitorFunction[parameters]];
	    stepped = True;
	    ]]];
    parameters]
       
SetAttributes[PowellHybrid, HoldAll]

(* Newton's method step.  Solve linearized equation J.x = -v. *)
NewtonStep[values_, jacobian_] :=
  Block[{nValues = Length[values], nParameters = Dimensions[jacobian][[2]]},
	If[nValues == nParameters, (* Simple case *)
	   LinearSolve[jacobian, -values],
	   Block[{u, w, v, wInverse},
	      (* SVD jacobian.  Number of nonsingular entries is (at most) the number of rows.
		 It does not seem to be necessary to put zero rows in the matrix first. *)
	      {u, w, v} = SingularValueDecomposition[jacobian, nValues];
	      wInverse = Table[If[i==j, 1/w[[i,j]], 0], {i, nValues}, {j, nValues}];
	      -v . wInverse . Transpose[u] . values
	   ]]]

(* Find the step in the gradient direction that will minimize the error, under the approximation that
   values are linear in params *)
gradientStep[jacobian_, gradient_, gnorm_] := gradient (gnorm / Norm[jacobian . gradient])^2/2;

(* Find a linear combination of the vectors gstep and nstep that has the length stepsize.  We only call
   this code when gstepsize < stepsize < nstepsize, so this is always possible.  We return
   t*nstep + (1-t)*gstep where t solves the quadratic equation
   t^2 |nstep|^2 + (1-t)^2 |gstep|^2 + 2t(1-t) gstep . nstep = stepsize^2 *)

combineSteps[gstep_, gstepsize_, nstep_, nstepsize_, stepsize_, verbosity_] :=
  Module[{a = nstepsize^2 + gstepsize^2 - 2 gstep . nstep,
   	  b = gstep . nstep - gstepsize^2,
   	  c = gstepsize^2 - stepsize^2},
   Module[{t = (Sqrt[b^2 - 4 a c] - b)/(2 a)},
     Module[{step =      t nstep + (1-t) gstep},
       If[verbosity>=3, Print["Combination of ", t, " Newton + ", 1-t, " gradient gives ", step];
	 Print["Newton norm is ", Norm[t nstep]/(Norm[(1-t) gstep]+Norm[t nstep]), " of sum of norms"]];
       step]]]

(* Simple test:
   PowellHybrid[{#[[1]]^2-2}&, {1.0}, {{2#[[1]]}}&, Verbosity->True,MaxIterations-> 4] *)


End[]

EndPackage[]
