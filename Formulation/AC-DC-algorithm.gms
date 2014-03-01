***** License ******************************************************************
* Optimal Design of Mixed AC-DC Distribution Systems For Commercial Buildings: *
* A Nonconvex Generalized Benders Decomposition Approach                       *
*                                                                              *
* Copyright (C) 2014  Stephen M. Frank (stephen.frank@ieee.org)                *
*                                                                              *
* This program is free software: you can redistribute it and/or modify         *
* it under the terms of the GNU General Public License as published by         *
* the Free Software Foundation, either version 3 of the License, or            *
* (at your option) any later version.                                          *
*                                                                              *
* This program is distributed in the hope that it will be useful,              *
* but WITHOUT ANY WARRANTY; without even the implied warranty of               *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
* GNU General Public License for more details.                                 *
*                                                                              *
* You should have received a copy of the GNU General Public License            *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.        *
********************************************************************************

$TITLE Solution algorithm for AC-DC Model
$ONTEXT

    Filename : AC-DC-algorithm.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1
    
    This file implements a solution algorithm for the AC-DC model based on
    nonconvex generalized Benders' decomposition (NGBD). This algorithm is
    described in detail in the article "Optimal Design of Mixed AC-DC
    Distribution Systems For Commercial Buildings: A Nonconvex Generalized
    Benders Decomposition Approach" by Stephen Frank and Steffen Rebennack.
    
    This script assumes that the model and data are already loaded. It is
    designed to be included (via $INCLUDE) in other .gms files, and relies on
    the following compile-time parameters to be set in the parent file (which
    are documented in AC-DC-defaults.gms):
        LBP_rel                 Relative tolerance for lower bounding problem
        LBP_abs                 Absolute tolerance for lower bounding problem
        UBP_rel                 Relative tolerance for upper bounding problem
        UBP_abs                 Absolute tolerance for upper bounding problem
        LBP_reslim              Resource limit for lower bounding problem
        UBP_reslim              Resource limit for upper bounding subproblems
        TimeLimit               Time limit for overall algorithm
        ConvexRelaxation        Specify type of convex relaxation
                                (determines the class of UBP solver)
        ExportEveryIteration    Whether to write a GDX after each NGBD iteration
        FindAllFeasibleDesigns  Whether to search for all feasible designs
        
    This script also requires the compile-time parameter %OutFile%, which
    specifies the root name of the output GDX files.

    NOTES:
    1. To attempt to find all feasible system designs (that is, all feasible
       integer solutions), set %FindAllFeasibleDesigns%==yes. This setting
       disables the convergence check; the algorithm will therefore run until
       the LBP becomes infeasible. Caution: there may be a very large number of
       feasible designs for large problems!
	2. If %ExportEveryIteration%==yes, then the algorithm will write a GDX file
	   after each NGBD iteration, numbered sequentially by iteration. This can
	   be useful for troubleshooting, but use caution: this option can generate
	   a very large number of files which use a lot of disk space.
    3. A number of conditions generate run-time warnings in the algorithm.
       Since there is not (to my knowledge) a method to trigger a run-time
       warning in GAMS itself, these warning messages are instead written to the
       .lst output with the prefix "WARNING:". Therefore, check the end of the 
       .lst file after executing the algorithm to determine if any run-time
       warnings were generated.
    
$OFFTEXT

* =============================================================================
*  Setup
* =============================================================================
* ------------------------------------------------------------------------------
*  Bookkeeping
* ------------------------------------------------------------------------------
* Storage of generated solutions
Parameters
    xB_LBP(r,a,i,k) 'xB LBP solution at iteration r'
    xS_LBP(r,a,s)   'xS LBP solution at iteration r'
    xL_LBP(r,a,l)   'xL LBP solution at iteration r'
    xC_LBP(r,c)     'xC LBP solution at iteration r' ;

* Storage of objective function values
Parameters
    zLBP_LB(r)  'Lower bound on LBP objective function value for iteration r'
    zLBP_UB(r)  'Best found LBP objective function value for iteration r'
    zUBP_LB(r)  'Lower bound on UBP objective function value for iteration r'
    zUBP_UB(r)  'Best found UBP objective function value for iteration r'
    zLB         'Tightest lower bound on LBP objective function (up to current iteration)'
    zBest       'Best feasible objective function value (tightest upper bound)' ;

* Loop control and bookkeeping
Parameters
    cvg                 'Convergence flag (yes/no)'
    isFeas(r)           'UBP feasibility flag for iteration r (yes/no)'
    rBest               'Ordinal of iteration r with best feasible solution'
    algStTime           'Algorithm starting time (for measuring time limits)'
    iterStTime          'Iteration starting time (for measuring iteration times)'
    algTiming(r,*,*)    'Report timing for LBP, UBP, etc. at each iteration r' ;

* Storage of binary variable bounds
* (These must be restored for the LBP after fixing them for the UBP)
Parameters
    xB_box(a,i,k,*) 'Storage of upper and lower bounds on xB'
    xS_box(a,s,*)   'Storage of upper and lower bounds on xS'
    xL_box(a,l,*)   'Storage of upper and lower bounds on xL'
    xC_box(c,*)     'Storage of upper and lower bounds on xC' ;

xB_box(BX(a,i,k),'lo') =    xB.lo(a,i,k);
xB_box(BX(a,i,k),'up') =    xB.up(a,i,k);

xS_box(SS(a,s),'lo') =      xS.lo(a,s);
xS_box(SS(a,s),'up') =      xS.up(a,s);

xL_box(LL(a,l),'lo') =      xL.lo(a,l);
xL_box(LL(a,l),'up') =      xL.up(a,l);

xC_box(c,'lo') =            xC.lo(c);
xC_box(c,'up') =            xC.up(c);


* ------------------------------------------------------------------------------
*  Solution elimination constraints
* ------------------------------------------------------------------------------
* Declaration
equation eq_SolElim(r)  'Solution elimination constraints for LBP' ;

* Definition
eq_SolElim(RR(r))..
    sum(BX(a,i,k)$(xB_LBP(r,a,i,k) = 0), xB(a,i,k) - xB_LBP(r,a,i,k) )
    + sum(BX(a,i,k)$(xB_LBP(r,a,i,k) = 1), xB_LBP(r,a,i,k) - xB(a,i,k) )
    + sum(SS(a,s)$(xS_LBP(r,a,s) = 0), xS(a,s) - xS_LBP(r,a,s) )
    + sum(SS(a,s)$(xS_LBP(r,a,s) = 1), xS_LBP(r,a,s) - xS(a,s) )
    + sum(LL(a,l)$(xL_LBP(r,a,l) = 0), xL(a,l) - xL_LBP(r,a,l) )
    + sum(LL(a,l)$(xL_LBP(r,a,l) = 1), xL_LBP(r,a,l) - xL(a,l) )
    + sum(c$(xC_LBP(r,c) = 0), xC(c) - xC_LBP(r,c) )
    + sum(c$(xC_LBP(r,c) = 1), xC_LBP(r,c) - xC(c) )
    =g= 1 ;


* ------------------------------------------------------------------------------
*  Lower bound constraint
* ------------------------------------------------------------------------------
* Constrains objective function value of LBP to be greater than or equal to the
* previous lower bound

* Declaration
equation eq_LBPLB   'Previous lower bound constraint for LBP' ;

* Definition
eq_LBPLB..
    z =g= zLB ;


* =============================================================================
*  Models
* =============================================================================
* Uses and augments the models from AC-DC-models.gms

** Models
* Lower Bounding Problem
*   -- Includes solution elimination constraint
*   -- Excludes obj. function value cut
model LBPalg        / LBP + eq_SolElim / ;

* Upper Bounding Problem
model UBPalg        / UBPr / ;

** Model attributes
* Convert fixed variables to constants (pre-process them out)
LBPalg.holdfixed = 1;
UBPalg.holdfixed = 1;

* Set relative optimality tolerances
LBPalg.optCR = %LBP_rel%;
UBPalg.optCR = %UBP_rel%;

* Set absolute optimality tolerances
LBPalg.optCA = %LBP_abs%;
UBPalg.optCA = %UBP_abs%;

* Set resource limits, if specified
$IF SET LBP_reslim          LBPalg.reslim = %LBP_reslim%;
$IF SET UBP_reslim          UBPalg.reslim = %UBP_reslim%;

* Allow branching priorities in LBP
LBPalg.prioropt = 5;


* =============================================================================
*  Algorithm
* =============================================================================
** Initialization
* Initialize iteration count to zero
RR(r) = no;
isFeas(r) = no;

* Set up for exporting GDX on each iteration
$IFTHEN %ExportEveryIteration% == yes
file theoutput;
put theoutput;
$ENDIF

* Initialize stored solutions to zero over an empty set
* (because without initialization the algorithm breaks)
xB_LBP(RR(r),a,i,k) = 0;
xS_LBP(RR(r),a,s)   = 0;
xL_LBP(RR(r),a,l)   = 0;
xC_LBP(RR(r),c)     = 0;

* Initialize convergence flag to 'no'
cvg = no;

* Initialize lower and upper bounds
zLB = -Inf;
zBest = +Inf;
rBest = 0;

** Loop
* Record starting time
algStTime = timeElapsed;

* Loops until convergence or max. iterations exhausted
loop(r2$(not cvg),
    // Abort if time limit exceeded
    if (timeElapsed - algStTime > %TimeLimit%,
        execute_unload '%OutFile%_on_abort';
        abort 'Algorithm total time limit exceeded';
    );

    // Record iteration start time
    iterStTime = timeElapsed;

    // ---Lower Bounding Problem / Setup---
    // Unfix x's
    xB.lo(BX(a,i,k)) =  xB_box(a,i,k,'lo');
    xB.up(BX(a,i,k)) =  xB_box(a,i,k,'up');

    xS.lo(SS(a,s)) =    xS_box(a,s,'lo');
    xS.up(SS(a,s)) =    xS_box(a,s,'up');

    xL.lo(LL(a,l)) =    xL_box(a,l,'lo');
    xL.up(LL(a,l)) =    xL_box(a,l,'up');

    xC.lo(c) =          xC_box(c,'lo');
    xC.up(c) =          xC_box(c,'up');

    // Begin by assuming that this iteration will yield a feasible solution
    isFeas(r2) = yes;

    // Turn on all time periods
    TT(t) = yes;

    // Load previous LBP solution (on iterations 2+)
    if( ord(r2) > 1,
        execute_loadpoint '%OutFile%_after_last_LBP',
            V, E, F,
            PB, QB, PS, QS, PL, QL, PCin, PCout, QCin, QCout,
            IReB, IImB, IReS, IImS, IReL, IImL,
            IReCout, IImCout, IReCin, IImCin,
            wVV, wEE, wEF, wFF, wC, wCC,
            xB, xS, xL, xC,
            z ;
    );

    // ---Lower Bounding Problem / Solve---
    // Solve LBP
$IFTHEN %ConvexRelaxation% == linear
    solve LBPalg using MIP minimizing z;      // Mixed Integer Linear
$ELSEIF %ConvexRelaxation% == quadratic
    solve LBPalg using MIQCP minimizing z;    // Convex Quadratically Constrained
$ELSEIF %ConvexRelaxation% == piecewise
    solve LBPalg using MIP minimizing z;      // Mixed Integer Linear
$ENDIF

    // Record time it took
    algTiming(r2, 'LBP', 'Total') = LBPalg.resUsd;

    // Check for abnormal termination
    if( LBPalg.solvestat > %SolveStat.ResourceInterrupt%,
        // Abnormal solver status; abort
        // NOTE: Reaching resource or iteration limits is acceptable
        execute_unload '%OutFile%_on_abort';
        abort 'LBP solver encountered a problem:',
            LBPalg.solvestat, LBPalg.modelstat;

        // Check/warn for exceeding resource limits
        if( LBPalg.solvestat > %SolveStat.NormalCompletion%,
            display 'WARNING: LBP solver terminated early; solvestat = ',
                LBPalg.solvestat;
        );

    // Check for feasible solution
    elseif LBPalg.modelstat = %ModelStat.Optimal% or
           LBPalg.modelstat = %ModelStat.LocallyOptimal% or
           LBPalg.modelstat = %ModelStat.IntermediateNonoptimal% or
           LBPalg.modelstat = %ModelStat.IntegerSolution%,
        // Store objective function values
        zLB = max(zLB, LBPalg.objEst);
        zLBP_LB(r2) = zLB;
        zLBP_UB(r2) = z.l;

    // Check for infeasible solution
    elseif LBPalg.modelstat = %ModelStat.Infeasible% or
           LBPalg.modelstat = %ModelStat.LocallyInfeasible% or
           LBPalg.modelstat = %ModelStat.IntermediateInfeasible% or
           LBPalg.modelstat = %ModelStat.IntegerInfeasible% or
           LBPalg.modelstat = %ModelStat.Infeasible-NoSolution%,

        // LBP infeasible -> meets convergence criteria
        zLB = max(zLB, LBPalg.objEst);
        zLBP_LB(r2) = zLB;
        zLBP_UB(r2) = +Inf;
        isFeas(r2) = no;
        cvg = yes;

        // If intermediate infeasible, display a warning message
        if( LBPalg.modelstat = %ModelStat.IntermediateInfeasible%,
            display 'WARNING: LBP solution is intermediate infeasible; algorithm cannot guarantee proper termination.';
        );

    else
        // Something is quite wrong! Terminate algorithm.
        execute_unload '%OutFile%_on_abort';
        abort 'Unexpected outcome when solving LBP:',
            LBPalg.solvestat, LBPalg.modelstat;
    );

    // LBP feasible -> store LBP solution
    xB_LBP(r2,a,i,k) =  xB.l(a,i,k);
    xS_LBP(r2,a,s) =    xS.l(a,s);
    xL_LBP(r2,a,l) =    xL.l(a,l);
    xC_LBP(r2,c) =      xC.l(c);

    // Save progress this far
    execute_unload '%OutFile%_after_last_LBP';

    // ---Upper Bounding Problem / Setup---
    // Fix x's
    // (Round to nearest integer to avoid numerical errors)
    xB.fx(BX(a,i,k)) =  round( xB.l(a,i,k) );
    xS.fx(SS(a,s)) =    round( xS.l(a,s) );
    xL.fx(LL(a,l)) =    round( xL.l(a,l) );
    xC.fx(c) =          round( xC.l(c) );

    // Initialize objective function values
    if( isFeas(r2),
         // LBP is feasible -> UBP may be feasible
         zUBP_LB(r2) = 0;
         zUBP_UB(r2) = 0;
    else
         // LBP is infeasible -> UBP is also infeasible
         zUBP_LB(r2) = +Inf;
         zUBP_UB(r2) = +Inf;
    );

    // ---Upper Bounding Problem / Solve---
    // Loop over each time period
    TT(t) = no;
    loop(t2$isFeas(r2),
        // Turn time period 't' on
        TT(t2) = yes;

        // Solve UBP for this time period
        solve UBPalg using MIQCP minimizing z;

        // Record time it took
        algTiming(r2, 'UBP', t2) = UBPalg.resUsd;

        // Check for abnormal termination
        if( UBPalg.solvestat > %SolveStat.ResourceInterrupt%,
            // Abnormal solver status; abort
            // NOTE: Reaching resource or iteration limits is acceptable
            execute_unload '%OutFile%_on_abort';
            abort 'UBP solver encountered a problem:',
                UBPalg.solvestat, UBPalg.modelstat;

            // Check/warn for exceeding resource limits
            if( UBPalg.solvestat > %SolveStat.NormalCompletion%,
                display 'WARNING: UBP solver terminated early; solvestat = ',
                    UBPalg.solvestat;
            );

        // Check for feasible solution
        elseif UBPalg.modelstat = %ModelStat.Optimal% or
               UBPalg.modelstat = %ModelStat.LocallyOptimal% or
               UBPalg.modelstat = %ModelStat.IntermediateNonoptimal% or
               UBPalg.modelstat = %ModelStat.IntegerSolution% or
               UBPalg.modelstat = %ModelStat.SolvedUnique% or
               UBPalg.modelstat = %ModelStat.Solved%,

            // UBP is feasible -> record obj. function value and best bound
            zUBP_LB(r2) = zUBP_LB(r2) + UBPalg.objEst;
            zUBP_UB(r2) = zUBP_UB(r2) + z.l;

        // Check for infeasible solution + resource interrupt
        elseif UBPalg.solvestat = %SolveStat.ResourceInterrupt% and
               (
               UBPalg.modelstat = %ModelStat.Infeasible% or
               UBPalg.modelstat = %ModelStat.LocallyInfeasible% or
               UBPalg.modelstat = %ModelStat.IntermediateInfeasible% or
               UBPalg.modelstat = %ModelStat.IntegerInfeasible% or
               UBPalg.modelstat = %ModelStat.Infeasible-NoSolution% or
               UBPalg.modelstat = %ModelStat.NoSolutionReturned%
               ),

            // UBP infeasible -> set UB to +Inf but keep going
            zUBP_LB(r2) = zUBP_LB(r2) + UBPalg.objEst;
            zUBP_UB(r2) = +Inf;

        // Check for infeasible solution
        elseif UBPalg.modelstat = %ModelStat.Infeasible% or
               UBPalg.modelstat = %ModelStat.LocallyInfeasible% or
               UBPalg.modelstat = %ModelStat.IntermediateInfeasible% or
               UBPalg.modelstat = %ModelStat.IntegerInfeasible% or
               UBPalg.modelstat = %ModelStat.Infeasible-NoSolution% or
               UBPalg.modelstat = %ModelStat.NoSolutionReturned%,

            // UBP infeasible -> meets convergence criteria
            zUBP_LB(r2) = +Inf;
            zUBP_UB(r2) = +Inf;
            isFeas(r2) = no;

            // If infeasibility is not guaranteed, display a warning message
            if( LBPalg.modelstat ne %ModelStat.IntermediateInfeasible%,
                display 'WARNING: No UBP solution found, but algorithm cannot guarantee infeasibility: solvestat = ',
                    UBPalg.solvestat;
            );

        else
            // Something is quite wrong! Terminate algorithm.
            execute_unload '%OutFile%_on_abort';
            abort 'Unexpected outcome when solving UBP:',
                UBPalg.solvestat, UBPalg.modelstat;
        );

        // Turn time period 't' back off
        TT(t2) = no;
    );

    // Tally total UBP time
    algTiming(r2, 'UBP', 'Total') = sum(t, algTiming(r2, 'UBP', t));

    // ---Convergence and Control---
    // Check for new incumbent
    if ( zUBP_UB(r2) < zBest,
        // Update incumbent
        zBest = zUBP_UB(r2);
        rBest = ord(r2);

        // Store new incumbent solution to file for later
        // (Should contain all harmonic data from last LBP solve as well)
        execute_unload '%OutFile%_best_UBP';
    );

    // Double-check for infeasibility
    // (Needed in case of unconfirmed infeasbility from UBP subproblems)
    if ( zUBP_UB(r2) = +Inf,
        isFeas(r2) = no;
    );


    // Activate new solution elimination constraint by updating set RR
    RR(r2) = yes;

    // Check for convergence
$IFTHEN NOT %FindAllFeasibleDesigns% == yes
    if( zLB >= zBest * (1 - %UBP_rel%) - %UBP_abs%, cvg = yes; );
$ENDIF

    // Export if in the 'export every iteration' verbose mode
$IFTHEN %ExportEveryIteration% == yes
    put_utility 'gdxout' / '%OutFile%_NGBD_iteration_' r2.tl:0;
    execute_unload;
$ENDIF

    // Warn if reached maximum iterations without convergence
    if (ord(r2) ge %MaxIter%,
        display 'WARNING: Maximum number of iterations (%MaxIter%) reached prior to algorithm convergence.';
    );

);

** Extract optimum
* Load the best incumbent (only relevant info)
*   This sets all the variable levels and marginals to the correct values at
*   the optimum. It no longer reads equations (for now).
execute_loadpoint '%OutFile%_best_UBP',
    // Variables
    V, E, F,
    PB, QB, PS, QS, PL, QL, PCin, PCout, QCin, QCout,
    IReB, IImB, IReS, IImS, IReL, IImL, IReCout, IImCout, IReCin, IImCin,
    wVV, wEE, wEF, wFF, wC, wCC,
    xB, xS, xL, xC,
    z ;
