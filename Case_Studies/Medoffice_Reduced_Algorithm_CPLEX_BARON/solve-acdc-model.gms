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

$TITLE Solve script for the mixed AC-DC distribution system design problem
$ONTEXT

    Filename : solve-acdc-model.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    This script executes a full solution to the mixed AC-DC distribution system
    design problem starting with a set of input data. The default steps are:
        1. Load data and initialize variables
        2. Compute big M values
        3. Execute bound tightening
        4. Execute NGBD solution algorithm
        5. Export solution to GDX
    However, by changing the %ExecMode% compile time flag, the user can
    alter the default execution mode to one of the following:
        %ExecMode% == OnlyBigM      Runs the big M computation procedure using
                                    the data in %DataFile%, then stops.
        %ExecMode% == OnlyBT        Runs the bound tightening procedure using
                                    the data in %DataFile%, then stops
        %ExecMode% == OnlyLBP       Performs a 1-shot solve of (LBP) using
                                    the data in %DataFile%
        %ExecMode% == OnlyUBP       Performs a 1-shot solve of (UBP) using
                                    the data in %DataFile%
        %ExecMode% == OnlyMonolith  Performs a 1-shot solve of the monolith (P)
                                    using the data in %DataFile%
        %ExecMode% == OnlyAlgorithm Executes the solution algorithm directly
                                    using the data in %DataFile%
        %ExecMode% == SkipBT        As the default mode, but skips bound
                                    tightening

    The output files differ by what problem is solved:
        %ExecMode%          GDX output file(s)
        Default             %OutFile%_presolve, %OutFile%_big_M,
                            %OutFile%_tightened, %OutFile%_solved
        OnlyAlgorithm       %OutFile%_presolve, %OutFile%_solved
        SkipBT              %OutFile%_presolve, %OutFile%_big_M,
                            %OutFile%_solved
        (Anything else)     %OutFile%
    IE suffixes are only appended when there is more than one output file.

    The script is designed to be called from the command line, either via
    another GAMS script or from MATLAB. Therefore, the following must be
    specified when the script is called:
        DataFile            Name of GDX file from which to import initial data

    In addition to the compile-time flags documented in 'AC-DC-defaults.gms',
    the script supports the following optional compile time settings:
        OutFile             Root name of output GDX file, to which various
                            suffixes will be appended. (Default = %DataFile%)
        FromMatlab          Set to 'yes' if calling from MATLAB via gdxmrw.
                            (Default = no)
        UseLocalSolverOpts  Set to 'yes' if the desired solver options are
                            specified in this file; 'no' if you want to pass
                            them at the command line. (Default = no)
        ExternalScript      If this is set, the contents of the file it names
                            are included just before running the algorithm.
        ExecMode            (Documented above)

    If %UseLocalSolverOpts% == no, then you will want to specify default MIP
    and MIQCP solvers at the command line to ensure you get what you want.
$OFFTEXT

** Enable C and C++ style comments
* End-of-line
$EOLCOM //
$ONEOLCOM

* In-line
$INLINECOM /* */
$ONINLINE

* =============================================================================
*  MATLAB Settings
* =============================================================================
* This section sets up the script for access from MATLAB. If calling from
* from MATLAB, set the command line parameter --FromMatlab=yes to enable
* this to work correctly.
$IF NOT SET FromMatlab          $SET FromMatlab no
$IFTHENi.ml %FromMatlab% == yes
$SET matout "'ToMatlab.GDX'"
$ENDIF.ml


* =============================================================================
*  User Settings
* =============================================================================
** Execution Mode
* Leave as 'default' to execute the full solution algorithm
$IF NOT SET ExecMode            $SET ExecMode default

** Compile-time Settings
* These may also be passed at the command line in lieu of setting them here.
* Anything not set here uses the values in 'AC-DC-defaults.gms'

* Relative path to AC-DC files
$IF NOT SET RelPath             $SET RelPath ../Formulation

* What convex relaxation for the LBP
*   Choose one of: linear, quadratic, piecewise
*   If 'piecewise', also choose piecewise type: mccormick, logtransform
$IF NOT SET ConvexRelaxation    $SET ConvexRelaxation linear
$IF NOT SET PiecewiseType       $SET PiecewiseType logtransform


** Solver Options
* Use these locally specified options?
$IF NOT SET UseLocalSolverOpts  $SET UseLocalSolverOpts no
$IFTHENi.solveopts %UseLocalSolverOpts% == yes

* Select solvers:
option MIP =    CPLEX;      // LBP w/ linear or piecewise linear relaxation
option RMIP =   CPLEX;      // Bound tightening w/ linear relaxation
option MIQCP =  GLOMIQO;    // LBP w/ quadratic relaxation, UBP, Monolith
option RMIQCP = GLOMIQO;    // Bound tightening w/ quadratic relaxation

* Set GAMS solver options to increase speed
option limrow = 0;
option limcol = 0;
option solprint = silent;
option solvelink = %solvelink.LoadLibrary%;

* Default resource limit per solve command (sec)
option reslim = 600;

$ENDIF.solveopts


* =============================================================================
*  Setup
* =============================================================================
** Timings
* Manually record the start and end time (s) for various activities
Parameter Timing(*,*) 'Algorithm execution times';

** File I/O
* GDX file to import from (pass on command line)
$IF NOT SET DataFile            $ABORT 'DataFile' must be specified.

* GDX file root name to write
$IF NOT SET OutFile             $SET OutFile %DataFile%

** Required Flags
* Defaults
$SET LoadSets       yes
$SET LoadParams     yes
$SET LoadVars       no
$SET LoadBigM       no
$SET InitVars       yes

* Overrides for special execution modes
$IFi %ExecMode% == OnlyBT           $SET LoadBigM       yes
$IFi %ExecMode% == OnlyLBP          $SET LoadBigM       yes
$IFi %ExecMode% == OnlyUBP          $SET LoadBigM       yes
$IFi %ExecMode% == OnlyMonolith     $SET LoadBigM       yes
$IFi %ExecMode% == OnlyAlgorithm    $SET LoadBigM       yes

$IFi %ExecMode% == OnlyLBP          $SET LoadVars       yes
$IFi %ExecMode% == OnlyUBP          $SET LoadVars       yes
$IFi %ExecMode% == OnlyMonolith     $SET LoadVars       yes
$IFi %ExecMode% == OnlyAlgorithm    $SET LoadVars       yes

$IFi %ExecMode% == OnlyLBP          $SET InitVars       no
$IFi %ExecMode% == OnlyUBP          $SET InitVars       no
$IFi %ExecMode% == OnlyMonolith     $SET InitVars       no
$IFi %ExecMode% == OnlyAlgorithm    $SET InitVars       no

** Default LBP optimality tolerances
* (For "OnlyLBP" execution mode)
$IFTHENi %ExecMode% == OnlyLBP

* Relative
$IF NOT SET LBP_rel                 $SET LBP_rel        1e-3

* Absolute
$IFTHENi.tols %ConvexRelaxation% == linear
$IF NOT SET LBP_abs                 $SET LBP_abs        0

$ELSEIFi.tols %ConvexRelaxation% == piecewise
$IF NOT SET LBP_abs                 $SET LBP_abs        0

$ELSEIFi.tols %ConvexRelaxation% == quadratic
$IF NOT SET LBP_abs                 $SET LBP_abs        1e-6

$ENDIF.tols

$ENDIF

** Default UBP optimality tolerances
* (For "OnlyUBP" and "OnlyMonolith" execution modes)
$IFTHENi %ExecMode% == OnlyUBP

* Relative
$IF NOT SET UBP_rel                 $SET UBP_rel        1e-3

* Absolute
$IF NOT SET UBP_abs                 $SET UBP_abs        1e-6

$ENDIF


* =============================================================================
*  Data Import and Model Initialization
* =============================================================================
* Record start time
Timing('initialization','start') = timeElapsed;

* Default compile time settings
$INCLUDE %RelPath%/AC-DC-defaults.gms

* Set, parameter, and variable Definitions
$INCLUDE %RelPath%/AC-DC-defs.gms

* Import & process sets & data
$INCLUDE %RelPath%/AC-DC-data.gms

* Equation definitions
$INCLUDE %RelPath%/AC-DC-equations-reform.gms
$INCLUDE %RelPath%/AC-DC-equations-relax.gms

* Models
$INCLUDE %RelPath%/AC-DC-models.gms

* Record end time and duration
Timing('initialization','end') = timeElapsed;
Timing('initialization','duration') =
    Timing('initialization','end') - Timing('initialization','start');

* Run external script, if set
$IF SET ExternalScript
$INCLUDE %ExternalScript%


* =============================================================================
*  Execution Modes
* =============================================================================
* Special modes come first. If none of them execute, then we use the default
* algorithm.

* -----------------------------------------------------------------------------
*  Big M's
* -----------------------------------------------------------------------------
$IFTHENi.execmode %ExecMode% == OnlyBigM
* Record start time
Timing('bigM','start') = timeElapsed;

* Run big M computation
$INCLUDE %RelPath%/AC-DC-big-M.gms

* Save results to GDX
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile% bigM . sets data variables

* Record end time and duration
Timing('bigM','end') = timeElapsed;
Timing('bigM','duration') = Timing('bigM','end') - Timing('bigM','start');


* -----------------------------------------------------------------------------
*  Bound Tightening
* -----------------------------------------------------------------------------
$ELSEIFi.execmode %ExecMode% == OnlyBT
* Record start time
Timing('boundTightening','start') = timeElapsed;

* Run bound tightening
$INCLUDE %RelPath%/AC-DC-bound-tightening.gms

* Save results to GDX
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile% bound_tightening . sets data variables

* Record end time and duration
Timing('boundTightening','end') = timeElapsed;
Timing('boundTightening','duration') =
    Timing('boundTightening','end') - Timing('boundTightening','start');


* -----------------------------------------------------------------------------
*  LBP
* -----------------------------------------------------------------------------
$ELSEIFi.execmode %ExecMode% == OnlyLBP
* Record start time
Timing('LBP','start') = timeElapsed;

** Setup for solve
* Map fixed variables -> constants
LBP.holdfixed = 1 ;

* Optimality tolerances
LBP.optCR = %LBP_rel%;
LBP.optCA = %LBP_abs%;

* Branching priorities
LBP.prioropt = 1;

* Set resource limit, if specified
$IF SET LBP_reslim          LBP.reslim = %LBP_reslim%;


** Solve command - based on type of relaxation
* Mixed Integer Linear
$IFTHENi.go %ConvexRelaxation% == linear
solve LBP using MIP minimizing z;

$ELSEIFi.go %ConvexRelaxation% == piecewise
solve LBP using MIP minimizing z;

* Convex Quadratically Constrained
$ELSEIFi.go %ConvexRelaxation% == quadratic
solve LBP using MIQCP minimizing z;

$ENDIF.go


** Export GDX
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile% . LBP sets data variables equations

* Record end time and duration
Timing('LBP','end') = timeElapsed;
Timing('LBP','duration') = Timing('LBP','end') - Timing('LBP','start');


* -----------------------------------------------------------------------------
*  UBP
* -----------------------------------------------------------------------------
$ELSEIFi.execmode %ExecMode% == OnlyUBP
* Record start time
Timing('UBP','start') = timeElapsed;

** Setup for solve
* Map fixed variables -> constants
UBP.holdfixed = 1 ;

* Fix the 'x' variables at their current levels
xB.fx(BX(a,i,k)) =  round( xB.l(a,i,k) );
xS.fx(SS(a,s)) =    round( xS.l(a,s) );
xL.fx(LL(a,l)) =    round( xL.l(a,l) );
xC.fx(c) =          round( xC.l(c) );

* Optimality tolerances
UBP.optCR = %UBP_rel%;
UBP.optCA = %UBP_abs%;

* Set resource limit, if specified
$IF SET UBP_reslim          UBP.reslim = %UBP_reslim%;


** Solve command
* Since holdfixed = 1, should be able to just use a QCP or NLP solver
* However, it must be specified MIQCP or MINLP b/c the x's are present,
* even though fixed...
solve UBP using MIQCP minimizing z;


** Export GDX
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile% . UBP sets data variables equations

* Record end time and duration
Timing('UBP','end') = timeElapsed;
Timing('UBP','duration') = Timing('UBP','end') - Timing('UBP','start');


* -----------------------------------------------------------------------------
*  Monolith
* -----------------------------------------------------------------------------
$ELSEIFi.execmode %ExecMode% == OnlyMonolith
* Record start time
Timing('Monolith','start') = timeElapsed;

** Setup for solve
* Map fixed variables -> constants
Monolith.holdfixed = 1 ;

* Optimality tolerances
Monolith.optCR = %UBP_rel%;
Monolith.optCA = %UBP_abs%;

* Set resource limit, if specified
$IF SET UBP_reslim          Monolith.reslim = %UBP_reslim%;


** Solve command
* Use MIQCP
solve Monolith using MIQCP minimizing z;


** Export GDX
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile% . Monolith sets data variables

* Record end time and duration
Timing('Monolith','end') = timeElapsed;
Timing('Monolith','duration') =
    Timing('Monolith','end') - Timing('Monolith','start');


* -----------------------------------------------------------------------------
*  Algorithm (skip bound tightening and big M's)
* -----------------------------------------------------------------------------
$ELSEIFi.execmode %ExecMode% == OnlyAlgorithm

** Presolve save
* Save a pre-solution version of everthing
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_presolve . All sets data variables


** Algorithm
* Record start time
Timing('algorithm','start') = timeElapsed;

* Include the algorithm
$INCLUDE %RelPath%/AC-DC-algorithm.gms

* Export GDX of results
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_solved main . sets data variables

* Record end time and duration
Timing('algorithm','end') = timeElapsed;
Timing('algorithm','duration') =
    Timing('algorithm','end') - Timing('algorithm','start');


* -----------------------------------------------------------------------------
*  Algorithm (skip bound tightening)
* -----------------------------------------------------------------------------
$ELSEIFi.execmode %ExecMode% == SkipBT

** Presolve save
* Save a pre-solution version of everthing
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_presolve . All sets data variables


** Big M's
* Record start time
Timing('bigM','start') = timeElapsed;

* Run big M computation
$INCLUDE %RelPath%/AC-DC-big-M.gms

* Save intermediate results to GDX
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_big_M bigM . sets data variables

* Record end time and duration
Timing('bigM','end') = timeElapsed;
Timing('bigM','duration') = Timing('bigM','end') - Timing('bigM','start');


** Algorithm
* Record start time
Timing('algorithm','start') = timeElapsed;

* Include the algorithm
$INCLUDE %RelPath%/AC-DC-algorithm.gms

* Export GDX of results
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_solved main . sets data variables

* Record end time and duration
Timing('algorithm','end') = timeElapsed;
Timing('algorithm','duration') =
    Timing('algorithm','end') - Timing('algorithm','start');


* -----------------------------------------------------------------------------
*  Algorithm (default execution mode)
* -----------------------------------------------------------------------------
$ELSE.execmode

** Presolve save
* Save a pre-solution version of everthing
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_presolve . All sets data variables


** Big M's
* Record start time
Timing('bigM','start') = timeElapsed;

* Run big M computation
$INCLUDE %RelPath%/AC-DC-big-M.gms

* Save intermediate results to GDX
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_big_M bigM . sets data variables

* Record end time and duration
Timing('bigM','end') = timeElapsed;
Timing('bigM','duration') = Timing('bigM','end') - Timing('bigM','start');


** Bound Tightening
* Record start time
Timing('boundTightening','start') = timeElapsed;

* Run bound tightening
$INCLUDE %RelPath%/AC-DC-bound-tightening.gms

* Save intermediate results to GDX
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_tightened bound_tightening . sets data variables

* Record end time and duration
Timing('boundTightening','end') = timeElapsed;
Timing('boundTightening','duration') =
    Timing('boundTightening','end') - Timing('boundTightening','start');


** Algorithm
* Record start time
Timing('algorithm','start') = timeElapsed;

* Include the algorithm
$INCLUDE %RelPath%/AC-DC-algorithm.gms

* Export GDX of results
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_solved main . sets data variables

* Record end time and duration
Timing('algorithm','end') = timeElapsed;
Timing('algorithm','duration') =
    Timing('algorithm','end') - Timing('algorithm','start');

$ENDIF.execmode


* =============================================================================
*  Finish
* =============================================================================
* Export timings
execute_unload '%OutFile%_timings' Timing;

* Send results to MATLAB if applicable
$IFTHENi.ml %FromMatlab% == yes
execute_unload %matout%;
$ENDIF.ml
