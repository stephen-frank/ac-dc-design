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

$TITLE Test Script: Solving the AC-DC Model
$ONTEXT

    Filename : test-solve.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    Provides a solve test for the AC-DC model or its relaxation. Valid data
    is required.

    Compile-time control settings:
        RootFilename            Root name of the GDX file for input and output
                                (Used to conveniently specify %DataFile% and
                                %OutFile%)
        DataFile                Name of the input data file
                                (defaults to %RootFilename%_big_M)
        OutFile                 Root name of output GDX file
                                (defaults to %RootFilename%)
        RelPath                 Relative path to AC-DC model files
        SolveLBP (yes/no)       Solve the lower bounding problem (LBP)?
        SolveUBP (yes/no)       Solve the upper bounding problem (UBP)?
        SolveMonolith (yes/no)  Solve the monolith directly?
        UBPOneShot (yes/no)     yes = single solve, no = split by time periods
        ExtraCuts (yes/no)      Generate extra valid inequalities for the LBP?
        IncludeHarmonics (yes/no) Include harmonic frequencies in the model?
        ConvexRelaxation        Specify one of:
                                linear = McCormick inequalities + valid cuts
                                quadratic = Convex quadratic underestimator,
                                    linear overestimator
                                piecewise = Piecewise linear
                                none/no/[Anything else] = Does nothing
        PiecewiseType           For bilinear terms; specify one of:
                                mccormick = Uses McCormick inequalities
                                logtransform = Uses a log transformation
        NumSeg                  Number of segments for linear and piecewise
                                linear relaxations

    To execute the test:
        1. Set the compile-time control settings to appropriate values
        2. Ensure that %DataFile% is in the current working directory
        3. Select what to solve:
            * The LBP (convex relaxation),
            * The UBP (nonconvex but with the binaries fixed), and/or
            * The monolith (everything)
        4. Select appropriate solvers for MINLP, MIP, etc., as appropriate for
           the settings you have chosen
        5. If desired, adjust the optimality criteria and resource limits under
           the "Other" subsection of the "Setup" section of the script
        6. If desired, uncomment and/or edit the appropriate code in the
           "Fix Data Here" section to limit the design space
        7. Run the script
        8. Examine the solution (if actually solving) in the files
            %OutFile%_LBP_solved.gdx
            %OutFile%_UBP_solved.gdx
            %OutFile%_monolith_solved.gdx

    NOTES:
    1. %DataFile% must include big M values as computed using the script
       'AC-DC-big-M.gms'. The outputs from 'test-big-M.gms' and
       'test-bound-tightening.gms' satisfy this requirement.
    2. The default configuration uses %RootFilename% to conveniently change the
       input and output files simultaneously. The input file name
       '%RootFilename%_big_M' is preconfigured to match the output from the
       'test-big-M.gms' script. To match the output from the bound tightening
       test script, set %DataFile% to %RootFilename%_tightened.
    3. Solves proceed in order LBP -> UBP -> Monolith, omitting any models
       that are disabled. The results of each solve are used to initialize the
       next solve.
    4. The monolith is pretty terrible even for small instances; you probably
       don't want to actually attempt a solve except for trivial cases or by
       first getting a good starting point from the LBP and/or UBP.
    5. You can limit the design space by fixing variables in the "Fix Data Here"
       section. This can be useful to examine the solver performance for various
       numbers of binary decision variables or for specific designs. Preset
       lists are included for fully fixing the three bus and medium office test
       cases; edit or comment portions as desired to limit the design space as
       you like.

$OFFTEXT

** Enable C and C++ style comments
* End-of-line
$EOLCOM //
$ONEOLCOM

* In-line
$INLINECOM /* */
$ONINLINE

* * * MODIFY SETTINGS HERE  * * * * * * * * * * * * * * * * * * * * * * * * * *

** Test Controls
* Root file name
$SET RootFilename threebus

* GDX file to import from
$SET DataFile %RootFilename%_big_M

* GDX file root name to write
$SET OutFile %RootFilename%

* Relative path to AC-DC files
$SET RelPath ../Formulation

* What to solve
$SET SolveLBP           yes
$SET SolveUBP           yes
$SET SolveMonolith      no

* If solving UBP, solve it in a single solve?
$SET UBPOneShot         no

* Extra cuts?
$SET ExtraCuts          yes

* Include harmonics?
$SET IncludeHarmonics   yes

* What convex relaxation for the LBP
*   Choose one of: linear, quadratic, piecewise
*   If 'piecewise', also choose piecewise type: mccormick, logtransform
$SET ConvexRelaxation linear
$SET PiecewiseType mccormick

* Number of segments to use in linear and piecewise linear relaxations:
$SET NumSeg 3

* Set solvers
option MIP =    CPLEX;      // LBP
option MIQCP =  GLOMIQO;    // UBP, Monolith

* NOTE: See also "Other" section under "Setup", which contains settings that
* must be set after the model is initialized.

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


* =============================================================================
*  Setup
* =============================================================================
** Required Flags
* What to load?
$SET LoadSets               yes
$SET LoadParams             yes
$SET LoadVars               yes
$SET InitVars               no
$SET LoadBigM               yes

** Include Files
* Default compile time settings
$INCLUDE %RelPath%/AC-DC-defaults.gms

* Set, Parameter, and Variable Definitions
$INCLUDE %RelPath%/AC-DC-defs.gms

* Set & Data Processing
$INCLUDE %RelPath%/AC-DC-data.gms

* Equations
$INCLUDE %RelPath%/AC-DC-equations-reform.gms
$IFTHEN NOT %ConvexRelaxation% == none
$INCLUDE %RelPath%/AC-DC-equations-relax.gms
$ENDIF

* Models
$INCLUDE %RelPath%/AC-DC-models.gms

** Other
* Set optimality criteria
option optCA =      1e-06;  // Absolute tolerance
option optCR =      1e-03;  // Relative tolerance

* Set time limit
option reslim =     1800;

** Export GDX
* Save a pre-solution version of everthing
* (both for reference and for pulling unfixed x limits back in)
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_presolve . All sets data variables equations

* =============================================================================
*  Fix Data Here
* =============================================================================
** For Three Bus Test Case
$ONTEXT
* If you want to fix any binary variables, do so here
xB.fx('DC','1','2') = 1;
xB.fx('AC','1','3') = 1;
xB.fx('DC','2','3') = 1;

xS.fx('DC','pv') = 1;

xL.fx('DC','light') = 1;

xC.fx('rectifier') = 1;
xC.fx('inverter') = 1;
$OFFTEXT

** For Medium Office Test Case
$ONTEXT
display xB.lo, xB.up;

* Fix branches
xB.fx('AC',' 1',' 2') = 1;
xB.fx('AC',' 2',' 3') = 1;
xB.fx('AC',' 2',' 4') = 1;
xB.fx('AC',' 2',' 5') = 1;
xB.fx('AC',' 2',' 6') = 1;
xB.fx('AC',' 2',' 7') = 1;
xB.fx('AC',' 2',' 8') = 1;
xB.fx('AC',' 2',' 9') = 1;
xB.fx('AC',' 2','10') = 0;
xB.fx('AC',' 3','14') = 1;
xB.fx('AC',' 3','21') = 1;
xB.fx('AC',' 4','15') = 1;
xB.fx('AC',' 5','16') = 1;
xB.fx('AC',' 6','17') = 1;
xB.fx('AC',' 6','22') = 1;
xB.fx('AC',' 7','18') = 1;
xB.fx('AC',' 8','19') = 1;
xB.fx('AC',' 9','20') = 1;
xB.fx('AC',' 9','23') = 1;
xB.fx('DC','10','11') = 1;
xB.fx('DC','10','12') = 1;
xB.fx('DC','10','13') = 1;
xB.fx('DC','11','24') = 1;
xB.fx('DC','12','25') = 1;
xB.fx('DC','13','26') = 1;

* Fix sources
xS.fx('AC','utility') = 1;      // Technically, already fixed by formulation
xS.fx('AC','pv') = 1;

* Fix loads
xL.fx('AC','lighting1') = 1;
xL.fx('AC','lighting2') = 1;
xL.fx('AC','lighting3') = 1;
xL.fx('AC','parkinglights') = 1;
xL.fx('AC','HVAC1') = 1;        // Technically, already fixed by formulation
xL.fx('AC','HVAC2') = 1;        // Technically, already fixed by formulation
xL.fx('AC','HVAC3') = 1;        // Technically, already fixed by formulation
xL.fx('AC','HVACRoof') = 1;
xL.fx('AC','elevators') = 1;
xL.fx('DC','receptacle1') = 1;
xL.fx('DC','receptacle2') = 1;
xL.fx('DC','receptacle3') = 1;

* Fix converters
xC.fx('R1') = 0;
xC.fx('R2') = 0;
xC.fx('R3') = 0;
xC.fx('R4') = 0;
xC.fx('R5') = 0;
xC.fx('R6') = 0;
xC.fx('R7') = 0;
xC.fx('R8') = 0;
xC.fx('I1') = 0;
xC.fx('I2') = 0;
xC.fx('I3') = 0;
xC.fx('D1') = 0;
$OFFTEXT

* =============================================================================
*  Solve LBP
* =============================================================================
$IFTHEN.lbp %SolveLBP% == yes
** Setup for solve
* Map fixed variables -> constants
LBP.holdfixed = 1 ;

** Solve command - based on type of relaxation
* Mixed Integer Linear
$IFTHEN %ConvexRelaxation% == linear
solve LBP using MIP minimizing z;

$ELSEIF %ConvexRelaxation% == piecewise
solve LBP using MIP minimizing z;

* Convex Quadratically Constrained
$ELSEIF %ConvexRelaxation% == quadratic
solve LBP using MIQCP minimizing z;

$ENDIF

** Export GDX
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_LBP_solved . LBP sets data variables equations

$ENDIF.lbp


* =============================================================================
*  Solve UBP
* =============================================================================
$IFTHEN.ubp %SolveUBP% == yes
** Setup for solve
* Map fixed variables -> constants
UBP.holdfixed = 1 ;

$ONTEXT
* Resets bounds to original data
$IFTHEN %SolveLBP% == yes
* Reset all variables
execute_load '%OutFile%_presolve' V, E, F, wVV, wEE, wEF, wFF, wC, wCC,
    PB, QB, PS, QS, PL, QL, PCout, QCout, PCin, QCin,
    IReB, IImB, IReS, IImS, IReL, IImL, IReCout, IImCout, IReCin, IImCin,
    xB, xS, xL, xC;

* Load LEVELS only from LBP solve
execute_loadpoint '%OutFile%_LBP_solved' V, E, F, wVV, wEE, wEF, wFF, wC, wCC,
    PB, QB, PS, QS, PL, QL, PCout, QCout, PCin, QCin,
    IReB, IImB, IReS, IImS, IReL, IImL, IReCout, IImCout, IReCin, IImCin,
    xB, xS, xL, xC;
$ENDIF
$OFFTEXT

* Fix the 'x' variables at their current levels
xB.fx(BX(a,i,k)) =  round( xB.l(a,i,k) );
xS.fx(SS(a,s)) =    round( xS.l(a,s) );
xL.fx(LL(a,l)) =    round( xL.l(a,l) );
xC.fx(c) =          round( xC.l(c) );

$IFTHEN.oneshot %UBPOneShot% == yes
** Solve command
* Since holdfixed = 1, should be able to just use a QCP or NLP solver
* However, it must be specified MIQCP or MINLP b/c the x's are present,
* even though fixed...
solve UBP using MIQCP minimizing z;

$ELSE.oneshot
** Alternative - Solve by individual time periods
* Emulates behavior of AC-DC-algorithm.gms
Parameter zz    'Aggregate objective function value';
zz = 0;

* Loop by time periods
TT(t) = no;
loop(t2,
    TT(t2) = yes;
    solve UBP using MIQCP minimizing z;
    zz = zz + z.l;
    TT(t2) = no;
)

* Display aggregate objective function value
display zz;
$ENDIF.oneshot

** Export GDX
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_UBP_solved . UBP sets data variables equations

$ENDIF.ubp


* =============================================================================
*  Solve Monolith
* =============================================================================
$IFTHEN.mono %SolveMonolith% == yes
** Setup for solve
* Map fixed variables -> constants
Monolith.holdfixed = 1 ;

* Restore proper bounds to x's
* (x variables were fixed for UBP; need to unfix them)
$IFTHEN %SolveUBP% == yes
execute_load '%OutFile%_presolve' xB, xS, xL, xC ;
execute_loadpoint '%OutFile%_UBP_solved' xB, xS, xL, xC ;
$ENDIF

** Solve command
* Use MIQCP
solve Monolith using MIQCP minimizing z;

** Export GDX
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_Monolith_solved . Monolith sets data variables
$ENDIF.mono

* =============================================================================
*  Finish
* =============================================================================
* Timing
Parameter Timing(*) 'Records time elapsed';
Timing('CompilationTime')  = timeComp;
Timing('ExecutionTime')    = timeExec;
Timing('TotalElapsedTime') = timeElapsed;

* Display timing
display Timing;
