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

$TITLE Test Script: Testing the AC-DC Model solution algorithm
$ONTEXT

    Filename : test-algorithm.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    Tests the NGBD algorithm for the AC-DC model. Valid data is required.
    
    Compile-time control settings:
        RootFilename            Root name of the GDX file for input and output
                                (Used to conveniently specify %DataFile% and
                                %OutFile%)
        DataFile                Name of the input data file
                                (defaults to %RootFilename%_big_M)
        OutFile                 Root name of output GDX file
                                (defaults to %RootFilename%)
        RelPath                 Relative path to AC-DC model files
        ExtraCuts (yes/no)      Generate extra valid inequalities for the LBP?
        IncludeHarmonics (yes/no) Include harmonic frequencies in the model?
        FindAllFeasibleDesigns (yes/no) Attempt to find all feasible designs?
        ExportEveryIteration (yes/no) Export a GDX after each NGBD iteration?
        ConvexRelaxation        Specify one of:
                                linear = McCormick inequalities + valid cuts
                                quadratic = Convex quadratic underestimator,
                                    linear overestimator
                                piecewise = Piecewise linear
        PiecewiseType           For bilinear terms; specify one of:
                                mccormick = Uses McCormick inequalities
                                logtransform = Uses a log transformation
        NumSeg                  Number of segments for linear and piecewise
                                linear relaxations
        MaxIter                 Maximum number of NGBD iterations
        LBP_rel                 Relative tolerance for lower bounding problem
        LBP_abs                 Absolute tolerance for lower bounding problem
        UBP_rel                 Relative tolerance for upper bounding problem
        UBP_abs                 Absolute tolerance for upper bounding problem

    To execute the test:
        1. Set the compile-time control settings to appropriate values
        2. Ensure that %DataFile% is in the current working directory
        3. Select appropriate LBP and UBP solvers as appropriate for the
           relaxation you have chosen
        4. Run the script
        5. Examine the solution in the output file %OutFile%_solved.gdx and/or
           the intermediate results in the other output files produced by
           'AC-DC-algorithm.gms'.
           
    NOTES:
    1. %DataFile% must include big M values as computed using the script
       'AC-DC-big-M.gms'. The outputs from 'test-big-M.gms' and
       'test-bound-tightening.gms' satisfy this requirement.
    2. The default configuration uses %RootFilename% to conveniently change the
       input and output files simultaneously. The input file name
       '%RootFilename%_big_M' is preconfigured to match the output from the
       'test-big-M.gms' script. To match the output from the bound tightening
       test script, set %DataFile% to %RootFilename%_tightened.
    3. The comments in 'AC-DC-algorithm.gms' provide more insight into the
       effects of the various compile-time settings on algorithm behavior and
       performance.

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

* Extra cuts?
$SET ExtraCuts yes

* Include harmonics?
$SET IncludeHarmonics yes

* Find all feasible designs?
* (Use for troubleshooting only)
$SET FindAllFeasibleDesigns no

* Export on every iteration?
* (Use for troubleshooting only)
$SET ExportEveryIteration no

* What convex relaxation for the LBP
*   Choose one of: linear, quadratic, piecewise
*   If 'piecewise', also choose piecewise type: mccormick, logtransform
$SET ConvexRelaxation linear
$SET PiecewiseType mccormick

* Number of segments to use in linear and piecewise linear relaxations:
$SET NumSeg 3

* Max. iterations for test
$SET MaxIter 100

* Set relative LBP and UBP tolerances
$SET LBP_rel            1.00e-3
$SET UBP_rel            1.00e-2

* Set absolute LBP and UBP tolerances
$SET LBP_abs            0
$SET UBP_abs            1e-6

** Solver Options
* Set solvers
option MIP =    CPLEX;          // LBP (linear or piecewise relaxation)
option MIQCP =  GLOMIQO;        // UBP, LBP (quadratic relaxation)

* Set time limit per solve to during testing
option reslim = 120;            // 2 minutes

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

** Set GAMS solve options to increase speed
option solprint = silent;
option solvelink = %solvelink.LoadLibrary%;

** Include Files
* Default compile time settings
$INCLUDE %RelPath%/AC-DC-defaults.gms

* Set, Parameter, and Variable Definitions
$INCLUDE %RelPath%/AC-DC-defs.gms

* Set & Data Processing
$INCLUDE %RelPath%/AC-DC-data.gms

* Equations
$INCLUDE %RelPath%/AC-DC-equations-reform.gms
$INCLUDE %RelPath%/AC-DC-equations-relax.gms

* Models
$INCLUDE %RelPath%/AC-DC-models.gms

** Export GDX
* Save a pre-solution version of everthing
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_presolve . All sets data variables


* =============================================================================
*  Algorithm
* =============================================================================
* Include the algorithm
$INCLUDE %RelPath%/AC-DC-algorithm.gms

* Export GDX of results
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_solved main All sets data variables
