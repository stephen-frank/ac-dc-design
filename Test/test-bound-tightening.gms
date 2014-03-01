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

$TITLE Test Script: Testing the AC-DC Model bound tightening algorithm
$ONTEXT

    Filename : test-bound-tightening.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1
    
    Provides a test of the bound tightening algorithm for the AC-DC model.
    Requires a valid data file with big M values computed.
    
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
        IncludeHarmonics (yes/no) Include harmonic frequencies during bound
                                tightening?
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
        MaxIter                 Maximum number of bound tightening iterations

    To execute the test:
        1. Set an appropriate file name for %DataFile% (directly or by using
           %RootFilename%)
        2. Ensure that %DataFile% is in the current working directory
        3. Set the GAMS solvers as appropriate for the relaxation type
           selected
        4. Run the script
        5. View the bound tightening variable definitions in the GDX
          '%OutFile%_tightened.gdx'
           
    NOTES:
    1. %DataFile% must include big M values as computed using the script
       'AC-DC-big-M.gms'. The output from 'test-big-M.gms' satisfies this
       requirement.
    2. The default configuration uses %RootFilename% to conveniently change the
       input and output files simultaneously. The input file name
       '%RootFilename%_big_M' is preconfigured to match the output from the
       'test-big-M.gms' script.
    3. The compile-time setting %PiecewiseType% has no effect unless
       %ConvexRelaxation%==piecewise

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
$SET ExtraCuts          yes

* Include harmonics?
$SET IncludeHarmonics   yes

* What convex relaxation for the LBP
*   Choose one of: linear, quadratic, piecewise
*   If 'piecewise', also choose piecewise type: mccormick, logtransform
$SET ConvexRelaxation linear
$SET PiecewiseType mccormick

* Number of segments to use in linear and piecewise linear relaxations:
$SET NumSeg 10

* Max. iterations for test
$SET MaxIter 10

* Set solvers:
option RMIP =   CPLEX;      // Linear relaxation
option MIP =    CPLEX;      // Piecewise linear relaxation
option RMIQCP = GLOMIQO;    // Quadratic relaxation

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

* =============================================================================
*  Setup
* =============================================================================
** Required Flags
* What to load?
$SET LoadSets               yes
$SET LoadParams             yes
$SET LoadVars               no
$SET InitVars               yes
$SET LoadBigM               yes

** Turn off excessive listing
* Equations
option Limrow = 0;

* Variables
option Limcol = 0;

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


** Other
* Set time limit to 5 min during testing
option reslim =     300;

** Export GDX
* Save a pre-solution version of everthing
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_presolve . All sets data variables


* =============================================================================
*  Algorithm
* =============================================================================
* Include the algorithm
$INCLUDE %RelPath%/AC-DC-bound-tightening.gms

* Export GDX of results
$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_tightened bound_tightening LBP sets data variables
*$BATINCLUDE %RelPath%/AC-DC-gdxout.gms %OutFile%_tightened test
