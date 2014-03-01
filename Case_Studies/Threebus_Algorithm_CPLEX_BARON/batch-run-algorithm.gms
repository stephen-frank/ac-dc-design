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

$TITLE Batch run script for the Mixed AC-DC model
$ONTEXT

    Filename : batch-run-algorithm.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1
    
    This script is configured to run the following test cases:
    
    System          LBP Solver      UBP Solver      Bound Tightening
    ------          ----------      ----------      ----------------
    Three Bus       CPLEX           BARON           No
    Three Bus       CPLEX           BARON           Yes

$OFFTEXT

* =============================================================================
*  Setup
* =============================================================================
* Enable end-of-line comments
$EOLCOM //
$ONEOLCOM

* Path to license file
$SET LicPath ../../license.lic

* Set up the move command (as per the shell)
*$SET MoveCmd move /Y               // Windows
$SET MoveCmd mv -f                  // Unix

* GAMS command (just use 'gams' for normal operation)
$SET GAMSCmd gams

* Solvers
$SET MIPsolver          CPLEX       // LBP Solver
$SET MIQCPsolver        BARON       // UBP Solver

* Algorithm settings
$SET ExtraCuts          yes         // Extra cuts in relaxation?
$SET NumSeg             3           // Number of breakpoints in each relaxation
$SET LBP_rel            5e-3        // Lower bounding problem relative tolerance
$SET UBP_rel            1e-2        // Upper bounding problem relative tolerance
$SET LBP_reslim         1800        // Lower bounding problem time limit (sec)
$SET UBP_reslim         600         // Upper bounding problem time limit (sec)
$SET RESLIM             9000        // Default resource limit (sec)
$SET TimeLimit          86400       // Time limit for BT and main algorithms (sec)
$SET MaxIter            1000        // Maximum number of NGBD iterations

* Local control settings
$SET RunBigM            yes         // Run big M computation?
$SET RunBT              yes         // Run bound tightening?
$SET RunWithBT          yes         // Run algorithm with bound tightening?
$SET RunWithoutBT       yes         // Run algorithm without bound tightening?

* Relative path to AC-DC model files
$SET RelPath            ../../Formulation

* Name of the test case to solve
$SET TestCase           threebus


* =============================================================================
*  Create Option Files
* =============================================================================
* Separate option files are now in use for various stages of the solution
* procedure:
*   1   Presolve/big M's/bound tightening
*   2   Algorithm

* -----------------------------------------------------------------------------
*  Optfile 1 - Presolve
* -----------------------------------------------------------------------------

** CPLEX
*   - No special settings
$ONECHO > cplex.opt
eprhs 1e-8
epint 1e-6
$OFFECHO

** BARON
*   - Decrease the feasibility tolerances
*   - Decrease isolation tolerance for distinct solutions
$ONECHO > baron.opt
AbsConFeasTol 1e-8
RelConFeasTol 1e-12
AbsIntFeasTol 1e-6
RelIntFeasTol 1e-12
BoxTol 1e-8
IsolTol 1e-6
$OFFECHO

* -----------------------------------------------------------------------------
*  Optfile 2 - Algorithm
* -----------------------------------------------------------------------------

** CPLEX
*   - Use barrier method at root node
*   - Default strategy
*   - Tailored cuts
$ONECHO > cplex.op2
eprhs 1e-8
epint 1e-6
startalg 4
mipemphasis 0
zerohalfcuts -1
cliques 0
covers 0
disjcuts -1
fraccuts 0
flowcovers 0
flowpaths -1
gubcovers -1
implbd 0
mcfcuts -1
mircuts 0
symmetry -1
$OFFECHO

** BARON
*   - Decrease the feasibility tolerances
*   - Decrease isolation tolerance for distinct solutions
$ONECHO > baron.op2
AbsConFeasTol 1e-8
RelConFeasTol 1e-12
AbsIntFeasTol 1e-6
RelIntFeasTol 1e-12
BoxTol 1e-8
IsolTol 1e-6
$OFFECHO


* =============================================================================
*  Execute Setup
* =============================================================================
** File Names
* Various output names for various processes
$SET OutFileBigM                 %TestCase%-%MIPsolver%-big-M
$SET OutFileBT                   %TestCase%-%MIPsolver%-bound-tightened
$SET OutFileAlgWithBT            %TestCase%-%MIPsolver%-%MIQCPsolver%-with-BT
$SET OutFileAlgWoutBT            %TestCase%-%MIPsolver%-%MIQCPsolver%-without-BT

* Various input names for various processes
$SET DataFileBigM                %TestCase%
$SET DataFileBT                  %OutFileBigM%
$SET DataFileAlgWithBT           %OutFileBT%
$SET DataFileAlgWoutBT           %OutFileBigM%

* -----------------------------------------------------------------------------
*  Big M's
* -----------------------------------------------------------------------------
** Compute Big M's
$IFTHENi %RunBigM% == yes
* Run GAMS
execute '%GAMSCmd% solve-acdc-model.gms --ExecMode=OnlyBigM --DataFile=%DataFileBigM% --OutFile=%OutFileBigM% --RelPath=%RelPath% --UseLocalSolverOpts=no --ExtraCuts=%ExtraCuts% --ConvexRelaxation=linear --NumSeg=%NumSeg% -LP=%MIPsolver% -MIP=%MIPsolver% -RMIP=%MIPsolver% -MIQCP=%MIQCPsolver% -RMIQCP=%MIQCPsolver% -reslim=%RESLIM% -limrow=0 -limcol=0 -solprint=0 -solvelink=%solvelink.LoadLibrary% -optfile=1 -license=%LicPath% -lo=2 -logfile="%OutFileBigM%.log"'

* Move the resulting listing file
execute '%MoveCmd% solve-acdc-model.lst "%OutFileBigM%.lst"'
$ENDIF

* -----------------------------------------------------------------------------
*  Bound Tightening
* -----------------------------------------------------------------------------
** Run Bound Tightening
$IFTHENi %RunBT% == yes
* Run GAMS
execute '%GAMSCmd% solve-acdc-model.gms --ExecMode=OnlyBT --DataFile=%DataFileBT% --OutFile=%OutFileBT% --RelPath=%RelPath% --UseLocalSolverOpts=no --ExtraCuts=%ExtraCuts% --ConvexRelaxation=linear --NumSeg=%NumSeg% -LP=%MIPsolver% -MIP=%MIPsolver% -RMIP=%MIPsolver% -MIQCP=%MIQCPsolver% -RMIQCP=%MIQCPsolver% -reslim=%RESLIM% -limrow=0 -limcol=0 -solprint=0 -solvelink=%solvelink.LoadLibrary% -optfile=1 -license=%LicPath% -lo=2 -logfile="%OutFileBT%.log"'

* Move the resulting listing file
execute '%MoveCmd% solve-acdc-model.lst "%OutFileBT%.lst"'
$ENDIF


* =============================================================================
*  With Bound Tightening
* =============================================================================
$IFTHENi %RunWithBT% == yes
** Execute Algorithm
* Run GAMS
execute '%GAMSCmd% solve-acdc-model.gms --ExecMode=OnlyAlgorithm --DataFile=%DataFileAlgWithBT% --OutFile=%OutFileAlgWithBT% --RelPath=%RelPath% --UseLocalSolverOpts=no --ExtraCuts=%ExtraCuts% --ConvexRelaxation=linear --NumSeg=%NumSeg% --LBP_rel=%LBP_rel% --UBP_rel=%UBP_rel% --LBP_reslim=%LBP_reslim% --UBP_reslim=%UBP_reslim% --TimeLimit=%TimeLimit% --MaxIter=%MaxIter% -MIP=%MIPsolver% -RMIP=%MIPsolver% -MIQCP=%MIQCPsolver% -RMIQCP=%MIQCPsolver% -reslim=%RESLIM% -limrow=0 -limcol=0 -solprint=0 -solvelink=%solvelink.LoadLibrary% -optfile=2 -license=%LicPath% -lo=2 -logfile="%OutFileAlgWithBT%.log"'

* Move the resulting listing file
execute '%MoveCmd% solve-acdc-model.lst "%OutFileAlgWithBT%.lst"'

* Check the power flow equations
execute '%GAMSCmd% check-acdc-model.gms --DataFile=%OutFileAlgWithBT%_solved --OutFile=%OutFileAlgWithBT% --RelPath=%RelPath% -lo=2 -logfile="%OutFileAlgWithBT%_checked.log"'

$ENDIF


* =============================================================================
*  Without Bound Tightening
* =============================================================================
$IFTHENi %RunWithoutBT% == yes
** Execute Algorithm
* Run GAMS
execute '%GAMSCmd% solve-acdc-model.gms --ExecMode=OnlyAlgorithm --DataFile=%DataFileAlgWoutBT% --OutFile=%OutFileAlgWoutBT% --RelPath=%RelPath% --UseLocalSolverOpts=no --ExtraCuts=%ExtraCuts% --ConvexRelaxation=linear --NumSeg=%NumSeg% --LBP_rel=%LBP_rel% --UBP_rel=%UBP_rel% --LBP_reslim=%LBP_reslim% --UBP_reslim=%UBP_reslim% --TimeLimit=%TimeLimit% --MaxIter=%MaxIter% -MIP=%MIPsolver% -RMIP=%MIPsolver% -MIQCP=%MIQCPsolver% -RMIQCP=%MIQCPsolver% -reslim=%RESLIM% -limrow=0 -limcol=0 -solprint=0 -solvelink=%solvelink.LoadLibrary% -optfile=2 -license=%LicPath% -lo=2 -logfile="%OutFileAlgWoutBT%.log"'

* Move the resulting listing file
execute '%MoveCmd% solve-acdc-model.lst "%OutFileAlgWoutBT%.lst"'

* Check the power flow equations
execute '%GAMSCmd% check-acdc-model.gms --DataFile=%OutFileAlgWoutBT%_solved --OutFile=%OutFileAlgWoutBT% --RelPath=%RelPath% -lo=2 -logfile="%OutFileAlgWoutBT%_checked.log"'

$ENDIF


