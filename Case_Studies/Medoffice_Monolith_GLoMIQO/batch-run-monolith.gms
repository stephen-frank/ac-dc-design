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

    Filename : batch-run-monolith.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    This script is configured to run the following test cases:
    
    System          Monolith Solver		Bound Tightening
    ------          ---------------     ----------------
    Medium Office   GLoMIQO             No
    Medium Office   GLoMIQO             Yes

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

* LBP Solver
$SET MIPsolver          CPLEX       // For Big M's / Bound Tightening
$SET MIQCPsolver        GLoMIQO     // For Monolith

* Algorithm settings
$SET ExtraCuts          yes         // Extra cuts in relaxation?
$SET NumSeg             3           // Number of breakpoints in each relaxation
$SET UBP_rel            1e-2        // Monolith relative tolerance
$SET UBP_reslim         1209600     // Monolith problem time limit (sec)
$SET RESLIM             9000        // Default resource limit (sec)

* Local control settings
$SET RunBigM            yes         // Run big M computation?
$SET RunBT              yes         // Run bound tightening?
$SET RunWithBT          yes         // Run monolith with bound tightening?
$SET RunWithoutBT       yes         // Run monolith without bound tightening?

* Relative path to AC-DC model files
$SET RelPath            ../../Formulation

* Name of the test case to solve
$SET TestCase           medoffice


* =============================================================================
*  Create Option Files
* =============================================================================
* Seperate option files are now in use for various stages of the solution
* procedure:
*   1   Presolve/big M's/bound tightening
*   2   Monolith

* -----------------------------------------------------------------------------
*  Optfile 1 - Presolve
* -----------------------------------------------------------------------------

** CPLEX
*   - No special settings
$ONECHO > cplex.opt
eprhs 1e-8
epint 1e-6
$OFFECHO

** GLoMIQO
*   - Decrease the feasibility tolerances
$ONECHO > glomiqo.opt
feas_tolerance 1e-8
$OFFECHO

* -----------------------------------------------------------------------------
*  Optfile 2 - Monolith
* -----------------------------------------------------------------------------

** GLoMIQO
*   - Decrease the feasibility tolerances
$ONECHO > glomiqo.op2
feas_tolerance 1e-8
$OFFECHO


* =============================================================================
*  Execute Setup
* =============================================================================
** File Names
* Various output names for various processes
$SET OutFileBigM                 %TestCase%-%MIPsolver%-big-M
$SET OutFileBT                   %TestCase%-%MIPsolver%-bound-tightened
$SET OutFileMonolithWithBT       %TestCase%-%MIQCPsolver%-Monolith-with-BT
$SET OutFileMonolithWoutBT       %TestCase%-%MIQCPsolver%-Monolith-without-BT

* Various input names for various processes
$SET DataFileBigM                %TestCase%
$SET DataFileBT                  %OutFileBigM%
$SET DataFileMonolithWithBT      %OutFileBT%
$SET DataFileMonolithWoutBT      %OutFileBigM%

* -----------------------------------------------------------------------------
*  Big M's
* -----------------------------------------------------------------------------
** Compute Big M's
$IFTHENi %RunBigM% == yes
* Run GAMS
execute '%GAMSCmd% solve-acdc-model.gms --ExecMode=OnlyBigM --DataFile=%DataFileBigM% --OutFile=%OutFileBigM% --RelPath=%RelPath% --UseLocalSolverOpts=no --ExtraCuts=%ExtraCuts% --ConvexRelaxation=linear --NumSeg=%NumSeg% -MIP=%MIPsolver% -RMIP=%MIPsolver% -MIQCP=%MIQCPsolver% -RMIQCP=%MIQCPsolver% -reslim=%RESLIM% -limrow=0 -limcol=0 -solprint=0 -solvelink=%solvelink.LoadLibrary% -optfile=1 -license=%LicPath% -lo=2 -logfile="%OutFileBigM%.log"'

* Move the resulting listing file
execute '%MoveCmd% solve-acdc-model.lst "%OutFileBigM%.lst"'
$ENDIF

* -----------------------------------------------------------------------------
*  Bound Tightening
* -----------------------------------------------------------------------------
** Run Bound Tightening
$IFTHENi %RunBT% == yes
* Run GAMS
execute '%GAMSCmd% solve-acdc-model.gms --ExecMode=OnlyBT --DataFile=%DataFileBT% --OutFile=%OutFileBT% --RelPath=%RelPath% --UseLocalSolverOpts=no --ExtraCuts=%ExtraCuts% --ConvexRelaxation=linear --NumSeg=%NumSeg% -MIP=%MIPsolver% -RMIP=%MIPsolver% -MIQCP=%MIQCPsolver% -RMIQCP=%MIQCPsolver% -reslim=%RESLIM% -limrow=0 -limcol=0 -solprint=0 -solvelink=%solvelink.LoadLibrary% -optfile=1 -license=%LicPath% -lo=2 -logfile="%OutFileBT%.log"'

* Move the resulting listing file
execute '%MoveCmd% solve-acdc-model.lst "%OutFileBT%.lst"'
$ENDIF


* =============================================================================
*  Monolith - With Bound Tightening
* =============================================================================
$IFTHENi %RunWithBT% == yes

** Execute Monolith
* Run GAMS
execute '%GAMSCmd% solve-acdc-model.gms --ExecMode=OnlyMonolith --DataFile=%DataFileMonolithWithBT% --OutFile=%OutFileMonolithWithBT% --RelPath=%RelPath% --UseLocalSolverOpts=no --UBP_rel=%UBP_rel% --UBP_reslim=%UBP_reslim% --TimeLimit=%TimeLimit% -MIQCP=%MIQCPsolver% -reslim=%RESLIM% -limrow=0 -limcol=0 -solprint=0 -solvelink=%solvelink.LoadLibrary% -optfile=2 -license=%LicPath% -lo=2 -logfile="%OutFileMonolithWithBT%.log"'

* Move the resulting listing file
execute '%MoveCmd% solve-acdc-model.lst "%OutFileMonolithWithBT%.lst"'

* Check the power flow equations
execute '%GAMSCmd% check-acdc-model.gms --DataFile=%OutFileMonolithWithBT% --OutFile=%OutFileMonolithWithBT% --RelPath=%RelPath% -lo=2 -logfile="%OutFileMonolithWithBT%_checked.log"'

$ENDIF


* =============================================================================
*  Monolith - Without Bound Tightening
* =============================================================================
$IFTHENi %RunWithoutBT% == yes

** Execute Monolith
* Run GAMS
execute '%GAMSCmd% solve-acdc-model.gms --ExecMode=OnlyMonolith --DataFile=%DataFileMonolithWoutBT% --OutFile=%OutFileMonolithWoutBT% --RelPath=%RelPath% --UseLocalSolverOpts=no --UBP_rel=%UBP_rel% --UBP_reslim=%UBP_reslim% --TimeLimit=%TimeLimit% -MIQCP=%MIQCPsolver% -reslim=%RESLIM% -limrow=0 -limcol=0 -solprint=0 -solvelink=%solvelink.LoadLibrary% -optfile=2 -license=%LicPath% -lo=2 -logfile="%OutFileMonolithWoutBT%.log"'

* Move the resulting listing file
execute '%MoveCmd% solve-acdc-model.lst "%OutFileMonolithWoutBT%.lst"'

* Check the power flow equations
execute '%GAMSCmd% check-acdc-model.gms --DataFile=%OutFileMonolithWoutBT% --OutFile=%OutFileMonolithWoutBT% --RelPath=%RelPath% -lo=2 -logfile="%OutFileMonolithWoutBT%_checked.log"'

$ENDIF
