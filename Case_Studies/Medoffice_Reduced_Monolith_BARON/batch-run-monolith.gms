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
    
    System                  Monolith Solver     Bound Tightening
    ------                  ---------------     ----------------
    Medium Office (Reduced) BARON               No
    Medium Office (Reduced) BARON               Yes
    
    NOTES:
    1. The test case reduction is done by fixing binary variables in the "Test
       Case Reduction" section under "Execute Setup" below.

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
$SET MIQCPsolver        BARON       // For Monolith

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
*  Optfile 2 - Monolith
* -----------------------------------------------------------------------------

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

* -----------------------------------------------------------------------------
*  Test Case Reduction
* -----------------------------------------------------------------------------
** Description
$ONTEXT
This scenario reduces the available design choices to (a) loads and converters
connected to the rooftop panel HMR, and (b) everything in the 120/208 V part of 
the distribution system. All other portions of the system are fixed to AC.
$OFFTEXT

** Create External Script
* Name to use for external script
$SET ExtScriptName      scenario.gms

* Set up external script for writing
file externalscript 'External script to execute prior to algorithm' / %ExtScriptName% /;
externalscript.ap = 0;

* Write to external script
put externalscript;
$ONPUT
* Fix branches
xB.fx('AC',' 1',' 2') = 1;       // Utility - MHDP
xB.fx('AC',' 2',' 3') = 1;       // MHDP - HL1
xB.fx('AC',' 2',' 4') = 1;       // MHDP - HL2
xB.fx('AC',' 2',' 5') = 1;       // MHDP - HL3
xB.fx('AC',' 2',' 6') = 1;       // MHDP - HM1
xB.fx('AC',' 2',' 7') = 1;       // MHDP - HM2
xB.fx('AC',' 2',' 8') = 1;       // MHDP - HM3
xB.fx('AC',' 2',' 9') = 1;       // MHDP - HMR
xB.fx('AC',' 3','14') = 1;       // HL1 - 1st Floor Lighting
xB.fx('AC',' 3','21') = 1;       // HL1 - Parking Lot Lighting
xB.fx('AC',' 4','15') = 1;       // HL2 - 2nd Floor Lighting
xB.fx('AC',' 5','16') = 1;       // HL3 - 3rd Floor Lighting
xB.fx('AC',' 6','17') = 1;       // HM1 - 1st Floor HVAC
xB.fx('AC',' 6','22') = 1;       // HM1 - Elevators
xB.fx('AC',' 7','18') = 1;       // HM2 - 2nd Floor HVAC
xB.fx('AC',' 8','19') = 1;       // HM3 - 3rd Floor HVAC
*xB.fx('AC',' 9','20') = 1;       // HMR - Rooftop HVAC
*xB.fx('AC',' 9','23') = 1;       // HMR - Rooftop PV Array
*xB.fx('AC',' 8','19') = 1;       // HM3 - 3rd Floor HVAC
*xB.fx('AC','10','11') = 1;       // LDP - LP1
*xB.fx('AC','10','12') = 1;       // LDP - LP2
*xB.fx('AC','10','13') = 1;       // LDP - LP3
*xB.fx('AC','11','24') = 1;       // LP1 - 1st Floor Receptacle Loads
*xB.fx('AC','12','25') = 1;       // LP2 - 2nd Floor Receptacle Loads
*xB.fx('AC','13','26') = 1;       // LP3 - 3rd Floor Receptacle Loads

* Fix sources
xS.fx('AC','utility') = 1;      // Technically, already fixed by formulation
*xS.fx('AC','pv') = 1;

* Fix loads
xL.fx('AC','lighting1') = 1;
xL.fx('AC','lighting2') = 1;
xL.fx('AC','lighting3') = 1;
xL.fx('AC','parkinglights') = 1;
xL.fx('AC','HVAC1') = 1;        // Technically, already fixed by formulation
xL.fx('AC','HVAC2') = 1;        // Technically, already fixed by formulation
xL.fx('AC','HVAC3') = 1;        // Technically, already fixed by formulation
*xL.fx('AC','HVACRoof') = 1;
xL.fx('AC','elevators') = 1;
*xL.fx('AC','receptacle1') = 1;
*xL.fx('AC','receptacle2') = 1;
*xL.fx('AC','receptacle3') = 1;

* Fix converters
xC.fx('R1') = 0;
*xC.fx('R2') = 0;
*xC.fx('R3') = 0;
xC.fx('R4') = 0;
xC.fx('R5') = 0;
xC.fx('R6') = 0;
xC.fx('R7') = 0;
*xC.fx('R8') = 0;
xC.fx('I1') = 0;
*xC.fx('I2') = 0;
xC.fx('I3') = 0;
xC.fx('D1') = 0;
$OFFPUT
putclose externalscript;


* =============================================================================
*  Monolith - With Bound Tightening
* =============================================================================
$IFTHENi %RunWithBT% == yes

** Execute Monolith
* Run GAMS
execute '%GAMSCmd% solve-acdc-model.gms --ExecMode=OnlyMonolith --DataFile=%DataFileMonolithWithBT% --OutFile=%OutFileMonolithWithBT% --RelPath=%RelPath% --ExternalScript=%ExtScriptName% --UseLocalSolverOpts=no --UBP_rel=%UBP_rel% --UBP_reslim=%UBP_reslim% --TimeLimit=%TimeLimit% -MIQCP=%MIQCPsolver% -reslim=%RESLIM% -limrow=0 -limcol=0 -solprint=0 -solvelink=%solvelink.LoadLibrary% -optfile=2 -license=%LicPath% -lo=2 -logfile="%OutFileMonolithWithBT%.log"'

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
execute '%GAMSCmd% solve-acdc-model.gms --ExecMode=OnlyMonolith --DataFile=%DataFileMonolithWoutBT% --OutFile=%OutFileMonolithWoutBT% --RelPath=%RelPath% --ExternalScript=%ExtScriptName% --UseLocalSolverOpts=no --UBP_rel=%UBP_rel% --UBP_reslim=%UBP_reslim% --TimeLimit=%TimeLimit% -MIQCP=%MIQCPsolver% -reslim=%RESLIM% -limrow=0 -limcol=0 -solprint=0 -solvelink=%solvelink.LoadLibrary% -optfile=2 -license=%LicPath% -lo=2 -logfile="%OutFileMonolithWoutBT%.log"'

* Move the resulting listing file
execute '%MoveCmd% solve-acdc-model.lst "%OutFileMonolithWoutBT%.lst"'

* Check the power flow equations
execute '%GAMSCmd% check-acdc-model.gms --DataFile=%OutFileMonolithWoutBT% --OutFile=%OutFileMonolithWoutBT% --RelPath=%RelPath% -lo=2 -logfile="%OutFileMonolithWoutBT%_checked.log"'

$ENDIF
