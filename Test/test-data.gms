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

$TITLE Test Script: AC-DC Data
$ONTEXT

    Filename : test-data.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    Provides a brief test script for testing the data import for the AC-DC
    model.

    Compile-time control settings:
        DataFile                Name of the input data file
        OutFile                 Root name of output GDX file
                                (defaults to %DataFile%)
        RelPath                 Relative path to AC-DC model files

    To execute the test:
        1. Set an appropriate file name in %DataFile%
        2. Ensure that %DataFile% is in the current working directory
        3. Run the script
        4. View the generated GDX '%OutFile%_data_import.gdx'

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
* GDX file to import from
$SET DataFile threebus

* GDX file root name to write
$SET OutFile %DataFile%

* Relative path to AC-DC files
$SET RelPath ../Formulation

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


* =============================================================================
*  Perform Test
* =============================================================================
** Setup
* Required Flags
$SET LoadSets       yes
$SET LoadParams     yes
$SET LoadVars       no
$SET LoadBigM       no
$SET InitVars       yes

** Include Files
* Default compile time settings
$INCLUDE %RelPath%/AC-DC-defaults.gms

* Set, Parameter, and Variable Definitions
$INCLUDE %RelPath%/AC-DC-defs.gms

* Set Data Processing
$INCLUDE %RelPath%/AC-DC-data.gms

** Re-export
* Write a full GDX of everything imported and computed, for manual review
execute_unload '%OutFile%_data_import.gdx';
