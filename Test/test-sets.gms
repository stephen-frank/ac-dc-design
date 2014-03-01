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

$TITLE Test Script: AC-DC Sets
$ONTEXT

    Filename : test-sets.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    Provides a brief test script for testing sets, subsets, and computed
    auxiliary sets for the AC-DC model.
    
    Compile-time control settings:
        CreateTestGDX (yes/no)  Set to "yes" to create the test sets, "no" to
                                read the created test sets
        TestGDXFilename         Name of output GDX file
        RelPath                 Relative path to AC-DC model files
    
    To execute the test:
        1. Set %CreateTestGDX%=yes
        2. Adjust the output filename %TestGDXFilename% as desired
        3. Execute the script to create the GDX with the sets for testing
        4. Set %CreateTestGDX%=no
        5. Re-execute the script to perform the test
        6. Examine the listing file for the test results to ensure they match
           the description given below

    Diagram of Test Network:

                1-------
                    |
                    | (1,2)
                    |
                2-------
                  |  \
            (2,3) |   \ (2,4)
                  |    \
            3-------  -------4
                         \
                          \ (4,5)
                           \
                         -------5

    Test Network Data:
        Buses: {1, 2, 3, 4, 5}

        Branches:   AC      DC
        (1,2)       Y
        (2,3)       Y
        (2,4)       Y       Y
        (4,5)               Y

        Sources:    AC      DC      bus
        utility     Y               1
        pv          Y       Y       2
        fuelcell            Y       5

        Loads:      AC      DC      bus
        1           Y       Y       3
        2           Y               4
        3           Y               2

        Converters: Type            bus
        1           DC->AC          2
        2           AC->DC          4
        3           DC->AC          4

        Note: Load 3 is tied to branch (1,2)

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
* Select mode: creating the GDX file or running the test
$SET CreateTestGDX yes

* Select name for test GDX file
$SET TestGDXFilename SetTest

* Relative path to AC-DC files
$SET RelPath ../Formulation

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


* =============================================================================
*  Create Text GDX File
* =============================================================================
$IFTHEN %CreateTestGDX%==yes

* Master sets
Sets
    A   /DC, AC/
    N   /1*5/
    S   /utility, pv, fuelcell/
    L   /1*3/
    C   /1*3/
    H   /0*1/
    T   /1*5/ ;

Alias (N, I, K);

* Declare subsets
Sets
    B(i,k)          / 1.2, 2.3, 2.4, 4.5 /
    BX(a,i,k)       / AC.1.2, AC.2.3, AC.2.4,
                                      DC.2.4, DC.4.5 /
    SX(a,i,s)       / AC.1.utility, AC.2.pv,
                                    DC.2.pv, DC.5.fuelcell /
    LX(a,i,l)       / AC.3.1, AC.4.2, AC.2.3
                      DC.3.1 /
    CXi(a,i,c)      / DC.2.1, AC.4.2, DC.4.3 /
    CXo(a,i,c)      / AC.2.1, DC.4.2, AC.4.3 /
    BLX(i,k,l)      / 1.2.3 / ;

* Export GDX
execute_unload '%TestGDXFilename%' ;

$ENDIF

$IF %CreateTestGDX%==yes $EXIT


* =============================================================================
*  Perform Test
* =============================================================================
** Setup
* Data File
$SET DataFile %TestGDXFilename%

* Required Flags
$SET LoadSets   yes
$SET LoadParams no
$SET LoadVars   no
$SET LoadBigM   no
$SET InitVars   no

** Include Files
* Default compile time settings
$INCLUDE %RelPath%/AC-DC-defaults.gms

* Set, Parameter, and Variable Definitions
$INCLUDE %RelPath%/AC-DC-defs.gms

* Set Data Processing
$INCLUDE %RelPath%/AC-DC-data.gms

** Display Sets
* Examine these manually to ensure everything looks right...
display A, N, S, L, C, B, H, T;
display BX, SX, LX, CXi, CXo, BLX;
display NN, SS, LL, CCi, CCo;
display ConnD, ConnU;
display HH, TT;

