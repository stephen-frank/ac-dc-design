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

    Filename : check-acdc-model.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    This script checks a solution to the mixed AC-DC distribution system to
    ensure that the nonlinear power flow equations are within tolerance.

    The script is designed to be called from the command line, either via
    another GAMS script or from MATLAB. Therefore, the following must be
    specified when the script is called:
        DataFile            Name of GDX file from which to acquire the solved
                            problem data.

    In addition to the compile-time flags documented in 'AC-DC-defaults.gms',
    the script supports the following optional compile time settings:
        OutFile             Root name of output GDX file, to which the suffix
                            'checked' will be appended. (Default = %DataFile%)

$OFFTEXT

** Enable C and C++ style comments
* End-of-line
$EOLCOM //
$ONEOLCOM

* In-line
$INLINECOM /* */
$ONINLINE

* =============================================================================
*  Setup
* =============================================================================
** File I/O
* Relative path to AC-DC files
$IF NOT SET RelPath             $SET RelPath ../Formulation

* GDX file to import from (pass on command line)
$IF NOT SET DataFile            $ABORT 'DataFile' must be specified.

* GDX file root name to write
$IF NOT SET OutFile             $SET OutFile %DataFile%

** Required Flags
* Load all data from file
$SET LoadSets       yes
$SET LoadParams     yes
$SET LoadVars       yes
$SET LoadBigM       yes
$SET InitVars       no

* =============================================================================
*  Check Data
* =============================================================================

* -----------------------------------------------------------------------------
*  Load Data
* -----------------------------------------------------------------------------
* Default compile time settings
$INCLUDE %RelPath%/AC-DC-defaults.gms

* Set, parameter, and variable Definitions
$INCLUDE %RelPath%/AC-DC-defs.gms

* Load sets & data
$INCLUDE %RelPath%/AC-DC-data.gms

* Fix the 'x' variables at their current levels
xB.fx(BX(a,i,k)) =  round( xB.l(a,i,k) );
xS.fx(SS(a,s)) =    round( xS.l(a,s) );
xL.fx(LL(a,l)) =    round( xL.l(a,l) );
xC.fx(c) =          round( xC.l(c) );

* -----------------------------------------------------------------------------
*  Compute Theoretical Power Injections
* -----------------------------------------------------------------------------
* Theoretical calculation of PB without reformulation variables
parameter PBtarg(i,k,h,t)   'Target values for PB';
PBtarg(i,k,h,t)$B(i,k) =
    // DC power flow
    (
        (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * V.l(i,h,t) * V.l(i,h,t)
        - (1/BranchData(i,k,h,t,'A'))
            * BranchData(i,k,h,t,'gSe')
            * V.l(i,h,t) * V.l(k,h,t)
    )$(xB.l('DC',i,k) = 1 and HH('DC',h))
    // AC power flow
    + (
        (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * V.l(i,h,t) * V.l(i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( E.l(i,h,t) * E.l(k,h,t) + F.l(i,h,t) * F.l(k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * (   BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( E.l(i,h,t) * F.l(k,h,t) - F.l(i,h,t) * E.l(k,h,t) )
    )$(xB.l('AC',i,k) = 1 and HH('AC',h))  ;
PBtarg(k,i,h,t)$B(i,k) =
    // DC power flow
    (
        ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * V.l(k,h,t) * V.l(k,h,t)
        - (1/BranchData(i,k,h,t,'A'))
            * BranchData(i,k,h,t,'gSe')
            * V.l(i,h,t) * V.l(k,h,t)
    )$(xB.l('DC',i,k) = 1 and HH('DC',h))
    // AC power flow
    + (
        ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * V.l(k,h,t) * V.l(k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( E.l(i,h,t) * E.l(k,h,t) + F.l(i,h,t) * F.l(k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * (   BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( E.l(i,h,t) * F.l(k,h,t) - F.l(i,h,t) * E.l(k,h,t) )
    )$(xB.l('AC',i,k) = 1 and HH('AC',h))  ;

* Difference from actual
parameter PBdiff(i,k,h,t)   'Differences in PB between target and actual';
PBdiff(i,k,h,t) = PBtarg(i,k,h,t) - PB.l(i,k,h,t);

* Theoretical calculation of QB without reformulation variables
parameter QBtarg(i,k,h,t)   'Target values for QB';
QBtarg(i,k,h,t)$(BX('AC',i,k) and HH('AC',h)) =
    // AC power flow
    xB.l('AC',i,k) * (
        - (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * V.l(i,h,t) * V.l(i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( + BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( E.l(i,h,t) * E.l(k,h,t) + F.l(i,h,t) * F.l(k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * (   BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( E.l(i,h,t) * F.l(k,h,t) - F.l(i,h,t) * E.l(k,h,t) )
    ) ;
QBtarg(k,i,h,t)$(BX('AC',i,k) and HH('AC',h)) =
    // AC power flow
    xB.l('AC',i,k) * (
        - ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * V.l(k,h,t) * V.l(k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( E.l(i,h,t) * E.l(k,h,t) + F.l(i,h,t) * F.l(k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( E.l(i,h,t) * F.l(k,h,t) - F.l(i,h,t) * E.l(k,h,t) )
    ) ;

* -----------------------------------------------------------------------------
*  Compute Mismatches
* -----------------------------------------------------------------------------
* Difference from actual
parameter QBdiff(i,k,h,t)   'Differences in QB between target and actual';
QBdiff(i,k,h,t) = QBtarg(i,k,h,t) - QB.l(i,k,h,t);

* =============================================================================
*  Save Result
* =============================================================================
* Save to GDX
execute_unload '%OutFile%_checked.gdx' PB, PBtarg, PBdiff, QB, QBtarg, QBdiff;
