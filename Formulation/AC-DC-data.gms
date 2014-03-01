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

$TITLE Data file for AC-DC model
$ONTEXT

    Filename : AC-DC-data.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    This file provides the data import and initialization for the mixed AC-DC
    distribution system design problem. It is designed to be included (via
    $INCLUDE) in other .gms files, and relies on the following external
    compile time flags and settings (which are documented in
    AC-DC-defaults.gms):
        DataFile           Name of the data file to load
        RelaxBinaries      Relax the binary variables to continuous?
        LoadSets           Load the set definitions?
        LoadParams         Load the parameter definitions?
        LoadVars           Load variable levels and bounds?
        LoadBigM           Load previously computed big M values?
        InitVars           Initialize variables (if not loading them)?

    NOTES:
    1. When calling this script, it is possible to select individually whether
       to load sets, parameters, and variables via the appropriate compile-time
       flags. See AC-DC-defaults.gms for the default settings.
    2. If the variables bounds and levels are not loaded, then instead the
       the script will initialize them if %InitVars%==yes. Note that setting
       %LoadVars%==yes implies %InitVars%==no. %InitVars%==yes is the default.
    3. Each load operation depends on having valid sets from either a previous
       load or other user input. So while it is possible to try to initialize
       variables before loading parameters, it will fail unless the parameters
       have been initialized elsewhere. Here are the dependencies:
            Load Sets:              (Nothing)
            Load Parameters:        Valid sets
            Load Variables:         Valid sets
            Load Big M values:      Valid sets
            Initialize Variables:   Valid sets and parameters
    4. In assigning bounds to the reformulation variables, it is necessary to
       consider the signs on the bounds of the underlying variables. For
       reformulation variables based on bilinear terms, this is done by
       checking the four cases:
            x(i).lo * x(k).lo
            x(i).lo * x(k).up
            x(i).up * x(k).lo
            x(i).up * x(k).up
    
$OFFTEXT

* =============================================================================
*  Setup
* =============================================================================
* Note: defaults for compile-time flags are set in 'AC-DC-defaults.gms'

* Initialize data file to import from
$GDXIN %DataFile%

* Don't recompute variable bounds/starting levels if loading them
* (Otherwise, all the values from the load get overwritten)
$IF %LoadVars% == yes           $SET InitVars               no

* =============================================================================
*  Data Import & Initialization
* =============================================================================

* -----------------------------------------------------------------------------
*  Sets & Indices
* -----------------------------------------------------------------------------
* BEGIN SETS
$IFTHEN %LoadSets% == yes

** Import
* Master Sets
$LOADDC N S L C B H T

* Subsets
$LOADDC BX SX LX CXi CXo BLX

** Compute Auxiliary Sets
* Sources that can connect to each network
SS(a,s)$(sum(i,SX(a,i,s))) = yes;

* Loads that can connect to each network
LL(a,l)$(sum(i,LX(a,i,l))) = yes;

* Converters whose input or output can connect to each network
CCi(a,c)$(sum(i,CXi(a,i,c))) = yes;
CCo(a,c)$(sum(i,CXo(a,i,c))) = yes;

* System buses which can be present in each network
NN(a,i)$( sum(k,BX(a,i,k)) + sum(k,BX(a,k,i))
    + sum(s,SX(a,i,s)) + sum(l,LX(a,i,l))
    + sum(c,CXi(a,i,c)) + sum(c,CXo(a,i,c))) = yes;

* Connectivity (directed)
ConnD(a,i,i)$NN(a,i) = yes;
ConnD(a,i,k)$BX(a,i,k) = yes;

* Connectivity (undirected)
ConnU(a,i,i)$NN(a,i) = yes;
ConnU(a,i,k)$BX(a,i,k) = yes;
ConnU(a,i,k)$BX(a,k,i) = yes;

* Mapping of harmonic indices to architectures
HH('DC',h)$sameas(h,'0') = yes;
HH('AC',h)$sameas(h,'1') = yes;

* Mapping of harmonic indices to categories
HX('DC',h)$sameas(h,'0') = yes;
HX('AC',h)$sameas(h,'1') = yes;
HX('Harm',h)$(not (sameas(h,'0') or sameas(h,'1'))) = yes;

* By default, all time periods are under consideration
*   (This set is used later in the decoupled solution algorithm.)
TT(t) = yes;

* END SETS
$ENDIF

* -----------------------------------------------------------------------------
*  Parameters
* -----------------------------------------------------------------------------
* BEGIN PARAMETERS
$IFTHEN %LoadParams% == yes

* All parameters
$LOADDC BusData BranchData SourceData LoadData ConvData ConvData2 Tau

* Correct zero [nominal] values of tap ratio 'A' to 1.0
BranchData(B(i,k),h,t,'A')$(BranchData(i,k,h,t,'A') = 0) = 1.0;

* SPECIAL: CONVERTER HARMONIC SET
* (Must be computed after parameters, if parameters are loaded)
$IFTHEN.B %LoadSets% == yes

* Converters capable of harmonic voltage control
CXh(i,c)$(CXi('AC',i,c) and sum(HX('Harm',h), abs(ConvData2(c,h,'Iimax'))) > 0) = yes;
CXh(i,c)$(CXo('AC',i,c) and sum(HX('Harm',h), abs(ConvData2(c,h,'Iomax'))) > 0) = yes;

$ENDIF.B

* END PARAMETERS
$ENDIF


* -----------------------------------------------------------------------------
*  Big M Values
* -----------------------------------------------------------------------------
* BEGIN BIG M
$IFTHEN %LoadBigM% == yes

* All parameters
$LOADDC MB MS ML

* END BIG M
$ENDIF

* -----------------------------------------------------------------------------
*  Variables
* -----------------------------------------------------------------------------
* BEGIN VARIABLES
$IFTHEN %LoadVars% == yes

* Voltages
$LOADDC V E F

* Powers
$LOADDC PB QB PS QS PL QL PCout QCout PCin QCin

* Currents
$LOADDC IReB IImB IReS IImS IReL IImL IReCout IImCout IReCin IImCin

* Reformulation variables
$LOADDC wVV wEE wEF wFF wC wCC

* Decision variables
$LOADDC xB xS xL xC

* END VARIABLES
$ENDIF


* BEGIN VARIABLE INITIALIZATION -- This controls the rest of the script
$IFTHEN.VARS %InitVars% == yes

* =============================================================================
*  Variable Bounds
* =============================================================================
* Voltage magnitude -- DC and AC fundamental
loop( (a),
    V.lo(i,h,t)$(NN(a,i) and (HX('DC',h) or HX('AC',h))) = BusData(i,h,'Vmin');
    V.up(i,h,t)$(NN(a,i) and (HX('DC',h) or HX('AC',h))) = BusData(i,h,'Vmax');
);

* Voltage in rectangular coordinates -- AC fundamental and AC harmonic
E.lo(i,h,t)$(NN('AC',i) and (HX('AC',h) or HX('Harm',h))) = BusData(i,h,'Emin');
E.up(i,h,t)$(NN('AC',i) and (HX('AC',h) or HX('Harm',h))) = BusData(i,h,'Emax');
F.lo(i,h,t)$(NN('AC',i) and (HX('AC',h) or HX('Harm',h))) = BusData(i,h,'Fmin');
F.up(i,h,t)$(NN('AC',i) and (HX('AC',h) or HX('Harm',h))) = BusData(i,h,'Fmax');

* Converter output real power
PCout.lo(c,h,t) = min( ConvData2(c,h,'Pomin'), 0);
PCout.up(c,h,t) = max( ConvData2(c,h,'Pomax'), 0);
PCout.fx(c,h,t)$(not (CCo('DC',c) or CCo('AC',c))) = 0;

* Converter input real power
PCin.lo(c,h,t) = min( ConvData2(c,h,'Pimin'), 0);
PCin.up(c,h,t) = max( ConvData2(c,h,'Pimax'), 0);
PCout.fx(c,h,t)$(not (CCi('DC',c) or CCi('AC',c))) = 0;

* Converter output reactive power
QCout.lo(c,h,t)$CCo('AC',c) = min( ConvData2(c,h,'Qomin'), 0);
QCout.up(c,h,t)$CCo('AC',c) = max( ConvData2(c,h,'Qomax'), 0);
QCout.fx(c,h,t)$(not CCo('AC',c)) = 0;

* Converter input reactive power
QCin.lo(c,h,t)$CCi('AC',c) = min( ConvData2(c,h,'Qimin'), 0);
QCin.up(c,h,t)$CCi('AC',c) = max( ConvData2(c,h,'Qimax'), 0);
QCin.fx(c,h,t)$(not CCi('AC',c)) = 0;

* Converter output harmonic current
IReCout.lo(c,h,t)$CCo('AC',c) = min( -ConvData2(c,h,'Iomax'), 0);
IReCout.up(c,h,t)$CCo('AC',c) = max(  ConvData2(c,h,'Iomax'), 0);
IReCout.fx(c,h,t)$(not CCo('AC',c)) = 0;

IImCout.lo(c,h,t)$CCo('AC',c) = min( -ConvData2(c,h,'Iomax'), 0);
IImCout.up(c,h,t)$CCo('AC',c) = max(  ConvData2(c,h,'Iomax'), 0);
IImCout.fx(c,h,t)$(not CCo('AC',c)) = 0;

* Converter input harmonic current
IReCin.lo(c,h,t)$CCi('AC',c) = min( -ConvData2(c,h,'Iimax'), 0);
IReCin.up(c,h,t)$CCi('AC',c) = max(  ConvData2(c,h,'Iimax'), 0);
IReCin.fx(c,h,t)$(not CCi('AC',c)) = 0;

IImCin.lo(c,h,t)$CCi('AC',c) = min( -ConvData2(c,h,'Iimax'), 0);
IImCin.up(c,h,t)$CCi('AC',c) = max(  ConvData2(c,h,'Iimax'), 0);
IImCin.fx(c,h,t)$(not CCi('AC',c)) = 0;

* Reformulation variables
$BATINCLUDE %RelPath%/AC-DC-bounds-breaks.gms wbounds


* =============================================================================
*  Fixed Variables
* =============================================================================

* -----------------------------------------------------------------------------
*  Fixed by the Formulation
* -----------------------------------------------------------------------------
** Branches
* Exclusions from each network architecture
xB.fx(a,i,k)$(not BX(a,i,k)) = 0;

** Sources
* Exclusions from each network architecture
xS.fx(a,s)$(not SS(a,s)) = 0;

* Inclusions in each network architecture (based on required exclusions)
xS.fx('AC',s)$(SS('AC',s) and not SS('DC',s)) = 1;
xS.fx('DC',s)$(SS('DC',s) and not SS('AC',s)) = 1;

** Loads
* Exclusions from each network architecture
xL.fx(a,l)$(not LL(a,l)) = 0;

* Inclusions in each network architecture (based on required exclusions)
xL.fx('AC',l)$(LL('AC',l) and not LL('DC',l) and not sum((i,k),BLX(i,k,l))) = 1;
xL.fx('DC',l)$(LL('DC',l) and not LL('AC',l) and not sum((i,k),BLX(i,k,l))) = 1;

* NOTE: Some other variables should be automatically eliminated from the
* formulation during complication; I have not bothered to fix their values.

* =============================================================================
*  Initial Variable Levels
* =============================================================================
** Continuous Variables
* Sets levels for voltages and reformulation voltage variables
$BATINCLUDE %RelPath%/AC-DC-bounds-breaks.gms vlevels

* Sets reformulation converter power variables to zero or their minimum value
wC.l(c,t) = max(0, wC.lo(c,t));
wCC.l(c,t) = max(0, wCC.lo(c,t));

** Binary Decision Variables
* Branches: Use all; use AC network over DC network
xB.l(BX('AC',i,k)) = 1;
xB.l(BX('DC',i,k))$(not(BX('AC',i,k))) = 1;

* Sources: Place in AC network over DC network
xS.l(SS('AC',s)) = 1;
xS.l(SS('DC',s))$(not(SS('AC',s))) = 1;

* Loads: Place in AC network over DC network
xL.l(LL('AC',l))$(not sum((i,k),BLX(i,k,l))) = 1;
xL.l(LL('DC',l))$(not(LL('AC',l) or sum((i,k),BLX(i,k,l)))) = 1;

* Converter assignments - start with all converters disabled
* (May not actually be feasible in some situations.)
xC.l(c) = 0;

* =============================================================================
*  Branching Priorties
* =============================================================================
* Use Converters -> Loads -> Sources -> Branches
xC.prior(c) = 1;
xL.prior(LL(a,l)) = 2;
xS.prior(SS(a,s)) = 3;
xB.prior(BX(a,i,k)) = 4;

* END VARIABLE INITIALIZATION
$ENDIF.VARS
