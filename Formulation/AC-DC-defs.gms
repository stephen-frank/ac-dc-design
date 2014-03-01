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

$TITLE Definitions file for AC-DC model
$ONTEXT

    Filename : AC-DC-defs.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    This file provides the parameter and variable definitions for the mixed
    AC-DC distribution system design problem. It is designed to be included
    (via $INCLUDE) in other .gms files.

    NOTES:
    1. Variables are declared over full sets, but can be defined in the model
       only for subsets using careful constraint construction. See:
       http://support.gams.com/doku.php?id=gams:define_variables_over_subsets

$OFFTEXT

* =============================================================================
*  Structure
* =============================================================================

* -----------------------------------------------------------------------------
*  Sets & Indices
* -----------------------------------------------------------------------------
** Network Design
* Architectures
Set A       'Set of network architectures'          /DC, AC/;

* Harmonic classifications
Set U       'Set of harmonic index classifications' /DC, AC, Harm/;

* Buses, sources, loads, converters
Sets
    N       'Set of all system buses'
    S       'Set of all system sources'
    L       'Set of all system loads'
    C       'Set of all system converters';

* Alias indices 'i' and 'k' for referencing buses and branches
Alias (N,I,K);

* Branches
Sets
    B(i,k)  'Set of all system branches';

** Other
Sets
    H       'Set of harmonics of interest'
    T       'Set of time periods';

* Aliases required for loop control
Alias (A, A2);
Alias (I, I2);
Alias (K, K2);
Alias (S, S2);
Alias (L, L2);
Alias (C, C2);
Alias (H, H2);
Alias (T, T2);

** Sets / Indices for Looping
* Iterations
Sets
    R       'Iterations' /1*%MaxIter%/
    RR(r)   'Set of completed iterations and solution elimination constraints' ;

* Aliases (required for loop control)
Alias (R, R2);

* -----------------------------------------------------------------------------
*  Subsets
* -----------------------------------------------------------------------------
** Subsets by network architecture
* Branches
Set BX(a,i,k)   'Set of all branches that can be assigned to each network (DC or AC)';

** Subsets by network architecture and bus interconnection
* Sources
Set SX(a,i,s)   'Set of sources that can connect to each bus in each network (DC or AC)';

* Loads
Set LX(a,i,l)   'Set of loads that can connect to each bus in each network (DC or AC)';

* Converters
Sets
    CXi(a,i,c)  'Set of converters with input terminals that can connect to each bus in each network (DC or AC)'
    CXo(a,i,c)  'Set of converters with output terminals that can connect to each bus in each network (DC or AC)';

** Subset of branches that are also loads
* This mainly is used for transformers
Set BLX(i,k,l)  'Maps a branch ik to a load l for branches that are also loads';

* -----------------------------------------------------------------------------
*  Aux. Sets
* -----------------------------------------------------------------------------
* Used to define domains for certain constraints
Sets
    NN(a,i)     'Set of all buses that can be present in each network (DC or AC)'
    SS(a,s)     'Set of all sources that can connect to each network (DC or AC)'
    LL(a,l)     'Set of all loads that can connect to each network (DC or AC)'
    CCi(a,c)    'Set of all converters whose input terminals can connect to each network (DC or AC)'
    CCo(a,c)    'Set of all converters whose output terminals can connect to each network (DC or AC)'
    CXh(i,c)    'Set of all converters capable of harmonic voltage control at each bus'
    HH(a,h)     'Matching set that maps DC to h = 0 and AC to h = 1'
    HX(u,h)     'Set of harmonics by classification (DC, AC Fundamental, or AC Harmonic)';

* NOTE:
* HH(a,h) and HX(u,h) have very similar purpose, and in fact could be the same
* except for GAMS' rules about set dependency. It is necessary in some places
* to cross-reference architecture (a) with harmonic (h), hence the need for set
* HH(a,h). In other places, it is necessary to restrict an equation or variable
* by it's category: DC, AC fundamental, or AC harmonic. Since HH(a,h) cannot
* include the 'AC harmonic' category, a seperate set is required.

* Used to define domain for reformulation variables
Set ConnD(a,i,k) 'Directed connectivity set for nodes i and k for each network architecture';
Set ConnU(a,i,k) 'Undirected connectivity set for nodes i and k for each network architecture';

* Used for decoupling by time period
Set TT(t)       'Set of currently active time periods';


* =============================================================================
*  Data
* =============================================================================

* -----------------------------------------------------------------------------
*  Parameters
* -----------------------------------------------------------------------------
** Field Specification
* Bus Data:
*   Vmin    Minimum allowable voltage magnitude at bus i at harmonic h
*   Vmax    Maximum allowable voltage magnitude at bus i at harmonic h
*   Emin    Minimum allowable real component of complex voltage at bus i at harmonic h
*   Emax    Maximum allowable real component of complex voltage at bus i at harmonic h
*   Fmin    Minimum allowable imaginary component of complex voltage at bus i at harmonic h
*   Fmax    Maximum allowable imaginary component of complex voltage at bus i at harmonic h
Set BusFields 'Valid fields for bus data' /Vmin, Vmax, Emin, Emax, Fmin, Fmax/ ;

* Branch Data:
*   gSe     Series conductance of branch (i,k) at harmonic h in time period t
*   bSe     Series susceptance of branch (i,k) at harmonic h in time period t
*   gSh     Shunt conductance of branch (i,k) at harmonic h in time period t
*   bSh     Shunt susceptance of branch (i,k) at harmonic h in time period t
*   A       Magnitude of complex turns ratio of branch (i,k) at harmonic h in time period t
*   Phi     Angle of complex turns ratio of branch (i,k) at harmonic h in time period t
Set BranchFields 'Valid fields for branch data' /gSe, bSe, gSh, bSh, A, phi/ ;

* Source Data:
*   P       Constant real power for source s at harmonic h in time period t
*   Q       Constant reactive power for source s at harmonic h in time period t
*   ID      Constant current, direct component, for source s at harmonic h in time period t
*   IQ      Constant current, quadrature component, for source s at harmonic h in time period t
*   IRe     Constant current, real component, for source s at harmonic h in time period t
*   IIm     Constant current, imaginary component, for source s at harmonic h in time period t
*   g       Constant conductance for source s at harmonic h in time period t
*   b       Constant susceptance for source s at harmonic h in time period t
Set SourceFields 'Valid fields for source data' /P, Q, ID, IQ, IRe, IIm, g, b/ ;

* NOTE: There must be an AC source s='Utility'. The formulation will minimize
* the AC power drawn from this source.

* Load Data:
*   P       Constant real power for load l at harmonic h in time period t
*   Q       Constant reactive power for load l at harmonic h in time period t
*   ID      Constant current, direct component, for load l at harmonic h in time period t
*   IQ      Constant current, quadrature component, for load l at harmonic h in time period t
*   IRe     Constant current, real component, for load l at harmonic h in time period t
*   IIm     Constant current, imaginary component, for load l at harmonic h in time period t
*   g       Constant conductance for load l at harmonic h in time period t
*   b       Constant susceptance for load l at harmonic h in time period t
Set LoadFields 'Valid fields for load data' /P, Q, ID, IQ, IRe, IIm, g, b/ ;

* Converter Data:
*   alpha   Constant power loss term for converter c
*   beta    Linear loss coefficient for converter c
*   gamma   Quadratic loss coefficient for converter c
*   Pomin   Minimum allowable real power output for converter c at harmonic h
*   Pomax   Maximum allowable real power output for converter c at harmonic h
*   Pimin   Minimum allowable real power input for converter c at harmonic h
*   Pimax   Maximum allowable real power input for converter c at harmonic h
*   Qomin   Minimum allowable reactive power output for converter c at harmonic h
*   Qomax   Maximum allowable reactive power output for converter c at harmonic h
*   Qimin   Minimum allowable reactive power input for converter c at harmonic h
*   Qimax   Maximum allowable reactive power input for converter c at harmonic h
*   Iomax   Maximum allowable harmonic current output for converter c at harmonic h
*   Iimax   Maximum allowable harmonic current input for converter c at harmonic h
Set ConvFields 'Valid fields for converter loss coefficients' /alpha, beta, gamma/ ;
Set ConvFields2 'Valid fields for converter power and current limits' /Pomin, Pomax, Pimin, Pimax, Qomin, Qomax, Qimin, Qimax, Iomax, Iimax/ ;

* NOTE: Fields Pimin and Pimax need not be specified in the input data; these
* fields can be automatically calculated later.

* Primary data
Parameters
    BusData(i,h,BusFields)              'Bus state variable data at each harmonic'
    BranchData(i,k,h,t,BranchFields)    'Branch admittance data at each harmonic'
    SourceData(s,h,t,SourceFields)      'Source ZIP parameters at each harmonic in each time period'
    LoadData(l,h,t,LoadFields)          'Load ZIP parameters at each harmonic in each time period'
    ConvData(c,ConvFields)              'Converter data (loss coefficients)'
    ConvData2(c,h,ConvFields2)          'Converter data (power and current limits)';

* Time Periods
Parameter
    Tau(t)          'Duration of time period t';

* Big M's
Set Mq              'Quantity of big M' /P, Q, IRe, IIm/;
Set Mtype           'Type (direction) of big M: positive or negative relaxation' /pos, neg/;

Parameter
    MB(i,k,h,t,mq,mtype)    'Big M values for branch flow constraints'
    MS(s,h,t,mq,mtype)      'Big M values for source power constraints'
    ML(l,h,t,mq,mtype)      'Big M values for load power constraints';


* =============================================================================
*  Variables
* =============================================================================
* Voltages
Positive variables
    V(i,h,t)        'Voltage magnitude at bus i at harmonic h in time period t';
Variables
    E(i,h,t)        'Real component of complex voltage at bus i at harmonic h in time period t'
    F(i,h,t)        'Imaginary component of complex voltage at bus i at harmonic h in time period t';

* Branch, source, and load powers
Variables
    PB(i,k,h,t)     'Branch real power flowing from bus i toward bus k at harmonic h in time period t'
    QB(i,k,h,t)     'Branch real power flowing from bus i toward bus k at harmonic h in time period t'
    PS(s,h,t)       'Net real power supplied by source s at harmonic h in time period t'
    QS(s,h,t)       'Net reactive power supplied by source s at harmonic h in time period t'
    PL(l,h,t)       'Net real power consumed by load l at harmonic h in time period t'
    QL(l,h,t)       'Net reactive power consumed by load l at harmonic h in time period t';

* Branch, source, and load currents
Variables
    IReB(i,k,h,t)   'Branch harmonic current, real component, flowing from bus i toward bus k at harmonic h in time period t'
    IImB(i,k,h,t)   'Branch harmonic current, imaginary component, flowing from bus i toward bus k at harmonic h in time period t'
    IReS(s,h,t)     'Net harmonic current, real component, supplied by source s at harmonic h in time period t'
    IImS(s,h,t)     'Net harmonic current, imaginary component, supplied by source s at harmonic h in time period t'
    IReL(l,h,t)     'Net harmonic current, real component, consumed by load l at harmonic h in time period t'
    IImL(l,h,t)     'Net harmonic current, imaginary component, consumed by load l at harmonic h in time period t';

* Converter powers
Variables
    PCout(c,h,t)    'Real power delivered at the output terminal of converter c at harmonic h in time period t'
    QCout(c,h,t)    'Reactive power delivered at the output terminal of converter c at harmonic h in time period t'
    PCin(c,h,t)     'Real power consumed at the input terminal of converter c at harmonic h in time period t'
    QCin(c,h,t)     'Reactive power consumed at the input terminal of converter c at harmonic h in time period t';

* Converter currents
Variables
    IReCout(c,h,t)  'Harmonic current, real component, delivered at the output terminal of converter c at harmonic h in time period t'
    IImCout(c,h,t)  'Harmonic current, imaginary component, delivered at the output terminal of converter c at harmonic h in time period t'
    IReCin(c,h,t)   'Harmonic current, real component, consumed at the input terminal of converter c at harmonic h in time period t'
    IImCin(c,h,t)   'Harmonic current, imaginary component, consumed at the input terminal of converter c at harmonic h in time period t';

* Auxiliary (reformulation) variables
Variables
    wVV(i,k,h,t)    'Reformulation variable: voltage magnitude * voltage magnitude'
    wEE(i,k,h,t)    'Reformulation variable: real voltage component * real voltage component'
    wEF(i,k,h,t)    'Reformulation variable: real voltage component * imaginary voltage component'
    wFF(i,k,h,t)    'Reformulation variable: imaginary voltage component * imaginary voltage component'
    wC(c,t)         'Reformulation variable: total converter output power'
    wCC(c,t)        'Reformulation variable: square of total converter output power';

Positive Variables wCC;

* Design decisions
Binary Variables
    xB(a,i,k)       'Assignment of branches (i,k) to networks a (DC, AC)'
    xS(a,s)         'Assignment of sources s to networks a (DC, AC)'
    xL(a,l)         'Assignment of loads l to networks a (DC, AC)'
    xC(c)           'Inclusions (or not) of converters in the system';

* Objective function value
Variable
    z               'Objective function value';
