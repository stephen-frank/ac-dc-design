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

$TITLE Bound tightening algorithm for AC-DC Model
$ONTEXT

    Filename : AC-DC-bound-tightening.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    Implements bound tightening using an iterative solution to the relaxed
    problem. It is designed to be included (via $INCLUDE) in other
    .gms files, and relies on the following external compile time flags
    (which are documented in AC-DC-defaults.gms):
        ConvexRelaxation    Specify type of convex relaxation
        IncludeHarmonics    Include equations for harmonic currents?
        MaxIter             Maximum number of iterations to perform

$OFFTEXT

* Assumption: data and equations are already loaded; all equations are available

* =============================================================================
*  Setup
* =============================================================================

* -----------------------------------------------------------------------------
*  Sets / Indices for Looping
* -----------------------------------------------------------------------------
** Sets
* Selection of a single bus, branch, source, etc.
Sets
    BusSelect(a,i,h,t)      'Select 1 bus + harmonic + time period'
    BranchSelect(a,i,k,h,t) 'Select 1 branch + harmonic + time period'
    SourceSelect(a,i,s,h,t) 'Select 1 source + harmonic + time period'
    LoadSelect(a,i,l,h,t)   'Select 1 load + harmonic + time period'
    ConvSelect(c,t)         'Select 1 converter + time period'
    ConvSelectH(c,h,t)      'Select 1 converter + harmonic + time period'
    xBSelect(a,i,k)         'Select 1 xB variable'
    xSSelect(a,s)           'Select 1 xS variable'
    xLSelect(a,l)           'Select 1 xL variable'
    xCSelect(c)             'Select 1 xC variable' ;

BusSelect(a,i,h,t) = no;
BranchSelect(a,i,k,h,t) = no;
SourceSelect(a,i,s,h,t) = no;
LoadSelect(a,i,l,h,t) = no;
ConvSelect(c,t) = no;
ConvSelectH(c,h,t) = no;
xBSelect(a,i,k) = no;
xSSelect(a,s) = no;
xLSelect(a,l) = no;
xCSelect(c) = no;


* -----------------------------------------------------------------------------
*  Bookkeeping
* -----------------------------------------------------------------------------
* Storage of past bound history
Parameters
    boundV(i,h,t,r,*)       'Bound tightening history for V'
    boundE(i,h,t,r,*)       'Bound tightening history for E'
    boundF(i,h,t,r,*)       'Bound tightening history for F'
//  boundPB(i,k,h,t,r,*)    'Bound tightening history for PB'
//  boundQB(i,k,h,t,r,*)    'Bound tightening history for QB'
//  boundIReB(i,k,h,t,r,*)  'Bound tightening history for IReB'
//  boundIImB(i,k,h,t,r,*)  'Bound tightening history for IImB'
//  boundPS(s,h,t,r,*)      'Bound tightening history for PS'
//  boundQS(s,h,t,r,*)      'Bound tightening history for QS'
//  boundIReS(s,h,t,r,*)    'Bound tightening history for IReS'
//  boundIImS(s,h,t,r,*)    'Bound tightening history for IImS'
//  boundPL(l,h,t,r,*)      'Bound tightening history for PL'
//  boundQL(l,h,t,r,*)      'Bound tightening history for QL'
//  boundIReL(l,h,t,r,*)    'Bound tightening history for IReL'
//  boundIImL(l,h,t,r,*)    'Bound tightening history for IImL'
//  boundPCout(c,h,t,r,*)   'Bound tightening history for PCout'
//  boundQCout(c,h,t,r,*)   'Bound tightening history for QCout'
//  boundIReCout(c,h,t,r,*) 'Bound tightening history for IReCout'
//  boundIImCout(c,h,t,r,*) 'Bound tightening history for IImCout'
//  boundPCin(c,h,t,r,*)    'Bound tightening history for PCin'
//  boundQCin(c,h,t,r,*)    'Bound tightening history for QCin'
//  boundIReCin(c,h,t,r,*)  'Bound tightening history for IReCin'
//  boundIImCin(c,h,t,r,*)  'Bound tightening history for IImCin'
    boundwVV(i,k,h,t,r,*)   'Bound tightening history for wVV'
    boundwEE(i,k,h,t,r,*)   'Bound tightening history for wEE'
    boundwEF(i,k,h,t,r,*)   'Bound tightening history for wEF'
    boundwFF(i,k,h,t,r,*)   'Bound tightening history for wFF'
    boundwC(c,t,r,*)        'Bound tightening history for wC'
    boundwCC(c,t,r,*)       'Bound tightening history for wCC' ;

* Relative change from iteration to iteration
Parameters
    changeV(i,h,t,r,*)          'Relative bound change history for V'
    changeE(i,h,t,r,*)          'Relative bound change history for E'
    changeF(i,h,t,r,*)          'Relative bound change history for F'
//  changePB(i,k,h,t,r,*)       'Relative bound change history for PB'
//  changeQB(i,k,h,t,r,*)       'Relative bound change history for QB'
//  changeIReB(i,k,h,t,r,*)     'Relative bound change history for IReB'
//  changeIImB(i,k,h,t,r,*)     'Relative bound change history for IImB'
//  changePS(s,h,t,r,*)         'Relative bound change history for PS'
//  changeQS(s,h,t,r,*)         'Relative bound change history for QS'
//  changeIReS(s,h,t,r,*)       'Relative bound change history for IReS'
//  changeIImS(s,h,t,r,*)       'Relative bound change history for IImS'
//  changePL(l,h,t,r,*)         'Relative bound change history for PL'
//  changeQL(l,h,t,r,*)         'Relative bound change history for QL'
//  changeIReL(l,h,t,r,*)       'Relative bound change for IReL'
//  changeIImL(l,h,t,r,*)       'Relative bound change for IImL'
//  changePCout(c,h,t,r,*)      'Relative bound change history for PCout'
//  changeQCout(c,h,t,r,*)      'Relative bound change history for QCout'
//  changeIReCout(c,h,t,r,*)    'Relative bound change history for IReCout'
//  changeIImCout(c,h,t,r,*)    'Relative bound change history for IImCout'
//  changePCin(c,h,t,r,*)       'Relative bound change history for PCin'
//  changeQCin(c,h,t,r,*)       'Relative bound change history for QCin'
//  changeIReCin(c,h,t,r,*)     'Relative bound change history for IReCin'
//  changeIImCin(c,h,t,r,*)     'Relative bound change history for IImCin'
    changewVV(i,k,h,t,r,*)      'Relative bound change history for wVV'
    changewEE(i,k,h,t,r,*)      'Relative bound change history for wEE'
    changewEF(i,k,h,t,r,*)      'Relative bound change history for wEF'
    changewFF(i,k,h,t,r,*)      'Relative bound change history for wFF'
    changewC(c,t,r,*)           'Relative bound change history for wC'
    changewCC(c,t,r,*)          'Relative bound change history for wCC' ;

* Loop control and bookkeeping
Parameters
    cvg                     'Convergence flag (yes/no)'
    end_iter                'Records the last iteration performed'
    btStTime                'Bound tightening starting time (for measuring time limits)' ;

* Parameters for saving existing data
Parameters
    save_xB(a,i,k)          'Save initial xB values'
    save_xS(a,s)            'Save initial xS values'
    save_xL(a,l)            'Save initial xL values'
    save_xC(c)              'Save initial xC values' ;

* Temporary storage
Parameter tmp               'Temporary storage for calculations' ;


* -----------------------------------------------------------------------------
*  Initialization
* -----------------------------------------------------------------------------
* Save the existing design decisions
save_xB(a,i,k)  = xB.l(a,i,k);
save_xS(a,s)    = xS.l(a,s);
save_xL(a,l)    = xL.l(a,l);
save_xC(c)      = xC.l(c);

* What type of solver to use
$IFTHEN %ConvexRelaxation% == linear    $SETLOCAL SolverType    RMIP
$ELSEIF %ConvexRelaxation% == quadratic $SETLOCAL SolverType    RMIQCP
$ELSEIF %ConvexRelaxation% == piecewise $SETLOCAL SolverType    MIP
$ENDIF

* Set optimality tolerances
* IMPORTANT: These MUST be as close to zero as possible, as the optimized bound
* is only tightened within these tolerances. (NOTE: Some solvers, e.g. GLoMIQO,
* require a nonzero optCA).
$IF %ConvexRelaxation% == linear        $SETLOCAL OptCR 0
$IF NOT %ConvexRelaxation% == linear    $SETLOCAL OptCR 1e-6
$IF %ConvexRelaxation% == quadratic     $SETLOCAL OptCA 1e-6
$IF NOT %ConvexRelaxation% == quadratic $SETLOCAL OptCA 0

option optCR = %OptCR%;
option optCA = %OptCA%;

** Macro for Troubleshooting
* Detects unboundedness or infeasibility which would imply an invalid solve.
* (Bound tightening solves should always yield feasible and bounded solutions.)
$MACRO assertValidSolve(mod) \
    abort$(&mod&.modelstat = %ModelStat.Unbounded%) \
        'Unboundedness detected in bound tightening calculation'; \
    abort$(&mod&.modelstat = %ModelStat.Infeasible%) \
        'Infeasibility detected in bound tightening calculation';


* =============================================================================
*  Equations
* =============================================================================

* -----------------------------------------------------------------------------
*  Voltages
* -----------------------------------------------------------------------------
** Declarations
Equations
    eq_bound_V(a,i,h,t)     'Isolate each voltage magnitude V(i,h,t)'
    eq_bound_E(a,i,h,t)     'Isolate each voltage real component E(i,h,t)'
    eq_bound_F(a,i,h,t)     'Isolate each voltage imaginary component F(i,h,t)' ;

** Definitions
* Voltage magnitude
eq_bound_V(BusSelect(a,i,h,t))..
    z =e= V(i,h,t) ;

* Voltage real component
eq_bound_E(BusSelect(a,i,h,t))..
    z =e= E(i,h,t) ;

* Voltage magnitude
eq_bound_F(BusSelect(a,i,h,t))..
    z =e= F(i,h,t) ;


* -----------------------------------------------------------------------------
*  Branch Power Flow
* -----------------------------------------------------------------------------
** Declarations
Equations
    eq_bound_PBik(a,i,k,h,t)    'Isolate each branch real power PB(i,k,h,t)'
    eq_bound_PBki(a,i,k,h,t)    'Isolate each branch real power PB(k,i,h,t)'
    eq_bound_QBik(a,i,k,h,t)    'Isolate each branch reactive power QB(i,k,h,t)'
    eq_bound_QBki(a,i,k,h,t)    'Isolate each branch reactive power QB(k,i,h,t)' ;

** Definitions
* Branch real power i -> k
eq_bound_PBik(BranchSelect(a,i,k,h,t))..
    z =e= PB(i,k,h,t) ;

* Branch real power k -> i
eq_bound_PBki(BranchSelect(a,i,k,h,t))..
    z =e= PB(k,i,h,t) ;

* Branch reactive power i -> k
eq_bound_QBik(BranchSelect(a,i,k,h,t))..
    z =e= QB(i,k,h,t) ;

* Branch real power k -> i
eq_bound_QBki(BranchSelect(a,i,k,h,t))..
    z =e= QB(k,i,h,t) ;


* -----------------------------------------------------------------------------
*  Branch Current
* -----------------------------------------------------------------------------
$IFTHEN.harm %IncludeHarmonics% == yes
** Declarations
Equations
    eq_bound_IReBik(a,i,k,h,t)  'Isolate each branch harmonic current, real component IReB(i,k,h,t)'
    eq_bound_IReBki(a,i,k,h,t)  'Isolate each branch harmonic current, real component IReB(k,i,h,t)'
    eq_bound_IImBik(a,i,k,h,t)  'Isolate each branch harmonic current, imaginary component IImB(i,k,h,t)'
    eq_bound_IImBki(a,i,k,h,t)  'Isolate each branch harmonic current, imaginary component IImB(k,i,h,t)' ;

** Definitions
* Branch harmonic current, real component i -> k
eq_bound_IReBik(BranchSelect(a,i,k,h,t))..
    z =e= IReB(i,k,h,t) ;

* Branch harmonic current, real component k -> i
eq_bound_IReBki(BranchSelect(a,i,k,h,t))..
    z =e= IReB(k,i,h,t) ;

* Branch harmonic current, imaginary component i -> k
eq_bound_IImBik(BranchSelect(a,i,k,h,t))..
    z =e= IImB(i,k,h,t) ;

* Branch harmonic current, imaginary component k -> i
eq_bound_IImBki(BranchSelect(a,i,k,h,t))..
    z =e= IImB(k,i,h,t) ;
$ENDIF.harm


* -----------------------------------------------------------------------------
*  Source & Load Power
* -----------------------------------------------------------------------------
** Declarations
Equations
    eq_bound_PS(a,i,s,h,t)    'Isolate each source real power PS(s,h,t)'
    eq_bound_QS(a,i,s,h,t)    'Isolate each source reactive power QS(s,h,t)'
    eq_bound_PL(a,i,l,h,t)    'Isolate each source real power PL(l,h,t)'
    eq_bound_QL(a,i,l,h,t)    'Isolate each source reactive power QL(l,h,t)' ;

** Definitions
* Source real power
eq_bound_PS(SourceSelect(a,i,s,h,t))..
    z =e= PS(s,h,t) ;

* Source reactive power
eq_bound_QS(SourceSelect(a,i,s,h,t))..
    z =e= QS(s,h,t) ;

* Load real power
eq_bound_PL(LoadSelect(a,i,l,h,t))..
    z =e= PL(l,h,t) ;

* Load reactive power
eq_bound_QL(LoadSelect(a,i,l,h,t))..
    z =e= QL(l,h,t) ;


* -----------------------------------------------------------------------------
*  Source & Load Current
* -----------------------------------------------------------------------------
$IFTHEN.harm %IncludeHarmonics% == yes
** Declarations
Equations
    eq_bound_IReS(a,i,s,h,t)    'Isolate each source harmonic current, real component IReS(s,h,t)'
    eq_bound_IImS(a,i,s,h,t)    'Isolate each source harmonic current, imaginary component IImS(s,h,t)'
    eq_bound_IReL(a,i,l,h,t)    'Isolate each source harmonic current, real component IReL(l,h,t)'
    eq_bound_IImL(a,i,l,h,t)    'Isolate each source harmonic current, imaginary component IImL(l,h,t)' ;

** Definitions
* Source harmonic current, real component
eq_bound_IReS(SourceSelect(a,i,s,h,t))..
    z =e= IReS(s,h,t) ;

* Source harmonic current, imaginary component
eq_bound_IImS(SourceSelect(a,i,s,h,t))..
    z =e= IImS(s,h,t) ;

* Load harmonic current, real component
eq_bound_IReL(LoadSelect(a,i,l,h,t))..
    z =e= IReL(l,h,t) ;

* Load harmonic current, imaginary component
eq_bound_IImL(LoadSelect(a,i,l,h,t))..
    z =e= IImL(l,h,t) ;
$ENDIF.harm


* -----------------------------------------------------------------------------
*  Converter Power
* -----------------------------------------------------------------------------
** Declarations
Equations
    eq_bound_PCout(c,h,t)   'Isolate each converter output real power PCout(c,t)'
    eq_bound_QCout(c,h,t)   'Isolate each converter output reactive power QCout(c,t)'
    eq_bound_PCin(c,h,t)    'Isolate each converter input power PCin(c,t)'
    eq_bound_QCin(c,h,t)    'Isolate each converter input reactive power QCin(c,t)'
    eq_bound_wC(c,t)        'Isolate each reformulation variable wC(c,t)'
    eq_bound_wCC(c,t)       'Isolate each reformulation variable wCC(c,t)' ;

** Definitions
* Converter output power
eq_bound_PCout(ConvSelectH(c,h,t))..
    z =e= PCout(c,h,t) ;
eq_bound_QCout(ConvSelectH(c,h,t))..
    z =e= QCout(c,h,t) ;

* Converter input power
eq_bound_PCin(ConvSelectH(c,h,t))..
    z =e= PCin(c,h,t) ;
eq_bound_QCin(ConvSelectH(c,h,t))..
    z =e= QCin(c,h,t) ;

* Converter output power reformulation
eq_bound_wC(ConvSelect(c,t))..
    z =e= wC(c,t) ;
eq_bound_wCC(ConvSelect(c,t))..
    z =e= wCC(c,t) ;


* -----------------------------------------------------------------------------
*  Converter Current
* -----------------------------------------------------------------------------
$IFTHEN.harm %IncludeHarmonics% == yes
** Declarations
Equations
    eq_bound_IReCout(c,h,t)   'Isolate each converter output harmonic current, real component IReCout(c,t)'
    eq_bound_IImCout(c,h,t)   'Isolate each converter output harmonic current, imaginary component IReCout(c,t)'
    eq_bound_IReCin(c,h,t)    'Isolate each converter input harmonic current, real component IReCin(c,t)'
    eq_bound_IImCin(c,h,t)    'Isolate each converter input harmonic current, imaginary component IReCin(c,t)' ;

** Definitions
* Converter output harmonic current
eq_bound_IReCout(ConvSelectH(c,h,t))..
    z =e= IReCout(c,h,t) ;
eq_bound_IImCout(ConvSelectH(c,h,t))..
    z =e= IImCout(c,h,t) ;

* Converter input harmonic current
eq_bound_IReCin(ConvSelectH(c,h,t))..
    z =e= IReCin(c,h,t) ;
eq_bound_IImCin(ConvSelectH(c,h,t))..
    z =e= IImCin(c,h,t) ;
$ENDIF.harm


* -----------------------------------------------------------------------------
*  Reformulation Variables (Voltages)
* -----------------------------------------------------------------------------
** Declarations
Equations
    eq_bound_wVV(a,i,k,h,t) 'Isolate each reformulation variable wVV(i,k,h,t)'
    eq_bound_wEE(a,i,k,h,t) 'Isolate each reformulation variable wEE(i,k,h,t)'
    eq_bound_wEF(a,i,k,h,t) 'Isolate each reformulation variable wEF(i,k,h,t)'
    eq_bound_wFF(a,i,k,h,t) 'Isolate each reformulation variable wFF(i,k,h,t)' ;

** Definitions
* wVV
eq_bound_wVV(BranchSelect(a,i,k,h,t))..
    z =e= wVV(i,k,h,t) ;

* wEE
eq_bound_wEE(BranchSelect(a,i,k,h,t))..
    z =e= wEE(i,k,h,t) ;

* wEF
eq_bound_wEF(BranchSelect(a,i,k,h,t))..
    z =e= wEF(i,k,h,t) ;

* wFF
eq_bound_wFF(BranchSelect(a,i,k,h,t))..
    z =e= wFF(i,k,h,t) ;


* -----------------------------------------------------------------------------
*  Binary Decisions
* -----------------------------------------------------------------------------
** Declarations
Equations
    eq_bound_xB(a,i,k)      'Isolate each branch decision xB'
    eq_bound_xS(a,s)        'Isolate each source decision xS'
    eq_bound_xL(a,l)        'Isolate each load decision xL'
    eq_bound_xC(c)          'Isolate each converter decision xC' ;

** Definitions
* xB
eq_bound_xB(xBSelect(a,i,k))..
    z =e= xB(a,i,k) ;

* xS
eq_bound_xS(xSSelect(a,s))..
    z =e= xS(a,s) ;

* xL
eq_bound_xL(xLSelect(a,l))..
    z =e= xL(a,l) ;

* xC
eq_bound_xC(xCSelect(c))..
    z =e= xC(c) ;


* =============================================================================
*  Models
* =============================================================================

* -----------------------------------------------------------------------------
*  Definitions
* -----------------------------------------------------------------------------
* Each model uses the relaxation with it's own objective

* Feasibility check
model mod_feas          / LBP - eq_obj + eq_dummy / ;

* Voltages
model mod_bound_V       / LBP - eq_obj + eq_bound_V / ;
model mod_bound_E       / LBP - eq_obj + eq_bound_E / ;
model mod_bound_F       / LBP - eq_obj + eq_bound_F / ;

* Reformulation voltages
model mod_bound_wVV     / LBP - eq_obj + eq_bound_wVV / ;
model mod_bound_wEE     / LBP - eq_obj + eq_bound_wEE / ;
model mod_bound_wEF     / LBP - eq_obj + eq_bound_wEF / ;
model mod_bound_wFF     / LBP - eq_obj + eq_bound_wFF / ;

* Branch powers
model mod_bound_PBik    / LBP - eq_obj + eq_bound_PBik / ;
model mod_bound_QBik    / LBP - eq_obj + eq_bound_QBik / ;
model mod_bound_PBki    / LBP - eq_obj + eq_bound_PBki / ;
model mod_bound_QBki    / LBP - eq_obj + eq_bound_QBki / ;

* Source powers
model mod_bound_PS      / LBP - eq_obj + eq_bound_PS / ;
model mod_bound_QS      / LBP - eq_obj + eq_bound_QS / ;

* Load powers
model mod_bound_PL      / LBP - eq_obj + eq_bound_PL / ;
model mod_bound_QL      / LBP - eq_obj + eq_bound_QL / ;

* Converter powers
model mod_bound_PCout   / LBP - eq_obj + eq_bound_PCout / ;
model mod_bound_QCout   / LBP - eq_obj + eq_bound_QCout / ;
model mod_bound_PCin    / LBP - eq_obj + eq_bound_PCin  / ;
model mod_bound_QCin    / LBP - eq_obj + eq_bound_QCin  / ;

$IFTHEN.harm %IncludeHarmonics% == yes
* Branch currents
model mod_bound_IReBik  / LBP - eq_obj + eq_bound_IReBik / ;
model mod_bound_IImBik  / LBP - eq_obj + eq_bound_IImBik / ;
model mod_bound_IReBki  / LBP - eq_obj + eq_bound_IReBki / ;
model mod_bound_IImBki  / LBP - eq_obj + eq_bound_IImBki / ;

* Source currents
model mod_bound_IReS    / LBP - eq_obj + eq_bound_IReS / ;
model mod_bound_IImS    / LBP - eq_obj + eq_bound_IImS / ;

* Load currents
model mod_bound_IReL    / LBP - eq_obj + eq_bound_IReL / ;
model mod_bound_IImL    / LBP - eq_obj + eq_bound_IImL / ;

* Converter currents
model mod_bound_IReCout / LBP - eq_obj + eq_bound_IReCout / ;
model mod_bound_IImCout / LBP - eq_obj + eq_bound_IImCout / ;
model mod_bound_IReCin  / LBP - eq_obj + eq_bound_IReCin  / ;
model mod_bound_IImCin  / LBP - eq_obj + eq_bound_IImCin  / ;
$ENDIF.harm

* Reformulation converter powers
model mod_bound_wC      / LBP - eq_obj + eq_bound_wC /  ;
model mod_bound_wCC     / LBP - eq_obj + eq_bound_wCC / ;

* Binary decisions
model mod_bound_xB      / LBP - eq_obj + eq_bound_xB /  ;
model mod_bound_xS      / LBP - eq_obj + eq_bound_xS /  ;
model mod_bound_xL      / LBP - eq_obj + eq_bound_xL /  ;
model mod_bound_xC      / LBP - eq_obj + eq_bound_xC /  ;


* -----------------------------------------------------------------------------
*  Model Attributes
* -----------------------------------------------------------------------------
* Convert fixed variables to constants (pre-process them out)
mod_feas.holdfixed         = 1;

mod_bound_V.holdfixed      = 1;
mod_bound_E.holdfixed      = 1;
mod_bound_F.holdfixed      = 1;

mod_bound_wVV.holdfixed    = 1;
mod_bound_wEE.holdfixed    = 1;
mod_bound_wEF.holdfixed    = 1;
mod_bound_wFF.holdfixed    = 1;

mod_bound_PBik.holdfixed   = 1;
mod_bound_PBki.holdfixed   = 1;
mod_bound_QBik.holdfixed   = 1;
mod_bound_QBki.holdfixed   = 1;

mod_bound_PS.holdfixed     = 1;
mod_bound_QS.holdfixed     = 1;

mod_bound_PL.holdfixed     = 1;
mod_bound_QL.holdfixed     = 1;

mod_bound_PCout.holdfixed  = 1;
mod_bound_QCout.holdfixed  = 1;
mod_bound_PCin.holdfixed   = 1;
mod_bound_QCin.holdfixed   = 1;

$IFTHEN.harm %IncludeHarmonics% == yes
mod_bound_IReBik.holdfixed   = 1;
mod_bound_IReBki.holdfixed   = 1;
mod_bound_IImBik.holdfixed   = 1;
mod_bound_IImBki.holdfixed   = 1;

mod_bound_IReS.holdfixed     = 1;
mod_bound_IImS.holdfixed     = 1;

mod_bound_IReL.holdfixed     = 1;
mod_bound_IImL.holdfixed     = 1;

mod_bound_IReCout.holdfixed  = 1;
mod_bound_IImCout.holdfixed  = 1;
mod_bound_IReCin.holdfixed   = 1;
mod_bound_IImCin.holdfixed   = 1;
$ENDIF.harm

mod_bound_wC.holdfixed     = 1;
mod_bound_wCC.holdfixed    = 1;

mod_bound_xB.holdfixed     = 1;
mod_bound_xS.holdfixed     = 1;
mod_bound_xL.holdfixed     = 1;
mod_bound_xC.holdfixed     = 1;


* =============================================================================
*  Algorithm
* =============================================================================
** Initialization
* Initialize convergence flag to 'no'
cvg = no;

* Initialize time periods to all off
TT(t) = no;

* Record starting time
btStTime = timeElapsed;


** Feasibility Check Loop
* Test each time period
loop(t2$(timeElapsed - btStTime <= %TimeLimit%),
    // Enable this time period
    TT(t2) = yes;

    // Solve LBP over each time period to verify feasibility
    solve mod_feas using %SolverType% minimizing z;

    // Check for infeasibility
    if( mod_feas.modelstat = %ModelStat.Infeasible% or
        mod_feas.modelstat = %ModelStat.IntegerInfeasible% or
        mod_feas.modelstat = %ModelStat.Infeasible-NoSolution%,

        // LBP infeasible -> abort
        execute_unload '%OutFile%_on_abort';
        abort 'LBP is infeasible; aborting bound tightening.',
            mod_feas.modelstat;
    );

    // Disable this time period
    TT(t2) = no;
);


** Bound Tightening Loop
* Loops until convergence or max. iterations exhausted
loop(r$(not cvg and (timeElapsed - btStTime <= %TimeLimit%)),
    // Start by assuming convergence
    cvg = yes;

    // ---Tighten Voltage Bounds---
    loop ( (NN(a2,i2),h2,t2)$HH(a2,h2),
        // Enable this time period
        TT(t2) = yes;

        // Select this bus
        BusSelect(a2,i2,h2,t2) = yes;

        // Minimize V; tighten bound (within tolerance)
        solve mod_bound_V using %SolverType% minimizing z;
        assertValidSolve(mod_bound_V)
        V.lo(i2,h2,t2) = max( V.lo(i2,h2,t2),
            z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
        tmp = min( sqr(V.lo(i2,h2,t2)), sqr(V.up(i2,h2,t2)) )$(
            sign(V.lo(i2,h2,t2)) = sign(V.up(i2,h2,t2)) );
        wVV.lo(i2,i2,h2,t2) = max( wVV.lo(i2,i2,h2,t2), tmp );

        // Maximize V; tighten bound
        solve mod_bound_V using %SolverType% maximizing z;
        assertValidSolve(mod_bound_V)
        V.up(i2,h2,t2) = min( V.up(i2,h2,t2),
            z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
        tmp = max( sqr(V.lo(i2,h2,t2)), sqr(V.up(i2,h2,t2)) );
        wVV.up(i2,i2,h2,t2) = min( wVV.up(i2,i2,h2,t2), tmp );

        // De-select this bus
        BusSelect(a2,i2,h2,t2) = no;

        // Disable this time period
        TT(t2) = no;
    );
    loop ( (NN('AC',i2),h2,t2)$HH('AC',h2),
        // Enable this time period
        TT(t2) = yes;

        // Select this bus
        BusSelect('AC',i2,h2,t2) = yes;

        // Minimize E; tighten bound
        solve mod_bound_E using %SolverType% minimizing z;
        assertValidSolve(mod_bound_E)
        E.lo(i2,h2,t2) = max( E.lo(i2,h2,t2),
            z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
        tmp = min( sqr(E.lo(i2,h2,t2)), sqr(E.up(i2,h2,t2)) )$(
            sign(E.lo(i2,h2,t2)) = sign(E.up(i2,h2,t2)) );
        wEE.lo(i2,i2,h2,t2) = max( wEE.lo(i2,i2,h2,t2), tmp );

        // Minimize F; tighten bound
        solve mod_bound_F using %SolverType% minimizing z;
        assertValidSolve(mod_bound_F)
        F.lo(i2,h2,t2) = max( F.lo(i2,h2,t2),
            z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
        tmp = min( sqr(F.lo(i2,h2,t2)), sqr(F.up(i2,h2,t2)) )$(
            sign(F.lo(i2,h2,t2)) = sign(F.up(i2,h2,t2)) );
        wFF.lo(i2,i2,h2,t2) = max( wFF.lo(i2,i2,h2,t2), tmp );

        // Maximize E; tighten bound
        solve mod_bound_E using %SolverType% maximizing z;
        assertValidSolve(mod_bound_E)
        E.up(i2,h2,t2) = min( E.up(i2,h2,t2),
            z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
        tmp = max( sqr(E.lo(i2,h2,t2)), sqr(E.up(i2,h2,t2)) );
        wEE.up(i2,i2,h2,t2) = min( wEE.up(i2,i2,h2,t2), tmp );

        // Maximize F; tighten bound
        solve mod_bound_F using %SolverType% maximizing z;
        assertValidSolve(mod_bound_F)
        F.up(i2,h2,t2) = min( F.up(i2,h2,t2),
            z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
        tmp = max( sqr(F.lo(i2,h2,t2)), sqr(F.up(i2,h2,t2)) );
        wFF.up(i2,i2,h2,t2) = min( wFF.up(i2,i2,h2,t2), tmp );

        // De-select this bus
        BusSelect('AC',i2,h2,t2) = no;

        // Disable this time period
        TT(t2) = no;
    );

    // ---Tighten Reformulation Variable (Voltage) Bounds---
    loop ( (ConnD(a2,i2,k2),h2,t2)$HH(a2,h2),
        // Enable this time period
        TT(t2) = yes;

        // Select this voltage pair
        BranchSelect(a2,i2,k2,h2,t2) = yes;

        // Minimize wVV; tighten bound (within tolerance)
        solve mod_bound_wVV using %SolverType% minimizing z;
        assertValidSolve(mod_bound_wVV)
        wVV.lo(i2,k2,h2,t2) = max( wVV.lo(i2,k2,h2,t2),
            z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
        if( sameas(i2,k2),
            // |V| >= sqrt(wVV.lo) -> V >= sqrt(wVV.lo) OR V <= -sqrt(wVV.lo)
            tmp = sqrt(wVV.lo(i2,k2,h2,t2));
            if ( V.lo(i2,h2,t2) > 0, // V > 0 -> V >= sqrt(wVV.lo)
                V.lo(i2,h2,t2) = max( V.lo(i2,h2,t2), tmp );
            );
            if ( V.up(i2,h2,t2) < 0, // V < 0 -> V <= -sqrt(wVV.lo)
                V.up(i2,h2,t2) = min( V.up(i2,h2,t2), -tmp );
            );
        );

        // Maximize wVV; tighten bound
        solve mod_bound_wVV using %SolverType% maximizing z;
        assertValidSolve(mod_bound_wVV)
        wVV.up(i2,k2,h2,t2) = min( wVV.up(i2,k2,h2,t2),
            z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
        if( sameas(i2,k2),
            // |V| <= sqrt(wVV.up) -> V <= sqrt(wVV.up), V >= -sqrt(wVV.up)
            tmp = sqrt(wVV.up(i2,k2,h2,t2));
            V.up(i2,h2,t2) = min( V.up(i2,h2,t2), tmp );
            V.lo(i2,h2,t2) = max( V.lo(i2,h2,t2), -tmp );
        );

        // De-select this voltage pair
        BranchSelect(a2,i2,k2,h2,t2) = no;

        // Disable this time period
        TT(t2) = no;
    );
    loop ( (ConnD('AC',i2,k2),h2,t2)$HH('AC',h2),
        // Enable this time period
        TT(t2) = yes;

        // Select this voltage pair
        BranchSelect(a2,i2,k2,h2,t2) = yes;

        // Minimize wEE; tighten bound (within tolerance)
        solve mod_bound_wEE using %SolverType% minimizing z;
        assertValidSolve(mod_bound_wEE)
        wEE.lo(i2,k2,h2,t2) = max( wEE.lo(i2,k2,h2,t2),
            z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
        if( sameas(i2,k2),
            // |E| >= sqrt(wEE.lo) -> E >= sqrt(wEE.lo) OR E <= -sqrt(wEE.lo)
            tmp = sqrt(wEE.lo(i2,k2,h2,t2));
            if ( E.lo(i2,h2,t2) > 0, // E > 0 -> E >= sqrt(wEE.lo)
                E.lo(i2,h2,t2) = max( E.lo(i2,h2,t2), tmp );
            );
            if ( E.up(i2,h2,t2) < 0, // E < 0 -> E <= -sqrt(wEE.lo)
                E.up(i2,h2,t2) = min( E.up(i2,h2,t2), -tmp );
            );
        );

        // Minimize wFF; tighten bound (within tolerance)
        solve mod_bound_wFF using %SolverType% minimizing z;
        assertValidSolve(mod_bound_wFF)
        wFF.lo(i2,k2,h2,t2) = max( wFF.lo(i2,k2,h2,t2),
            z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
        if( sameas(i2,k2),
            // |F| >= sqrt(wFF.lo) -> F >= sqrt(wFF.lo) OR F <= -sqrt(wFF.lo)
            tmp = sqrt(wFF.lo(i2,k2,h2,t2));
            if ( F.lo(i2,h2,t2) > 0, // F > 0 -> F >= sqrt(wFF.lo)
                F.lo(i2,h2,t2) = max( F.lo(i2,h2,t2), tmp );
            );
            if ( F.up(i2,h2,t2) < 0, // F < 0 -> F <= -sqrt(wFF.lo)
                F.up(i2,h2,t2) = min( F.up(i2,h2,t2), -tmp );
            );
        );

        // Maximize wEE; tighten bound
        solve mod_bound_wEE using %SolverType% maximizing z;
        assertValidSolve(mod_bound_wEE)
        wEE.up(i2,k2,h2,t2) = min( wEE.up(i2,k2,h2,t2),
            z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
        if( sameas(i2,k2),
            // |E| <= sqrt(wEE.up) -> E <= sqrt(wEE.up), E >= -sqrt(wEE.up)
            tmp = sqrt(wEE.up(i2,k2,h2,t2));
            E.up(i2,h2,t2) = min( E.up(i2,h2,t2), tmp );
            E.lo(i2,h2,t2) = max( E.lo(i2,h2,t2), -tmp );
        );

        // Maximize wFF; tighten bound
        solve mod_bound_wFF using %SolverType% maximizing z;
        assertValidSolve(mod_bound_wFF)
        wFF.up(i2,k2,h2,t2) = min( wFF.up(i2,k2,h2,t2),
            z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
        if( sameas(i2,k2),
            // |F| <= sqrt(wFF.up) -> F <= sqrt(wFF.up), F >= -sqrt(wFF.up)
            tmp = sqrt(wFF.up(i2,k2,h2,t2));
            F.up(i2,h2,t2) = min( F.up(i2,h2,t2), tmp );
            F.lo(i2,h2,t2) = max( F.lo(i2,h2,t2), -tmp );
        );

        // De-select this voltage pair
        BranchSelect(a2,i2,k2,h2,t2) = no;

        // Disable this time period
        TT(t2) = no;
    );
    loop ( (ConnU('AC',i2,k2),h2,t2)$HH('AC',h2),
        // Enable this time period
        TT(t2) = yes;

        // Select this voltage pair
        BranchSelect(a2,i2,k2,h2,t2) = yes;

        // Minimize wEF; tighten bound (within tolerance)
        solve mod_bound_wEF using %SolverType% minimizing z;
        assertValidSolve(mod_bound_wEF)
        wEF.lo(i2,k2,h2,t2) = max( wEF.lo(i2,k2,h2,t2),
            z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );

        // Maximize wEF; tighten bound
        solve mod_bound_wEF using %SolverType% maximizing z;
        assertValidSolve(mod_bound_wEF)
        wEF.up(i2,k2,h2,t2) = min( wEF.up(i2,k2,h2,t2),
            z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );

        // De-select this voltage pair
        BranchSelect(a2,i2,k2,h2,t2) = no;

        // Disable this time period
        TT(t2) = no;
    );

    // ---Tighten Converter Power Bounds---
    loop ( (c2,t2),
        // Enable this time period
        TT(t2) = yes;

        // Select this converter
        ConvSelect(c2,t2) = yes;

        // Minimize wC; tighten bound (within tolerance)
        solve mod_bound_wC using %SolverType% minimizing z;
        assertValidSolve(mod_bound_wC)
        wC.lo(c2,t2) = max( wC.lo(c2,t2),
            z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
        tmp = min( sqr(wC.lo(c2,t2)), sqr(wC.up(c2,t2)) )$(
            sign(wC.lo(c2,t2)) = sign(wC.up(c2,t2)) );
        wCC.lo(c2,t2) = max( wCC.lo(c2,t2), tmp );

        // Minimize wCC; tighten bound (within tolerance)
        solve mod_bound_wCC using %SolverType% minimizing z;
        assertValidSolve(mod_bound_wCC)
        wCC.lo(c2,t2) = max( wCC.lo(c2,t2),
            z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
        // |wC| >= sqrt(wCC.lo) -> wC >= sqrt(wCC.lo) OR wC <= -sqrt(wCC.lo)
            tmp = sqrt(wCC.lo(c2,t2));
            if ( wC.lo(c2,t2) > 0, // wC > 0 -> wC >= sqrt(wCC.lo)
                wC.lo(c2,t2) = max( wC.lo(c2,t2), tmp );
            );
            if ( wC.up(c2,t2) < 0, // wC < 0 -> wC <= -sqrt(wCC.lo)
                wC.up(c2,t2) = min( wC.up(c2,t2), -tmp );
            );

        // Maximize wC; tighten bound (within tolerance)
        solve mod_bound_wC using %SolverType% maximizing z;
        assertValidSolve(mod_bound_wC)
        wC.up(c2,t2) = min( wC.up(c2,t2),
            z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
        tmp = max( sqr(wC.lo(c2,t2)), sqr(wC.up(c2,t2)) );
        wCC.up(c2,t2) = min( wCC.up(c2,t2), tmp );

        // Maximize wCC; tighten bound (within tolerance)
        solve mod_bound_wCC using %SolverType% maximizing z;
        assertValidSolve(mod_bound_wCC)
        wCC.up(c2,t2) = min( wCC.up(c2,t2),
            z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
        // |wC| <= sqrt(wCC.up) -> wC <= sqrt(wCC.up), wC >= -sqrt(wCC.up)
            tmp = sqrt(wCC.up(c2,t2));
            wC.up(c2,t2) = min( wC.up(c2,t2), tmp );
            wC.lo(c2,t2) = max( wC.lo(c2,t2), -tmp );

        // De-select this converter
        ConvSelect(c2,t2) = no;

        // Disable this time period
        TT(t2) = no;
    );

    // ---Housekeeping---
    // Update bounds on w variables and linear/piecewise linear breakpoints
$BATINCLUDE %RelPath%/AC-DC-bounds-breaks.gms wbounds2 breakpoints

    // Update bound history
    boundV(i,h,t,r,'lo') = V.lo(i,h,t);
    boundV(i,h,t,r,'up') = V.up(i,h,t);

    boundE(i,h,t,r,'lo') = E.lo(i,h,t);
    boundE(i,h,t,r,'up') = E.up(i,h,t);

    boundF(i,h,t,r,'lo') = F.lo(i,h,t);
    boundF(i,h,t,r,'up') = F.up(i,h,t);

    boundwVV(i,k,h,t,r,'lo') = wVV.lo(i,k,h,t);
    boundwVV(i,k,h,t,r,'up') = wVV.up(i,k,h,t);

    boundwEE(i,k,h,t,r,'lo') = wEE.lo(i,k,h,t);
    boundwEE(i,k,h,t,r,'up') = wEE.up(i,k,h,t);

    boundwEF(i,k,h,t,r,'lo') = wEF.lo(i,k,h,t);
    boundwEF(i,k,h,t,r,'up') = wEF.up(i,k,h,t);

    boundwFF(i,k,h,t,r,'lo') = wFF.lo(i,k,h,t);
    boundwFF(i,k,h,t,r,'up') = wFF.up(i,k,h,t);

    boundwC(c,t,r,'lo') = wC.lo(c,t);
    boundwC(c,t,r,'up') = wC.up(c,t);

    boundwCC(c,t,r,'lo') = wCC.lo(c,t);
    boundwCC(c,t,r,'up') = wCC.up(c,t);

    // Update relative change in bounds
    changeV(i,h,t,r,'lo') =
        0$(boundV(i,h,t,r-1,'lo') = 0 and boundV(i,h,t,r,'lo') = 0)
        + Inf$(boundV(i,h,t,r-1,'lo') = 0 and boundV(i,h,t,r,'lo') <> 0)
        + abs(
            ( boundV(i,h,t,r-1,'lo') - boundV(i,h,t,r,'lo') )
            / boundV(i,h,t,r-1,'lo')
        )$(boundV(i,h,t,r-1,'lo') <> 0) ;
    changeV(i,h,t,r,'up') =
        0$(boundV(i,h,t,r-1,'up') = 0 and boundV(i,h,t,r,'up') = 0)
        + Inf$(boundV(i,h,t,r-1,'up') = 0 and boundV(i,h,t,r,'up') <> 0)
        + abs(
            ( boundV(i,h,t,r-1,'up') - boundV(i,h,t,r,'up') )
            / boundV(i,h,t,r-1,'up')
        )$(boundV(i,h,t,r-1,'up') <> 0) ;

    changeE(i,h,t,r,'lo') =
        0$(boundE(i,h,t,r-1,'lo') = 0 and boundE(i,h,t,r,'lo') = 0)
        + Inf$(boundE(i,h,t,r-1,'lo') = 0 and boundE(i,h,t,r,'lo') <> 0)
        + abs(
            ( boundE(i,h,t,r-1,'lo') - boundE(i,h,t,r,'lo') )
            / boundE(i,h,t,r-1,'lo')
        )$(boundE(i,h,t,r-1,'lo') <> 0) ;
    changeE(i,h,t,r,'up') =
        0$(boundE(i,h,t,r-1,'up') = 0 and boundE(i,h,t,r,'up') = 0)
        + Inf$(boundE(i,h,t,r-1,'up') = 0 and boundE(i,h,t,r,'up') <> 0)
        + abs(
            ( boundE(i,h,t,r-1,'up') - boundE(i,h,t,r,'up') )
            / boundE(i,h,t,r-1,'up')
        )$(boundE(i,h,t,r-1,'up') <> 0) ;

    changeF(i,h,t,r,'lo') =
        0$(boundF(i,h,t,r-1,'lo') = 0 and boundF(i,h,t,r,'lo') = 0)
        + Inf$(boundF(i,h,t,r-1,'lo') = 0 and boundF(i,h,t,r,'lo') <> 0)
        + abs(
            ( boundF(i,h,t,r-1,'lo') - boundF(i,h,t,r,'lo') )
            / boundF(i,h,t,r-1,'lo')
        )$(boundF(i,h,t,r-1,'lo') <> 0) ;
    changeF(i,h,t,r,'up') =
        0$(boundF(i,h,t,r-1,'up') = 0 and boundF(i,h,t,r,'up') = 0)
        + Inf$(boundF(i,h,t,r-1,'up') = 0 and boundF(i,h,t,r,'up') <> 0)
        + abs(
            ( boundF(i,h,t,r-1,'up') - boundF(i,h,t,r,'up') )
            / boundF(i,h,t,r-1,'up')
        )$(boundF(i,h,t,r-1,'up') <> 0) ;

    changewVV(i,k,h,t,r,'lo') =
        0$(boundwVV(i,k,h,t,r-1,'lo') = 0 and boundwVV(i,k,h,t,r,'lo') = 0)
        + Inf$(boundwVV(i,k,h,t,r-1,'lo') = 0 and boundwVV(i,k,h,t,r,'lo') <> 0)
        + abs(
            ( boundwVV(i,k,h,t,r-1,'lo') - boundwVV(i,k,h,t,r,'lo') )
            / boundwVV(i,k,h,t,r-1,'lo')
        )$(boundwVV(i,k,h,t,r-1,'lo') <> 0) ;
    changewVV(i,k,h,t,r,'up') =
        0$(boundwVV(i,k,h,t,r-1,'up') = 0 and boundwVV(i,k,h,t,r,'up') = 0)
        + Inf$(boundwVV(i,k,h,t,r-1,'up') = 0 and boundwVV(i,k,h,t,r,'up') <> 0)
        + abs(
            ( boundwVV(i,k,h,t,r-1,'up') - boundwVV(i,k,h,t,r,'up') )
            / boundwVV(i,k,h,t,r-1,'up')
        )$(boundwVV(i,k,h,t,r-1,'up') <> 0) ;

    changewEE(i,k,h,t,r,'lo') =
        0$(boundwEE(i,k,h,t,r-1,'lo') = 0 and boundwEE(i,k,h,t,r,'lo') = 0)
        + Inf$(boundwEE(i,k,h,t,r-1,'lo') = 0 and boundwEE(i,k,h,t,r,'lo') <> 0)
        + abs(
            ( boundwEE(i,k,h,t,r-1,'lo') - boundwEE(i,k,h,t,r,'lo') )
            / boundwEE(i,k,h,t,r-1,'lo')
        )$(boundwEE(i,k,h,t,r-1,'lo') <> 0) ;
    changewEE(i,k,h,t,r,'up') =
        0$(boundwEE(i,k,h,t,r-1,'up') = 0 and boundwEE(i,k,h,t,r,'up') = 0)
        + Inf$(boundwEE(i,k,h,t,r-1,'up') = 0 and boundwEE(i,k,h,t,r,'up') <> 0)
        + abs(
            ( boundwEE(i,k,h,t,r-1,'up') - boundwEE(i,k,h,t,r,'up') )
            / boundwEE(i,k,h,t,r-1,'up')
        )$(boundwEE(i,k,h,t,r-1,'up') <> 0) ;

    changewEF(i,k,h,t,r,'lo') =
        0$(boundwEF(i,k,h,t,r-1,'lo') = 0 and boundwEF(i,k,h,t,r,'lo') = 0)
        + Inf$(boundwEF(i,k,h,t,r-1,'lo') = 0 and boundwEF(i,k,h,t,r,'lo') <> 0)
        + abs(
            ( boundwEF(i,k,h,t,r-1,'lo') - boundwEF(i,k,h,t,r,'lo') )
            / boundwEF(i,k,h,t,r-1,'lo')
        )$(boundwEF(i,k,h,t,r-1,'lo') <> 0) ;
    changewEF(i,k,h,t,r,'up') =
        0$(boundwEF(i,k,h,t,r-1,'up') = 0 and boundwEF(i,k,h,t,r,'up') = 0)
        + Inf$(boundwEF(i,k,h,t,r-1,'up') = 0 and boundwEF(i,k,h,t,r,'up') <> 0)
        + abs(
            ( boundwEF(i,k,h,t,r-1,'up') - boundwEF(i,k,h,t,r,'up') )
            / boundwEF(i,k,h,t,r-1,'up')
        )$(boundwEF(i,k,h,t,r-1,'up') <> 0) ;

    changewFF(i,k,h,t,r,'lo') =
        0$(boundwFF(i,k,h,t,r-1,'lo') = 0 and boundwFF(i,k,h,t,r,'lo') = 0)
        + Inf$(boundwFF(i,k,h,t,r-1,'lo') = 0 and boundwFF(i,k,h,t,r,'lo') <> 0)
        + abs(
            ( boundwFF(i,k,h,t,r-1,'lo') - boundwFF(i,k,h,t,r,'lo') )
            / boundwFF(i,k,h,t,r-1,'lo')
        )$(boundwFF(i,k,h,t,r-1,'lo') <> 0) ;
    changewFF(i,k,h,t,r,'up') =
        0$(boundwFF(i,k,h,t,r-1,'up') = 0 and boundwFF(i,k,h,t,r,'up') = 0)
        + Inf$(boundwFF(i,k,h,t,r-1,'up') = 0 and boundwFF(i,k,h,t,r,'up') <> 0)
        + abs(
            ( boundwFF(i,k,h,t,r-1,'up') - boundwFF(i,k,h,t,r,'up') )
            / boundwFF(i,k,h,t,r-1,'up')
        )$(boundwFF(i,k,h,t,r-1,'up') <> 0) ;

    changewC(c,t,r,'lo') =
        0$(boundwC(c,t,r-1,'lo') = 0 and boundwC(c,t,r,'lo') = 0)
        + Inf$(boundwC(c,t,r-1,'lo') = 0 and boundwC(c,t,r,'lo') <> 0)
        + abs(
            ( boundwC(c,t,r-1,'lo') - boundwC(c,t,r,'lo') )
            / boundwC(c,t,r-1,'lo')
        )$(boundwC(c,t,r-1,'lo') <> 0) ;
    changewC(c,t,r,'up') =
        0$(boundwC(c,t,r-1,'up') = 0 and boundwC(c,t,r,'up') = 0)
        + Inf$(boundwC(c,t,r-1,'up') = 0 and boundwC(c,t,r,'up') <> 0)
        + abs(
            ( boundwC(c,t,r-1,'up') - boundwC(c,t,r,'up') )
            / boundwC(c,t,r-1,'up')
        )$(boundwC(c,t,r-1,'up') <> 0) ;

    changewCC(c,t,r,'lo') =
        0$(boundwCC(c,t,r-1,'lo') = 0 and boundwCC(c,t,r,'lo') = 0)
        + Inf$(boundwCC(c,t,r-1,'lo') = 0 and boundwCC(c,t,r,'lo') <> 0)
        + abs(
            ( boundwCC(c,t,r-1,'lo') - boundwCC(c,t,r,'lo') )
            / boundwCC(c,t,r-1,'lo')
        )$(boundwCC(c,t,r-1,'lo') <> 0) ;
    changewCC(c,t,r,'up') =
        0$(boundwCC(c,t,r-1,'up') = 0 and boundwCC(c,t,r,'up') = 0)
        + Inf$(boundwCC(c,t,r-1,'up') = 0 and boundwCC(c,t,r,'up') <> 0)
        + abs(
            ( boundwCC(c,t,r-1,'up') - boundwCC(c,t,r,'up') )
            / boundwCC(c,t,r-1,'up')
        )$(boundwCC(c,t,r-1,'up') <> 0) ;

    // ---Convergence Check---
    // Check for changes in bound history
    if( sum( (NN(a,i),h,t)$HH(a,h), changeV(i,h,t,r,'lo') > %CvgTol% ) > 0,
        cvg = no;
    );
    if( sum( (NN(a,i),h,t)$HH(a,h), changeV(i,h,t,r,'up') > %CvgTol% ) > 0,
        cvg = no;
    );

    if( sum( (NN('AC',i),h,t)$HH('AC',h), changeE(i,h,t,r,'lo') > %CvgTol% ) > 0,
        cvg = no;
    );
    if( sum( (NN('AC',i),h,t)$HH('AC',h), changeE(i,h,t,r,'up') > %CvgTol% ) > 0,
        cvg = no;
    );

    if( sum( (NN('AC',i),h,t)$HH('AC',h), changeF(i,h,t,r,'lo') > %CvgTol% ) > 0,
        cvg = no;
    );
    if( sum( (NN('AC',i),h,t)$HH('AC',h), changeF(i,h,t,r,'up') > %CvgTol% ) > 0,
        cvg = no;
    );

    if( sum( (ConnD(a,i,k),h,t)$HH(a,h), changewVV(i,k,h,t,r,'lo') > %CvgTol% ) > 0,
        cvg = no;
    );
    if( sum( (ConnD(a,i,k),h,t)$HH(a,h), changewVV(i,k,h,t,r,'up') > %CvgTol% ) > 0,
        cvg = no;
    );

    if( sum( (ConnD('AC',i,k),h,t)$HH('AC',h), changewEE(i,k,h,t,r,'lo') > %CvgTol% ) > 0,
        cvg = no;
    );
    if( sum( (ConnD('AC',i,k),h,t)$HH('AC',h), changewEE(i,k,h,t,r,'up') > %CvgTol% ) > 0,
        cvg = no;
    );

    if( sum( (ConnD('AC',i,k),h,t)$HH('AC',h), changewFF(i,k,h,t,r,'lo') > %CvgTol% ) > 0,
        cvg = no;
    );
    if( sum( (ConnD('AC',i,k),h,t)$HH('AC',h), changewFF(i,k,h,t,r,'up') > %CvgTol% ) > 0,
        cvg = no;
    );

    if( sum( (ConnU('AC',i,k),h,t)$HH('AC',h), changewEF(i,k,h,t,r,'lo') > %CvgTol% ) > 0,
        cvg = no;
    );
    if( sum( (ConnU('AC',i,k),h,t)$HH('AC',h), changewEF(i,k,h,t,r,'up') > %CvgTol% ) > 0,
        cvg = no;
    );

    if( sum( (c,t), changewC(c,t,r,'lo') > %CvgTol% ) > 0,
        cvg = no;
    );
    if( sum( (c,t), changewC(c,t,r,'up') > %CvgTol% ) > 0,
        cvg = no;
    );

    if( sum( (c,t), changewCC(c,t,r,'lo') > %CvgTol% ) > 0,
        cvg = no;
    );
    if( sum( (c,t), changewCC(c,t,r,'up') > %CvgTol% ) > 0,
        cvg = no;
    );

    // Record this as the highest iteration performed (so far)
    end_iter = ord(r);
);


* =============================================================================
*  Tighten Harmonic Voltage Bounds
* =============================================================================
$IFTHEN.harm %IncludeHarmonics% == yes
* These are not performed in the main loop because they have no associated
* reformulation variables. Therefore, it is unlikely that tightening their
* bounds iteratively would be useful.

* Complex Bus Voltages
loop ( (NN('AC',i2),h2,t2)$(HX('Harm',h2)
    and (timeElapsed - btStTime <= %TimeLimit%)),
    // Enable this time period
    TT(t2) = yes;

    // Select this bus
    BusSelect('AC',i2,h2,t2) = yes;

    // Minimize E; tighten bound
    solve mod_bound_E using %SolverType% minimizing z;
    assertValidSolve(mod_bound_E)
    E.lo(i2,h2,t2) = max( E.lo(i2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );

    // Minimize F; tighten bound
    solve mod_bound_F using %SolverType% minimizing z;
    assertValidSolve(mod_bound_F)
    F.lo(i2,h2,t2) = max( F.lo(i2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );

    // Maximize E; tighten bound
    solve mod_bound_E using %SolverType% maximizing z;
    assertValidSolve(mod_bound_E)
    E.up(i2,h2,t2) = min( E.up(i2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );

    // Maximize F; tighten bound
    solve mod_bound_F using %SolverType% maximizing z;
    assertValidSolve(mod_bound_F)
    F.up(i2,h2,t2) = min( F.up(i2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );

    // De-select this bus
    BusSelect('AC',i2,h2,t2) = no;

    // Disable this time period
    TT(t2) = no;
);
$ENDIF.harm


* =============================================================================
*  Tighten Power Bounds / Big M Values
* =============================================================================
* These are not performed in the main loop because they have no associated
* reformulation variables. Therefore, it is unlikely that tightening their
* bounds iteratively would be useful.

* Branches
loop ( (BX(a2,i2,k2),h2,t2)$(HH(a2,h2) and
    (timeElapsed - btStTime <= %TimeLimit%)),
    // Enable this time period
    TT(t2) = yes;

    // Select this branch
    BranchSelect(a2,i2,k2,h2,t2) = yes;

    // Minimize branch real power i -> k; tighten bound (within tolerance)
    solve mod_bound_PBik using %SolverType% minimizing z;
    assertValidSolve(mod_bound_PBik)
    PB.lo(i2,k2,h2,t2) = max( PB.lo(i2,k2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
    MB(i2,k2,h2,t2,'P','neg') = max( -PB.lo(i2,k2,h2,t2) + %BigMSlack%, 0);

    // Maximize branch real power k -> i; tighten bound (within tolerance)
    solve mod_bound_PBki using %SolverType% maximizing z;
    assertValidSolve(mod_bound_PBki)
    PB.up(k2,i2,h2,t2) = min( PB.up(k2,i2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
    MB(k2,i2,h2,t2,'P','pos') = max( PB.up(k2,i2,h2,t2) + %BigMSlack%, 0);

    // Minimize branch real power k -> i; tighten bound (within tolerance)
    solve mod_bound_PBki using %SolverType% minimizing z;
    assertValidSolve(mod_bound_PBki)
    PB.lo(k2,i2,h2,t2) = max( PB.lo(k2,i2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
    MB(k2,i2,h2,t2,'P','neg') = max( -PB.lo(k2,i2,h2,t2) + %BigMSlack%, 0);

    // Maximize branch real power i -> k; tighten bound (within tolerance)
    solve mod_bound_PBik using %SolverType% maximizing z;
    assertValidSolve(mod_bound_PBik)
    PB.up(i2,k2,h2,t2) = min( PB.up(i2,k2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
    MB(i2,k2,h2,t2,'P','pos') = max( PB.up(i2,k2,h2,t2) + %BigMSlack%, 0);

    // Branch reactive power
    if ( sameas(a2, 'AC'),
        // Minimize branch reactive power i -> k; tighten bound (within tolerance)
        solve mod_bound_QBik using %SolverType% minimizing z;
        assertValidSolve(mod_bound_QBik)
        QB.lo(i2,k2,h2,t2) = max( QB.lo(i2,k2,h2,t2),
            z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
        MB(i2,k2,h2,t2,'Q','neg') = max( -QB.lo(i2,k2,h2,t2) + %BigMSlack%, 0);

        // Maximize branch reactive power k -> i; tighten bound (within tolerance)
        solve mod_bound_QBki using %SolverType% maximizing z;
        assertValidSolve(mod_bound_QBki)
        QB.up(k2,i2,h2,t2) = min( QB.up(k2,i2,h2,t2),
            z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
        MB(k2,i2,h2,t2,'Q','pos') = max( QB.up(k2,i2,h2,t2) + %BigMSlack%, 0);

        // Minimize branch reactive power k -> i; tighten bound (within tolerance)
        solve mod_bound_QBki using %SolverType% minimizing z;
        assertValidSolve(mod_bound_QBki)
        QB.lo(k2,i2,h2,t2) = max( QB.lo(k2,i2,h2,t2),
            z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
        MB(k2,i2,h2,t2,'Q','neg') = max( -QB.lo(k2,i2,h2,t2) + %BigMSlack%, 0);

        // Maximize branch reactive power i -> k; tighten bound (within tolerance)
        solve mod_bound_QBik using %SolverType% maximizing z;
        assertValidSolve(mod_bound_QBik)
        QB.up(i2,k2,h2,t2) = min( QB.up(i2,k2,h2,t2),
            z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
        MB(i2,k2,h2,t2,'Q','pos') = max( QB.up(i2,k2,h2,t2) + %BigMSlack%, 0);
    );

    // De-select this branch
    BranchSelect(a2,i2,k2,h2,t2) = no;

    // Disable this time period
    TT(t2) = no;
);

* Sources
loop ( (SX(a2,i2,s2),h2,t2)$(HH(a2,h2) and
    (timeElapsed - btStTime <= %TimeLimit%)),
    // Enable this time period
    TT(t2) = yes;

    // Select this source
    SourceSelect(a2,i2,s2,h2,t2) = yes;

    // Minimize source real power; tighten bound (within tolerance)
    solve mod_bound_PS using %SolverType% minimizing z;
    assertValidSolve(mod_bound_PS)
    PS.lo(s2,h2,t2) = max( PS.lo(s2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
    MS(s2,h2,t2,'P','neg') = max( -PS.lo(s2,h2,t2) + %BigMSlack%, 0);

    // Maximize source real power; tighten bound (within tolerance)
    solve mod_bound_PS using %SolverType% maximizing z;
    assertValidSolve(mod_bound_PS)
    PS.up(s2,h2,t2) = min( PS.up(s2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
    MS(s2,h2,t2,'P','pos') = max( PS.up(s2,h2,t2) + %BigMSlack%, 0);

    // Source reactive power
    if ( sameas(a2, 'AC'),
        // Minimize source reactive power; tighten bound (within tolerance)
        solve mod_bound_QS using %SolverType% minimizing z;
        assertValidSolve(mod_bound_QS)
        QS.lo(s2,h2,t2) = max( QS.lo(s2,h2,t2),
            z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
        MS(s2,h2,t2,'Q','neg') = max( -QS.lo(s2,h2,t2) + %BigMSlack%, 0);

        // Maximize source reactive power; tighten bound (within tolerance)
        solve mod_bound_QS using %SolverType% maximizing z;
        assertValidSolve(mod_bound_QS)
        QS.up(s2,h2,t2) = min( QS.up(s2,h2,t2),
            z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
        MS(s2,h2,t2,'Q','pos') = max( QS.up(s2,h2,t2) + %BigMSlack%, 0);
    );

    // De-select this source
    SourceSelect(a2,i2,s2,h2,t2) = no;

    // Disable this time period
    TT(t2) = no;
);

* Loads
loop ( (LX(a2,i2,l2),h2,t2)$(HH(a2,h2) and
    (timeElapsed - btStTime <= %TimeLimit%)),
    // Enable this time period
    TT(t2) = yes;

    // Select this load
    LoadSelect(a2,i2,l2,h2,t2) = yes;

    // Minimize load real power; tighten bound (within tolerance)
    solve mod_bound_PL using %SolverType% minimizing z;
    assertValidSolve(mod_bound_PL)
    PL.lo(l2,h2,t2) = max( PL.lo(l2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
    ML(l2,h2,t2,'P','neg') = max( -PL.lo(l2,h2,t2) + %BigMSlack%, 0);

    // Maximize load real power; tighten bound (within tolerance)
    solve mod_bound_PL using %SolverType% maximizing z;
    assertValidSolve(mod_bound_PL)
    PL.up(l2,h2,t2) = min( PL.up(l2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
    ML(l2,h2,t2,'P','pos') = max( PL.up(l2,h2,t2) + %BigMSlack%, 0);

    // Load reactive power
    if ( sameas(a2, 'AC'),
        // Minimize load reactive power; tighten bound (within tolerance)
        solve mod_bound_QL using %SolverType% minimizing z;
        assertValidSolve(mod_bound_QL)
        QL.lo(l2,h2,t2) = max( QL.lo(l2,h2,t2),
            z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
        ML(l2,h2,t2,'Q','neg') = max( -QL.lo(l2,h2,t2) + %BigMSlack%, 0);

        // Maximize load reactive power; tighten bound (within tolerance)
        solve mod_bound_QL using %SolverType% maximizing z;
        assertValidSolve(mod_bound_QL)
        QL.up(l2,h2,t2) = min( QL.up(l2,h2,t2),
            z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
        ML(l2,h2,t2,'Q','pos') = max( QL.up(l2,h2,t2) + %BigMSlack%, 0);
    );

    // De-select this load
    LoadSelect(a2,i2,l2,h2,t2) = no;

    // Disable this time period
    TT(t2) = no;
);

* Converters
loop ( (c2,h2,t2)$(timeElapsed - btStTime <= %TimeLimit%),
    // Enable this time period
    TT(t2) = yes;

    // Select this converter + harmonic
    ConvSelectH(c2,h2,t2) = yes;

    // Minimize PCout; tighten bound (within tolerance)
    solve mod_bound_PCout using %SolverType% minimizing z;
    assertValidSolve(mod_bound_PCout)
    PCout.lo(c2,h2,t2) = max( PCout.lo(c2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );

    // Minimize PCin; tighten bound (within tolerance)
    solve mod_bound_PCin using %SolverType% minimizing z;
    assertValidSolve(mod_bound_PCin)
    PCin.lo(c2,h2,t2) = max( PCin.lo(c2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );

    // Maximize PCout; tighten bound (within tolerance)
    solve mod_bound_PCout using %SolverType% maximizing z;
    assertValidSolve(mod_bound_PCout)
    PCout.up(c2,h2,t2) = min( PCout.up(c2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );

    // Maximize PCin; tighten bound (within tolerance)
    solve mod_bound_PCin using %SolverType% maximizing z;
    assertValidSolve(mod_bound_PCin)
    PCin.up(c2,h2,t2) = min( PCin.up(c2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );

    // Input reactive power
    if ( CCi('AC',c2),
        // Minimize QCin; tighten bound (within tolerance)
        solve mod_bound_QCin using %SolverType% minimizing z;
        assertValidSolve(mod_bound_QCin)
        QCin.lo(c2,h2,t2) = max( QCin.lo(c2,h2,t2),
            z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );

        // Maximize QCin; tighten bound (within tolerance)
        solve mod_bound_QCin using %SolverType% maximizing z;
        assertValidSolve(mod_bound_QCin)
        QCin.up(c2,h2,t2) = min( QCin.up(c2,h2,t2),
            z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
    );

    // Output reactive power
    if ( CCo('AC',c2),
        // Minimize QCout; tighten bound (within tolerance)
        solve mod_bound_QCout using %SolverType% minimizing z;
        assertValidSolve(mod_bound_QCout)
        QCout.lo(c2,h2,t2) = max( QCout.lo(c2,h2,t2),
            z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );

        // Maximize QCout; tighten bound (within tolerance)
        solve mod_bound_QCout using %SolverType% maximizing z;
        assertValidSolve(mod_bound_QCout)
        QCout.up(c2,h2,t2) = min( QCout.up(c2,h2,t2),
            z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
    );

    // De-select this converter + harmonic
    ConvSelectH(c2,h2,t2) = no;

    // Disable this time period
    TT(t2) = no;
);


* =============================================================================
*  Tighten Current Bounds / Big M Values
* =============================================================================
$IFTHEN.harm %IncludeHarmonics% == yes
* These are not performed in the main loop because they have no associated
* reformulation variables. Therefore, it is unlikely that tightening their
* bounds iteratively would be useful.

* Branches
loop ( (BX('AC',i2,k2),h2,t2)$(HX('Harm',h2) and
    (timeElapsed - btStTime <= %TimeLimit%)),
    // Enable this time period
    TT(t2) = yes;

    // Select this branch
    BranchSelect('AC',i2,k2,h2,t2) = yes;

    // Minimize branch harmonic current, real component i -> k; tighten bound (within tolerance)
    solve mod_bound_IReBik using %SolverType% minimizing z;
    assertValidSolve(mod_bound_IReBik)
    IReB.lo(i2,k2,h2,t2) = max( IReB.lo(i2,k2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
    MB(i2,k2,h2,t2,'IRe','neg') = max( -IReB.lo(i2,k2,h2,t2) + %BigMSlack%, 0);

    // Maximize branch harmonic current, real component k -> i; tighten bound (within tolerance)
    solve mod_bound_IReBki using %SolverType% maximizing z;
    assertValidSolve(mod_bound_IReBki)
    IReB.up(k2,i2,h2,t2) = min( IReB.up(k2,i2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
    MB(k2,i2,h2,t2,'IRe','pos') = max( IReB.up(k2,i2,h2,t2) + %BigMSlack%, 0);

    // Minimize branch harmonic current, real component k -> i; tighten bound (within tolerance)
    solve mod_bound_IReBki using %SolverType% minimizing z;
    assertValidSolve(mod_bound_IReBki)
    IReB.lo(k2,i2,h2,t2) = max( IReB.lo(k2,i2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
    MB(k2,i2,h2,t2,'IRe','neg') = max( -IReB.lo(k2,i2,h2,t2) + %BigMSlack%, 0);

    // Maximize branch harmonic current, real component i -> k; tighten bound (within tolerance)
    solve mod_bound_IReBik using %SolverType% maximizing z;
    assertValidSolve(mod_bound_IReBik)
    IReB.up(i2,k2,h2,t2) = min( IReB.up(i2,k2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
    MB(i2,k2,h2,t2,'IRe','pos') = max( IReB.up(i2,k2,h2,t2) + %BigMSlack%, 0);

    // Minimize branch harmonic current, imaginary component i -> k; tighten bound (within tolerance)
    solve mod_bound_IImBik using %SolverType% minimizing z;
    assertValidSolve(mod_bound_IImBik)
    IImB.lo(i2,k2,h2,t2) = max( IImB.lo(i2,k2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
    MB(i2,k2,h2,t2,'IIm','neg') = max( -IImB.lo(i2,k2,h2,t2) + %BigMSlack%, 0);

    // Maximize branch harmonic current, imaginary component k -> i; tighten bound (within tolerance)
    solve mod_bound_IImBki using %SolverType% maximizing z;
    assertValidSolve(mod_bound_IImBki)
    IImB.up(k2,i2,h2,t2) = min( IImB.up(k2,i2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
    MB(k2,i2,h2,t2,'IIm','pos') = max( IImB.up(k2,i2,h2,t2) + %BigMSlack%, 0);

    // Minimize branch harmonic current, imaginary component k -> i; tighten bound (within tolerance)
    solve mod_bound_IImBki using %SolverType% minimizing z;
    assertValidSolve(mod_bound_IImBki)
    IImB.lo(k2,i2,h2,t2) = max( IImB.lo(k2,i2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
    MB(k2,i2,h2,t2,'IIm','neg') = max( -IImB.lo(k2,i2,h2,t2) + %BigMSlack%, 0);

    // Maximize branch harmonic current, imaginary component i -> k; tighten bound (within tolerance)
    solve mod_bound_IImBik using %SolverType% maximizing z;
    assertValidSolve(mod_bound_IImBik)
    IImB.up(i2,k2,h2,t2) = min( IImB.up(i2,k2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
    MB(i2,k2,h2,t2,'IIm','pos') = max( IImB.up(i2,k2,h2,t2) + %BigMSlack%, 0);

    // De-select this branch
    BranchSelect('AC',i2,k2,h2,t2) = no;

    // Disable this time period
    TT(t2) = no;
);

* Sources
loop ( (SX('AC',i2,s2),h2,t2)$(HX('Harm',h2) and
    (timeElapsed - btStTime <= %TimeLimit%)),
    // Enable this time period
    TT(t2) = yes;

    // Select this source
    SourceSelect('AC',i2,s2,h2,t2) = yes;

    // Minimize source harmonic current, real component; tighten bound (within tolerance)
    solve mod_bound_IReS using %SolverType% minimizing z;
    assertValidSolve(mod_bound_IReS)
    IReS.lo(s2,h2,t2) = max( IReS.lo(s2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
    MS(s2,h2,t2,'IRe','neg') = max( -IReS.lo(s2,h2,t2) + %BigMSlack%, 0);

    // Maximize source harmonic current, real component; tighten bound (within tolerance)
    solve mod_bound_IReS using %SolverType% maximizing z;
    assertValidSolve(mod_bound_IReS)
    IReS.up(s2,h2,t2) = min( IReS.up(s2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
    MS(s2,h2,t2,'IRe','pos') = max( IReS.up(s2,h2,t2) + %BigMSlack%, 0);

    // Minimize source harmonic current, imaginary component; tighten bound (within tolerance)
    solve mod_bound_IImS using %SolverType% minimizing z;
    assertValidSolve(mod_bound_IImS)
    IImS.lo(s2,h2,t2) = max( IImS.lo(s2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
    MS(s2,h2,t2,'IIm','neg') = max( -IImS.lo(s2,h2,t2) + %BigMSlack%, 0);

    // Maximize source harmonic current, imaginary component; tighten bound (within tolerance)
    solve mod_bound_IImS using %SolverType% maximizing z;
    assertValidSolve(mod_bound_IImS)
    IImS.up(s2,h2,t2) = min( IImS.up(s2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
    MS(s2,h2,t2,'IIm','pos') = max( IImS.up(s2,h2,t2) + %BigMSlack%, 0);

    // De-select this source
    SourceSelect('AC',i2,s2,h2,t2) = no;

    // Disable this time period
    TT(t2) = no;
);

* Loads
loop ( (LX('AC',i2,l2),h2,t2)$(HX('Harm',h2) and
    (timeElapsed - btStTime <= %TimeLimit%)),
    // Enable this time period
    TT(t2) = yes;

    // Select this load
    LoadSelect('AC',i2,l2,h2,t2) = yes;

    // Minimize load harmonic current, real component; tighten bound (within tolerance)
    solve mod_bound_IReL using %SolverType% minimizing z;
    assertValidSolve(mod_bound_IReL)
    IReL.lo(l2,h2,t2) = max( IReL.lo(l2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
    ML(l2,h2,t2,'IRe','neg') = max( -IReL.lo(l2,h2,t2) + %BigMSlack%, 0);

    // Maximize load harmonic current, real component; tighten bound (within tolerance)
    solve mod_bound_IReL using %SolverType% maximizing z;
    assertValidSolve(mod_bound_IReL)
    IReL.up(l2,h2,t2) = min( IReL.up(l2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
    ML(l2,h2,t2,'IRe','pos') = max( IReL.up(l2,h2,t2) + %BigMSlack%, 0);

    // Minimize load harmonic current, imaginary component; tighten bound (within tolerance)
    solve mod_bound_IImL using %SolverType% minimizing z;
    assertValidSolve(mod_bound_IImL)
    IImL.lo(l2,h2,t2) = max( IImL.lo(l2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );
    ML(l2,h2,t2,'IIm','neg') = max( -IImL.lo(l2,h2,t2) + %BigMSlack%, 0);

    // Maximize load harmonic current, imaginary component; tighten bound (within tolerance)
    solve mod_bound_IImL using %SolverType% maximizing z;
    assertValidSolve(mod_bound_IImL)
    IImL.up(l2,h2,t2) = min( IImL.up(l2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );
    ML(l2,h2,t2,'IIm','pos') = max( IImL.up(l2,h2,t2) + %BigMSlack%, 0);

    // De-select this load
    LoadSelect('AC',i2,l2,h2,t2) = no;

    // Disable this time period
    TT(t2) = no;
);

* Converters
loop ( (c2,h2,t2)$(HX('Harm',h2) and (timeElapsed - btStTime <= %TimeLimit%)),
    // Enable this time period
    TT(t2) = yes;

    // Select this converter + harmonic
    ConvSelectH(c2,h2,t2) = yes;

    // Minimize IReCout; tighten bound (within tolerance)
    solve mod_bound_IReCout using %SolverType% minimizing z;
    assertValidSolve(mod_bound_IReCout)
    IReCout.lo(c2,h2,t2) = max( IReCout.lo(c2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );

    // Minimize IReCin; tighten bound (within tolerance)
    solve mod_bound_IReCin using %SolverType% minimizing z;
    assertValidSolve(mod_bound_IReCin)
    IReCin.lo(c2,h2,t2) = max( IReCin.lo(c2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );

    // Maximize IReCout; tighten bound (within tolerance)
    solve mod_bound_IReCout using %SolverType% maximizing z;
    assertValidSolve(mod_bound_IReCout)
    IReCout.up(c2,h2,t2) = min( IReCout.up(c2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );

    // Maximize IReCin; tighten bound (within tolerance)
    solve mod_bound_IReCin using %SolverType% maximizing z;
    assertValidSolve(mod_bound_IReCin)
    IReCin.up(c2,h2,t2) = min( IReCin.up(c2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );

    // Minimize IImCout; tighten bound (within tolerance)
    solve mod_bound_IImCout using %SolverType% minimizing z;
    assertValidSolve(mod_bound_IImCout)
    IImCout.lo(c2,h2,t2) = max( IImCout.lo(c2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );

    // Minimize IImCin; tighten bound (within tolerance)
    solve mod_bound_IImCin using %SolverType% minimizing z;
    assertValidSolve(mod_bound_IImCin)
    IImCin.lo(c2,h2,t2) = max( IImCin.lo(c2,h2,t2),
        z.l * (1 - sign(z.l)*%OptCR%) - %OptCA% - %BoundSlack% );

    // Maximize IImCout; tighten bound (within tolerance)
    solve mod_bound_IImCout using %SolverType% maximizing z;
    assertValidSolve(mod_bound_IImCout)
    IImCout.up(c2,h2,t2) = min( IImCout.up(c2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );

    // Maximize IImCin; tighten bound (within tolerance)
    solve mod_bound_IImCin using %SolverType% maximizing z;
    assertValidSolve(mod_bound_IImCin)
    IImCin.up(c2,h2,t2) = min( IImCin.up(c2,h2,t2),
        z.l * (1 + sign(z.l)*%OptCR%) + %OptCA% + %BoundSlack% );

    // De-select this converter + harmonic
    ConvSelectH(c2,h2,t2) = no;

    // Disable this time period
    TT(t2) = no;
);
$ENDIF.harm


* =============================================================================
*  Tighten Binary Decision Variables
* =============================================================================
* These are not performed in the main loop because it would be time consuming
* to do so and the gains would likely be minimal. Instead, they are tightened
* after the fact. For each loop, only unfixed decisions are considered.

* Branches
loop ( (a2,i2,k2,t2)$(xB.lo(a2,i2,k2) <> xB.up(a2,i2,k2) and
    (timeElapsed - btStTime <= %TimeLimit%)),
    // Enable this time period
    TT(t2) = yes;

    // Select this branch
    xBSelect(a2,i2,k2) = yes;

    // Minimize xB; tighten bound
    solve mod_bound_xB using %SolverType% minimizing z;
    assertValidSolve(mod_bound_xB)
    xB.lo(a2,i2,k2) = ceil(z.l - %IntegerTol%);

    // Maximize xB; tighten bound
    solve mod_bound_xB using %SolverType% maximizing z;
    assertValidSolve(mod_bound_xB)
    xB.up(a2,i2,k2) = floor(z.l + %IntegerTol%);

    // De-select this branch
    xBSelect(a2,i2,k2) = no;

    // Disable this time period
    TT(t2) = no;
);

* Sources
loop ( (a2,s2,t2)$(xS.lo(a2,s2) <> xS.up(a2,s2) and
    (timeElapsed - btStTime <= %TimeLimit%)),
    // Enable this time period
    TT(t2) = yes;

    // Select this source
    xSSelect(a2,s2) = yes;

    // Minimize xS; tighten bound
    solve mod_bound_xS using %SolverType% minimizing z;
    assertValidSolve(mod_bound_xS)
    xS.lo(a2,s2) = ceil(z.l - %IntegerTol%);

    // Maximize xS; tighten bound
    solve mod_bound_xS using %SolverType% maximizing z;
    assertValidSolve(mod_bound_xS)
    xS.up(a2,s2) = floor(z.l + %IntegerTol%);

    // De-select this source
    xSSelect(a2,s2) = no;

    // Disable this time period
    TT(t2) = no;
);

* Loads
loop ( (a2,l2,t2)$(xL.lo(a2,l2) <> xL.up(a2,l2) and
    (timeElapsed - btStTime <= %TimeLimit%)),
    // Enable this time period
    TT(t2) = yes;

    // Select this load
    xLSelect(a2,l2) = yes;

    // Minimize xL; tighten bound
    solve mod_bound_xL using %SolverType% minimizing z;
    assertValidSolve(mod_bound_xL)
    xL.lo(a2,l2) = ceil(z.l - %IntegerTol%);

    // Maximize xL; tighten bound
    solve mod_bound_xL using %SolverType% maximizing z;
    assertValidSolve(mod_bound_xL)
    xL.up(a2,l2) = floor(z.l + %IntegerTol%);

    // De-select this load
    xLSelect(a2,l2) = no;

    // Disable this time period
    TT(t2) = no;
);

* Converters
loop ( (c2,t2)$(xC.lo(c2) <> xC.up(c2) and
    (timeElapsed - btStTime <= %TimeLimit%)),
    // Enable this time period
    TT(t2) = yes;

    // Select this converter
    xCSelect(c2) = yes;

    // Minimize xC; tighten bound
    solve mod_bound_xC using %SolverType% minimizing z;
    assertValidSolve(mod_bound_xC)
    xC.lo(c2) = ceil(z.l - %IntegerTol%);

    // Maximize xC; tighten bound
    solve mod_bound_xC using %SolverType% maximizing z;
    assertValidSolve(mod_bound_xC)
    xC.up(c2) = floor(z.l + %IntegerTol%);

    // De-select this converter
    xCSelect(c2) = no;

    // Disable this time period
    TT(t2) = no;
);


* =============================================================================
*  Cleanup
* =============================================================================
* Print log warning if time exceeded
if (timeElapsed - btStTime > %TimeLimit%,
    display 'Warning: Bound tightening time limit exceeded';
);

* Restore original design decisions
xB.l(a,i,k) = save_xB(a,i,k);
xS.l(a,s)   = save_xS(a,s);
xL.l(a,l)   = save_xL(a,l);
xC.l(c)     = save_xC(c);

* Correct any design decisions that have been fixed by the bound tightening
xB.l(a,i,k) = min( xB.l(a,i,k), xB.up(a,i,k) );
xB.l(a,i,k) = max( xB.l(a,i,k), xB.lo(a,i,k) );

xS.l(a,s) = min( xS.l(a,s), xS.up(a,s) );
xS.l(a,s) = max( xS.l(a,s), xS.lo(a,s) );

xL.l(a,l) = min( xL.l(a,l), xL.up(a,l) );
xL.l(a,l) = max( xL.l(a,l), xL.lo(a,l) );

xC.l(c) = min( xC.l(c), xC.up(c) );
xC.l(c) = max( xC.l(c), xC.lo(c) );

* Reset all powers to zero, or to the closest bound to zero
PB.l(i,k,h,t) = 0
    + PB.lo(i,k,h,t)$(PB.lo(i,k,h,t) > 0)
    + PB.up(i,k,h,t)$(PB.up(i,k,h,t) < 0);
QB.l(i,k,h,t) = 0
    + QB.lo(i,k,h,t)$(QB.lo(i,k,h,t) > 0)
    + QB.up(i,k,h,t)$(QB.up(i,k,h,t) < 0);

PS.l(s,h,t) = 0
    + PS.lo(s,h,t)$(PS.lo(s,h,t) > 0)
    + PS.up(s,h,t)$(PS.up(s,h,t) < 0);
QS.l(s,h,t) = 0
    + QS.lo(s,h,t)$(QS.lo(s,h,t) > 0)
    + QS.up(s,h,t)$(QS.up(s,h,t) < 0);

PL.l(l,h,t) = 0
    + PL.lo(l,h,t)$(PL.lo(l,h,t) > 0)
    + PL.up(l,h,t)$(PL.up(l,h,t) < 0);
QL.l(l,h,t) = 0
    + QL.lo(l,h,t)$(QL.lo(l,h,t) > 0)
    + QL.up(l,h,t)$(QL.up(l,h,t) < 0);

PCout.l(c,h,t) = 0
    + PCout.lo(c,h,t)$(PCout.lo(c,h,t) > 0)
    + PCout.up(c,h,t)$(PCout.up(c,h,t) < 0);
QCout.l(c,h,t) = 0
    + QCout.lo(c,h,t)$(QCout.lo(c,h,t) > 0)
    + QCout.up(c,h,t)$(QCout.up(c,h,t) < 0);

PCin.l(c,h,t) = 0
    + PCin.lo(c,h,t)$(PCin.lo(c,h,t) > 0)
    + PCin.up(c,h,t)$(PCin.up(c,h,t) < 0);
QCin.l(c,h,t) = 0
    + QCin.lo(c,h,t)$(QCin.lo(c,h,t) > 0)
    + QCin.up(c,h,t)$(QCin.up(c,h,t) < 0);

* Reset all currents to zero
$IFTHEN.harm %IncludeHarmonics% == yes
IReB.l(i,k,h,t) = 0
    + IReB.lo(i,k,h,t)$(IReB.lo(i,k,h,t) > 0)
    + IReB.up(i,k,h,t)$(IReB.up(i,k,h,t) < 0);
IImB.l(i,k,h,t) = 0
    + IImB.lo(i,k,h,t)$(IImB.lo(i,k,h,t) > 0)
    + IImB.up(i,k,h,t)$(IImB.up(i,k,h,t) < 0);

IReS.l(s,h,t) = 0
    + IReS.lo(s,h,t)$(IReS.lo(s,h,t) > 0)
    + IReS.up(s,h,t)$(IReS.up(s,h,t) < 0);
IImS.l(s,h,t) = 0
    + IImS.lo(s,h,t)$(IImS.lo(s,h,t) > 0)
    + IImS.up(s,h,t)$(IImS.up(s,h,t) < 0);

IReL.l(l,h,t) = 0
    + IReL.lo(l,h,t)$(IReL.lo(l,h,t) > 0)
    + IReL.up(l,h,t)$(IReL.up(l,h,t) < 0);
IImL.l(l,h,t) = 0
    + IImL.lo(l,h,t)$(IImL.lo(l,h,t) > 0)
    + IImL.up(l,h,t)$(IImL.up(l,h,t) < 0);

IReCout.l(c,h,t) = 0
    + IReCout.lo(c,h,t)$(IReCout.lo(c,h,t) > 0)
    + IReCout.up(c,h,t)$(IReCout.up(c,h,t) < 0);
IImCout.l(c,h,t) = 0
    + IImCout.lo(c,h,t)$(IImCout.lo(c,h,t) > 0)
    + IImCout.up(c,h,t)$(IImCout.up(c,h,t) < 0);

IReCin.l(c,h,t) = 0
    + IReCin.lo(c,h,t)$(IReCin.lo(c,h,t) > 0)
    + IReCin.up(c,h,t)$(IReCin.up(c,h,t) < 0);
IImCin.l(c,h,t) = 0
    + IImCin.lo(c,h,t)$(IImCin.lo(c,h,t) > 0)
    + IImCin.up(c,h,t)$(IImCin.up(c,h,t) < 0);
$ENDIF.harm

* Reset voltages to a flat start state
$BATINCLUDE %RelPath%/AC-DC-bounds-breaks.gms vlevels

