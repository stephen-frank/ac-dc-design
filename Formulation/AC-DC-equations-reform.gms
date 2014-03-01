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

$TITLE Reformulated equations for AC-DC model
$ONTEXT

    Filename : AC-DC-equations-reform.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    This file implements the reformulation for the mixed AC-DC distribution
    system design problem. It is designed to be included (via $INCLUDE) in other
    .gms files, and relies on the following external compile time flags
    (which are documented in AC-DC-defaults.gms):
        BilinearEq          Generate bilinear equations for u = y f(x)?
        BigM                Generate big M equations for u = y f(x)?
        ExtraCuts           Generate valid feasiblity cuts?
        ExactReform         Generate exact definitions for reformulation vars?
        IncludeHarmonics    Include equations for harmonic currents?

    NOTES:
    1. Some of the equations defined here are identical to those in
       'AC-DC-model.gms' but they have been included here as well for
       completeness (and to prevent needing to mix and match the two files).
    2. Setting BilinearEq == yes generates functions for the powers of
       the form u = y f(x). These are nonconvex, but become convex in x if y
       is fixed.
    3. Setting BigM == yes generates reformulations for functions of the form
       u = y f(x) using big M's. This isn't needed if fixing the binaries (y's),
       and if used in such a case just creates extra constraints.
    4. Setting ExtraCuts == yes generates valid inequalities for the voltages,
       branch powers, etc. These aren't needed in the exact representation
       because they are redundant, but are useful for the relaxations.
    5. Setting ExactReform == yes generates w12 = x1 * x2 type reformulation
       constraints for the continuous-continuous bilinear terms. These are
       exact, but nonconvex.
    7. None of the yes/no flags are mutually exclusive; you can use any or all
       of them and specify the constraint names in the model later to get the
       subset you want (as opposed to using /all/). Enabling everything is the
       default mode.

$OFFTEXT

* ==============================================================================
*  Equations
* ==============================================================================
* ------------------------------------------------------------------------------
*  Equation Declaration
* ------------------------------------------------------------------------------
* Main constraints
Equations
    // Objective Functions
    eq_obj                  'Objective function'
    eq_dummy                'Placeholder objective for feasibility test'
    // Power Balance
    eq_PBal(a,i,h,t)        'Real power balance'
    eq_QBal(a,i,h,t)        'Reactive power balance'
    // Current Balance
$IFTHEN.harm %IncludeHarmonics% == yes
    eq_IReBal(a,i,h,t)      'Harmonic current balance, real component'
    eq_IImBal(a,i,h,t)      'Harmonic current balance, imaginary component'
$ENDIF.harm
    // Branch Power Flow
$IFTHEN %BilinearEq% == yes
    eq_Pik(a,i,k,h,t)       'Branch ik (from -> to) real power flow'
    eq_Pki(a,i,k,h,t)       'Branch ki (to -> from) real power flow'
    eq_Qik(a,i,k,h,t)       'Branch ik (from -> to) reactive power flow'
    eq_Qki(a,i,k,h,t)       'Branch ki (to -> from) reactive power flow'
$ENDIF
$IFTHEN %BigM% == yes
    eq_Pik_1(a,i,k,h,t)     'Branch ik (from -> to) real power flow - 1 of 4'
    eq_Pik_2(a,i,k,h,t)     'Branch ik (from -> to) real power flow - 2 of 4'
    eq_Pik_3(a,i,k,h,t)     'Branch ik (from -> to) real power flow - 3 of 4'
    eq_Pik_4(a,i,k,h,t)     'Branch ik (from -> to) real power flow - 4 of 4'
    eq_Pki_1(a,i,k,h,t)     'Branch ki (to -> from) real power flow - 1 of 4'
    eq_Pki_2(a,i,k,h,t)     'Branch ki (to -> from) real power flow - 2 of 4'
    eq_Pki_3(a,i,k,h,t)     'Branch ki (to -> from) real power flow - 3 of 4'
    eq_Pki_4(a,i,k,h,t)     'Branch ki (to -> from) real power flow - 4 of 4'
    eq_Qik_1(a,i,k,h,t)     'Branch ik (from -> to) reactive power flow - 1 of 4'
    eq_Qik_2(a,i,k,h,t)     'Branch ik (from -> to) reactive power flow - 2 of 4'
    eq_Qik_3(a,i,k,h,t)     'Branch ik (from -> to) reactive power flow - 3 of 4'
    eq_Qik_4(a,i,k,h,t)     'Branch ik (from -> to) reactive power flow - 4 of 4'
    eq_Qki_1(a,i,k,h,t)     'Branch ki (to -> from) reactive power flow - 1 of 4'
    eq_Qki_2(a,i,k,h,t)     'Branch ki (to -> from) reactive power flow - 2 of 4'
    eq_Qki_3(a,i,k,h,t)     'Branch ki (to -> from) reactive power flow - 3 of 4'
    eq_Qki_4(a,i,k,h,t)     'Branch ki (to -> from) reactive power flow - 4 of 4'
$ENDIF
$IFTHEN %ExtraCuts% == yes
    eq_PBLoss_L(a,i,k,h,t)  'Lower bound on branch ik real power loss'
    eq_QBLoss_L(a,i,k,h,t)  'Lower bound on branch ik reactive power loss'
    eq_QBLoss_U(a,i,k,h,t)  'Upper bound on branch ik reactive power loss'
$ENDIF
    // Branch Current
$IFTHEN.harm %IncludeHarmonics% == yes
$IFTHEN %BilinearEq% == yes
    eq_IReik(a,i,k,h,t)     'Branch ik (from -> to) harmonic current, real component'
    eq_IReki(a,i,k,h,t)     'Branch ki (to -> from) harmonic current, real component'
    eq_IImik(a,i,k,h,t)     'Branch ik (from -> to) harmonic current, imaginary component'
    eq_IImki(a,i,k,h,t)     'Branch ki (to -> from) harmonic current, imaginary component'
$ENDIF
$IFTHEN %BigM% == yes
    eq_IReik_1(a,i,k,h,t)   'Branch ik (from -> to) harmonic current, real component - 1 of 4'
    eq_IReik_2(a,i,k,h,t)   'Branch ik (from -> to) harmonic current, real component - 2 of 4'
    eq_IReik_3(a,i,k,h,t)   'Branch ik (from -> to) harmonic current, real component - 3 of 4'
    eq_IReik_4(a,i,k,h,t)   'Branch ik (from -> to) harmonic current, real component - 4 of 4'
    eq_IReki_1(a,i,k,h,t)   'Branch ki (to -> from) harmonic current, real component - 1 of 4'
    eq_IReki_2(a,i,k,h,t)   'Branch ki (to -> from) harmonic current, real component - 2 of 4'
    eq_IReki_3(a,i,k,h,t)   'Branch ki (to -> from) harmonic current, real component - 3 of 4'
    eq_IReki_4(a,i,k,h,t)   'Branch ki (to -> from) harmonic current, real component - 4 of 4'
    eq_IImik_1(a,i,k,h,t)   'Branch ik (from -> to) harmonic current, imaginary component - 1 of 4'
    eq_IImik_2(a,i,k,h,t)   'Branch ik (from -> to) harmonic current, imaginary component - 2 of 4'
    eq_IImik_3(a,i,k,h,t)   'Branch ik (from -> to) harmonic current, imaginary component - 3 of 4'
    eq_IImik_4(a,i,k,h,t)   'Branch ik (from -> to) harmonic current, imaginary component - 4 of 4'
    eq_IImki_1(a,i,k,h,t)   'Branch ki (to -> from) harmonic current, imaginary component - 1 of 4'
    eq_IImki_2(a,i,k,h,t)   'Branch ki (to -> from) harmonic current, imaginary component - 2 of 4'
    eq_IImki_3(a,i,k,h,t)   'Branch ki (to -> from) harmonic current, imaginary component - 3 of 4'
    eq_IImki_4(a,i,k,h,t)   'Branch ki (to -> from) harmonic current, imaginary component - 4 of 4'
$ENDIF
$ENDIF.harm
    // Source Power
$IFTHEN %BilinearEq% == yes
    eq_PSrc(a,i,s,h,t)      'Source real power'
    eq_QSrc(a,i,s,h,t)      'Source reactive power'
$ENDIF
$IFTHEN %BigM% == yes
    eq_PSrc_1(a,i,s,h,t)    'Source real power - 1 of 4'
    eq_PSrc_2(a,i,s,h,t)    'Source real power - 2 of 4'
    eq_PSrc_3(a,i,s,h,t)    'Source real power - 3 of 4'
    eq_PSrc_4(a,i,s,h,t)    'Source real power - 4 of 4'
    eq_QSrc_1(a,i,s,h,t)    'Source reactive power - 1 of 4'
    eq_QSrc_2(a,i,s,h,t)    'Source reactive power - 2 of 4'
    eq_QSrc_3(a,i,s,h,t)    'Source reactive power - 3 of 4'
    eq_QSrc_4(a,i,s,h,t)    'Source reactive power - 4 of 4'
$ENDIF
    // Source Current
$IFTHEN.harm %IncludeHarmonics% == yes
$IFTHEN %BilinearEq% == yes
    eq_IReSrc(a,i,s,h,t)    'Source harmonic current, real component'
    eq_IImSrc(a,i,s,h,t)    'Source harmonic current, imaginary component'
$ENDIF
$IFTHEN %BigM% == yes
    eq_IReSrc_1(a,i,s,h,t)  'Source harmonic current, real component - 1 of 4'
    eq_IReSrc_2(a,i,s,h,t)  'Source harmonic current, real component - 2 of 4'
    eq_IReSrc_3(a,i,s,h,t)  'Source harmonic current, real component - 3 of 4'
    eq_IReSrc_4(a,i,s,h,t)  'Source harmonic current, real component - 4 of 4'
    eq_IImSrc_1(a,i,s,h,t)  'Source harmonic current, imaginary component - 1 of 4'
    eq_IImSrc_2(a,i,s,h,t)  'Source harmonic current, imaginary component - 2 of 4'
    eq_IImSrc_3(a,i,s,h,t)  'Source harmonic current, imaginary component - 3 of 4'
    eq_IImSrc_4(a,i,s,h,t)  'Source harmonic current, imaginary component - 4 of 4'
$ENDIF
$ENDIF.harm
    // Load Power
$IFTHEN %BilinearEq% == yes
    eq_PLoad(a,i,l,h,t)     'Load real power'
    eq_QLoad(a,i,l,h,t)     'Load reactive power'
$ENDIF
$IFTHEN %BigM% == yes
    eq_PLoad_1(a,i,l,h,t)   'Load real power - 1 of 4'
    eq_PLoad_2(a,i,l,h,t)   'Load real power - 2 of 4'
    eq_PLoad_3(a,i,l,h,t)   'Load real power - 3 of 4'
    eq_PLoad_4(a,i,l,h,t)   'Load real power - 4 of 4'
    eq_QLoad_1(a,i,l,h,t)   'Load reactive power - 1 of 4'
    eq_QLoad_2(a,i,l,h,t)   'Load reactive power - 2 of 4'
    eq_QLoad_3(a,i,l,h,t)   'Load reactive power - 3 of 4'
    eq_QLoad_4(a,i,l,h,t)   'Load reactive power - 4 of 4'
$ENDIF
    // Load Current
$IFTHEN.harm %IncludeHarmonics% == yes
$IFTHEN %BilinearEq% == yes
    eq_IReLoad(a,i,l,h,t)   'Load harmonic current, real component'
    eq_IImLoad(a,i,l,h,t)   'Load harmonic current, imaginary component'
$ENDIF
$IFTHEN %BigM% == yes
    eq_IReLoad_1(a,i,l,h,t) 'Load harmonic current, real component - 1 of 4'
    eq_IReLoad_2(a,i,l,h,t) 'Load harmonic current, real component - 2 of 4'
    eq_IReLoad_3(a,i,l,h,t) 'Load harmonic current, real component - 3 of 4'
    eq_IReLoad_4(a,i,l,h,t) 'Load harmonic current, real component - 4 of 4'
    eq_IImLoad_1(a,i,l,h,t) 'Load harmonic current, imaginary component - 1 of 4'
    eq_IImLoad_2(a,i,l,h,t) 'Load harmonic current, imaginary component - 2 of 4'
    eq_IImLoad_3(a,i,l,h,t) 'Load harmonic current, imaginary component - 3 of 4'
    eq_IImLoad_4(a,i,l,h,t) 'Load harmonic current, imaginary component - 4 of 4'
$ENDIF
$ENDIF.harm
    // Converter Power
    eq_PConvOutMin(c,h,t)   'Logic switch for min. converter real power output'
    eq_PConvOutMax(c,h,t)   'Logic switch for max. converter real power output'
    eq_PConvInMin(c,h,t)    'Logic switch for min. converter real power input'
    eq_PConvInMax(c,h,t)    'Logic switch for max. converter real power input'
    eq_PConvInOut(c,t)      'Converter input power as function of output power'
    eq_QConvOutMin(c,h,t)   'Logic switch for min. converter reactive power output'
    eq_QConvOutMax(c,h,t)   'Logic switch for max. converter reactive power output'
    eq_QConvInMin(c,h,t)    'Logic switch for min. converter reactive power input'
    eq_QConvInMax(c,h,t)    'Logic switch for max. converter reactive power input'
    // Converter Current
$IFTHEN.harm %IncludeHarmonics% == yes
    eq_IReConvOutMin(c,h,t) 'Logic switch for min. converter output harmonic current, real component'
    eq_IReConvOutMax(c,h,t) 'Logic switch for max. converter output harmonic current, real component'
    eq_IReConvInMin(c,h,t)  'Logic switch for min. converter input harmonic current, real component'
    eq_IReConvInMax(c,h,t)  'Logic switch for max. converter input harmonic current, real component'
    eq_IImConvOutMin(c,h,t) 'Logic switch for min. converter output harmonic current, imaginary component'
    eq_IImConvOutMax(c,h,t) 'Logic switch for max. converter output harmonic current, imaginary component'
    eq_IImConvInMin(c,h,t)  'Logic switch for min. converter input harmonic current, imaginary component'
    eq_IImConvInMax(c,h,t)  'Logic switch for max. converter input harmonic current, imaginary component'
$ENDIF.harm
    // State Variable Relationships
    eq_VoltMag(i,h,t)       'Computation of voltage magnitude'
    eq_wC(c,t)              'Definition of reformulation variable wC'
$IFTHEN.harm %IncludeHarmonics% == yes
    eq_EHarmMin(i,c,h,t)    'Logic switch for min. real voltage at converter terminal'
    eq_EHarmMax(i,c,h,t)    'Logic switch for max. real voltage at converter terminal'
    eq_FHarmMin(i,c,h,t)    'Logic switch for min. imaginary voltage at converter terminal'
    eq_FHarmMax(i,c,h,t)    'Logic switch for max. imaginary voltage at converter terminal'
$ENDIF.harm
$IFTHEN %ExtraCuts% == yes
    eq_VEF_1(i,h,t)         'Voltage magnitude geometry valid inequality - 1 of 3'
    eq_VEF_2(i,h,t)         'Voltage magnitude geometry valid inequality - 2 of 3'
    eq_VEF_3(i,h,t)         'Voltage magnitude geometry valid inequality - 3 of 3'
$ENDIF
    // Binary Logic
    eq_BranchAsgn(i,k)      'Assignment of branches to AC or DC network'
    eq_SrcAsgn(s)           'Assignment of sources to AC or DC network'
    eq_LoadAsgn(l)          'Assignment of loads to AC or DC network'
    eq_BranchLoad(a,i,k,l)  'Maps branches that are also loads'
    eq_SrcConn(a,i,s)       'Ensure feasible connection of sources'
    eq_LoadConn(a,i,l)      'Ensure feasible connection of loads'
    eq_ConvIn(a,i,c)        'Ensure feasible connection of converter inputs'
    eq_ConvOut(a,i,c)       'Exclude suboptimal connection of converter outputs' ;

* Reformulation Constraints
$IFTHEN %ExactReform% == yes
Equations
    eq_wVV(a,i,k,h,t)       'Exact definition of reformulation variable wVV'
    eq_wEE(a,i,k,h,t)       'Exact definition of reformulation variable wEE'
    eq_wEF(a,i,k,h,t)       'Exact definition of reformulation variable wEF'
    eq_wFF(a,i,k,h,t)       'Exact definition of reformulation variable wFF'
    eq_wCC(c,t)             'Exact definition of reformulation variable wCC' ;
$ENDIF

* ------------------------------------------------------------------------------
*  Equation Definitions
* ------------------------------------------------------------------------------
** Objective Functions
* Objective = minimize utility energy at fundamental AC (h=1)
eq_obj..
    z =e= sum((HH('AC',h),TT(t)), Tau(t) * PS('utility',h,t) ) ;

* Dummy objective for feasibility tests
eq_dummy..
    z =e= 0 ;


** Power Balance
* Real power balance
eq_PBal(NN(a,i),h,TT(t))$HH(a,h)..
    sum(SX(a,i,s), PS(s,h,t)) - sum(LX(a,i,l), PL(l,h,t))
        + sum(CXo(a,i,c), PCout(c,h,t)) - sum(CXi(a,i,c), PCin(c,h,t))
        =e=
        sum(BX(a,i,k), PB(i,k,h,t)) + sum(BX(a,k,i), PB(i,k,h,t)) ;

* AC reactive power balance
eq_QBal(NN('AC',i),h,TT(t))$HH('AC',h)..
    sum(SX('AC',i,s), QS(s,h,t)) - sum(LX('AC',i,l), QL(l,h,t))
        + sum(CXo('AC',i,c), QCout(c,h,t)) - sum(CXi('AC',i,c), QCin(c,h,t))
        =e=
        sum(BX('AC',i,k), QB(i,k,h,t)) + sum(BX('AC',k,i), QB(i,k,h,t)) ;


** Current Balance
$IFTHEN.harm %IncludeHarmonics% == yes
* Harmonic current balance, real component
eq_IReBal(NN('AC',i),h,TT(t))$HX('Harm',h)..
    sum(SX('AC',i,s), IReS(s,h,t)) - sum(LX('AC',i,l), IReL(l,h,t))
        + sum(CXo('AC',i,c), IReCout(c,h,t)) - sum(CXi('AC',i,c), IReCin(c,h,t))
        =e=
        sum(BX('AC',i,k), IReB(i,k,h,t)) + sum(BX('AC',k,i), IReB(i,k,h,t)) ;

* Harmonic current balance, imaginary component
eq_IImBal(NN('AC',i),h,TT(t))$HX('Harm',h)..
    sum(SX('AC',i,s), IImS(s,h,t)) - sum(LX('AC',i,l), IImL(l,h,t))
        + sum(CXo('AC',i,c), IImCout(c,h,t)) - sum(CXi('AC',i,c), IImCin(c,h,t))
        =e=
        sum(BX('AC',i,k), IImB(i,k,h,t)) + sum(BX('AC',k,i), IImB(i,k,h,t)) ;
$ENDIF.harm


** Branch Power Flow
* Branch real power flow (from -> to)
$IFTHEN %BilinearEq% == yes
eq_Pik(BX(a,i,k),h,TT(t))$HH(a,h)..
    PB(i,k,h,t) =e=
    // DC power flow
    ( xB(a,i,k) * (
        (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * wVV(i,i,h,t)
        - (1/BranchData(i,k,h,t,'A'))
            * BranchData(i,k,h,t,'gSe')
            * wVV(i,k,h,t)
    ))$sameas(a,'DC')
    // AC power flow
    + ( xB(a,i,k) * (
        (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * wVV(i,i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( wEE(i,k,h,t) + wFF(i,k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * (   BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( wEF(i,k,h,t) - wEF(k,i,h,t) )
    ))$sameas(a,'AC') ;
$ENDIF
$IFTHEN %BigM% == yes
eq_Pik_1(BX(a,i,k),h,TT(t))$HH(a,h)..
    PB(i,k,h,t) =l=
    // DC power flow
    (
        (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * wVV(i,i,h,t)
        - (1/BranchData(i,k,h,t,'A'))
            * BranchData(i,k,h,t,'gSe')
            * wVV(i,k,h,t)
    )$sameas(a,'DC')
    // AC power flow
    + (
        (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * wVV(i,i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( wEE(i,k,h,t) + wFF(i,k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * (   BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( wEF(i,k,h,t) - wEF(k,i,h,t) )
    )$sameas(a,'AC')
    // Big M
    + MB(i,k,h,t,'P','neg') * (1 - xB(a,i,k)) ;

eq_Pik_2(BX(a,i,k),h,TT(t))$HH(a,h)..
    PB(i,k,h,t) =g=
    // DC power flow
    (
        (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * wVV(i,i,h,t)
        - (1/BranchData(i,k,h,t,'A'))
            * BranchData(i,k,h,t,'gSe')
            * wVV(i,k,h,t)
    )$sameas(a,'DC')
    // AC power flow
    + (
        (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * wVV(i,i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( wEE(i,k,h,t) + wFF(i,k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * (   BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( wEF(i,k,h,t) - wEF(k,i,h,t) )
    )$sameas(a,'AC')
    // Big M
    - MB(i,k,h,t,'P','pos') * (1 - xB(a,i,k)) ;

eq_Pik_3(BX(a,i,k),h,TT(t))$HH(a,h)..
    PB(i,k,h,t) =l=  MB(i,k,h,t,'P','pos') * xB(a,i,k) ;

eq_Pik_4(BX(a,i,k),h,TT(t))$HH(a,h)..
    PB(i,k,h,t) =g= -MB(i,k,h,t,'P','neg') * xB(a,i,k) ;
$ENDIF

* Branch real power flow (to -> from)
$IFTHEN %BilinearEq% == yes
eq_Pki(BX(a,i,k),h,TT(t))$HH(a,h)..
    PB(k,i,h,t) =e=
    // DC power flow
    ( xB(a,i,k) * (
        ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * wVV(k,k,h,t)
        - (1/BranchData(i,k,h,t,'A'))
            * BranchData(i,k,h,t,'gSe')
            * wVV(i,k,h,t)
    ))$sameas(a,'DC')
    // AC power flow
    + ( xB(a,i,k) * (
        ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * wVV(k,k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( wEE(i,k,h,t) + wFF(i,k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * (   BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( wEF(i,k,h,t) - wEF(k,i,h,t) )
    ))$sameas(a,'AC') ;
$ENDIF
$IFTHEN %BigM% == yes
eq_Pki_1(BX(a,i,k),h,TT(t))$HH(a,h)..
    PB(k,i,h,t) =l=
    // DC power flow
    (
        ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * wVV(k,k,h,t)
        - (1/BranchData(i,k,h,t,'A'))
            * BranchData(i,k,h,t,'gSe')
            * wVV(i,k,h,t)
    )$sameas(a,'DC')
    // AC power flow
    + (
        ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * wVV(k,k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( wEE(i,k,h,t) + wFF(i,k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * (   BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( wEF(i,k,h,t) - wEF(k,i,h,t) )
    )$sameas(a,'AC')
    // Big M
    + MB(k,i,h,t,'P','neg') * (1 - xB(a,i,k)) ;

eq_Pki_2(BX(a,i,k),h,TT(t))$HH(a,h)..
    PB(k,i,h,t) =g=
    // DC power flow
    (
        ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * wVV(k,k,h,t)
        - (1/BranchData(i,k,h,t,'A'))
            * BranchData(i,k,h,t,'gSe')
            * wVV(i,k,h,t)
    )$sameas(a,'DC')
    // AC power flow
    + (
        ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * wVV(k,k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( wEE(i,k,h,t) + wFF(i,k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * (   BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( wEF(i,k,h,t) - wEF(k,i,h,t) )
    )$sameas(a,'AC')
    // Big M
    - MB(k,i,h,t,'P','pos') * (1 - xB(a,i,k)) ;

eq_Pki_3(BX(a,i,k),h,TT(t))$HH(a,h)..
    PB(k,i,h,t) =l=  MB(k,i,h,t,'P','pos') * xB(a,i,k) ;

eq_Pki_4(BX(a,i,k),h,TT(t))$HH(a,h)..
    PB(k,i,h,t) =g= -MB(k,i,h,t,'P','neg') * xB(a,i,k) ;
$ENDIF

* Branch reactive power flow (from -> to)
$IFTHEN %BilinearEq% == yes
eq_Qik(BX('AC',i,k),h,TT(t))$HH('AC',h)..
    QB(i,k,h,t) =e=
    // AC power flow
    xB('AC',i,k) * (
        - (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * wVV(i,i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( + BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( wEE(i,k,h,t) + wFF(i,k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * (   BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( wEF(i,k,h,t) - wEF(k,i,h,t) )
    ) ;
$ENDIF
$IFTHEN %BigM% == yes
eq_Qik_1(BX('AC',i,k),h,TT(t))$HH('AC',h)..
    QB(i,k,h,t) =l=
    // AC power flow
    (
        - (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * wVV(i,i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( + BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( wEE(i,k,h,t) + wFF(i,k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * (   BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( wEF(i,k,h,t) - wEF(k,i,h,t) )
    )
    // Big M
    + MB(i,k,h,t,'Q','neg') * (1 - xB('AC',i,k)) ;

eq_Qik_2(BX('AC',i,k),h,TT(t))$HH('AC',h)..
    QB(i,k,h,t) =g=
    // AC power flow
    (
        - (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * wVV(i,i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( + BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( wEE(i,k,h,t) + wFF(i,k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * (   BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( wEF(i,k,h,t) - wEF(k,i,h,t) )
    )
    // Big M
    - MB(i,k,h,t,'Q','pos') * (1 - xB('AC',i,k)) ;

eq_Qik_3(BX('AC',i,k),h,TT(t))$HH('AC',h)..
    QB(i,k,h,t) =l=  MB(i,k,h,t,'Q','pos') * xB('AC',i,k) ;

eq_Qik_4(BX('AC',i,k),h,TT(t))$HH('AC',h)..
    QB(i,k,h,t) =g= -MB(i,k,h,t,'Q','neg') * xB('AC',i,k) ;
$ENDIF

* Branch reactive power flow (to -> from)
$IFTHEN %BilinearEq% == yes
eq_Qki(BX('AC',i,k),h,TT(t))$HH('AC',h)..
    QB(k,i,h,t) =e=
    // AC power flow
    xB('AC',i,k) * (
        - ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * wVV(k,k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( wEE(i,k,h,t) + wFF(i,k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( wEF(i,k,h,t) - wEF(k,i,h,t) )
    ) ;
$ENDIF
$IFTHEN %BigM% == yes
eq_Qki_1(BX('AC',i,k),h,TT(t))$HH('AC',h)..
    QB(k,i,h,t) =l=
    // AC power flow
    (
        - ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * wVV(k,k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( wEE(i,k,h,t) + wFF(i,k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( wEF(i,k,h,t) - wEF(k,i,h,t) )
    )
    // Big M
    + MB(k,i,h,t,'Q','neg') * (1 - xB('AC',i,k)) ;

eq_Qki_2(BX('AC',i,k),h,TT(t))$HH('AC',h)..
    QB(k,i,h,t) =g=
    // AC power flow
    (
        - ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * wVV(k,k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( wEE(i,k,h,t) + wFF(i,k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( wEF(i,k,h,t) - wEF(k,i,h,t) )
    )
    // Big M
    - MB(k,i,h,t,'Q','pos') * (1 - xB('AC',i,k)) ;

eq_Qki_3(BX('AC',i,k),h,TT(t))$HH('AC',h)..
    QB(k,i,h,t) =l=  MB(k,i,h,t,'Q','pos') * xB('AC',i,k) ;

eq_Qki_4(BX('AC',i,k),h,TT(t))$HH('AC',h)..
    QB(k,i,h,t) =g= -MB(k,i,h,t,'Q','neg') * xB('AC',i,k) ;
$ENDIF

$IFTHEN %ExtraCuts% == yes
* Lower bound on real power loss
eq_PBLoss_L(BX(a,i,k),h,TT(t))$(HH(a,h) and
    BranchData(i,k,h,t,'gSe') >= 0 and BranchData(i,k,h,t,'gSh') >= 0 )..
    PB(i,k,h,t) + PB(k,i,h,t) =g= 0 ;

* Lower bound on reactive power loss
eq_QBLoss_L(BX('AC',i,k),h,TT(t))$(HH('AC',h) and
    BranchData(i,k,h,t,'bSe') <= 0 and BranchData(i,k,h,t,'bSh') <= 0)..
    QB(i,k,h,t) + QB(k,i,h,t) =g= 0 ;

* Upper bound on reactive power loss
eq_QBLoss_U(BX('AC',i,k),h,TT(t))$(HH('AC',h) and
    BranchData(i,k,h,t,'bSe') >= 0 and BranchData(i,k,h,t,'bSh') >= 0)..
    QB(i,k,h,t) + QB(k,i,h,t) =l= 0 ;
$ENDIF


** Branch Current
$IFTHEN.harm %IncludeHarmonics% == yes
* Branch harmonic current, real component (from -> to)
$IFTHEN %BilinearEq% == yes
eq_IReik(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IReB(i,k,h,t) =e=
    xB('AC',i,k) * (
        (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * E(i,h,t)
        - (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * F(i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * E(k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( + BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * F(k,h,t)
    ) ;
$ENDIF
$IFTHEN %BigM% == yes
eq_IReik_1(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IReB(i,k,h,t) =l=
    // Current
    (
        (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * E(i,h,t)
        - (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * F(i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * E(k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( + BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * F(k,h,t)
    )
    // Big M
    + MB(i,k,h,t,'IRe','neg') * (1 - xB('AC',i,k)) ;

eq_IReik_2(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IReB(i,k,h,t) =g=
    // Current
    (
        (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * E(i,h,t)
        - (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * F(i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * E(k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( + BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * F(k,h,t)
    )
    // Big M
    - MB(i,k,h,t,'IRe','pos') * (1 - xB('AC',i,k)) ;

eq_IReik_3(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IReB(i,k,h,t) =l=  MB(i,k,h,t,'IRe','pos') * xB('AC',i,k) ;

eq_IReik_4(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IReB(i,k,h,t) =g= -MB(i,k,h,t,'IRe','neg') * xB('AC',i,k) ;
$ENDIF

* Branch harmonic current, real component (to -> from)
$IFTHEN %BilinearEq% == yes
eq_IReki(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IReB(k,i,h,t) =e=
    xB('AC',i,k) * (
        ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * E(k,h,t)
        - ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * F(k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * E(i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * F(i,h,t)
    ) ;
$ENDIF
$IFTHEN %BigM% == yes
eq_IReki_1(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IReB(k,i,h,t) =l=
    // Current
    (
        ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * E(k,h,t)
        - ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * F(k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * E(i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * F(i,h,t)
    )
    // Big M
    + MB(k,i,h,t,'IRe','neg') * (1 - xB('AC',i,k)) ;

eq_IReki_2(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IReB(k,i,h,t) =g=
    // Current
    (
        ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * E(k,h,t)
        - ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * F(k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * E(i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * F(i,h,t)
    )
    // Big M
    - MB(k,i,h,t,'IRe','pos') * (1 - xB('AC',i,k)) ;

eq_IReki_3(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IReB(k,i,h,t) =l=  MB(k,i,h,t,'IRe','pos') * xB('AC',i,k) ;

eq_IReki_4(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IReB(k,i,h,t) =g= -MB(k,i,h,t,'IRe','neg') * xB('AC',i,k) ;
$ENDIF

* Branch harmonic current, imaginary component (from -> to)
$IFTHEN %BilinearEq% == yes
eq_IImik(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IImB(i,k,h,t) =e=
    xB('AC',i,k) * (
        (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * E(i,h,t)
        + (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * F(i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * E(k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * F(k,h,t)
    ) ;
$ENDIF
$IFTHEN %BigM% == yes
eq_IImik_1(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IImB(i,k,h,t) =l=
    // Current
    (
        (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * E(i,h,t)
        + (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * F(i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * E(k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * F(k,h,t)
    )
    // Big M
    + MB(i,k,h,t,'IIm','neg') * (1 - xB('AC',i,k)) ;

eq_IImik_2(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IImB(i,k,h,t) =g=
    // Current
    (
        (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * E(i,h,t)
        + (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * F(i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * E(k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * F(k,h,t)
    )
    // Big M
    - MB(i,k,h,t,'IIm','pos') * (1 - xB('AC',i,k)) ;

eq_IImik_3(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IImB(i,k,h,t) =l=  MB(i,k,h,t,'IIm','pos') * xB('AC',i,k) ;

eq_IImik_4(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IImB(i,k,h,t) =g= -MB(i,k,h,t,'IIm','neg') * xB('AC',i,k) ;
$ENDIF

* Branch harmonic current, imaginary component (to -> from)
$IFTHEN %BilinearEq% == yes
eq_IImki(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IImB(k,i,h,t) =e=
    xB('AC',i,k) * (
        ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * E(k,h,t)
        + ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * F(k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( + BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * E(i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * F(i,h,t)
    ) ;
$ENDIF
$IFTHEN %BigM% == yes
eq_IImki_1(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IImB(k,i,h,t) =l=
    // Current
    (
        ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * E(k,h,t)
        + ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * F(k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( + BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * E(i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * F(i,h,t)
    )
    // Big M
    + MB(k,i,h,t,'IIm','neg') * (1 - xB('AC',i,k)) ;

eq_IImki_2(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IImB(k,i,h,t) =g=
    // Current
    (
        ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * E(k,h,t)
        + ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * F(k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( + BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * E(i,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * F(i,h,t)
    )
    // Big M
    - MB(k,i,h,t,'IIm','pos') * (1 - xB('AC',i,k)) ;

eq_IImki_3(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IImB(k,i,h,t) =l=  MB(k,i,h,t,'IIm','pos') * xB('AC',i,k) ;

eq_IImki_4(BX('AC',i,k),h,TT(t))$HX('Harm',h)..
    IImB(k,i,h,t) =g= -MB(k,i,h,t,'IIm','neg') * xB('AC',i,k) ;
$ENDIF
$ENDIF.harm


** Source Power
* Source real power
$IFTHEN %BilinearEq% == yes
eq_PSrc(SX(a,i,s),h,TT(t))$HH(a,h)..
    PS(s,h,t) =e=
    // Source
    xS(a,s) * (
        SourceData(s,h,t,'P')
        + V(i,h,t) * SourceData(s,h,t,'ID')
        + ( E(i,h,t) * SourceData(s,h,t,'IRe') )$sameas(a,'AC')
        + ( F(i,h,t) * SourceData(s,h,t,'IIm') )$sameas(a,'AC')
        - wVV(i,i,h,t) * SourceData(s,h,t,'g')
    ) ;
$ENDIF
$IFTHEN %BigM% == yes
eq_PSrc_1(SX(a,i,s),h,TT(t))$HH(a,h)..
    PS(s,h,t) =l=
    // Source
    (
        SourceData(s,h,t,'P')
        + V(i,h,t) * SourceData(s,h,t,'ID')
        + ( E(i,h,t) * SourceData(s,h,t,'IRe') )$sameas(a,'AC')
        + ( F(i,h,t) * SourceData(s,h,t,'IIm') )$sameas(a,'AC')
        - wVV(i,i,h,t) * SourceData(s,h,t,'g')
    )
    // Big M
    + MS(s,h,t,'P','neg') * (1 - xS(a,s)) ;

eq_PSrc_2(SX(a,i,s),h,TT(t))$HH(a,h)..
    PS(s,h,t) =g=
    // Source
    (
        SourceData(s,h,t,'P')
        + V(i,h,t) * SourceData(s,h,t,'ID')
        + ( E(i,h,t) * SourceData(s,h,t,'IRe') )$sameas(a,'AC')
        + ( F(i,h,t) * SourceData(s,h,t,'IIm') )$sameas(a,'AC')
        - wVV(i,i,h,t) * SourceData(s,h,t,'g')
    )
    // Big M
    - MS(s,h,t,'P','pos') * (1 - xS(a,s)) ;

eq_PSrc_3(SX(a,i,s),h,TT(t))$HH(a,h)..
    PS(s,h,t) =l=  MS(s,h,t,'P','pos') * xS(a,s) ;

eq_PSrc_4(SX(a,i,s),h,TT(t))$HH(a,h)..
    PS(s,h,t) =g= -MS(s,h,t,'P','neg') * xS(a,s) ;
$ENDIF

* Source reactive power
$IFTHEN %BilinearEq% == yes
eq_QSrc(SX('AC',i,s),h,TT(t))$HH('AC',h)..
    QS(s,h,t) =e=
    // Source
    xS('AC',s) * (
        SourceData(s,h,t,'Q')
        - V(i,h,t) * SourceData(s,h,t,'IQ')
        - E(i,h,t) * SourceData(s,h,t,'IIm')
        + F(i,h,t) * SourceData(s,h,t,'IRe')
        + wVV(i,i,h,t) * SourceData(s,h,t,'b')
    ) ;
$ENDIF
$IFTHEN %BigM% == yes
eq_QSrc_1(SX('AC',i,s),h,TT(t))$HH('AC',h)..
    QS(s,h,t) =l=
    // Source
    (
        SourceData(s,h,t,'Q')
        - V(i,h,t) * SourceData(s,h,t,'IQ')
        - E(i,h,t) * SourceData(s,h,t,'IIm')
        + F(i,h,t) * SourceData(s,h,t,'IRe')
        + wVV(i,i,h,t) * SourceData(s,h,t,'b')
    )
    // Big M
    + MS(s,h,t,'Q','neg') * (1 - xS('AC',s)) ;

eq_QSrc_2(SX('AC',i,s),h,TT(t))$HH('AC',h)..
    QS(s,h,t) =g=
    // Source
    (
        SourceData(s,h,t,'Q')
        - V(i,h,t) * SourceData(s,h,t,'IQ')
        - E(i,h,t) * SourceData(s,h,t,'IIm')
        + F(i,h,t) * SourceData(s,h,t,'IRe')
        + wVV(i,i,h,t) * SourceData(s,h,t,'b')
    )
    // Big M
    - MS(s,h,t,'Q','pos') * (1 - xS('AC',s)) ;

eq_QSrc_3(SX('AC',i,s),h,TT(t))$HH('AC',h)..
    QS(s,h,t) =l=  MS(s,h,t,'Q','pos') * xS('AC',s) ;

eq_QSrc_4(SX('AC',i,s),h,TT(t))$HH('AC',h)..
    QS(s,h,t) =g= -MS(s,h,t,'Q','neg') * xS('AC',s) ;
$ENDIF


** Source Current
$IFTHEN.harm %IncludeHarmonics% == yes
* Source harmonic current, real component
$IFTHEN %BilinearEq% == yes
eq_IReSrc(SX('AC',i,s),h,TT(t))$HX('Harm',h)..
    IReS(s,h,t) =e=
    xS('AC',s) * (
        SourceData(s,h,t,'IRe')
        - SourceData(s,h,t,'g') * E(i,h,t)
        + SourceData(s,h,t,'b') * F(i,h,t)
    ) ;
$ENDIF
$IFTHEN %BigM% == yes
eq_IReSrc_1(SX('AC',i,s),h,TT(t))$HX('Harm',h)..
    IReS(s,h,t) =l=
    // Current
    (
        SourceData(s,h,t,'IRe')
        - SourceData(s,h,t,'g') * E(i,h,t)
        + SourceData(s,h,t,'b') * F(i,h,t)
    )
    // Big M
    + MS(s,h,t,'IRe','neg') * (1 - xS('AC',s)) ;

eq_IReSrc_2(SX('AC',i,s),h,TT(t))$HX('Harm',h)..
    IReS(s,h,t) =g=
    // Current
    (
        SourceData(s,h,t,'IRe')
        - SourceData(s,h,t,'g') * E(i,h,t)
        + SourceData(s,h,t,'b') * F(i,h,t)
    )
    // Big M
    - MS(s,h,t,'IRe','pos') * (1 - xS('AC',s)) ;

eq_IReSrc_3(SX('AC',i,s),h,TT(t))$HX('Harm',h)..
    IReS(s,h,t) =l=  MS(s,h,t,'IRe','pos') * xS('AC',s) ;

eq_IReSrc_4(SX('AC',i,s),h,TT(t))$HX('Harm',h)..
    IReS(s,h,t) =g= -MS(s,h,t,'IRe','neg') * xS('AC',s) ;
$ENDIF

* Source harmonic current, imaginary component
$IFTHEN %BilinearEq% == yes
eq_IImSrc(SX('AC',i,s),h,TT(t))$HX('Harm',h)..
    IImS(s,h,t) =e=
    xS('AC',s) * (
        SourceData(s,h,t,'IIm')
        - SourceData(s,h,t,'b') * E(i,h,t)
        - SourceData(s,h,t,'g') * F(i,h,t)
    ) ;
$ENDIF
$IFTHEN %BigM% == yes
eq_IImSrc_1(SX('AC',i,s),h,TT(t))$HX('Harm',h)..
    IImS(s,h,t) =l=
    // Current
    (
        SourceData(s,h,t,'IIm')
        - SourceData(s,h,t,'b') * E(i,h,t)
        - SourceData(s,h,t,'g') * F(i,h,t)
    )
    // Big M
    + MS(s,h,t,'IIm','neg') * (1 - xS('AC',s)) ;

eq_IImSrc_2(SX('AC',i,s),h,TT(t))$HX('Harm',h)..
    IImS(s,h,t) =g=
    // Current
    (
        SourceData(s,h,t,'IIm')
        - SourceData(s,h,t,'b') * E(i,h,t)
        - SourceData(s,h,t,'g') * F(i,h,t)
    )
    // Big M
    - MS(s,h,t,'IIm','pos') * (1 - xS('AC',s)) ;

eq_IImSrc_3(SX('AC',i,s),h,TT(t))$HX('Harm',h)..
    IImS(s,h,t) =l=  MS(s,h,t,'IIm','pos') * xS('AC',s) ;

eq_IImSrc_4(SX('AC',i,s),h,TT(t))$HX('Harm',h)..
    IImS(s,h,t) =g= -MS(s,h,t,'IIm','neg') * xS('AC',s) ;
$ENDIF
$ENDIF.harm


** Load Power
* Load real power
$IFTHEN %BilinearEq% == yes
eq_PLoad(LX(a,i,l),h,TT(t))$HH(a,h)..
    PL(l,h,t) =e=
    // Load
    xL(a,l) * (
        LoadData(l,h,t,'P')
        + V(i,h,t) * LoadData(l,h,t,'ID')
        + ( E(i,h,t) * LoadData(l,h,t,'IRe') )$sameas(a,'AC')
        + ( F(i,h,t) * LoadData(l,h,t,'IIm') )$sameas(a,'AC')
        + wVV(i,i,h,t) * LoadData(l,h,t,'g')
    ) ;
$ENDIF
$IFTHEN %BigM% == yes
eq_PLoad_1(LX(a,i,l),h,TT(t))$HH(a,h)..
    PL(l,h,t) =l=
    // Load
    (
        LoadData(l,h,t,'P')
        + V(i,h,t) * LoadData(l,h,t,'ID')
        + ( E(i,h,t) * LoadData(l,h,t,'IRe') )$sameas(a,'AC')
        + ( F(i,h,t) * LoadData(l,h,t,'IIm') )$sameas(a,'AC')
        + wVV(i,i,h,t) * LoadData(l,h,t,'g')
    )
    // Big M
    + ML(l,h,t,'P','neg') * (1 - xL(a,l)) ;

eq_PLoad_2(LX(a,i,l),h,TT(t))$HH(a,h)..
    PL(l,h,t) =g=
    // Load
    (
        LoadData(l,h,t,'P')
        + V(i,h,t) * LoadData(l,h,t,'ID')
        + ( E(i,h,t) * LoadData(l,h,t,'IRe') )$sameas(a,'AC')
        + ( F(i,h,t) * LoadData(l,h,t,'IIm') )$sameas(a,'AC')
        + wVV(i,i,h,t) * LoadData(l,h,t,'g')
    )
    // Big M
    - ML(l,h,t,'P','pos') * (1 - xL(a,l)) ;

eq_PLoad_3(LX(a,i,l),h,TT(t))$HH(a,h)..
    PL(l,h,t) =l=  ML(l,h,t,'P','pos') * xL(a,l) ;

eq_PLoad_4(LX(a,i,l),h,TT(t))$HH(a,h)..
    PL(l,h,t) =g= -ML(l,h,t,'P','neg') * xL(a,l) ;
$ENDIF

* Load reactive power
$IFTHEN %BilinearEq% == yes
eq_QLoad(LX('AC',i,l),h,TT(t))$HH('AC',h)..
    QL(l,h,t) =e=
    // Load
    xL('AC',l) * (
        LoadData(l,h,t,'Q')
        - V(i,h,t) * LoadData(l,h,t,'IQ')
        - E(i,h,t) * LoadData(l,h,t,'IIm')
        + F(i,h,t) * LoadData(l,h,t,'IRe')
        - wVV(i,i,h,t) * LoadData(l,h,t,'b')
    ) ;
$ENDIF
$IFTHEN %BigM% == yes
eq_QLoad_1(LX('AC',i,l),h,TT(t))$HH('AC',h)..
    QL(l,h,t) =l=
    // Load
    (
        LoadData(l,h,t,'Q')
        - V(i,h,t) * LoadData(l,h,t,'IQ')
        - E(i,h,t) * LoadData(l,h,t,'IIm')
        + F(i,h,t) * LoadData(l,h,t,'IRe')
        - wVV(i,i,h,t) * LoadData(l,h,t,'b')
    )
    // Big M
    + ML(l,h,t,'Q','neg') * (1 - xL('AC',l)) ;

eq_QLoad_2(LX('AC',i,l),h,TT(t))$HH('AC',h)..
    QL(l,h,t) =g=
    // Load
    (
        LoadData(l,h,t,'Q')
        - V(i,h,t) * LoadData(l,h,t,'IQ')
        - E(i,h,t) * LoadData(l,h,t,'IIm')
        + F(i,h,t) * LoadData(l,h,t,'IRe')
        - wVV(i,i,h,t) * LoadData(l,h,t,'b')
    )
    // Big M
    - ML(l,h,t,'Q','pos') * (1 - xL('AC',l)) ;

eq_QLoad_3(LX('AC',i,l),h,TT(t))$HH('AC',h)..
    QL(l,h,t) =l=  ML(l,h,t,'Q','pos') * xL('AC',l) ;

eq_QLoad_4(LX('AC',i,l),h,TT(t))$HH('AC',h)..
    QL(l,h,t) =g= -ML(l,h,t,'Q','neg') * xL('AC',l) ;
$ENDIF


** Load Current
$IFTHEN.harm %IncludeHarmonics% == yes
* Load harmonic current, real component
$IFTHEN %BilinearEq% == yes
eq_IReLoad(LX('AC',i,l),h,TT(t))$HX('Harm',h)..
    IReL(l,h,t) =e=
    xL('AC',l) * (
        LoadData(l,h,t,'IRe')
        + LoadData(l,h,t,'g') * E(i,h,t)
        - LoadData(l,h,t,'b') * F(i,h,t)
    ) ;
$ENDIF
$IFTHEN %BigM% == yes
eq_IReLoad_1(LX('AC',i,l),h,TT(t))$HX('Harm',h)..
    IReL(l,h,t) =l=
    // Load
    (
        LoadData(l,h,t,'IRe')
        + LoadData(l,h,t,'g') * E(i,h,t)
        - LoadData(l,h,t,'b') * F(i,h,t)
    )
    // Big M
    + ML(l,h,t,'IRe','neg') * (1 - xL('AC',l)) ;

eq_IReLoad_2(LX('AC',i,l),h,TT(t))$HX('Harm',h)..
    IReL(l,h,t) =g=
    // Load
    (
        LoadData(l,h,t,'IRe')
        + LoadData(l,h,t,'g') * E(i,h,t)
        - LoadData(l,h,t,'b') * F(i,h,t)
    )
    // Big M
    - ML(l,h,t,'IRe','pos') * (1 - xL('AC',l)) ;

eq_IReLoad_3(LX('AC',i,l),h,TT(t))$HX('Harm',h)..
    IReL(l,h,t) =l=  ML(l,h,t,'IRe','pos') * xL('AC',l) ;

eq_IReLoad_4(LX('AC',i,l),h,TT(t))$HX('Harm',h)..
    IReL(l,h,t) =g= -ML(l,h,t,'IRe','neg') * xL('AC',l) ;
$ENDIF

* Load harmonic current, imaginary component
$IFTHEN %BilinearEq% == yes
eq_IImLoad(LX('AC',i,l),h,TT(t))$HX('Harm',h)..
    IImL(l,h,t) =e=
    xL('AC',l) * (
        LoadData(l,h,t,'IIm')
        + LoadData(l,h,t,'b') * E(i,h,t)
        + LoadData(l,h,t,'g') * F(i,h,t)
    ) ;
$ENDIF
$IFTHEN %BigM% == yes
eq_IImLoad_1(LX('AC',i,l),h,TT(t))$HX('Harm',h)..
    IImL(l,h,t) =l=
    // Load
    (
        LoadData(l,h,t,'IIm')
        + LoadData(l,h,t,'b') * E(i,h,t)
        + LoadData(l,h,t,'g') * F(i,h,t)
    )
    // Big M
    + ML(l,h,t,'IIm','neg') * (1 - xL('AC',l)) ;

eq_IImLoad_2(LX('AC',i,l),h,TT(t))$HX('Harm',h)..
    IImL(l,h,t) =g=
    // Load
    (
        LoadData(l,h,t,'IIm')
        + LoadData(l,h,t,'b') * E(i,h,t)
        + LoadData(l,h,t,'g') * F(i,h,t)
    )
    // Big M
    - ML(l,h,t,'IIm','pos') * (1 - xL('AC',l)) ;

eq_IImLoad_3(LX('AC',i,l),h,TT(t))$HX('Harm',h)..
    IImL(l,h,t) =l=  ML(l,h,t,'IIm','pos') * xL('AC',l) ;

eq_IImLoad_4(LX('AC',i,l),h,TT(t))$HX('Harm',h)..
    IImL(l,h,t) =g= -ML(l,h,t,'IIm','neg') * xL('AC',l) ;
$ENDIF
$ENDIF.harm


** Converter Power
* Converter output real power
eq_PConvOutMin(c,h,TT(t))..
    PCout(c,h,t) =g= xC(c) * ConvData2(c,h,'Pomin') ;
eq_PConvOutMax(c,h,TT(t))..
    PCout(c,h,t) =l= xC(c) * ConvData2(c,h,'Pomax') ;

* Converter input real power
eq_PConvInMin(c,h,TT(t))..
    PCin(c,h,t) =g= xC(c) * ConvData2(c,h,'Pimin') ;
eq_PConvInMax(c,h,TT(t))..
    PCin(c,h,t) =l= xC(c) * ConvData2(c,h,'Pimax') ;

* Converter input-output power relationship
eq_PConvInOut(c,TT(t))..
    sum(h, PCin(c,h,t)) =e=
    wC(c,t)
        + xC(c) * ConvData(c,'alpha')
        + ConvData(c, 'beta') * wC(c,t)
        + ConvData(c,'gamma') * wCC(c,t) ;

* Converter output reactive power
eq_QConvOutMin(c,h,TT(t))$(CCo('AC',c) and HH('AC',h))..
    QCout(c,h,t) =g= xC(c) * ConvData2(c,h,'Qomin') ;
eq_QConvOutMax(c,h,TT(t))$(CCo('AC',c) and HH('AC',h))..
    QCout(c,h,t) =l= xC(c) * ConvData2(c,h,'Qomax') ;

* Converter input reactive power
eq_QConvInMin(c,h,TT(t))$(CCi('AC',c) and HH('AC',h))..
    QCin(c,h,t) =g= xC(c) * ConvData2(c,h,'Qimin') ;
eq_QConvInMax(c,h,TT(t))$(CCi('AC',c) and HH('AC',h))..
    QCin(c,h,t) =l= xC(c) * ConvData2(c,h,'Qimax') ;


** Converter Current
$IFTHEN.harm %IncludeHarmonics% == yes
* Converter output harmonic current, real component
eq_IReConvOutMin(c,h,TT(t))$(CCo('AC',c) and HX('Harm',h))..
    IReCout(c,h,t) =g= -xC(c) * ConvData2(c,h,'Iomax') ;
eq_IReConvOutMax(c,h,TT(t))$(CCo('AC',c) and HX('Harm',h))..
    IReCout(c,h,t) =l= xC(c) * ConvData2(c,h,'Iomax') ;

* Converter output harmonic current, imaginary component
eq_IImConvOutMin(c,h,TT(t))$(CCo('AC',c) and HX('Harm',h))..
    IImCout(c,h,t) =g= -xC(c) * ConvData2(c,h,'Iomax') ;
eq_IImConvOutMax(c,h,TT(t))$(CCo('AC',c) and HX('Harm',h))..
    IImCout(c,h,t) =l= xC(c) * ConvData2(c,h,'Iomax') ;

* Converter input harmonic current, real component
eq_IReConvInMin(c,h,TT(t))$(CCo('AC',c) and HX('Harm',h))..
    IReCin(c,h,t) =g= -xC(c) * ConvData2(c,h,'Iimax') ;
eq_IReConvInMax(c,h,TT(t))$(CCo('AC',c) and HX('Harm',h))..
    IReCin(c,h,t) =l= xC(c) * ConvData2(c,h,'Iimax') ;

* Converter input harmonic current, imaginary component
eq_IImConvInMin(c,h,TT(t))$(CCo('AC',c) and HX('Harm',h))..
    IImCin(c,h,t) =g= -xC(c) * ConvData2(c,h,'Iimax') ;
eq_IImConvInMax(c,h,TT(t))$(CCo('AC',c) and HX('Harm',h))..
    IImCin(c,h,t) =l= xC(c) * ConvData2(c,h,'Iimax') ;
$ENDIF.harm


** State Variable Relationships
* Voltage magnitude
eq_VoltMag(i,h,TT(t))$(NN('AC',i) and HH('AC',h))..
    wVV(i,i,h,t) =e= wEE(i,i,h,t) + wFF(i,i,h,t) ;

* Total converter power output
eq_wC(c,TT(t))..
    wC(c,t) =e= sum(h, PCout(c,h,t)) ;

* Extra voltage magnitude cuts (based on geometry)
$IFTHEN %ExtraCuts% == yes
eq_VEF_1(i,h,TT(t))$(NN('AC',i) and HH('AC',h) and
    E.lo(i,h,t) >= 0 and F.lo(i,h,t) >= 0 )..
    V(i,h,t) =l= E(i,h,t) + F(i,h,t) ;

eq_VEF_2(i,h,TT(t))$(NN('AC',i) and HH('AC',h) and
    V.lo(i,h,t) >= 0 )..
    E(i,h,t) =l= V(i,h,t) ;

eq_VEF_3(i,h,TT(t))$(NN('AC',i) and HH('AC',h) and
    V.lo(i,h,t) >= 0 )..
    F(i,h,t) =l= V(i,h,t) ;
$ENDIF

* Converter-voltage relationships: E
$IFTHEN.harm %IncludeHarmonics% == yes
eq_EHarmMin(CXh(i,c),h,TT(t))$(NN('AC',i) and HX('Harm',h))..
    E(i,h,t) =g= (1 - xC(c)) * E.lo(i,h,t) ;
eq_EHarmMax(CXh(i,c),h,TT(t))$(NN('AC',i) and HX('Harm',h))..
    E(i,h,t) =l= (1 - xC(c)) * E.up(i,h,t) ;
$ENDIF.harm

* Converter-voltage relationships: F
$IFTHEN.harm %IncludeHarmonics% == yes
eq_FHarmMin(CXh(i,c),h,TT(t))$(NN('AC',i) and HX('Harm',h))..
    F(i,h,t) =g= (1 - xC(c)) * F.lo(i,h,t) ;
eq_FHarmMax(CXh(i,c),h,TT(t))$(NN('AC',i) and HX('Harm',h))..
    F(i,h,t) =l= (1 - xC(c)) * F.up(i,h,t) ;
$ENDIF.harm


** Binary Logic
* Branch assignment
eq_BranchAsgn(B(i,k))..
    sum(a, xB(a,i,k)) =l= 1;

* Source assignment
eq_SrcAsgn(s)..
    sum(a, xS(a,s)) =e= 1;

* Load assignment
eq_LoadAsgn(l)$(not sum((i,k),BLX(i,k,l)))..
    sum(a, xL(a,l)) =e= 1;

* Loads that are also branchs
eq_BranchLoad(a,BLX(i,k,l))..
    xB(a,i,k) =e= xL(a,l);

* Source connection feasibility
eq_SrcConn(SX(a,i,s))..
    xS(a,s) =l= sum(LX(a,i,l), xL(a,l))
        + sum(BX(a,i,k), xB(a,i,k))
        + sum(BX(a,k,i), xB(a,k,i))
        + sum(CXi(a,i,c), xC(c)) ;

* Load connection feasibility
eq_LoadConn(LX(a,i,l))..
    xL(a,l) =l= sum(SX(a,i,s), xS(a,s))
        + sum(BX(a,i,k), xB(a,i,k))
        + sum(BX(a,k,i), xB(a,k,i))
        + sum(CXo(a,i,c), xC(c)) ;

* Converter input connection feasibility
eq_ConvIn(CXi(a,i,c))..
    xC(c) =l= sum(SX(a,i,s), xS(a,s))
        + sum(BX(a,i,k), xB(a,i,k))
        + sum(BX(a,k,i), xB(a,k,i)) ;

* Converter output connection optimality
eq_ConvOut(CXo(a,i,c))..
    xC(c) =l= sum(LX(a,i,l), xL(a,l))
        + sum(BX(a,i,k), xB(a,i,k))
        + sum(BX(a,k,i), xB(a,k,i)) ;


** Reformulation equations (exact definition)
$IFTHEN %ExactReform% == yes
* wVV
eq_wVV(ConnD(a,i,k),h,TT(t))$HH(a,h)..
    wVV(i,k,h,t) =e= V(i,h,t) * V(k,h,t) ;

* wEE
eq_wEE(ConnD('AC',i,k),h,TT(t))$HH('AC',h)..
    wEE(i,k,h,t) =e= E(i,h,t) * E(k,h,t) ;

* wEF
eq_wEF(ConnU('AC',i,k),h,TT(t))$HH('AC',h)..
    wEF(i,k,h,t) =e= E(i,h,t) * F(k,h,t) ;

* wFF
eq_wFF(ConnD('AC',i,k),h,TT(t))$HH('AC',h)..
    wFF(i,k,h,t) =e= F(i,h,t) * F(k,h,t) ;

* wCC
eq_wCC(c,TT(t))..
    wCC(c,t) =e= wC(c,t) * wC(c,t) ;
$ENDIF
