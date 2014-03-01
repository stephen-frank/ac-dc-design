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

$TITLE Model file for AC-DC model
$ONTEXT

    Filename : AC-DC-equations-original.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    This file implements the formulation for the (monolith of the) mixed AC-DC
    distribution system design problem. It is designed to be included (via
    $INCLUDE) in other .gms files.

    NOTES:
    1. Equations are declared over full sets, but are defined only for subsets.
       See:
       http://support.gams.com/doku.php?id=gams:define_variables_over_subsets

$OFFTEXT

* ==============================================================================
*  Equations
* ==============================================================================
* ------------------------------------------------------------------------------
*  Equation Declaration
* ------------------------------------------------------------------------------
Equations
    // Objective Function
    eq_obj                  'Objective function'
    // Power Balance
    eq_PBal(a,i,h,t)        'Real power balance'
    eq_QBal(a,i,h,t)        'Reactive power balance'
    // Current Balance
    eq_IReBal(a,i,h,t)      'Harmonic current balance, real component'
    eq_IImBal(a,i,h,t)      'Harmonic current balance, imaginary component'
    // Branch Power Flow
    eq_Pik(a,i,k,h,t)       'Branch ik (from -> to) real power flow'
    eq_Pki(a,i,k,h,t)       'Branch ki (to -> from) real power flow'
    eq_Qik(a,i,k,h,t)       'Branch ik (from -> to) reactive power flow'
    eq_Qki(a,i,k,h,t)       'Branch ki (to -> from) reactive power flow'
    // Branch Current
    eq_IReik(a,i,k,h,t)     'Branch ik (from -> to) harmonic current, real component'
    eq_IReki(a,i,k,h,t)     'Branch ki (to -> from) harmonic current, real component'
    eq_IImik(a,i,k,h,t)     'Branch ik (from -> to) harmonic current, imaginary component'
    eq_IImki(a,i,k,h,t)     'Branch ki (to -> from) harmonic current, imaginary component'
    // Source Power
    eq_PSrc(a,i,s,h,t)      'Source real power'
    eq_QSrc(a,i,s,h,t)      'Source reactive power'
    // Source Current
    eq_IReSrc(a,i,s,h,t)    'Source harmonic current, real component'
    eq_IImSrc(a,i,s,h,t)    'Source harmonic current, imaginary component'
    // Load Power
    eq_PLoad(a,i,l,h,t)     'Load real power'
    eq_QLoad(a,i,l,h,t)     'Load reactive power'
    // Load Current
    eq_IReLoad(a,i,l,h,t)   'Load harmonic current, real component'
    eq_IImLoad(a,i,l,h,t)   'Load harmonic current, imaginary component'
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
    eq_IReConvOutMin(c,h,t) 'Logic switch for min. converter output harmonic current, real component'
    eq_IReConvOutMax(c,h,t) 'Logic switch for max. converter output harmonic current, real component'
    eq_IReConvInMin(c,h,t)  'Logic switch for min. converter input harmonic current, real component'
    eq_IReConvInMax(c,h,t)  'Logic switch for max. converter input harmonic current, real component'
    eq_IImConvOutMin(c,h,t) 'Logic switch for min. converter output harmonic current, imaginary component'
    eq_IImConvOutMax(c,h,t) 'Logic switch for max. converter output harmonic current, imaginary component'
    eq_IImConvInMin(c,h,t)  'Logic switch for min. converter input harmonic current, imaginary component'
    eq_IImConvInMax(c,h,t)  'Logic switch for max. converter input harmonic current, imaginary component'
    // State Variable Relationships
    eq_VoltMag(i,h,t)       'Computation of voltage magnitude'
    eq_EHarmMin(i,c,h,t)    'Logic switch for min. real voltage at converter terminal'
    eq_EHarmMax(i,c,h,t)    'Logic switch for max. real voltage at converter terminal'
    eq_FHarmMin(i,c,h,t)    'Logic switch for min. imaginary voltage at converter terminal'
    eq_FHarmMax(i,c,h,t)    'Logic switch for max. imaginary voltage at converter terminal'
    // Binary Logic
    eq_BranchAsgn(i,k)      'Assignment of branches to AC or DC network'
    eq_SrcAsgn(s)           'Assignment of sources to AC or DC network'
    eq_LoadAsgn(l)          'Assignment of loads to AC or DC network'
    eq_BranchLoad(a,i,k,l)  'Maps branches that are also loads'
    eq_SrcConn(a,i,s)       'Ensure feasible connection of sources'
    eq_LoadConn(a,i,l)      'Ensure feasible connection of loads'
    eq_ConvIn(a,i,c)        'Ensure feasible connection of converter inputs'
    eq_ConvOut(a,i,c)       'Exclude suboptimal connection of converter outputs' ;

* ------------------------------------------------------------------------------
*  Equation Definitions
* ------------------------------------------------------------------------------
** Objective Function
* Objective = minimize utility energy at fundamental AC (h=1)
eq_obj..
    z =e= sum((HH('AC',h),TT(t)), Tau(t) * PS('utility',h,t) ) ;


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


** Branch Power Flow
* Branch real power flow (from -> to)
eq_Pik(BX(a,i,k),h,TT(t))$HH(a,h)..
    PB(i,k,h,t) =e=
    // DC power flow
    ( xB(a,i,k) * (
        (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * V(i,h,t) * V(i,h,t)
        - (1/BranchData(i,k,h,t,'A'))
            * BranchData(i,k,h,t,'gSe')
            * V(i,h,t) * V(k,h,t)
    ))$sameas(a,'DC')
    // AC power flow
    + ( xB(a,i,k) * (
        (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * ( E(i,h,t) * E(i,h,t) + F(i,h,t) * F(i,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( E(i,h,t) * E(k,h,t) + F(i,h,t) * F(k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * (   BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( E(i,h,t) * F(k,h,t) - F(i,h,t) * E(k,h,t) )
    ))$sameas(a,'AC') ;

* Branch real power flow (to -> from)
eq_Pki(BX(a,i,k),h,TT(t))$HH(a,h)..
    PB(k,i,h,t) =e=
    // DC power flow
    ( xB(a,i,k) * (
        ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * V(k,h,t) * V(k,h,t)
        - (1/BranchData(i,k,h,t,'A'))
            * BranchData(i,k,h,t,'gSe')
            * V(i,h,t) * V(k,h,t)
    ))$sameas(a,'DC')
    // AC power flow
    + ( xB(a,i,k) * (
        ( BranchData(i,k,h,t,'gSe') + BranchData(i,k,h,t,'gSh')/2 )
            * ( E(k,h,t) * E(k,h,t) + F(k,h,t) * F(k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( E(i,h,t) * E(k,h,t) + F(i,h,t) * F(k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * (   BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( E(i,h,t) * F(k,h,t) - F(i,h,t) * E(k,h,t) )
    ))$sameas(a,'AC') ;

* Branch reactive power flow (from -> to)
eq_Qik(BX('AC',i,k),h,TT(t))$HH('AC',h)..
    QB(i,k,h,t) =e=
    // AC power flow
    xB('AC',i,k) * (
        - (1/sqr(BranchData(i,k,h,t,'A')))
            * ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * ( E(i,h,t) * E(i,h,t) + F(i,h,t) * F(i,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * ( + BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( E(i,h,t) * E(k,h,t) + F(i,h,t) * F(k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * (   BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( E(i,h,t) * F(k,h,t) - F(i,h,t) * E(k,h,t) )
    ) ;

* Branch reactive power flow (to -> from)
eq_Qki(BX('AC',i,k),h,TT(t))$HH('AC',h)..
    QB(k,i,h,t) =e=
    // AC power flow
    xB('AC',i,k) * (
        - ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * ( E(k,h,t) * E(k,h,t) + F(k,h,t) * F(k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( E(i,h,t) * E(k,h,t) + F(i,h,t) * F(k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( E(i,h,t) * F(k,h,t) - F(i,h,t) * E(k,h,t) )
    ) ;


** Branch Current
* Branch harmonic current, real component (from -> to)
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

* Branch harmonic current, real component (to -> from)
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

* Branch harmonic current, imaginary component (from -> to)
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

* Branch harmonic current, imaginary component (to -> from)
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


** Source Power
* Source real power
eq_PSrc(SX(a,i,s),h,TT(t))$HH(a,h)..
    PS(s,h,t) =e=
    xS(a,s) * (
        SourceData(s,h,t,'P')
        + V(i,h,t) * SourceData(s,h,t,'ID')
        + ( E(i,h,t) * SourceData(s,h,t,'IRe') )$sameas(a,'AC')
        + ( F(i,h,t) * SourceData(s,h,t,'IIm') )$sameas(a,'AC')
        - V(i,h,t) * V(i,h,t) * SourceData(s,h,t,'g')
    ) ;

* Source reactive power
eq_QSrc(SX('AC',i,s),h,TT(t))$HH('AC',h)..
    QS(s,h,t) =e=
    xS('AC',s) * (
        SourceData(s,h,t,'Q')
        - V(i,h,t) * SourceData(s,h,t,'IQ')
        - E(i,h,t) * SourceData(s,h,t,'IIm')
        + F(i,h,t) * SourceData(s,h,t,'IRe')
        + V(i,h,t) * V(i,h,t) * SourceData(s,h,t,'b')
    ) ;


** Source Current
* Source harmonic current, real component
eq_IReSrc(SX('AC',i,s),h,TT(t))$HX('Harm',h)..
    IReS(s,h,t) =e=
    xS('AC',s) * (
        SourceData(s,h,t,'IRe')
        - SourceData(s,h,t,'g') * E(i,h,t)
        + SourceData(s,h,t,'b') * F(i,h,t)
    ) ;

* Source harmonic current, imaginary component
eq_IImSrc(SX('AC',i,s),h,TT(t))$HX('Harm',h)..
    IImS(s,h,t) =e=
    xS('AC',s) * (
        SourceData(s,h,t,'IIm')
        - SourceData(s,h,t,'b') * E(i,h,t)
        - SourceData(s,h,t,'g') * F(i,h,t)
    ) ;


** Load Power
* Load real power
eq_PLoad(LX(a,i,l),h,TT(t))$HH(a,h)..
    PL(l,h,t) =e=
    xL(a,l) * (
        LoadData(l,h,t,'P')
        + V(i,h,t) * LoadData(l,h,t,'ID')
        + ( E(i,h,t) * LoadData(l,h,t,'IRe') )$sameas(a,'AC')
        + ( F(i,h,t) * LoadData(l,h,t,'IIm') )$sameas(a,'AC')
        + V(i,h,t) * V(i,h,t) * LoadData(l,h,t,'g')
    ) ;

* Load reactive power
eq_QLoad(LX('AC',i,l),h,TT(t))$HH('AC',h)..
    QL(l,h,t) =e=
    xL('AC',l) * (
        LoadData(l,h,t,'Q')
        - V(i,h,t) * LoadData(l,h,t,'IQ')
        - E(i,h,t) * LoadData(l,h,t,'IIm')
        + F(i,h,t) * LoadData(l,h,t,'IRe')
        - V(i,h,t) * V(i,h,t) * LoadData(l,h,t,'b')
    ) ;


** Load Current
* Load harmonic current, real component
eq_IReLoad(LX('AC',i,l),h,TT(t))$HX('Harm',h)..
    IReL(l,h,t) =e=
    xL('AC',l) * (
        LoadData(l,h,t,'IRe')
        + LoadData(l,h,t,'g') * E(i,h,t)
        - LoadData(l,h,t,'b') * F(i,h,t)
    ) ;

* Load harmonic current, imaginary component
eq_IImLoad(LX('AC',i,l),h,TT(t))$HX('Harm',h)..
    IImL(l,h,t) =e=
    xL('AC',l) * (
        LoadData(l,h,t,'IIm')
        + LoadData(l,h,t,'b') * E(i,h,t)
        + LoadData(l,h,t,'g') * F(i,h,t)
    ) ;


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
    xC(c) * (
        sum(h, PCout(c,h,t))
        + ConvData(c,'alpha')
        + ConvData(c, 'beta') * sum(h, PCout(c,h,t))
        + ConvData(c,'gamma') * sqr( sum(h, PCout(c,h,t)) )
    ) ;

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


** State Variable Relationships
* Voltage magnitude
eq_VoltMag(i,h,TT(t))$(NN('AC',i) and HH('AC',h))..
    V(i,h,t) * V(i,h,t) =e= E(i,h,t) * E(i,h,t) + F(i,h,t) * F(i,h,t) ;

* Converter-voltage relationships: E
eq_EHarmMin(CXh(i,c),h,TT(t))$(NN('AC',i) and HX('Harm',h))..
    E(i,h,t) =g= (1 - xC(c)) * E.lo(i,h,t) ;
eq_EHarmMax(CXh(i,c),h,TT(t))$(NN('AC',i) and HX('Harm',h))..
    E(i,h,t) =l= (1 - xC(c)) * E.up(i,h,t) ;

* Converter-voltage relationships: F
eq_FHarmMin(CXh(i,c),h,TT(t))$(NN('AC',i) and HX('Harm',h))..
    F(i,h,t) =g= (1 - xC(c)) * F.lo(i,h,t) ;
eq_FHarmMax(CXh(i,c),h,TT(t))$(NN('AC',i) and HX('Harm',h))..
    F(i,h,t) =l= (1 - xC(c)) * F.up(i,h,t) ;


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
