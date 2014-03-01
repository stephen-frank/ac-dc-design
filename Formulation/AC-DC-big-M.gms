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

$TITLE Compute big M's for AC-DC model
$ONTEXT

    Filename : AC-DC-big-M.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    This file computes valid big M values for the reformulation constraints in
    the AC-DC model. It is designed to be included (via $INCLUDE) in other
    .gms files, and relies on the following external compile time flags
    (which are documented in AC-DC-defaults.gms):
        IncludeHarmonics    Include equations for harmonic currents?

    NOTES:
    1. Big M's are computed via a series of small optimization problems which,
       for the current bounds on the voltage variables, find the minimum and
       maximum power for the associated constraint. Therefore, the tighter the
       voltage bounds, the better the big M values.
    2. In addition to setting the big M values, this script also indirectly
       defines the bounds on the associated real and reactive powers.

$OFFTEXT

* ==============================================================================
*  Setup
* ==============================================================================
** Setup
* Parameters for saving existing data
Parameters
    saveV(i,h,t,*)      'Save voltages V, E, F'
    saveW(i,k,h,t,*)    'Save reformulation variables wVV, wEE, wEF, wFF' ;

** Save the existing voltage state
* Voltages
saveV(i,h,t,'V') = V.l(i,h,t);
saveV(i,h,t,'E') = E.l(i,h,t);
saveV(i,h,t,'F') = F.l(i,h,t);

* Reformulation variables
saveW(i,k,h,t,'wVV') = wVV.l(i,k,h,t);
saveW(i,k,h,t,'wEE') = wEE.l(i,k,h,t);
saveW(i,k,h,t,'wEF') = wEF.l(i,k,h,t);
saveW(i,k,h,t,'wFF') = wFF.l(i,k,h,t);


** Macro for Troubleshooting
* Detects unboundedness, infeasibility, or suboptimality, which would imply
* an invalid M value
$MACRO assertValidM(mod) \
    abort$(&mod&.modelstat = %ModelStat.Unbounded%) \
        'Unboundedness detected in big M calculation'; \
    abort$(&mod&.modelstat = %ModelStat.Infeasible%) \
        'Infeasibility detected in big M calculation'; \
    abort$(&mod&.modelstat ne %ModelStat.Optimal%) \
        'Suboptimality detected in big M calculation';

** Set definitions
* Branch selection
Set BranchSelect(a,i,k,h,t);
BranchSelect(a,i,k,h,t) = no;

* Source selection
Set SourceSelect(a,i,s,h,t);
SourceSelect(a,i,s,h,t) = no;

* Load selection
Set LoadSelect(a,i,l,h,t);
LoadSelect(a,i,l,h,t) = no;


* ==============================================================================
*  Branch Power Flow
* ==============================================================================
** Local Model
* Set of equations that define optimization problems for 'MB' (branch big M's)
Equations
    // Branch Power Flow
    eq_Pik_bigM(a,i,k,h,t)  'Branch ik (from -> to) real power flow'
    eq_Pki_bigM(a,i,k,h,t)  'Branch ki (to -> from) real power flow'
    eq_Qik_bigM(a,i,k,h,t)  'Branch ik (from -> to) reactive power flow'
    eq_Qki_bigM(a,i,k,h,t)  'Branch ki (to -> from) reactive power flow' ;

* Branch real power flow (from -> to)
eq_Pik_bigM(BranchSelect(a,i,k,h,t))..
    z =e=
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
    )$sameas(a,'AC') ;

* Branch real power flow (to -> from)
eq_Pki_bigM(BranchSelect(a,i,k,h,t))..
    z =e=
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
    )$sameas(a,'AC') ;

* Branch reactive power flow (from -> to)
eq_Qik_bigM(BranchSelect(a,i,k,h,t))..
    z =e=
    // AC power flow
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
            *  ( wEF(i,k,h,t) - wEF(k,i,h,t) ) ;

* Branch reactive power flow (to -> from)
eq_Qki_bigM(BranchSelect(a,i,k,h,t))..
    z =e=
    // AC power flow
        - ( BranchData(i,k,h,t,'bSe') + BranchData(i,k,h,t,'bSh')/2 )
            * wVV(k,k,h,t)
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * sin(BranchData(i,k,h,t,'phi'))
                + BranchData(i,k,h,t,'bSe') * cos(BranchData(i,k,h,t,'phi')) )
            * ( wEE(i,k,h,t) + wFF(i,k,h,t) )
        + (1/BranchData(i,k,h,t,'A'))
            * ( - BranchData(i,k,h,t,'gSe') * cos(BranchData(i,k,h,t,'phi'))
                - BranchData(i,k,h,t,'bSe') * sin(BranchData(i,k,h,t,'phi')) )
            * ( wEF(i,k,h,t) - wEF(k,i,h,t) ) ;

* Define models
model bigM_branch_Pik /eq_Pik_bigM/ ;
model bigM_branch_Pki /eq_Pki_bigM/ ;
model bigM_branch_Qik /eq_Qik_bigM/ ;
model bigM_branch_Qki /eq_Qki_bigM/ ;


** Loop
* Loop through every branch: DC
loop ( (BX('DC',i2,k2),h2,t2)$HX('DC',h2),
    // Select this combination for the optimization problem
    BranchSelect('DC',i2,k2,h2,t2) = yes;

    // Branch real power
    solve bigM_branch_Pik using LP max z;
    assertValidM(bigM_branch_Pik);
    MB(i2,k2,h2,t2,'P','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_branch_Pik using LP min z;
    assertValidM(bigM_branch_Pik);
    MB(i2,k2,h2,t2,'P','neg') = max( -z.l + %BigMSlack%, 0);

    solve bigM_branch_Pki using LP max z;
    assertValidM(bigM_branch_Pki);
    MB(k2,i2,h2,t2,'P','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_branch_Pki using LP min z;
    assertValidM(bigM_branch_Pki);
    MB(k2,i2,h2,t2,'P','neg') = max( -z.l + %BigMSlack%, 0);

    // Deselect this combination
    BranchSelect('DC',i2,k2,h2,t2) = no;
);

* Loop through every branch: AC fundamental
loop ( (BX('AC',i2,k2),h2,t2)$HX('AC',h2),
    // Select this combination for the optimization problem
    BranchSelect('AC',i2,k2,h2,t2) = yes;

    // Branch real power
    solve bigM_branch_Pik using LP max z;
    assertValidM(bigM_branch_Pik);
    MB(i2,k2,h2,t2,'P','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_branch_Pik using LP min z;
    assertValidM(bigM_branch_Pik);
    MB(i2,k2,h2,t2,'P','neg') = max( -z.l + %BigMSlack%, 0);

    solve bigM_branch_Pki using LP max z;
    assertValidM(bigM_branch_Pki);
    MB(k2,i2,h2,t2,'P','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_branch_Pki using LP min z;
    assertValidM(bigM_branch_Pki);
    MB(k2,i2,h2,t2,'P','neg') = max( -z.l + %BigMSlack%, 0);

    // Branch reactive power
    solve bigM_branch_Qik using LP max z;
    assertValidM(bigM_branch_Qik);
    MB(i2,k2,h2,t2,'Q','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_branch_Qik using LP min z;
    assertValidM(bigM_branch_Qik);
    MB(i2,k2,h2,t2,'Q','neg') = max( -z.l + %BigMSlack%, 0);

    solve bigM_branch_Qki using LP max z;
    assertValidM(bigM_branch_Qki);
    MB(k2,i2,h2,t2,'Q','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_branch_Qki using LP min z;
    assertValidM(bigM_branch_Qki);
    MB(k2,i2,h2,t2,'Q','neg') = max( -z.l + %BigMSlack%, 0);

    // Deselect this combination
    BranchSelect('AC',i2,k2,h2,t2) = no;
);


* ==============================================================================
*  Branch Current
* ==============================================================================
$IFTHEN %IncludeHarmonics% == yes
** Local Model
* Set of equations that define optimization problems for 'MB' (branch big M's)
Equations
    // Branch Current
    eq_IReik_bigM(a,i,k,h,t)    'Branch ik (from -> to) harmonic current, real component'
    eq_IReki_bigM(a,i,k,h,t)    'Branch ki (to -> from) harmonic current, imaginary component'
    eq_IImik_bigM(a,i,k,h,t)    'Branch ik (from -> to) harmonic current, real component'
    eq_IImki_bigM(a,i,k,h,t)    'Branch ki (to -> from) harmonic current, imaginary component' ;

* Branch harmonic current, real component (from -> to)
eq_IReik_bigM(BranchSelect(a,i,k,h,t))..
    z =e=
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
            * F(k,h,t) ;

* Branch harmonic current, real component (to -> from)
eq_IReki_bigM(BranchSelect(a,i,k,h,t))..
    z =e=
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
            * F(i,h,t) ;

* Branch harmonic current, imaginary component (from -> to)
eq_IImik_bigM(BranchSelect(a,i,k,h,t))..
    z =e=
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
            * F(k,h,t) ;

* Branch harmonic current, imaginary component (to -> from)
eq_IImki_bigM(BranchSelect(a,i,k,h,t))..
    z =e=
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
            * F(i,h,t) ;

* Define models
model bigM_branch_IReik /eq_IReik_bigM/ ;
model bigM_branch_IReki /eq_IReki_bigM/ ;
model bigM_branch_IImik /eq_IImik_bigM/ ;
model bigM_branch_IImki /eq_IImki_bigM/ ;


** Loop
* Loop through every branch: AC harmonics
loop ( (BX('AC',i2,k2),h2,t2)$HX('Harm',h2),
    // Select this combination for the optimization problem
    BranchSelect('AC',i2,k2,h2,t2) = yes;

    // Branch current, real component
    solve bigM_branch_IReik using LP max z;
    assertValidM(bigM_branch_IReik);
    MB(i2,k2,h2,t2,'IRe','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_branch_IReik using LP min z;
    assertValidM(bigM_branch_IReik);
    MB(i2,k2,h2,t2,'IRe','neg') = max( -z.l + %BigMSlack%, 0);

    solve bigM_branch_IReki using LP max z;
    assertValidM(bigM_branch_IReki);
    MB(k2,i2,h2,t2,'IRe','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_branch_IReki using LP min z;
    assertValidM(bigM_branch_IReki);
    MB(k2,i2,h2,t2,'IRe','neg') = max( -z.l + %BigMSlack%, 0);

    // Branch current, imaginary component
    solve bigM_branch_IImik using LP max z;
    assertValidM(bigM_branch_IImik);
    MB(i2,k2,h2,t2,'IIm','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_branch_IImik using LP min z;
    assertValidM(bigM_branch_IImik);
    MB(i2,k2,h2,t2,'IIm','neg') = max( -z.l + %BigMSlack%, 0);

    solve bigM_branch_IImki using LP max z;
    assertValidM(bigM_branch_IImki);
    MB(k2,i2,h2,t2,'IIm','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_branch_IImki using LP min z;
    assertValidM(bigM_branch_IImki);
    MB(k2,i2,h2,t2,'IIm','neg') = max( -z.l + %BigMSlack%, 0);

    // Deselect this combination
    BranchSelect('AC',i2,k2,h2,t2) = no;
);
$ENDIF


* ==============================================================================
*  Source Power
* ==============================================================================
** Local Model
* Set of equations that define optimization problems for 'MS' (source big M's)
Equations
    // Source Power
    eq_PSrc_bigM(a,i,s,h,t)     'Source real power'
    eq_QSrc_bigM(a,i,s,h,t)     'Source reactive power' ;

* Source real power
eq_PSrc_bigM(SourceSelect(a,i,s,h,t))..
    z =e=
        SourceData(s,h,t,'P')
        + V(i,h,t) * SourceData(s,h,t,'ID')
        + ( E(i,h,t) * SourceData(s,h,t,'IRe') )$sameas(a,'AC')
        + ( F(i,h,t) * SourceData(s,h,t,'IIm') )$sameas(a,'AC')
        - wVV(i,i,h,t) * SourceData(s,h,t,'g') ;

* Source reactive power
eq_QSrc_bigM(SourceSelect(a,i,s,h,t))..
    z =e=
        SourceData(s,h,t,'Q')
        - V(i,h,t) * SourceData(s,h,t,'IQ')
        - E(i,h,t) * SourceData(s,h,t,'IIm')
        + F(i,h,t) * SourceData(s,h,t,'IRe')
        + wVV(i,i,h,t) * SourceData(s,h,t,'b') ;

* Define models
model bigM_src_P /eq_PSrc_bigM/ ;
model bigM_src_Q /eq_QSrc_bigM/ ;


** Loop
* Loop through every source: DC
loop ( (SX('DC',i2,s2),h2,t2)$HX('DC',h2),
    // Select this combination for the optimization problem
    SourceSelect('DC',i2,s2,h2,t2) = yes;

    // Source real power
    solve bigM_src_P using LP max z;
    assertValidM(bigM_src_P);
    MS(s2,h2,t2,'P','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_src_P using LP min z;
    assertValidM(bigM_src_P);
    MS(s2,h2,t2,'P','neg') = max( -z.l + %BigMSlack%, 0);

    // Deselect this combination
    SourceSelect('DC',i2,s2,h2,t2) = no;
);

* Loop through every source: AC fundamental
loop ( (SX('AC',i2,s2),h2,t2)$HX('AC',h2),
    // Select this combination for the optimization problem
    SourceSelect('AC',i2,s2,h2,t2) = yes;

    // Source real power
    solve bigM_src_P using LP max z;
    assertValidM(bigM_src_P);
    MS(s2,h2,t2,'P','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_src_P using LP min z;
    assertValidM(bigM_src_P);
    MS(s2,h2,t2,'P','neg') = max( -z.l + %BigMSlack%, 0);

    // Source reactive power
    solve bigM_src_Q using LP max z;
    assertValidM(bigM_src_Q);
    MS(s2,h2,t2,'Q','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_src_Q using LP min z;
    assertValidM(bigM_src_Q);
    MS(s2,h2,t2,'Q','neg') = max( -z.l + %BigMSlack%, 0);

    // Deselect this combination
    SourceSelect('AC',i2,s2,h2,t2) = no;
);


* ==============================================================================
*  Source Current
* ==============================================================================
$IFTHEN %IncludeHarmonics% == yes
** Local Model
* Set of equations that define optimization problems for 'MS' (source big M's)
Equations
    // Source Current
    eq_IReSrc_bigM(a,i,s,h,t)   'Source harmonic current, real component'
    eq_IImSrc_bigM(a,i,s,h,t)   'Source harmonic current, imaginary component' ;

* Source harmonic current, real component
eq_IReSrc_bigM(SourceSelect(a,i,s,h,t))..
    z =e=
        SourceData(s,h,t,'IRe')
        - SourceData(s,h,t,'g') * E(i,h,t)
        + SourceData(s,h,t,'b') * F(i,h,t) ;

* Source harmonic current, imaginary component
eq_IImSrc_bigM(SourceSelect(a,i,s,h,t))..
    z =e=
        SourceData(s,h,t,'IIm')
        - SourceData(s,h,t,'b') * E(i,h,t)
        - SourceData(s,h,t,'g') * F(i,h,t) ;

* Define models
model bigM_src_IRe  /eq_IReSrc_bigM/ ;
model bigM_src_IIm  /eq_IImSrc_bigM/ ;


** Loop
* Loop through every source: AC harmonic
loop ( (SX('AC',i2,s2),h2,t2)$HX('Harm',h2),
    // Select this combination for the optimization problem
    SourceSelect('AC',i2,s2,h2,t2) = yes;

    // Source harmonic current, real component
    solve bigM_src_IRe using LP max z;
    assertValidM(bigM_src_IRe);
    MS(s2,h2,t2,'IRe','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_src_IRe using LP min z;
    assertValidM(bigM_src_IRe);
    MS(s2,h2,t2,'IRe','neg') = max( -z.l + %BigMSlack%, 0);

    // Source harmonic current, imaginary component
    solve bigM_src_IIm using LP max z;
    assertValidM(bigM_src_IIm);
    MS(s2,h2,t2,'IIm','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_src_IIm using LP min z;
    assertValidM(bigM_src_IIm);
    MS(s2,h2,t2,'IIm','neg') = max( -z.l + %BigMSlack%, 0);

    // Deselect this combination
    SourceSelect('AC',i2,s2,h2,t2) = no;
);
$ENDIF


* ==============================================================================
*  Load Power
* ==============================================================================
** Local Model
* Set of equations that define optimization problems for 'ML' (load big M's)
Equations
    // Source Power
    eq_PLoad_bigM(a,i,l,h,t)    'Load real power'
    eq_QLoad_bigM(a,i,l,h,t)    'Load reactive power' ;

* Source real power
eq_PLoad_bigM(LoadSelect(a,i,l,h,t))..
    z =e=
        LoadData(l,h,t,'P')
        + V(i,h,t) * LoadData(l,h,t,'ID')
        + ( E(i,h,t) * LoadData(l,h,t,'IRe') )$sameas(a,'AC')
        + ( F(i,h,t) * LoadData(l,h,t,'IIm') )$sameas(a,'AC')
        + wVV(i,i,h,t) * LoadData(l,h,t,'g') ;

* Source reactive power
eq_QLoad_bigM(LoadSelect(a,i,l,h,t))..
    z =e=
        LoadData(l,h,t,'Q')
        - V(i,h,t) * LoadData(l,h,t,'IQ')
        - E(i,h,t) * LoadData(l,h,t,'IIm')
        + F(i,h,t) * LoadData(l,h,t,'IRe')
        - wVV(i,i,h,t) * LoadData(l,h,t,'b') ;

* Define models
model bigM_load_P /eq_PLoad_bigM/ ;
model bigM_load_Q /eq_QLoad_bigM/ ;


** Loop
* Loop through every load: DC
loop ( (LX('DC',i2,l2),h2,t2)$HX('DC',h2),
    // Select this combination for the optimization problem
    LoadSelect('DC',i2,l2,h2,t2) = yes;

    // Load real power
    solve bigM_load_P using LP max z;
    assertValidM(bigM_load_P);
    ML(l2,h2,t2,'P','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_load_P using LP min z;
    assertValidM(bigM_load_P);
    ML(l2,h2,t2,'P','neg') = max( -z.l + %BigMSlack%, 0);

    // Deselect this combination
    LoadSelect('DC',i2,l2,h2,t2) = no;
);

* Loop through every load: AC fundamental
loop ( (LX('AC',i2,l2),h2,t2)$HX('AC',h2),
    // Select this combination for the optimization problem
    LoadSelect('AC',i2,l2,h2,t2) = yes;

    // Load real power
    solve bigM_load_P using LP max z;
    assertValidM(bigM_load_P);
    ML(l2,h2,t2,'P','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_load_P using LP min z;
    assertValidM(bigM_load_P);
    ML(l2,h2,t2,'P','neg') = max( -z.l + %BigMSlack%, 0);

    // Load reactive power
    solve bigM_load_Q using LP max z;
    assertValidM(bigM_load_Q);
    ML(l2,h2,t2,'Q','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_load_Q using LP min z;
    assertValidM(bigM_load_Q);
    ML(l2,h2,t2,'Q','neg') = max( -z.l + %BigMSlack%, 0);

    // Deselect this combination
    LoadSelect('AC',i2,l2,h2,t2) = no;
);


* ==============================================================================
*  Load Current
* ==============================================================================
$IFTHEN %IncludeHarmonics% == yes
** Local Model
* Set of equations that define optimization problems for 'MS' (source big M's)
Equations
    // Source Current
    eq_IReLoad_bigM(a,i,l,h,t)  'Load harmonic current, real component'
    eq_IImLoad_bigM(a,i,l,h,t)  'Load harmonic current, imaginary component' ;

* Load harmonic current, real component
eq_IReLoad_bigM(LoadSelect(a,i,l,h,t))..
    z =e=
        LoadData(l,h,t,'IRe')
        + LoadData(l,h,t,'g') * E(i,h,t)
        - LoadData(l,h,t,'b') * F(i,h,t) ;

* Load harmonic current, imaginary component
eq_IImLoad_bigM(LoadSelect(a,i,l,h,t))..
    z =e=
        LoadData(l,h,t,'IIm')
        + LoadData(l,h,t,'b') * E(i,h,t)
        + LoadData(l,h,t,'g') * F(i,h,t) ;

* Define models
model bigM_load_IRe /eq_IReLoad_bigM/ ;
model bigM_load_IIm /eq_IImLoad_bigM/ ;


** Loop
* Loop through every load: AC harmonic
loop ( (LX('AC',i2,l2),h2,t2)$HX('Harm',h2),
    // Select this combination for the optimization problem
    LoadSelect('AC',i2,l2,h2,t2) = yes;

    // Load harmonic current, real component
    solve bigM_load_IRe using LP max z;
    assertValidM(bigM_load_IRe);
    ML(l2,h2,t2,'IRe','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_load_IRe using LP min z;
    assertValidM(bigM_load_IRe);
    ML(l2,h2,t2,'IRe','neg') = max( -z.l + %BigMSlack%, 0);

    // Load harmonic current, imaginary component
    solve bigM_load_IIm using LP max z;
    assertValidM(bigM_load_IIm);
    ML(l2,h2,t2,'IIm','pos') = max( z.l + %BigMSlack%, 0);

    solve bigM_load_IIm using LP min z;
    assertValidM(bigM_load_IIm);
    ML(l2,h2,t2,'IIm','neg') = max( -z.l + %BigMSlack%, 0);

    // Deselect this combination
    LoadSelect('AC',i2,l2,h2,t2) = no;
);
$ENDIF


* ==============================================================================
*  Cleanup
* ==============================================================================
** Restore the initial voltage state
* Voltages
V.l(i,h,t) = saveV(i,h,t,'V');
E.l(i,h,t) = saveV(i,h,t,'E');
F.l(i,h,t) = saveV(i,h,t,'F');

* Reformulation variables
wVV.l(i,k,h,t) = saveW(i,k,h,t,'wVV');
wEE.l(i,k,h,t) = saveW(i,k,h,t,'wEE');
wEF.l(i,k,h,t) = saveW(i,k,h,t,'wEF');
wFF.l(i,k,h,t) = saveW(i,k,h,t,'wFF');
