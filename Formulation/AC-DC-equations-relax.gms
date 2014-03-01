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

$TITLE Relaxation equations for AC-DC model
$ONTEXT

    Filename : AC-DC-equations-relax.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    This file implements the relaxation for the mixed AC-DC distribution system
    design problem. It is designed to be included (via $INCLUDE) in other .gms
    files, and relies on the following external compile time settings (which are
    documented in AC-DC-defaults.gms):
        ConvexRelaxation    Specify type of convex relaxation
        PiecewiseType       Specify piecewise relaxation subtype
        NumSeg              Specify number of segments to use in relaxations

    NOTES:
    1. ConvexRelaxation selects among a number of possible convex relaxations
       for the model. Note that some of the relaxations require certain bounds
       on the variables (e.g. the log transformation requires strictly positive
       lower voltage bounds).
    2. You can only generate one convex relaxation out of the specified options
       'ConvexRelaxation' and 'PiecewiseType'.

$OFFTEXT

* ==============================================================================
*  Equations
* ==============================================================================
* ------------------------------------------------------------------------------
*  Sanity Checks
* ------------------------------------------------------------------------------
** Sanity checks
* Quadratic convex relaxation requires nonnegative voltages:
$IFTHEN %ConvexRelaxation% == quadratic
if ( sum((HH(a,h),i,t)$NN(a,i), V.lo(i,h,t) < 0) or
     sum((HH('AC',h),i,t)$NN('AC',i), E.lo(i,h,t) < 0) or
     sum((HH('AC',h),i,t)$NN('AC',i), F.lo(i,h,t) < 0),
    abort 'For quadratic relaxation, lower bounds on voltages V, E, and F must be nonnegative';
);

* Log tranforming the bilinear terms requires strictly positive voltages:
$ELSEIF %ConvexRelaxation% == piecewise
$IFTHEN.PW %PiecewiseType% == logtransform
if ( sum((HH(a,h),i,t)$NN(a,i), V.lo(i,h,t) <= 0) or
     sum((HH('AC',h),i,t)$NN('AC',i), E.lo(i,h,t) <= 0) or
     sum((HH('AC',h),i,t)$NN('AC',i), F.lo(i,h,t) <= 0),
    abort 'For piecewise linear relaxation by log transformation, lower bounds on voltages V, E, and F must be strictly positive';
);
$ENDIF.PW

$ENDIF


* ------------------------------------------------------------------------------
*  Aux. Variables, Parameters, and Sets (For Piecewise Relaxations)
* ------------------------------------------------------------------------------
$IFTHEN NOT %ConvexRelaxation% == quadratic
** Setup
* Set of breakpoints - same for all variables
Set BP  'Breakpoints of (piecewise) linear relaxation'  /b0*b%NumSeg%/ ;

$IFTHEN.PW %ConvexRelaxation% == piecewise
* Set of segments/divisions - same for all variables
Set D   'Segments of piecewise linear relaxation'       /d1*d%NumSeg%/ ;
Alias (D, DD);
$ENDIF.PW

$ENDIF


$IFTHEN %ConvexRelaxation% == piecewise
** Univariate Quadratic Terms
* Variables for univariate piecewise linear relaxations
SOS2 Variables
    lambdaV(i,h,t,bp)   'SOS2 segment selection variables for V'
    lambdaE(i,h,t,bp)   'SOS2 segment selection variables for E'
    lambdaF(i,h,t,bp)   'SOS2 segment selection variables for F'
    lambdaC(c,t,bp)     'SOS2 segment selection variables for wC' ;

* Bounds on lambda variables
lambdaV.lo(i,h,t,bp) = 0;
lambdaE.lo(i,h,t,bp) = 0;
lambdaF.lo(i,h,t,bp) = 0;
lambdaC.lo(c,t,bp)   = 0;

lambdaV.up(i,h,t,bp) = 1;
lambdaE.up(i,h,t,bp) = 1;
lambdaF.up(i,h,t,bp) = 1;
lambdaC.up(c,t,bp)   = 1;
$ENDIF

$IFTHEN NOT %ConvexRelaxation% == quadratic
* Parameters for univariate linear and piecewise linear relaxations
Parameters
    breakV(i,h,t,bp)    'Breakpoints on V segments'
    breakE(i,h,t,bp)    'Breakpoints on E segments'
    breakF(i,h,t,bp)    'Breakpoints on F segments'
    breakC(c,t,bp)      'Breakpoints on wC segments' ;
$ENDIF


$IFTHEN %ConvexRelaxation% == piecewise
$IFTHEN.pw %PiecewiseType% == mccormick
** Bilinear Terms - Using McCormick Inequalities
* Variables for bilinear piecewise linear relaxations
Variables
    segVVi(i,k,h,t,d,dd)    'Piecewise linear variable for Vi in Vi-Vk pair'
    segVVk(i,k,h,t,d,dd)    'Piecewise linear variable for Vk in Vi-Vk pair'
    segEEi(i,k,h,t,d,dd)    'Piecewise linear variable for Ei in Ei-Ek pair'
    segEEk(i,k,h,t,d,dd)    'Piecewise linear variable for Ek in Ei-Ek pair'
    segEFi(i,k,h,t,d,dd)    'Piecewise linear variable for Ei in Ei-Fk pair'
    segEFk(i,k,h,t,d,dd)    'Piecewise linear variable for Fk in Ei-Fk pair'
    segFFi(i,k,h,t,d,dd)    'Piecewise linear variable for Fi in Fi-Fk pair'
    segFFk(i,k,h,t,d,dd)    'Piecewise linear variable for Fk in Fi-Fk pair' ;

Binary variables
    yVV(i,k,h,t,d,dd)       'Segment selection variable for V-V pairs'
    yEE(i,k,h,t,d,dd)       'Segment selection variable for E-E pairs'
    yEF(i,k,h,t,d,dd)       'Segment selection variable for E-F pairs'
    yFF(i,k,h,t,d,dd)       'Segment selection variable for F-F pairs' ;

* Parameters for bilinear piecewise linear relaxations
Parameters
    seglimV(i,h,t,d,*)      'Upper and lower limits on V segments'
    seglimE(i,h,t,d,*)      'Upper and lower limits on E segments'
    seglimF(i,h,t,d,*)      'Upper and lower limits on F segments' ;


$ELSEIF.pw %PiecewiseType% == logtransform
** Bilinear Terms - Using Logarithmic Transform
* Variables for bilinear piecewise linear relaxations
SOS2 Variables
    lambdaVV(i,k,h,t,bp)   'SOS2 segment selection variables for wVV'
    lambdaEE(i,k,h,t,bp)   'SOS2 segment selection variables for wEE'
    lambdaEF(i,k,h,t,bp)   'SOS2 segment selection variables for wEF'
    lambdaFF(i,k,h,t,bp)   'SOS2 segment selection variables for wFF' ;

Variables
    lnV(i,h,t)              'Log transformation variable for ln(V)'
    lnE(i,h,t)              'Log transformation variable for ln(E)'
    lnF(i,h,t)              'Log transformation variable for ln(F)'
    lnwVV(i,k,h,t)          'Log transformation variable for ln(wVV)'
    lnwEE(i,k,h,t)          'Log transformation variable for ln(wEE)'
    lnwEF(i,k,h,t)          'Log transformation variable for ln(wEF)'
    lnwFF(i,k,h,t)          'Log transformation variable for ln(wFF)' ;

* Bounds on lambda variables
lambdaVV.lo(i,k,h,t,bp) = 0;
lambdaEE.lo(i,k,h,t,bp) = 0;
lambdaEF.lo(i,k,h,t,bp) = 0;
lambdaFF.lo(i,k,h,t,bp) = 0;

lambdaVV.up(i,k,h,t,bp) = 1;
lambdaEE.up(i,k,h,t,bp) = 1;
lambdaEF.up(i,k,h,t,bp) = 1;
lambdaFF.up(i,k,h,t,bp) = 1;

* Parameters for log transformation
Parameters
    breakVV(i,k,h,t,bp) 'Breakpoints on wVV segments'
    breakEE(i,k,h,t,bp) 'Breakpoints on wEE segments'
    breakEF(i,k,h,t,bp) 'Breakpoints on wEF segments'
    breakFF(i,k,h,t,bp) 'Breakpoints on wFF segments' ;

$ENDIF.pw
$ENDIF

* Compute initial bounds and breakpoints
$BATINCLUDE %RelPath%/AC-DC-bounds-breaks.gms breakpoints

* ------------------------------------------------------------------------------
*  Relaxation Equation Declaration
* ------------------------------------------------------------------------------
* Note: Some equations are reused between different relaxation strategies.
* However, I have simply declared all equations for each relaxation strategy as
* a single group for clarity. However, in the equation definitions below
* the compile-time logic gets a bit more complicated.

$IFTHEN %ConvexRelaxation% == linear
Equations
    // Univariate terms
    eq_wVVii_over(a,i,h,t)          'Linear overestimator  for reformulation variable wVV for i = k'
    eq_wVVii_under(a,i,h,t,bp)      'Linear underestimator for reformulation variable wVV for i = k'

    eq_wEEii_over(a,i,h,t)          'Linear overestimator  for reformulation variable wEE for i = k'
    eq_wEEii_under(a,i,h,t,bp)      'Linear underestimator for reformulation variable wEE for i = k'

    eq_wFFii_over(a,i,h,t)          'Linear overestimator  for reformulation variable wFF for i = k'
    eq_wFFii_under(a,i,h,t,bp)      'Linear underestimator for reformulation variable wFF for i = k'

    eq_wCC_over(c,t)                'Linear overestimator  for reformulation variable wCC'
    eq_wCC_under(c,t,bp)            'Linear underestimator for reformulation variable wCC'

    // Bilinear terms
    eq_wVVik_over_1(a,i,k,h,t)      'Linear overestimator  for reformulation variable wVV - 1 of 2'
    eq_wVVik_over_2(a,i,k,h,t)      'Linear overestimator  for reformulation variable wVV - 2 of 2'
    eq_wVVik_under_1(a,i,k,h,t)     'Linear underestimator for reformulation variable wVV - 1 of 2'
    eq_wVVik_under_2(a,i,k,h,t)     'Linear underestimator for reformulation variable wVV - 2 of 2'

    eq_wEEik_over_1(a,i,k,h,t)      'Linear overestimator  for reformulation variable wEE - 1 of 2'
    eq_wEEik_over_2(a,i,k,h,t)      'Linear overestimator  for reformulation variable wEE - 2 of 2'
    eq_wEEik_under_1(a,i,k,h,t)     'Linear underestimator for reformulation variable wEE - 1 of 2'
    eq_wEEik_under_2(a,i,k,h,t)     'Linear underestimator for reformulation variable wEE - 2 of 2'

    eq_wEFik_over_1(a,i,k,h,t)      'Linear overestimator  for reformulation variable wEF - 1 of 2'
    eq_wEFik_over_2(a,i,k,h,t)      'Linear overestimator  for reformulation variable wEF - 2 of 2'
    eq_wEFik_under_1(a,i,k,h,t)     'Linear underestimator for reformulation variable wEF - 1 of 2'
    eq_wEFik_under_2(a,i,k,h,t)     'Linear underestimator for reformulation variable wEF - 2 of 2'

    eq_wFFik_over_1(a,i,k,h,t)      'Linear overestimator  for reformulation variable wFF - 1 of 2'
    eq_wFFik_over_2(a,i,k,h,t)      'Linear overestimator  for reformulation variable wFF - 2 of 2'
    eq_wFFik_under_1(a,i,k,h,t)     'Linear underestimator for reformulation variable wFF - 1 of 2'
    eq_wFFik_under_2(a,i,k,h,t)     'Linear underestimator for reformulation variable wFF - 2 of 2' ;

$ELSEIF %ConvexRelaxation% == quadratic
Equations
    // Univariate terms
    eq_wVVii_over(a,i,h,t)          'Linear overestimator     for reformulation variable wVV for i = k'
    eq_wVVii_under(a,i,h,t)         'Quadratic underestimator for reformulation variable wVV for i = k'

    eq_wEEii_over(a,i,h,t)          'Linear overestimator     for reformulation variable wEE for i = k'
    eq_wEEii_under(a,i,h,t)         'Quadratic underestimator for reformulation variable wEE for i = k'

    eq_wFFii_over(a,i,h,t)          'Linear overestimator     for reformulation variable wFF for i = k'
    eq_wFFii_under(a,i,h,t)         'Quadratic underestimator for reformulation variable wFF for i = k'

    eq_wCC_over(c,t)                'Linear overestimator     for reformulation variable wCC'
    eq_wCC_under(c,t)               'Quadratic underestimator for reformulation variable wCC'

    // Bilinear terms
    eq_wVVik_over_1(a,i,k,h,t)      'Linear overestimator     for reformulation variable wVV - 1 of 2'
    eq_wVVik_over_2(a,i,k,h,t)      'Linear overestimator     for reformulation variable wVV - 2 of 2'
    eq_wVVik_under(a,i,k,h,t)       'Quadratic underestimator for reformulation variable wVV'

    eq_wEEik_over_1(a,i,k,h,t)      'Linear overestimator     for reformulation variable wEE - 1 of 2'
    eq_wEEik_over_2(a,i,k,h,t)      'Linear overestimator     for reformulation variable wEE - 2 of 2'
    eq_wEEik_under(a,i,k,h,t)       'Quadratic underestimator for reformulation variable wEE'

    eq_wEFik_over_1(a,i,k,h,t)      'Linear overestimator     for reformulation variable wEF - 1 of 2'
    eq_wEFik_over_2(a,i,k,h,t)      'Linear overestimator     for reformulation variable wEF - 2 of 2'
    eq_wEFik_under(a,i,k,h,t)       'Quadratic underestimator for reformulation variable wEF'

    eq_wFFik_over_1(a,i,k,h,t)      'Linear overestimator     for reformulation variable wFF - 1 of 2'
    eq_wFFik_over_2(a,i,k,h,t)      'Linear overestimator     for reformulation variable wFF - 2 of 2'
    eq_wFFik_under(a,i,k,h,t)       'Quadratic underestimator for reformulation variable wFF' ;


$ELSEIF %ConvexRelaxation% == piecewise
Equations
    // Select piecewise segments for univariate terms
    eq_lambdaV(a,i,h,t)             'Restriction on segment selection for V'
    eq_breakV(a,i,h,t)              'Relate V to its piecewise linear breakpoints'

    eq_lambdaE(a,i,h,t)             'Restriction on segment selection for E'
    eq_breakE(a,i,h,t)              'Relate E to its piecewise linear breakpoints'

    eq_lambdaF(a,i,h,t)             'Restriction on segment selection for F'
    eq_breakF(a,i,h,t)              'Relate F to its piecewise linear breakpoints'

    eq_lambdaC(c,t)                 'Restriction on segment selection for wC'
    eq_breakC(c,t)                  'Relate wC to its piecewise linear breakpoints'

    // Univariate terms
    eq_wVVii_over(a,i,h,t)          'Piecewise linear overestimator for reformulation variable wVV for i = k'
    eq_wVVii_under(a,i,h,t,bp)      'Linear underestimator          for reformulation variable wVV for i = k'

    eq_wEEii_over(a,i,h,t)          'Piecewise linear overestimator for reformulation variable wEE for i = k'
    eq_wEEii_under(a,i,h,t,bp)      'Linear underestimator          for reformulation variable wEE for i = k'

    eq_wFFii_over(a,i,h,t)          'Piecewise linear overestimator for reformulation variable wFF for i = k'
    eq_wFFii_under(a,i,h,t,bp)      'Linear underestimator          for reformulation variable wFF for i = k'

    eq_wCC_over(c,t)                'Piecewise linear overestimator for reformulation variable wCC'
    eq_wCC_under(c,t,bp)            'Linear underestimator          for reformulation variable wCC' ;

$IFTHEN.PW %PiecewiseType% == mccormick
Equations
    // Select piecewise segments for bilinear terms
    eq_segsel_VV(a,i,k,h,t)         'Select exactly one segment for each Vi-Vk pair'
    eq_seg_VVi(a,i,k,h,t)           'Equate Vi with its segment variables for Vi-Vk pair'
    eq_seg_VVk(a,i,k,h,t)           'Equate Vk with its segment variables for Vi-Vk pair'
    eq_seg_VVi_lo(a,i,k,h,t,d,dd)   'Lower bound on Vi segment variables for each Vi-Vk pair'
    eq_seg_VVi_up(a,i,k,h,t,d,dd)   'Upper bound on Vi segment variables for each Vi-Vk pair'
    eq_seg_VVk_lo(a,i,k,h,t,d,dd)   'Lower bound on Vk segment variables for each Vi-Vk pair'
    eq_seg_VVk_up(a,i,k,h,t,d,dd)   'Upper bound on Vk segment variables for each Vi-Vk pair'

    eq_segsel_EE(a,i,k,h,t)         'Select exactly one segment for each Ei-Ek pair'
    eq_seg_EEi(a,i,k,h,t)           'Equate Ei with its segment variables for Ei-Ek pair'
    eq_seg_EEk(a,i,k,h,t)           'Equate Ek with its segment variables for Ei-Ek pair'
    eq_seg_EEi_lo(a,i,k,h,t,d,dd)   'Lower bound on Ei segment variables for each Ei-Ek pair'
    eq_seg_EEi_up(a,i,k,h,t,d,dd)   'Upper bound on Ei segment variables for each Ei-Ek pair'
    eq_seg_EEk_lo(a,i,k,h,t,d,dd)   'Lower bound on Ek segment variables for each Ei-Ek pair'
    eq_seg_EEk_up(a,i,k,h,t,d,dd)   'Upper bound on Ek segment variables for each Ei-Ek pair'

    eq_segsel_EF(a,i,k,h,t)         'Select exactly one segment for each Ei-Fk pair'
    eq_seg_EFi(a,i,k,h,t)           'Equate Ei with its segment variables for Ei-Fk pair'
    eq_seg_EFk(a,i,k,h,t)           'Equate Fk with its segment variables for Ei-Fk pair'
    eq_seg_EFi_lo(a,i,k,h,t,d,dd)   'Lower bound on Ei segment variables for each Ei-Fk pair'
    eq_seg_EFi_up(a,i,k,h,t,d,dd)   'Upper bound on Ei segment variables for each Ei-Fk pair'
    eq_seg_EFk_lo(a,i,k,h,t,d,dd)   'Lower bound on Fk segment variables for each Ei-Fk pair'
    eq_seg_EFk_up(a,i,k,h,t,d,dd)   'Upper bound on Fk segment variables for each Ei-Fk pair'

    eq_segsel_FF(a,i,k,h,t)         'Select exactly one segment for each Fi-Fk pair'
    eq_seg_FFi(a,i,k,h,t)           'Equate Fi with its segment variables for Fi-Fk pair'
    eq_seg_FFk(a,i,k,h,t)           'Equate Fk with its segment variables for Fi-Fk pair'
    eq_seg_FFi_lo(a,i,k,h,t,d,dd)   'Lower bound on Fi segment variables for each Fi-Fk pair'
    eq_seg_FFi_up(a,i,k,h,t,d,dd)   'Upper bound on Fi segment variables for each Fi-Fk pair'
    eq_seg_FFk_lo(a,i,k,h,t,d,dd)   'Lower bound on Fk segment variables for each Fi-Fk pair'
    eq_seg_FFk_up(a,i,k,h,t,d,dd)   'Upper bound on Fk segment variables for each Fi-Fk pair'

    // Bilinear terms
    eq_wVVik_over_1(a,i,k,h,t)      'Linear overestimator  for reformulation variable wVV - 1 of 2'
    eq_wVVik_over_2(a,i,k,h,t)      'Linear overestimator  for reformulation variable wVV - 2 of 2'
    eq_wVVik_under_1(a,i,k,h,t)     'Linear underestimator for reformulation variable wVV - 1 of 2'
    eq_wVVik_under_2(a,i,k,h,t)     'Linear underestimator for reformulation variable wVV - 2 of 2'

    eq_wEEik_over_1(a,i,k,h,t)      'Linear overestimator  for reformulation variable wEE - 1 of 2'
    eq_wEEik_over_2(a,i,k,h,t)      'Linear overestimator  for reformulation variable wEE - 2 of 2'
    eq_wEEik_under_1(a,i,k,h,t)     'Linear underestimator for reformulation variable wEE - 1 of 2'
    eq_wEEik_under_2(a,i,k,h,t)     'Linear underestimator for reformulation variable wEE - 2 of 2'

    eq_wEFik_over_1(a,i,k,h,t)      'Linear overestimator  for reformulation variable wEF - 1 of 2'
    eq_wEFik_over_2(a,i,k,h,t)      'Linear overestimator  for reformulation variable wEF - 2 of 2'
    eq_wEFik_under_1(a,i,k,h,t)     'Linear underestimator for reformulation variable wEF - 1 of 2'
    eq_wEFik_under_2(a,i,k,h,t)     'Linear underestimator for reformulation variable wEF - 2 of 2'

    eq_wFFik_over_1(a,i,k,h,t)      'Linear overestimator  for reformulation variable wFF - 1 of 2'
    eq_wFFik_over_2(a,i,k,h,t)      'Linear overestimator  for reformulation variable wFF - 2 of 2'
    eq_wFFik_under_1(a,i,k,h,t)     'Linear underestimator for reformulation variable wFF - 1 of 2'
    eq_wFFik_under_2(a,i,k,h,t)     'Linear underestimator for reformulation variable wFF - 2 of 2' ;

$ELSEIF.PW %PiecewiseType% == logtransform
Equations
    // Select piecewise segments for bilinear terms
    eq_lambdaVV(a,i,k,h,t)          'Restriction on segment selection for wVV'
    eq_breakVV(a,i,k,h,t)           'Relate wVV to its piecewise linear breakpoints'

    eq_lambdaEE(a,i,k,h,t)          'Restriction on segment selection for wEE'
    eq_breakEE(a,i,k,h,t)           'Relate wEE to its piecewise linear breakpoints'

    eq_lambdaEF(a,i,k,h,t)          'Restriction on segment selection for wEF'
    eq_breakEF(a,i,k,h,t)           'Relate wEF to its piecewise linear breakpoints'

    eq_lambdaFF(a,i,k,h,t)          'Restriction on segment selection for wFF'
    eq_breakFF(a,i,k,h,t)           'Relate wFF to its piecewise linear breakpoints'

    // Univariate log transformations
    eq_lnV_over(a,i,h,t,bp)         'Linear overestimator            for log transformation variable ln(V)'
    eq_lnV_under(a,i,h,t)           'Piecewise linear underestimator for log transformation variable ln(V)'

    eq_lnE_over(a,i,h,t,bp)         'Linear overestimator            for log transformation variable ln(E)'
    eq_lnE_under(a,i,h,t)           'Piecewise linear underestimator for log transformation variable ln(E)'

    eq_lnF_over(a,i,h,t,bp)         'Linear overestimator            for log transformation variable ln(F)'
    eq_lnF_under(a,i,h,t)           'Piecewise linear underestimator for log transformation variable ln(F)'

    eq_lnwVV_over(a,i,k,h,t,bp)     'Linear overestimator            for log transformation variable ln(wVV)'
    eq_lnwVV_under(a,i,k,h,t)       'Piecewise linear underestimator for log transformation variable ln(wVV)'

    eq_lnwEE_over(a,i,k,h,t,bp)     'Linear overestimator            for log transformation variable ln(wEE)'
    eq_lnwEE_under(a,i,k,h,t)       'Piecewise linear underestimator for log transformation variable ln(wEE)'

    eq_lnwEF_over(a,i,k,h,t,bp)     'Linear overestimator            for log transformation variable ln(wEF)'
    eq_lnwEF_under(a,i,k,h,t)       'Piecewise linear underestimator for log transformation variable ln(wEF)'

    eq_lnwFF_over(a,i,k,h,t,bp)     'Linear overestimator            for log transformation variable ln(wFF)'
    eq_lnwFF_under(a,i,k,h,t)       'Piecewise linear underestimator for log transformation variable ln(wFF)'

    // Bilinear terms
    eq_lnwVV(a,i,k,h,t)             'Define the linear relationship ln(wVV_ik) = ln(V_i) + ln(V_k)'
    eq_lnwEE(a,i,k,h,t)             'Define the linear relationship ln(wEE_ik) = ln(E_i) + ln(E_k)'
    eq_lnwEF(a,i,k,h,t)             'Define the linear relationship ln(wEF_ik) = ln(E_i) + ln(F_k)'
    eq_lnwFF(a,i,k,h,t)             'Define the linear relationship ln(wFF_ik) = ln(F_i) + ln(F_k)' ;

$ENDIF.PW
$ENDIF

* ------------------------------------------------------------------------------
*  Relaxation Equation Definitions - Univariate Terms
* ------------------------------------------------------------------------------
$IFTHEN %ConvexRelaxation% == piecewise
** Selection of Piecewise Segments (Univariate Terms)
* V
eq_lambdaV(NN(a,i),h,TT(t))$HH(a,h)..
    sum(bp, lambdaV(i,h,t,bp)) =e= 1 ;

eq_breakV(NN(a,i),h,TT(t))$HH(a,h)..
    V(i,h,t) =e= sum(bp, lambdaV(i,h,t,bp) * breakV(i,h,t,bp) ) ;

* E
eq_lambdaE(NN('AC',i),h,TT(t))$HH('AC',h)..
    sum(bp, lambdaE(i,h,t,bp)) =e= 1 ;

eq_breakE(NN('AC',i),h,TT(t))$HH('AC',h)..
    E(i,h,t) =e= sum(bp, lambdaE(i,h,t,bp) * breakE(i,h,t,bp) ) ;

* F
eq_lambdaF(NN('AC',i),h,TT(t))$HH('AC',h)..
    sum(bp, lambdaF(i,h,t,bp)) =e= 1 ;

eq_breakF(NN('AC',i),h,TT(t))$HH('AC',h)..
    F(i,h,t) =e= sum(bp, lambdaF(i,h,t,bp) * breakF(i,h,t,bp) ) ;

* wC
eq_lambdaC(c,TT(t))..
    sum(bp, lambdaC(c,t,bp)) =e= 1 ;

eq_breakC(c,TT(t))..
    wC(c,t) =e= sum(bp, lambdaC(c,t,bp) * breakC(c,t,bp) ) ;
$ENDIF


$IFTHEN NOT %ConvexRelaxation% == piecewise
** Overestimators - Linear or Quadratic Relaxation
* wVV
eq_wVVii_over(NN(a,i),h,TT(t))$HH(a,h)..
    wVV(i,i,h,t) =l=
          V.lo(i,h,t) * V(i,h,t)
        + V.up(i,h,t) * V(i,h,t)
        - V.lo(i,h,t) * V.up(i,h,t) ;

* wEE
eq_wEEii_over(NN('AC',i),h,TT(t))$HH('AC',h)..
    wEE(i,i,h,t) =l=
          E.lo(i,h,t) * E(i,h,t)
        + E.up(i,h,t) * E(i,h,t)
        - E.lo(i,h,t) * E.up(i,h,t) ;

* wFF
eq_wFFii_over(NN('AC',i),h,TT(t))$HH('AC',h)..
    wFF(i,i,h,t) =l=
          F.lo(i,h,t) * F(i,h,t)
        + F.up(i,h,t) * F(i,h,t)
        - F.lo(i,h,t) * F.up(i,h,t) ;

* wCC
eq_wCC_over(c,TT(t))..
    wCC(c,t) =l=
          wC.lo(c,t) * wC(c,t)
        + wC.up(c,t) * wC(c,t)
        - wC.lo(c,t) * wC.up(c,t) ;

$ELSE
** Overestimators - Piecewise Linear Relaxation
* wVV
eq_wVVii_over(NN(a,i),h,TT(t))$HH(a,h)..
    wVV(i,i,h,t) =l= sum(bp, lambdaV(i,h,t,bp) * sqr( breakV(i,h,t,bp) ) ) ;

* wEE
eq_wEEii_over(NN('AC',i),h,TT(t))$HH('AC',h)..
    wEE(i,i,h,t) =l= sum(bp, lambdaE(i,h,t,bp) * sqr( breakE(i,h,t,bp) ) ) ;

* wFF
eq_wFFii_over(NN('AC',i),h,TT(t))$HH('AC',h)..
    wFF(i,i,h,t) =l= sum(bp, lambdaF(i,h,t,bp) * sqr( breakF(i,h,t,bp) ) ) ;

* wCC
eq_wCC_over(c,TT(t))..
    wCC(c,t) =l= sum(bp, lambdaC(c,t,bp) * sqr( breakC(c,t,bp) ) ) ;

$ENDIF


$IFTHEN NOT %ConvexRelaxation% == quadratic
** Underestimators - Linear or Piecewise Linear Relaxation
* wVV
eq_wVVii_under(NN(a,i),h,TT(t),bp)$HH(a,h)..
    wVV(i,i,h,t) =g=
        2 * breakV(i,h,t,bp) * V(i,h,t)
        - sqr( breakV(i,h,t,bp) ) ;

* wEE
eq_wEEii_under(NN('AC',i),h,TT(t),bp)$HH('AC',h)..
    wEE(i,i,h,t) =g=
        2 * breakE(i,h,t,bp) * E(i,h,t)
        - sqr( breakE(i,h,t,bp) ) ;

* wFF
eq_wFFii_under(NN('AC',i),h,TT(t),bp)$HH('AC',h)..
    wFF(i,i,h,t) =g=
        2 * breakF(i,h,t,bp) * F(i,h,t)
        - sqr( breakF(i,h,t,bp) ) ;

* wCC
eq_wCC_under(c,TT(t),bp)..
    wCC(c,t) =g=
        2 * breakC(c,t,bp) * wC(c,t)
        - sqr( breakC(c,t,bp) ) ;

$ELSE
** Underestimators - Quadratic Relaxation
* wVV
eq_wVVii_under(NN(a,i),h,TT(t))$HH(a,h)..
    wVV(i,i,h,t) =g= sqr( V(i,h,t) ) ;

* wEE
eq_wEEii_under(NN('AC',i),h,TT(t))$HH('AC',h)..
    wEE(i,i,h,t) =g= sqr( E(i,h,t) ) ;

* wFF
eq_wFFii_under(NN('AC',i),h,TT(t))$HH('AC',h)..
    wFF(i,i,h,t) =g= sqr( F(i,h,t) ) ;

* wCC
eq_wCC_under(c,TT(t))..
    wCC(c,t) =g=  sqr( wC(c,t) ) ;

$ENDIF


* ------------------------------------------------------------------------------
*  Relaxation Equation Definitions - Bilinear Terms
* ------------------------------------------------------------------------------
$IFTHEN %ConvexRelaxation% == piecewise
** Selection of Piecewise Segments (Bilinear Terms)
$IFTHEN.PW %PiecewiseType% == mccormick
* wVV
eq_segsel_VV(ConnD(a,i,k),h,TT(t))$(HH(a,h) and not sameas(i,k))..
    sum((d,dd), yVV(i,k,h,t,d,dd)) =e= 1 ;

eq_seg_VVi(ConnD(a,i,k),h,TT(t))$(HH(a,h) and not sameas(i,k))..
    V(i,h,t) =e= sum((d,dd), segVVi(i,k,h,t,d,dd)) ;

eq_seg_VVk(ConnD(a,i,k),h,TT(t))$(HH(a,h) and not sameas(i,k))..
    V(k,h,t) =e= sum((d,dd), segVVk(i,k,h,t,d,dd)) ;

eq_seg_VVi_lo(ConnD(a,i,k),h,TT(t),d,dd)$(HH(a,h) and not sameas(i,k))..
    segVVi(i,k,h,t,d,dd) =g= seglimV(i,h,t,d,'lo') * yVV(i,k,h,t,d,dd) ;

eq_seg_VVi_up(ConnD(a,i,k),h,TT(t),d,dd)$(HH(a,h) and not sameas(i,k))..
    segVVi(i,k,h,t,d,dd) =l= seglimV(i,h,t,d,'up') * yVV(i,k,h,t,d,dd) ;

eq_seg_VVk_lo(ConnD(a,i,k),h,TT(t),d,dd)$(HH(a,h) and not sameas(i,k))..
    segVVk(i,k,h,t,d,dd) =g= seglimV(k,h,t,dd,'lo') * yVV(i,k,h,t,d,dd) ;

eq_seg_VVk_up(ConnD(a,i,k),h,TT(t),d,dd)$(HH(a,h) and not sameas(i,k))..
    segVVk(i,k,h,t,d,dd) =l= seglimV(k,h,t,dd,'up') * yVV(i,k,h,t,d,dd) ;

* wEE
eq_segsel_EE(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    sum((d,dd), yEE(i,k,h,t,d,dd)) =e= 1 ;

eq_seg_EEi(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    E(i,h,t) =e= sum((d,dd), segEEi(i,k,h,t,d,dd)) ;

eq_seg_EEk(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    E(k,h,t) =e= sum((d,dd), segEEk(i,k,h,t,d,dd)) ;

eq_seg_EEi_lo(ConnD('AC',i,k),h,TT(t),d,dd)$(HH('AC',h) and not sameas(i,k))..
    segEEi(i,k,h,t,d,dd) =g= seglimE(i,h,t,d,'lo') * yEE(i,k,h,t,d,dd) ;

eq_seg_EEi_up(ConnD('AC',i,k),h,TT(t),d,dd)$(HH('AC',h) and not sameas(i,k))..
    segEEi(i,k,h,t,d,dd) =l= seglimE(i,h,t,d,'up') * yEE(i,k,h,t,d,dd) ;

eq_seg_EEk_lo(ConnD('AC',i,k),h,TT(t),d,dd)$(HH('AC',h) and not sameas(i,k))..
    segEEk(i,k,h,t,d,dd) =g= seglimE(k,h,t,dd,'lo') * yEE(i,k,h,t,d,dd) ;

eq_seg_EEk_up(ConnD('AC',i,k),h,TT(t),d,dd)$(HH('AC',h) and not sameas(i,k))..
    segEEk(i,k,h,t,d,dd) =l= seglimE(k,h,t,dd,'up') * yEE(i,k,h,t,d,dd) ;

* wEF
eq_segsel_EF(ConnU('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    sum((d,dd), yEF(i,k,h,t,d,dd)) =e= 1 ;

eq_seg_EFi(ConnU('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    E(i,h,t) =e= sum((d,dd), segEFi(i,k,h,t,d,dd)) ;

eq_seg_EFk(ConnU('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    F(k,h,t) =e= sum((d,dd), segEFk(i,k,h,t,d,dd)) ;

eq_seg_EFi_lo(ConnU('AC',i,k),h,TT(t),d,dd)$(HH('AC',h) and not sameas(i,k))..
    segEFi(i,k,h,t,d,dd) =g= seglimE(i,h,t,d,'lo') * yEF(i,k,h,t,d,dd) ;

eq_seg_EFi_up(ConnU('AC',i,k),h,TT(t),d,dd)$(HH('AC',h) and not sameas(i,k))..
    segEFi(i,k,h,t,d,dd) =l= seglimE(i,h,t,d,'up') * yEF(i,k,h,t,d,dd) ;

eq_seg_EFk_lo(ConnU('AC',i,k),h,TT(t),d,dd)$(HH('AC',h) and not sameas(i,k))..
    segEFk(i,k,h,t,d,dd) =g= seglimF(k,h,t,dd,'lo') * yEF(i,k,h,t,d,dd) ;

eq_seg_EFk_up(ConnU('AC',i,k),h,TT(t),d,dd)$(HH('AC',h) and not sameas(i,k))..
    segEFk(i,k,h,t,d,dd) =l= seglimF(k,h,t,dd,'up') * yEF(i,k,h,t,d,dd) ;

* wFF
eq_segsel_FF(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    sum((d,dd), yFF(i,k,h,t,d,dd)) =e= 1 ;

eq_seg_FFi(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    F(i,h,t) =e= sum((d,dd), segFFi(i,k,h,t,d,dd)) ;

eq_seg_FFk(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    F(k,h,t) =e= sum((d,dd), segFFk(i,k,h,t,d,dd)) ;

eq_seg_FFi_lo(ConnD('AC',i,k),h,TT(t),d,dd)$(HH('AC',h) and not sameas(i,k))..
    segFFi(i,k,h,t,d,dd) =g= seglimF(i,h,t,d,'lo') * yFF(i,k,h,t,d,dd) ;

eq_seg_FFi_up(ConnD('AC',i,k),h,TT(t),d,dd)$(HH('AC',h) and not sameas(i,k))..
    segFFi(i,k,h,t,d,dd) =l= seglimF(i,h,t,d,'up') * yFF(i,k,h,t,d,dd) ;

eq_seg_FFk_lo(ConnD('AC',i,k),h,TT(t),d,dd)$(HH('AC',h) and not sameas(i,k))..
    segFFk(i,k,h,t,d,dd) =g= seglimF(k,h,t,dd,'lo') * yFF(i,k,h,t,d,dd) ;

eq_seg_FFk_up(ConnD('AC',i,k),h,TT(t),d,dd)$(HH('AC',h) and not sameas(i,k))..
    segFFk(i,k,h,t,d,dd) =l= seglimF(k,h,t,dd,'up') * yFF(i,k,h,t,d,dd) ;


$ELSEIF.PW %PiecewiseType% == logtransform
* wVV
eq_lambdaVV(ConnD(a,i,k),h,TT(t))$(HH(a,h) and not sameas(i,k))..
    sum(bp, lambdaVV(i,k,h,t,bp)) =e= 1 ;

eq_breakVV(ConnD(a,i,k),h,TT(t))$(HH(a,h) and not sameas(i,k))..
    wVV(i,k,h,t) =e= sum(bp, lambdaVV(i,k,h,t,bp) * breakVV(i,k,h,t,bp) ) ;

* wEE
eq_lambdaEE(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    sum(bp, lambdaEE(i,k,h,t,bp)) =e= 1 ;

eq_breakEE(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEE(i,k,h,t) =e= sum(bp, lambdaEE(i,k,h,t,bp) * breakEE(i,k,h,t,bp) ) ;

* wEF
eq_lambdaEF(ConnU('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    sum(bp, lambdaEF(i,k,h,t,bp)) =e= 1 ;

eq_breakEF(ConnU('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEF(i,k,h,t) =e= sum(bp, lambdaEF(i,k,h,t,bp) * breakEF(i,k,h,t,bp) ) ;

* wFF
eq_lambdaFF(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    sum(bp, lambdaFF(i,k,h,t,bp)) =e= 1 ;

eq_breakFF(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wFF(i,k,h,t) =e= sum(bp, lambdaFF(i,k,h,t,bp) * breakFF(i,k,h,t,bp) ) ;

** Log Transformation Relationships
* ln(wVV_ik) = ln(V_i) + ln(V_k)
eq_lnwVV(ConnD(a,i,k),h,TT(t))$(HH(a,h) and not sameas(i,k))..
    lnwVV(i,k,h,t) =e= lnV(i,h,t) + lnV(k,h,t) ;

* ln(wEE_ik) = ln(E_i) + ln(E_k)
eq_lnwEE(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    lnwEE(i,k,h,t) =e= lnE(i,h,t) + lnE(k,h,t) ;

* ln(wEF_ik) = ln(E_i) + ln(F_k)
eq_lnwEF(ConnU('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    lnwEF(i,k,h,t) =e= lnE(i,h,t) + lnF(k,h,t) ;

* ln(wFF_ik) = ln(F_i) + ln(F_k)
eq_lnwFF(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    lnwFF(i,k,h,t) =e= lnF(i,h,t) + lnF(k,h,t) ;

$ENDIF.PW
$ENDIF


$IFTHEN NOT %ConvexRelaxation% == piecewise
** Overestimators - Linear or Quadratic Relaxation
* wVV
eq_wVVik_over_1(ConnD(a,i,k),h,TT(t))$(HH(a,h) and not sameas(i,k))..
    wVV(i,k,h,t) =l=
          V.up(i,h,t) * V(k,h,t)
        + V(i,h,t)    * V.lo(k,h,t)
        - V.up(i,h,t) * V.lo(k,h,t) ;

eq_wVVik_over_2(ConnD(a,i,k),h,TT(t))$(HH(a,h) and not sameas(i,k))..
    wVV(i,k,h,t) =l=
          V.lo(i,h,t) * V(k,h,t)
        + V(i,h,t)    * V.up(k,h,t)
        - V.lo(i,h,t) * V.up(k,h,t) ;

* wEE
eq_wEEik_over_1(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEE(i,k,h,t) =l=
          E.up(i,h,t) * E(k,h,t)
        + E(i,h,t)    * E.lo(k,h,t)
        - E.up(i,h,t) * E.lo(k,h,t) ;

eq_wEEik_over_2(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEE(i,k,h,t) =l=
          E.lo(i,h,t) * E(k,h,t)
        + E(i,h,t)    * E.up(k,h,t)
        - E.lo(i,h,t) * E.up(k,h,t) ;

* wEF
eq_wEFik_over_1(ConnU('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEF(i,k,h,t) =l=
          E.up(i,h,t) * F(k,h,t)
        + E(i,h,t)    * F.lo(k,h,t)
        - E.up(i,h,t) * F.lo(k,h,t) ;

eq_wEFik_over_2(ConnU('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEF(i,k,h,t) =l=
          E.lo(i,h,t) * F(k,h,t)
        + E(i,h,t)    * F.up(k,h,t)
        - E.lo(i,h,t) * F.up(k,h,t) ;

* wFF
eq_wFFik_over_1(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wFF(i,k,h,t) =l=
          F.up(i,h,t) * F(k,h,t)
        + F(i,h,t)    * F.lo(k,h,t)
        - F.up(i,h,t) * F.lo(k,h,t) ;

eq_wFFik_over_2(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wFF(i,k,h,t) =l=
          F.lo(i,h,t) * F(k,h,t)
        + F(i,h,t)    * F.up(k,h,t)
        - F.lo(i,h,t) * F.up(k,h,t) ;


$ELSEIF %ConvexRelaxation% == piecewise
$IFTHEN.PW %PiecewiseType% == mccormick
** Overestimators - Piecewise Linear Relaxation w/ McCormick Inequalities
* wVV
eq_wVVik_over_1(ConnD(a,i,k),h,TT(t))$(HH(a,h) and not sameas(i,k))..
    wVV(i,k,h,t) =l= sum((d,dd),
          seglimV(i,h,t,d,'up') * segVVk(i,k,h,t,d,dd)
        + segVVi(i,k,h,t,d,dd)  * seglimV(k,h,t,dd,'lo')
        - seglimV(i,h,t,d,'up') * seglimV(k,h,t,dd,'lo') * yVV(i,k,h,t,d,dd)
    ) ;

eq_wVVik_over_2(ConnD(a,i,k),h,TT(t))$(HH(a,h) and not sameas(i,k))..
    wVV(i,k,h,t) =l= sum((d,dd),
          seglimV(i,h,t,d,'lo') * segVVk(i,k,h,t,d,dd)
        + segVVi(i,k,h,t,d,dd)  * seglimV(k,h,t,dd,'up')
        - seglimV(i,h,t,d,'lo') * seglimV(k,h,t,dd,'up') * yVV(i,k,h,t,d,dd)
    ) ;

* wEE
eq_wEEik_over_1(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEE(i,k,h,t) =l= sum((d,dd),
          seglimE(i,h,t,d,'up') * segEEk(i,k,h,t,d,dd)
        + segEEi(i,k,h,t,d,dd)  * seglimE(k,h,t,dd,'lo')
        - seglimE(i,h,t,d,'up') * seglimE(k,h,t,dd,'lo') * yEE(i,k,h,t,d,dd)
    ) ;

eq_wEEik_over_2(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEE(i,k,h,t) =l= sum((d,dd),
          seglimE(i,h,t,d,'up') * segEEk(i,k,h,t,d,dd)
        + segEEi(i,k,h,t,d,dd)  * seglimE(k,h,t,dd,'lo')
        - seglimE(i,h,t,d,'up') * seglimE(k,h,t,dd,'lo') * yEE(i,k,h,t,d,dd)
    ) ;

* wEF
eq_wEFik_over_1(ConnU('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEF(i,k,h,t) =l= sum((d,dd),
          seglimE(i,h,t,d,'up') * segEFk(i,k,h,t,d,dd)
        + segEFi(i,k,h,t,d,dd)  * seglimF(k,h,t,dd,'lo')
        - seglimE(i,h,t,d,'up') * seglimF(k,h,t,dd,'lo') * yEF(i,k,h,t,d,dd)
    ) ;

eq_wEFik_over_2(ConnU('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEF(i,k,h,t) =l= sum((d,dd),
          seglimE(i,h,t,d,'lo') * segEFk(i,k,h,t,d,dd)
        + segEFi(i,k,h,t,d,dd)  * seglimF(k,h,t,dd,'up')
        - seglimE(i,h,t,d,'lo') * seglimF(k,h,t,dd,'up') * yEF(i,k,h,t,d,dd)
    ) ;

* wFF
eq_wFFik_over_1(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wFF(i,k,h,t) =l= sum((d,dd),
          seglimF(i,h,t,d,'up') * segFFk(i,k,h,t,d,dd)
        + segFFi(i,k,h,t,d,dd)  * seglimF(k,h,t,dd,'lo')
        - seglimF(i,h,t,d,'up') * seglimF(k,h,t,dd,'lo') * yFF(i,k,h,t,d,dd)
    ) ;

eq_wFFik_over_2(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wFF(i,k,h,t) =l= sum((d,dd),
          seglimF(i,h,t,d,'lo') * segFFk(i,k,h,t,d,dd)
        + segFFi(i,k,h,t,d,dd)  * seglimF(k,h,t,dd,'up')
        - seglimF(i,h,t,d,'lo') * seglimF(k,h,t,dd,'up') * yFF(i,k,h,t,d,dd)
    ) ;


$ELSEIF.PW %PiecewiseType% == logtransform
** Overestimators - Piecewise Linear Relaxation w/ Log Transformation
* ln(V)
eq_lnV_over(NN(a,i),h,TT(t),bp)$HH(a,h)..
    lnV(i,h,t) =l=
        log( breakV(i,h,t,bp) ) + V(i,h,t) / breakV(i,h,t,bp) - 1 ;

* ln(E)
eq_lnE_over(NN('AC',i),h,TT(t),bp)$HH('AC',h)..
    lnE(i,h,t) =l=
        log( breakE(i,h,t,bp) ) + E(i,h,t) / breakE(i,h,t,bp) - 1 ;

* ln(F)
eq_lnF_over(NN('AC',i),h,TT(t),bp)$HH('AC',h)..
    lnF(i,h,t) =l=
        log( breakF(i,h,t,bp) ) + F(i,h,t) / breakF(i,h,t,bp) - 1 ;

* ln(wVV)
eq_lnwVV_over(ConnD(a,i,k),h,TT(t),bp)$(HH(a,h) and not sameas(i,k))..
    lnwVV(i,k,h,t) =l=
        log( breakVV(i,k,h,t,bp) ) + wVV(i,k,h,t) / breakVV(i,k,h,t,bp) - 1 ;

* ln(wEE)
eq_lnwEE_over(ConnD('AC',i,k),h,TT(t),bp)$(HH('AC',h) and not sameas(i,k))..
    lnwEE(i,k,h,t) =l=
        log( breakEE(i,k,h,t,bp) ) + wEE(i,k,h,t) / breakEE(i,k,h,t,bp) - 1 ;

* ln(wEF)
eq_lnwEF_over(ConnU('AC',i,k),h,TT(t),bp)$(HH('AC',h) and not sameas(i,k))..
    lnwEF(i,k,h,t) =l=
        log( breakEF(i,k,h,t,bp) ) + wEF(i,k,h,t) / breakEF(i,k,h,t,bp) - 1 ;

* ln(wFF)
eq_lnwFF_over(ConnD('AC',i,k),h,TT(t),bp)$(HH('AC',h) and not sameas(i,k))..
    lnwFF(i,k,h,t) =l=
        log( breakFF(i,k,h,t,bp) ) + wFF(i,k,h,t) / breakFF(i,k,h,t,bp) - 1 ;

$ENDIF.PW
$ENDIF


$IFTHEN %ConvexRelaxation% == linear
** Underestimators - Linear Relaxation
* wVV
eq_wVVik_under_1(ConnD(a,i,k),h,TT(t))$(HH(a,h) and not sameas(i,k))..
    wVV(i,k,h,t) =g=
          V.lo(i,h,t) * V(k,h,t)
        + V(i,h,t)    * V.lo(k,h,t)
        - V.lo(i,h,t) * V.lo(k,h,t) ;

eq_wVVik_under_2(ConnD(a,i,k),h,TT(t))$(HH(a,h) and not sameas(i,k))..
    wVV(i,k,h,t) =g=
          V.up(i,h,t) * V(k,h,t)
        + V(i,h,t)    * V.up(k,h,t)
        - V.up(i,h,t) * V.up(k,h,t) ;

* wEE
eq_wEEik_under_1(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEE(i,k,h,t) =g=
          E.lo(i,h,t) * E(k,h,t)
        + E(i,h,t)    * E.lo(k,h,t)
        - E.lo(i,h,t) * E.lo(k,h,t) ;

eq_wEEik_under_2(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEE(i,k,h,t) =g=
          E.up(i,h,t) * E(k,h,t)
        + E(i,h,t)    * E.up(k,h,t)
        - E.up(i,h,t) * E.up(k,h,t) ;

* wEF
eq_wEFik_under_1(ConnU('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEF(i,k,h,t) =g=
          E.lo(i,h,t) * F(k,h,t)
        + E(i,h,t)    * F.lo(k,h,t)
        - E.lo(i,h,t) * F.lo(k,h,t) ;

eq_wEFik_under_2(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEF(i,k,h,t) =g=
          E.up(i,h,t) * F(k,h,t)
        + E(i,h,t)    * F.up(k,h,t)
        - E.up(i,h,t) * F.up(k,h,t) ;

* wFF
eq_wFFik_under_1(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wFF(i,k,h,t) =g=
          F.lo(i,h,t) * F(k,h,t)
        + F(i,h,t)    * F.lo(k,h,t)
        - F.lo(i,h,t) * F.lo(k,h,t) ;

eq_wFFik_under_2(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wFF(i,k,h,t) =g=
          F.up(i,h,t) * F(k,h,t)
        + F(i,h,t)    * F.up(k,h,t)
        - F.up(i,h,t) * F.up(k,h,t) ;


$ELSEIF %ConvexRelaxation% == quadratic
** Underestimators - Quadratic Relaxation
* wVV
eq_wVVik_under(ConnD(a,i,k),h,TT(t))$(HH(a,h) and not sameas(i,k))..
    wVV(i,k,h,t) =g= V(i,h,t) * V(k,h,t) ;

* wEE
eq_wEEik_under(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEE(i,k,h,t) =g= E(i,h,t) * E(k,h,t) ;

* wEF
eq_wEFik_under(ConnU('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEF(i,k,h,t) =g= E(i,h,t) * F(k,h,t) ;

* wFF
eq_wFFik_under(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wFF(i,k,h,t) =g= F(i,h,t) * F(k,h,t) ;


$ELSEIF %ConvexRelaxation% == piecewise
$IFTHEN.PW %PiecewiseType% == mccormick
** Underestimators - Piecewise Linear Relaxation w/ McCormick Inequalities
* wVV
eq_wVVik_under_1(ConnD(a,i,k),h,TT(t))$(HH(a,h) and not sameas(i,k))..
    wVV(i,k,h,t) =g= sum((d,dd),
          seglimV(i,h,t,d,'lo') * segVVk(i,k,h,t,d,dd)
        + segVVi(i,k,h,t,d,dd)  * seglimV(k,h,t,dd,'lo')
        - seglimV(i,h,t,d,'lo') * seglimV(k,h,t,dd,'lo') * yVV(i,k,h,t,d,dd)
    ) ;

eq_wVVik_under_2(ConnD(a,i,k),h,TT(t))$(HH(a,h) and not sameas(i,k))..
    wVV(i,k,h,t) =g= sum((d,dd),
          seglimV(i,h,t,d,'up') * segVVk(i,k,h,t,d,dd)
        + segVVi(i,k,h,t,d,dd)  * seglimV(k,h,t,dd,'up')
        - seglimV(i,h,t,d,'up') * seglimV(k,h,t,dd,'up') * yVV(i,k,h,t,d,dd)
    ) ;

* wEE
eq_wEEik_under_1(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEE(i,k,h,t) =g= sum((d,dd),
          seglimE(i,h,t,d,'lo') * segEEk(i,k,h,t,d,dd)
        + segEEi(i,k,h,t,d,dd)  * seglimE(k,h,t,dd,'lo')
        - seglimE(i,h,t,d,'lo') * seglimE(k,h,t,dd,'lo') * yEE(i,k,h,t,d,dd)
    ) ;

eq_wEEik_under_2(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEE(i,k,h,t) =g= sum((d,dd),
          seglimE(i,h,t,d,'up') * segEEk(i,k,h,t,d,dd)
        + segEEi(i,k,h,t,d,dd)  * seglimE(k,h,t,dd,'up')
        - seglimE(i,h,t,d,'up') * seglimE(k,h,t,dd,'up') * yEE(i,k,h,t,d,dd)
    ) ;

* wEF
eq_wEFik_under_1(ConnU('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEF(i,k,h,t) =g= sum((d,dd),
          seglimE(i,h,t,d,'lo') * segEFk(i,k,h,t,d,dd)
        + segEFi(i,k,h,t,d,dd)  * seglimF(k,h,t,dd,'lo')
        - seglimE(i,h,t,d,'lo') * seglimF(k,h,t,dd,'lo') * yEF(i,k,h,t,d,dd)
    ) ;

eq_wEFik_under_2(ConnU('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wEF(i,k,h,t) =g= sum((d,dd),
          seglimE(i,h,t,d,'up') * segEFk(i,k,h,t,d,dd)
        + segEFi(i,k,h,t,d,dd)  * seglimF(k,h,t,dd,'up')
        - seglimE(i,h,t,d,'up') * seglimF(k,h,t,dd,'up') * yEF(i,k,h,t,d,dd)
    ) ;

* wFF
eq_wFFik_under_1(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wFF(i,k,h,t) =g= sum((d,dd),
          seglimF(i,h,t,d,'lo') * segFFk(i,k,h,t,d,dd)
        + segFFi(i,k,h,t,d,dd)  * seglimF(k,h,t,dd,'lo')
        - seglimF(i,h,t,d,'lo') * seglimF(k,h,t,dd,'lo') * yFF(i,k,h,t,d,dd)
    ) ;

eq_wFFik_under_2(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    wFF(i,k,h,t) =g= sum((d,dd),
          seglimF(i,h,t,d,'up') * segFFk(i,k,h,t,d,dd)
        + segFFi(i,k,h,t,d,dd)  * seglimF(k,h,t,dd,'up')
        - seglimF(i,h,t,d,'up') * seglimF(k,h,t,dd,'up') * yFF(i,k,h,t,d,dd)
    ) ;


$ELSEIF.PW %PiecewiseType% == logtransform
** Underestimators - Piecewise Linear Relaxation w/ Log Transformation
* ln(V)
eq_lnV_under(NN(a,i),h,TT(t))$HH(a,h)..
    lnV(i,h,t) =g= sum(bp, lambdaV(i,h,t,bp) * log( breakV(i,h,t,bp) ) ) ;

* ln(E)
eq_lnE_under(NN('AC',i),h,TT(t))$HH('AC',h)..
    lnE(i,h,t) =g= sum(bp, lambdaE(i,h,t,bp) * log( breakE(i,h,t,bp) ) ) ;

* ln(F)
eq_lnF_under(NN('AC',i),h,TT(t))$HH('AC',h)..
    lnF(i,h,t) =g= sum(bp, lambdaF(i,h,t,bp) * log( breakF(i,h,t,bp) ) ) ;

* ln(wVV)
eq_lnwVV_under(ConnD(a,i,k),h,TT(t))$(HH(a,h) and not sameas(i,k))..
    lnwVV(i,k,h,t) =g=
        sum(bp, lambdaVV(i,k,h,t,bp) * log( breakVV(i,k,h,t,bp) ) ) ;

* ln(wEE)
eq_lnwEE_under(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    lnwEE(i,k,h,t) =g=
        sum(bp, lambdaEE(i,k,h,t,bp) * log( breakEE(i,k,h,t,bp) ) ) ;

* ln(wEF)
eq_lnwEF_under(ConnU('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    lnwEF(i,k,h,t) =g=
        sum(bp, lambdaEF(i,k,h,t,bp) * log( breakEF(i,k,h,t,bp) ) ) ;

* ln(wFF)
eq_lnwFF_under(ConnD('AC',i,k),h,TT(t))$(HH('AC',h) and not sameas(i,k))..
    lnwFF(i,k,h,t) =g=
        sum(bp, lambdaFF(i,k,h,t,bp) * log( breakFF(i,k,h,t,bp) ) ) ;

$ENDIF.PW
$ENDIF
