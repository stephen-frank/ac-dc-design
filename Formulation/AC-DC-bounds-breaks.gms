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

$TITLE Reusable code for computing bounds and breakpoints for AC-DC model
$ONTEXT

    Filename : AC-DC-bounds-breaks.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    This file computes bounds on reformulation variables and/or breakpoints for
    linear and piecewise linear relaxations, as well as provides some other
    low-level reusable utility code. It is designed to be included using
    $BATINCLUDE; the procedures performed are controlled by the $BATINCLUDE
    parameters and compile-time flags.

    Batch Arguments:
    %1-%4   Keywords 'wbounds', 'wbounds2', 'breakpoints', and/or 'vlevels'
            (Any or all keywords may be included, in any order)

    In addition, this script relies on the following compile-time parameters to
    be set in the parent file (which are documented in AC-DC-defaults.gms):
        ConvexRelaxation    Specify type of convex relaxation
        PiecewiseType       Specify piecewise relaxation subtype

    Notes:
    1.  Passing 'wbounds' recomputes bounds on the reformulation variables 'w'
        using the bounds on their parent variables.
    2.  Passing 'wbounds2' recomputes bounds on the reformulation variables 'w'
        using the bounds on their parent variables, but does not update the
        bounds on 'w' if the existing bounds are tighter.
    3.  Passing 'breakpoints' recomputes breakpoints for linear and piecewise
        linear relaxations based on the bounds on the variables to be relaxed.
    4.  Passing 'vlevels' resets voltage levels to the midpoint between their
        bounds and the 'w' levels to match.

$OFFTEXT

* =============================================================================
*  Setup
* =============================================================================
** Process Arugments
* Arguments 1-2 = select what to do
$IFi %1 == wbounds              $SETLOCAL wbounds           yes
$IFi %2 == wbounds              $SETLOCAL wbounds           yes
$IFi %3 == wbounds              $SETLOCAL wbounds           yes
$IFi %4 == wbounds              $SETLOCAL wbounds           yes

$IFi %1 == wbounds2             $SETLOCAL wbounds2          yes
$IFi %2 == wbounds2             $SETLOCAL wbounds2          yes
$IFi %3 == wbounds2             $SETLOCAL wbounds2          yes
$IFi %4 == wbounds2             $SETLOCAL wbounds2          yes

$IFi %1 == breakpoints          $SETLOCAL breakpoints       yes
$IFi %2 == breakpoints          $SETLOCAL breakpoints       yes
$IFi %3 == breakpoints          $SETLOCAL breakpoints       yes
$IFi %4 == breakpoints          $SETLOCAL breakpoints       yes

$IFi %1 == vlevels              $SETLOCAL vlevels           yes
$IFi %2 == vlevels              $SETLOCAL vlevels           yes
$IFi %3 == vlevels              $SETLOCAL vlevels           yes
$IFi %4 == vlevels              $SETLOCAL vlevels           yes

** Defaults
* For things that the batch arguments didn't set
$IF NOT SET wbounds             $SETLOCAL wbounds           no
$IF NOT SET wbounds2            $SETLOCAL wbounds2          no
$IF NOT SET breakpoints         $SETLOCAL breakpoints       no
$IF NOT SET vlevels             $SETLOCAL vlevels           no


* =============================================================================
*  Reformulation Variable Bounds
* =============================================================================
$IFTHEN.bounds %wbounds% == yes
* Voltage magnitude * Voltage magnitude
loop( (a),
    wVV.lo(i,k,h,t)$(ConnD(a,i,k) and (HX('DC',h) or HX('AC',h))) =
        min( V.lo(i,h,t) * V.lo(k,h,t), V.lo(i,h,t) * V.up(k,h,t),
             V.up(i,h,t) * V.lo(k,h,t), V.up(i,h,t) * V.up(k,h,t)  );
    wVV.up(i,k,h,t)$(ConnD(a,i,k) and (HX('DC',h) or HX('AC',h))) =
        max( V.lo(i,h,t) * V.lo(k,h,t), V.lo(i,h,t) * V.up(k,h,t),
             V.up(i,h,t) * V.lo(k,h,t), V.up(i,h,t) * V.up(k,h,t)  );
);

* Voltage real component * Voltage real component
wEE.lo(i,k,h,t)$(ConnD('AC',i,k) and HX('AC',h)) =
    min( E.lo(i,h,t) * E.lo(k,h,t), E.lo(i,h,t) * E.up(k,h,t),
         E.up(i,h,t) * E.lo(k,h,t), E.up(i,h,t) * E.up(k,h,t)  );
wEE.up(i,k,h,t)$(ConnD('AC',i,k) and HX('AC',h)) =
    max( E.lo(i,h,t) * E.lo(k,h,t), E.lo(i,h,t) * E.up(k,h,t),
         E.up(i,h,t) * E.lo(k,h,t), E.up(i,h,t) * E.up(k,h,t)  );

* Voltage real component * Voltage imaginary component
wEF.lo(i,k,h,t)$(ConnU('AC',i,k) and HX('AC',h)) =
    min( E.lo(i,h,t) * F.lo(k,h,t), E.lo(i,h,t) * F.up(k,h,t),
         E.up(i,h,t) * F.lo(k,h,t), E.up(i,h,t) * F.up(k,h,t)  );
wEF.up(i,k,h,t)$(ConnU('AC',i,k) and HX('AC',h)) =
    max( E.lo(i,h,t) * F.lo(k,h,t), E.lo(i,h,t) * F.up(k,h,t),
         E.up(i,h,t) * F.lo(k,h,t), E.up(i,h,t) * F.up(k,h,t)  );

* Voltage imaginary component * Voltage imaginary component
wFF.lo(i,k,h,t)$(ConnD('AC',i,k) and HX('AC',h)) =
    min( F.lo(i,h,t) * F.lo(k,h,t), F.lo(i,h,t) * F.up(k,h,t),
         F.up(i,h,t) * F.lo(k,h,t), F.up(i,h,t) * F.up(k,h,t)  );
wFF.up(i,k,h,t)$(ConnD('AC',i,k) and HX('AC',h)) =
    max( F.lo(i,h,t) * F.lo(k,h,t), F.lo(i,h,t) * F.up(k,h,t),
         F.up(i,h,t) * F.lo(k,h,t), F.up(i,h,t) * F.up(k,h,t)  );

* Converter output power
wC.lo(c,t) = sum(h, ConvData2(c, h, 'Pomin'));
wC.up(c,t) = sum(h, ConvData2(c, h, 'Pomax'));

* Converter output power squared
wCC.lo(c,t) = min( sqr(wC.lo(c,t)), sqr(wC.up(c,t)) )$(
    sign(wC.lo(c,t)) = sign(wC.up(c,t)) );
wCC.up(c,t) = max( sqr(wC.lo(c,t)), sqr(wC.up(c,t)) );
$ENDIF.bounds

$IFTHEN.bounds %wbounds2% == yes
* Voltage magnitude * Voltage magnitude
loop( (a),
    wVV.lo(i,k,h,t)$(ConnD(a,i,k) and (HX('DC',h) or HX('AC',h))) =
        max( wVV.lo(i,k,h,t),
           min( V.lo(i,h,t) * V.lo(k,h,t), V.lo(i,h,t) * V.up(k,h,t),
                V.up(i,h,t) * V.lo(k,h,t), V.up(i,h,t) * V.up(k,h,t)  )
        );
    wVV.up(i,k,h,t)$(ConnD(a,i,k) and (HX('DC',h) or HX('AC',h))) =
        min( wVV.up(i,k,h,t),
            max( V.lo(i,h,t) * V.lo(k,h,t), V.lo(i,h,t) * V.up(k,h,t),
                 V.up(i,h,t) * V.lo(k,h,t), V.up(i,h,t) * V.up(k,h,t)  )
        );
);

* Voltage real component * Voltage real component
wEE.lo(i,k,h,t)$(ConnD('AC',i,k) and HX('AC',h)) =
    max( wEE.lo(i,k,h,t),
        min( E.lo(i,h,t) * E.lo(k,h,t), E.lo(i,h,t) * E.up(k,h,t),
             E.up(i,h,t) * E.lo(k,h,t), E.up(i,h,t) * E.up(k,h,t)  )
    );
wEE.up(i,k,h,t)$(ConnD('AC',i,k) and HX('AC',h)) =
    min( wEE.up(i,k,h,t),
        max( E.lo(i,h,t) * E.lo(k,h,t), E.lo(i,h,t) * E.up(k,h,t),
             E.up(i,h,t) * E.lo(k,h,t), E.up(i,h,t) * E.up(k,h,t)  )
    );

* Voltage real component * Voltage imaginary component
wEF.lo(i,k,h,t)$(ConnU('AC',i,k) and HX('AC',h)) =
    max( wEF.lo(i,k,h,t),
        min( E.lo(i,h,t) * F.lo(k,h,t), E.lo(i,h,t) * F.up(k,h,t),
             E.up(i,h,t) * F.lo(k,h,t), E.up(i,h,t) * F.up(k,h,t)  )
    );
wEF.up(i,k,h,t)$(ConnU('AC',i,k) and HX('AC',h)) =
    min( wEF.up(i,k,h,t),
        max( E.lo(i,h,t) * F.lo(k,h,t), E.lo(i,h,t) * F.up(k,h,t),
             E.up(i,h,t) * F.lo(k,h,t), E.up(i,h,t) * F.up(k,h,t)  )
    );

* Voltage imaginary component * Voltage imaginary component
wFF.lo(i,k,h,t)$(ConnD('AC',i,k) and HX('AC',h)) =
    max( wFF.lo(i,k,h,t),
        min( F.lo(i,h,t) * F.lo(k,h,t), F.lo(i,h,t) * F.up(k,h,t),
             F.up(i,h,t) * F.lo(k,h,t), F.up(i,h,t) * F.up(k,h,t)  )
    );
wFF.up(i,k,h,t)$(ConnD('AC',i,k) and HX('AC',h)) =
    min( wFF.up(i,k,h,t),
        max( F.lo(i,h,t) * F.lo(k,h,t), F.lo(i,h,t) * F.up(k,h,t),
             F.up(i,h,t) * F.lo(k,h,t), F.up(i,h,t) * F.up(k,h,t)  )
    );

* Converter output power
wC.lo(c,t) = max( wC.lo(c,t), sum(h, ConvData2(c, h, 'Pomin')) );
wC.up(c,t) = min( wC.up(c,t), sum(h, ConvData2(c, h, 'Pomax')) );

* Converter output power squared
wCC.lo(c,t) = max(
    wCC.lo(c,t),
    min( sqr(wC.lo(c,t)), sqr(wC.up(c,t)) )$(
        sign(wC.lo(c,t)) = sign(wC.up(c,t)) )
    );
wCC.up(c,t) = min(
    wCC.up(c,t),
    max( sqr(wC.lo(c,t)), sqr(wC.up(c,t)) )
    );
$ENDIF.bounds

* =============================================================================
*  Linear/Piecewise Linear Breakpoints
* =============================================================================
$IFTHEN.breaks %breakpoints% == yes

$IFTHEN NOT %ConvexRelaxation% == quadratic
** Univariate Breakpoints
loop( (a),
    breakV(i,h,t,bp)$(NN(a,i) and (HX('DC',h) or HX('AC',h))) =
        V.lo(i,h,t) + (V.up(i,h,t) - V.lo(i,h,t))*(ord(bp)-1)/(card(BP)-1) ;
);

breakE(i,h,t,bp)$(NN('AC',i) and HX('AC',h)) =
    E.lo(i,h,t) + (E.up(i,h,t) - E.lo(i,h,t))*(ord(bp)-1)/(card(BP)-1) ;

breakF(i,h,t,bp)$(NN('AC',i) and HX('AC',h)) =
    F.lo(i,h,t) + (F.up(i,h,t) - F.lo(i,h,t))*(ord(bp)-1)/(card(BP)-1) ;

breakC(c,t,bp) =
    wC.lo(c,t) + (wC.up(c,t) - wC.lo(c,t))*(ord(bp)-1)/(card(BP)-1) ;
$ENDIF


$IFTHEN %ConvexRelaxation% == piecewise
$IFTHEN.pw %PiecewiseType% == mccormick
** Bilinear Breakpoints - Using McCormick Inequalities
* Derive segment limits from previously computed breakpoints
seglimV(i,h,t,d,'lo') = sum(bp$(ord(bp) = ord(d)), breakV(i,h,t,bp)) ;
seglimV(i,h,t,d,'up') = sum(bp$(ord(bp) = ord(d)+1), breakV(i,h,t,bp)) ;

seglimE(i,h,t,d,'lo') = sum(bp$(ord(bp) = ord(d)), breakE(i,h,t,bp)) ;
seglimE(i,h,t,d,'up') = sum(bp$(ord(bp) = ord(d)+1), breakE(i,h,t,bp)) ;

seglimF(i,h,t,d,'lo') = sum(bp$(ord(bp) = ord(d)), breakF(i,h,t,bp)) ;
seglimF(i,h,t,d,'up') = sum(bp$(ord(bp) = ord(d)+1), breakF(i,h,t,bp)) ;

* Set bounds on piecewise linear segment variables
segVVi.lo(i,k,h,t,d,dd) = min( seglimV(i,h,t,d,'lo'), 0 );
segVVi.up(i,k,h,t,d,dd) = max( seglimV(i,h,t,d,'up'), 0 );
segVVk.lo(i,k,h,t,d,dd) = min( seglimV(k,h,t,dd,'lo'), 0 );
segVVk.up(i,k,h,t,d,dd) = max( seglimV(k,h,t,dd,'up'), 0 );

segEEi.lo(i,k,h,t,d,dd) = min( seglimE(i,h,t,d,'lo'), 0 );
segEEi.up(i,k,h,t,d,dd) = max( seglimE(i,h,t,d,'up'), 0 );
segEEk.lo(i,k,h,t,d,dd) = min( seglimE(k,h,t,dd,'lo'), 0 );
segEEk.up(i,k,h,t,d,dd) = max( seglimE(k,h,t,dd,'up'), 0 );

segEFi.lo(i,k,h,t,d,dd) = min( seglimE(i,h,t,d,'lo'), 0 );
segEFi.up(i,k,h,t,d,dd) = max( seglimE(i,h,t,d,'up'), 0 );
segEFk.lo(i,k,h,t,d,dd) = min( seglimF(k,h,t,dd,'lo'), 0 );
segEFk.up(i,k,h,t,d,dd) = max( seglimF(k,h,t,dd,'up'), 0 );

segFFi.lo(i,k,h,t,d,dd) = min( seglimF(i,h,t,d,'lo'), 0 );
segFFi.up(i,k,h,t,d,dd) = max( seglimF(i,h,t,d,'up'), 0 );
segFFk.lo(i,k,h,t,d,dd) = min( seglimF(k,h,t,dd,'lo'), 0 );
segFFk.up(i,k,h,t,d,dd) = max( seglimF(k,h,t,dd,'up'), 0 );


$ELSEIF.pw %PiecewiseType% == logtransform
** Bilinear Breakpoints - Using Logarithmic Transform
* Bounds on log transformation variables
loop( (a),
    lnV.lo(i,h,t)$(NN(a,i) and (HX('DC',h) or HX('AC',h))) = log( V.lo(i,h,t) ) ;
    lnV.up(i,h,t)$(NN(a,i) and (HX('DC',h) or HX('AC',h))) = log( V.up(i,h,t) ) ;
);

lnE.lo(i,h,t)$(NN('AC',i) and HX('AC',h)) = log( E.lo(i,h,t) ) ;
lnE.up(i,h,t)$(NN('AC',i) and HX('AC',h)) = log( E.up(i,h,t) ) ;

lnF.lo(i,h,t)$(NN('AC',i) and HX('AC',h)) = log( F.lo(i,h,t) ) ;
lnF.up(i,h,t)$(NN('AC',i) and HX('AC',h)) = log( F.up(i,h,t) ) ;

loop( (a),
    lnwVV.lo(i,k,h,t)$(ConnD(a,i,k) and (HX('DC',h) or HX('AC',h))) = log( wVV.lo(i,k,h,t) ) ;
    lnwVV.up(i,k,h,t)$(ConnD(a,i,k) and (HX('DC',h) or HX('AC',h))) = log( wVV.up(i,k,h,t) ) ;
);

lnwEE.lo(i,k,h,t)$(ConnD('AC',i,k) and HX('AC',h)) = log( wEE.lo(i,k,h,t) ) ;
lnwEE.up(i,k,h,t)$(ConnD('AC',i,k) and HX('AC',h)) = log( wEE.up(i,k,h,t) ) ;

lnwEF.lo(i,k,h,t)$(ConnU('AC',i,k) and HX('AC',h)) = log( wEF.lo(i,k,h,t) ) ;
lnwEF.up(i,k,h,t)$(ConnU('AC',i,k) and HX('AC',h)) = log( wEF.up(i,k,h,t) ) ;

lnwFF.lo(i,k,h,t)$(ConnD('AC',i,k) and HX('AC',h)) = log( wFF.lo(i,k,h,t) ) ;
lnwFF.up(i,k,h,t)$(ConnD('AC',i,k) and HX('AC',h)) = log( wFF.up(i,k,h,t) ) ;

* Compute breakpoints
loop( (a),
    breakVV(i,k,h,t,bp)$(ConnD(a,i,k) and (HX('DC',h) or HX('AC',h))) =
        wVV.lo(i,k,h,t)
        + (wVV.up(i,k,h,t) - wVV.lo(i,k,h,t))*(ord(bp)-1)/(card(BP)-1) ;
);

breakEE(i,k,h,t,bp)$(ConnD('AC',i,k) and HX('AC',h)) =
    wEE.lo(i,k,h,t)
    + (wEE.up(i,k,h,t) - wEE.lo(i,k,h,t))*(ord(bp)-1)/(card(BP)-1) ;

breakEF(i,k,h,t,bp)$(ConnU('AC',i,k) and HX('AC',h)) =
    wEF.lo(i,k,h,t)
    + (wEF.up(i,k,h,t) - wEF.lo(i,k,h,t))*(ord(bp)-1)/(card(BP)-1) ;

breakFF(i,k,h,t,bp)$(ConnD('AC',i,k) and HX('AC',h)) =
    wFF.lo(i,k,h,t)
    + (wFF.up(i,k,h,t) - wFF.lo(i,k,h,t))*(ord(bp)-1)/(card(BP)-1) ;

$ENDIF.pw
$ENDIF

$ENDIF.breaks

* =============================================================================
*  Voltage Variable Levels
* =============================================================================
$IFTHEN.levs %vlevels% == yes
* Voltages (flat start)
loop( (a),
    V.l(i,h,t)$(NN(a,i) and (HX('DC',h) or HX('AC',h))) = ( V.lo(i,h,t) + V.up(i,h,t) ) / 2;
);

E.l(i,h,t)$(NN('AC',i) and HX('AC',h)) = V.lo(i,h,t);
F.l(i,h,t)$(NN('AC',i) and HX('AC',h)) = 0;

* Reformulation Variables
loop( (a),
    wVV.l(i,k,h,t)$(ConnD(a,i,k) and (HX('DC',h) or HX('AC',h))) = V.l(i,h,t) * V.l(k,h,t);
);

wEE.l(i,k,h,t)$(ConnD('AC',i,k) and HX('AC',h)) = E.l(i,h,t) * E.l(k,h,t);
wEF.l(i,k,h,t)$(ConnU('AC',i,k) and HX('AC',h)) = E.l(i,h,t) * F.l(k,h,t);
wFF.l(i,k,h,t)$(ConnD('AC',i,k) and HX('AC',h)) = F.l(i,h,t) * F.l(k,h,t);
$ENDIF.levs
