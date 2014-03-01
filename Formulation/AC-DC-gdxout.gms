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

$TITLE Flexible GDX output script for AC-DC model and algorithm
$ONTEXT

    Filename : AC-DC-gdxout.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    This file exports sets, parameters, equations, etc. for the AC-DC model to
    a user-specified GDX file. It is designed to be included using $BATINCLUDE;
    what exactly is exported is controlled by the $BATINCLUDE parameters.

    Batch Arguments:
    %1      Name of output file
    %2      Which algorithm is being used: bigM, bound_tightening, main, test
            (Exports some history and tracking info specific to the algorithm)
    %3      Which problem is being worked with: LBP, UBP, UBPr, Monolith, or All
            (Used to select equation set + certain data and variables)
    %4-%7   Keywords 'sets', 'data', 'variables', 'equations'
            (Any number of keywords may be included, in any order)

    In addition, this script relies on the following compile-time parameters to
    be set in the parent file (which are documented in AC-DC-defaults.gms):
        ConvexRelaxation    Specify type of convex relaxation
        ExtraCuts           Whether additional feasibility cuts are in use
        PiecewiseType       Specify piecewise relaxation subtype

    Notes:
    1.  The keywords in arguments %4 up to %7 determine what gets exported:
        sets & aliases, parameters (data), variables, and/or equations. Any or
        none of these keywords may be specified; the script sorts them out and
        sets appropriate flags.
    2.  Argument %2 controls whether to export some algorithm-specific
        information, such as convergence data, iteration history, etc. The four
        options are:
            bigM: Exports stuff related to initial big-M calculations
            bound_tightening: Exports bound tightening history/convergence data
            main: Exports bound and convergence data for LBP and UBP from the
                main algorithm
            test: Special; overrides everything else and produces a full dump of
                everything in active GAMS memory to the GDX file.
            {anything else}: skips any algorithm-specific processing
    3.  The set of equations to export is controlled by argument %3. (This
        avoids having unused equations cluttering up the GDX output.)

$OFFTEXT

* ==============================================================================
*  Setup
* ==============================================================================
** Process Arugments
* Argument 1 = Output file name
$SETLOCAL OutFile %1

* Argument 2 = Algorithm
$SETLOCAL TheAlgorithm %2

* Argument 3 = Which model
$IFi %3 == LBP                  $SETLOCAL ExportLBP         yes
$IFi %3 == UBP                  $SETLOCAL ExportUBP         yes
$IFi %3 == UBPr                 $SETLOCAL ExportUBPr        yes
$IFi %3 == Monolith             $SETLOCAL ExportMonolith    yes
$IFTHENi %3 == All
$SETLOCAL                                 ExportLBP         yes
$SETLOCAL                                 ExportUBP         yes
$SETLOCAL                                 ExportUBPr        yes
$SETLOCAL                                 ExportMonolith    yes
$ENDIF

* Arguments 4-7 = select what else to export
$IFi %4 == sets                 $SETLOCAL ExportSets        yes
$IFi %5 == sets                 $SETLOCAL ExportSets        yes
$IFi %6 == sets                 $SETLOCAL ExportSets        yes
$IFi %7 == sets                 $SETLOCAL ExportSets        yes

$IFi %4 == data                 $SETLOCAL ExportData        yes
$IFi %5 == data                 $SETLOCAL ExportData        yes
$IFi %6 == data                 $SETLOCAL ExportData        yes
$IFi %7 == data                 $SETLOCAL ExportData        yes

$IFi %4 == variables            $SETLOCAL ExportVars        yes
$IFi %5 == variables            $SETLOCAL ExportVars        yes
$IFi %6 == variables            $SETLOCAL ExportVars        yes
$IFi %7 == variables            $SETLOCAL ExportVars        yes

$IFi %4 == equations            $SETLOCAL ExportEquations   yes
$IFi %5 == equations            $SETLOCAL ExportEquations   yes
$IFi %6 == equations            $SETLOCAL ExportEquations   yes
$IFi %7 == equations            $SETLOCAL ExportEquations   yes

** Defaults
* For things that the batch arguments didn't set
$IF NOT SET ExportLBP           $SETLOCAL ExportLBP         no
$IF NOT SET ExportUBP           $SETLOCAL ExportUBP         no
$IF NOT SET ExportUBPr          $SETLOCAL ExportUBPr        no
$IF NOT SET ExportMonolith      $SETLOCAL ExportMonolith    no
$IF NOT SET ExportSets          $SETLOCAL ExportSets        no
$IF NOT SET ExportData          $SETLOCAL ExportData        no
$IF NOT SET ExportVars          $SETLOCAL ExportVars        no
$IF NOT SET ExportEquations     $SETLOCAL ExportEquations   no


* ==============================================================================
*  GDX Export
* ==============================================================================
** Everything
* If %TheAlgorithm% == test, then perform a "core dump", so to speak
$IFTHENi.coredump %TheAlgorithm% == test
execute_unload '%OutFile%' ;

** Some Things
* Otherwise, be more selective
$ELSE.coredump
execute_unload '%OutFile%'
* ------------------------------------------------------------------------------
*  Sets
* ------------------------------------------------------------------------------
$IFTHEN.sets %ExportSets% == yes
* Primary sets
    A, N, S, L, C, I, K, B, H, T
* Subsets
    BX, SX, LX, CXi, CXo, BLX
* Auxiliary sets
    NN, SS, LL, CCi, CCo, CXh, HH, HX, TT, ConnD, ConnU
$ENDIF.sets
* ------------------------------------------------------------------------------
*  Data
* ------------------------------------------------------------------------------
$IFTHEN.data %ExportData% == yes
* Model data
    BusData, BranchData, SourceData, LoadData, ConvData, ConvData2, Tau
* Computed big M values
    MB, MS, ML
* Computed relaxation data
$IFTHENi.mod %ExportLBP% == yes
$IFTHENi.algorithm NOT %TheAlgorithm% == bigM
$IFTHEN.relax NOT %ConvexRelaxation% == quadratic
    breakV, breakE, breakF, breakC
$ENDIF.relax
$IFTHEN.relax %ConvexRelaxation% == piecewise
$IFTHEN.pw %PiecewiseType% == mccormick
    seglimV, seglimE, seglimF
$ELSEIF.pw %PiecewiseType% == logtransform
    breakVV, breakEE, breakEF, breakFF
$ENDIF.pw
$ENDIF.relax
$ENDIF.algorithm
$ENDIF.mod
$ENDIF.data
* ------------------------------------------------------------------------------
*  Variables
* ------------------------------------------------------------------------------
$IFTHEN.vars %ExportVars% == yes
* Voltages
    V, E, F
* Powers
    PB, QB, PS, QS, PL, QL, PCout, QCout, PCin, QCin
* Currents
    IReB, IImB, IReS, IImS, IReL, IImL, IReCout, IImCout, IReCin, IImCin
* Auxiliary (reformulation) variables
    wVV, wEE, wEF, wFF, wC, wCC
* Design decisions
    xB, xS, xL, xC
* Objective function value
    z
* Relaxation variables
$IFTHENi.mod %ExportLBP% == yes
$IFTHENi.algorithm NOT %TheAlgorithm% == bigM
$IFTHEN.relax %ConvexRelaxation% == piecewise
    lambdaV, lambdaE, lambdaF, lambdaC
$IFTHEN.pw %PiecewiseType% == mccormick
    segVVi, segVVk, segEEi, segEEk, segEFi, segEFk, segFFi, segFFk
    yVV, yEE, yEF, yFF
$ELSEIF.pw %PiecewiseType% == logtransform
    lambdaVV, lambdaEE, lambdaEF, lambdaFF
    lnV, lnE, lnF, lnwVV, lnwEE, lnwEF, lnwFF
$ENDIF.pw
$ENDIF.relax
$ENDIF.algorithm
$ENDIF.mod
$ENDIF.vars
* ------------------------------------------------------------------------------
*  Equations
* ------------------------------------------------------------------------------
$IFTHEN.equ %ExportEquations% == yes
* Lower Bounding Problem
$IFTHEN.mod %ExportLBP% == yes
    eq_obj
    eq_PBal, eq_QBal
$IFTHEN %IncludeHarmonics% == yes
    eq_IReBal, eq_IImBal
$ENDIF
    eq_Pik_1, eq_Pik_2, eq_Pik_3, eq_Pik_4
    eq_Pki_1, eq_Pki_2, eq_Pki_3, eq_Pki_4
    eq_Qik_1, eq_Qik_2, eq_Qik_3, eq_Qik_4
    eq_Qki_1, eq_Qki_2, eq_Qki_3, eq_Qki_4
    eq_PSrc_1, eq_PSrc_2, eq_PSrc_3, eq_PSrc_4
    eq_QSrc_1, eq_QSrc_2, eq_QSrc_3, eq_QSrc_4
    eq_PLoad_1, eq_PLoad_2, eq_PLoad_3, eq_PLoad_4
    eq_QLoad_1, eq_QLoad_2, eq_QLoad_3, eq_QLoad_4
$IFTHEN %IncludeHarmonics% == yes
    eq_IReik_1, eq_IReik_2, eq_IReik_3, eq_IReik_4
    eq_IReki_1, eq_IReki_2, eq_IReki_3, eq_IReki_4
    eq_IImik_1, eq_IImik_2, eq_IImik_3, eq_IImik_4
    eq_IImki_1, eq_IImki_2, eq_IImki_3, eq_IImki_4
    eq_IReSrc_1, eq_IReSrc_2, eq_IReSrc_3, eq_IReSrc_4
    eq_IImSrc_1, eq_IImSrc_2, eq_IImSrc_3, eq_IImSrc_4
    eq_IReLoad_1, eq_IReLoad_2, eq_IReLoad_3, eq_IReLoad_4
    eq_IImLoad_1, eq_IImLoad_2, eq_IImLoad_3, eq_IImLoad_4
$ENDIF
    eq_PConvOutMin, eq_PConvOutMax, eq_PConvInMin, eq_PConvInMax
    eq_PConvInOut
    eq_QConvOutMin, eq_QConvOutMax, eq_QConvInMin, eq_QConvInMax
$IFTHEN %IncludeHarmonics% == yes
    eq_IReConvOutMin, eq_IReConvOutMax
    eq_IReConvInMin, eq_IReConvInMax
    eq_IImConvOutMin, eq_IImConvOutMax
    eq_IImConvInMin, eq_IImConvInMax
$ENDIF
    eq_VoltMag, eq_wC,
$IFTHEN %IncludeHarmonics% == yes
    eq_EHarmMin, eq_EHarmMax, eq_FHarmMin, eq_FHarmMax
$ENDIF
    eq_BranchAsgn, eq_SrcAsgn, eq_LoadAsgn
    eq_SrcConn, eq_LoadConn, eq_ConvIn, eq_ConvOut
$IFTHEN %ExtraCuts% == yes
    eq_PBLoss_L, eq_QBLoss_L, eq_QBLoss_U
    eq_VEF_1, eq_VEF_2, eq_VEF_3
$ENDIF
$IFTHEN %ConvexRelaxation% == linear
    eq_wVVii_over, eq_wVVii_under
    eq_wEEii_over, eq_wEEii_under
    eq_wFFii_over, eq_wFFii_under
    eq_wCC_over, eq_wCC_under
    eq_wVVik_over_1, eq_wVVik_over_2, eq_wVVik_under_1, eq_wVVik_under_2
    eq_wEEik_over_1, eq_wEEik_over_2, eq_wEEik_under_1, eq_wEEik_under_2
    eq_wEFik_over_1, eq_wEFik_over_2, eq_wEFik_under_1, eq_wEFik_under_2
    eq_wFFik_over_1, eq_wFFik_over_2, eq_wFFik_under_1, eq_wFFik_under_2
$ELSEIF %ConvexRelaxation% == quadratic
    eq_wVVii_over, eq_wVVii_under
    eq_wEEii_over, eq_wEEii_under
    eq_wFFii_over, eq_wFFii_under
    eq_wCC_over, eq_wCC_under
    eq_wVVik_over_1, eq_wVVik_over_2, eq_wVVik_under
    eq_wEEik_over_1, eq_wEEik_over_2, eq_wEEik_under
    eq_wEFik_over_1, eq_wEFik_over_2, eq_wEFik_under
    eq_wFFik_over_1, eq_wFFik_over_2, eq_wFFik_under
$ELSEIF %ConvexRelaxation% == piecewise
    eq_lambdaV, eq_breakV, eq_wVVii_over, eq_wVVii_under
    eq_lambdaE, eq_breakE, eq_wEEii_over, eq_wEEii_under
    eq_lambdaF, eq_breakF, eq_wFFii_over, eq_wFFii_under
    eq_lambdaC, eq_breakC, eq_wCC_over, eq_wCC_under
$IFTHEN.pw %PiecewiseType% == mccormick
    eq_segsel_VV, eq_segsel_EE, eq_segsel_EF, eq_segsel_FF
    eq_seg_VVi, eq_seg_VVk, eq_seg_EEi, eq_seg_EEk
    eq_seg_EFi, eq_seg_EFk, eq_seg_FFi, eq_seg_FFk
    eq_seg_VVi_lo, eq_seg_VVi_up, eq_seg_VVk_lo, eq_seg_VVk_up
    eq_seg_EEi_lo, eq_seg_EEi_up, eq_seg_EEk_lo, eq_seg_EEk_up
    eq_seg_EFi_lo, eq_seg_EFi_up, eq_seg_EFk_lo, eq_seg_EFk_up
    eq_seg_FFi_lo, eq_seg_FFi_up, eq_seg_FFk_lo, eq_seg_FFk_up
    eq_wVVik_over_1, eq_wVVik_over_2, eq_wVVik_under_1, eq_wVVik_under_2
    eq_wEEik_over_1, eq_wEEik_over_2, eq_wEEik_under_1, eq_wEEik_under_2
    eq_wEFik_over_1, eq_wEFik_over_2, eq_wEFik_under_1, eq_wEFik_under_2
    eq_wFFik_over_1, eq_wFFik_over_2, eq_wFFik_under_1, eq_wFFik_under_2
$ELSEIF.pw %PiecewiseType% == logtransform
    eq_lambdaVV, eq_breakVV
    eq_lambdaEE, eq_breakEE
    eq_lambdaEF, eq_breakEF
    eq_lambdaFF, eq_breakFF
    eq_lnV_over, eq_lnV_under
    eq_lnE_over, eq_lnE_under
    eq_lnF_over, eq_lnF_under
    eq_lnwVV, eq_lnwVV_over, eq_lnwVV_under
    eq_lnwEE, eq_lnwEE_over, eq_lnwEE_under
    eq_lnwEF, eq_lnwEF_over, eq_lnwEF_under
    eq_lnwFF, eq_lnwFF_over, eq_lnwFF_under
$ENDIF.pw
$ENDIF
$IFTHENi.algorithm %TheAlgorithm% == main
    eq_SolElim
$ENDIF.algorithm
$ENDIF.mod
* Upper Bounding Problem
$IFTHEN.mod %ExportUBP% == yes
    eq_obj
    eq_PBal, eq_QBal
$IFTHEN %IncludeHarmonics% == yes
    eq_IReBal, eq_IImBal
$ENDIF
    eq_Pik, eq_Pki, eq_Qik, eq_Qki
    eq_PSrc, eq_QSrc
    eq_PLoad, eq_QLoad
$IFTHEN %IncludeHarmonics% == yes
    eq_IReik, eq_IReki, eq_IImik, eq_IImki
    eq_IReSrc, eq_IImSrc
    eq_IReLoad, eq_IImLoad
$ENDIF
    eq_PConvOutMin, eq_PConvOutMax, eq_PConvInMin, eq_PConvInMax
    eq_PConvInOut,
    eq_QConvOutMin, eq_QConvOutMax, eq_QConvInMin, eq_QConvInMax
$IFTHEN %IncludeHarmonics% == yes
    eq_IReConvOutMin, eq_IReConvOutMax
    eq_IReConvInMin, eq_IReConvInMax
    eq_IImConvOutMin, eq_IImConvOutMax
    eq_IImConvInMin, eq_IImConvInMax
$ENDIF
    eq_VoltMag, eq_wC
$IFTHEN %IncludeHarmonics% == yes
    eq_EHarmMin, eq_EHarmMax, eq_FHarmMin, eq_FHarmMax
$ENDIF
    eq_BranchAsgn, eq_SrcAsgn, eq_LoadAsgn
    eq_SrcConn, eq_LoadConn, eq_ConvIn, eq_ConvOut
$IFTHEN %ExtraCuts% == yes
    eq_PBLoss_L, eq_QBLoss_L, eq_QBLoss_U
    eq_VEF_1, eq_VEF_2, eq_VEF_3
$ENDIF
    eq_wVV, eq_wEE, eq_wEF, eq_wFF, eq_wCC
$ENDIF.mod
* Reduced Upper Bounding Problem
$IFTHEN.mod %ExportUBPr% == yes
    eq_obj
    eq_PBal, eq_QBal
    eq_Pik, eq_Pki, eq_Qik, eq_Qki
    eq_PSrc, eq_QSrc
    eq_PLoad, eq_QLoad
    eq_PConvOutMin, eq_PConvOutMax, eq_PConvInMin, eq_PConvInMax
    eq_PConvInOut,
    eq_QConvOutMin, eq_QConvOutMax, eq_QConvInMin, eq_QConvInMax
    eq_VoltMag, eq_wC
    eq_BranchAsgn, eq_SrcAsgn, eq_LoadAsgn
    eq_SrcConn, eq_LoadConn, eq_ConvIn, eq_ConvOut
$IFTHEN %ExtraCuts% == yes
    eq_PBLoss_L, eq_QBLoss_L, eq_QBLoss_U
    eq_VEF_1, eq_VEF_2, eq_VEF_3
$ENDIF
    eq_wVV, eq_wEE, eq_wEF, eq_wFF, eq_wCC
$ENDIF.mod
* Monolith
$IFTHEN.mod %ExportMonolith% == yes
    eq_obj
    eq_PBal, eq_QBal
$IFTHEN %IncludeHarmonics% == yes
    eq_IReBal, eq_IImBal
$ENDIF
    eq_Pik_1, eq_Pik_2, eq_Pik_3, eq_Pik_4
    eq_Pki_1, eq_Pki_2, eq_Pki_3, eq_Pki_4
    eq_Qik_1, eq_Qik_2, eq_Qik_3, eq_Qik_4
    eq_Qki_1, eq_Qki_2, eq_Qki_3, eq_Qki_4
    eq_PSrc_1, eq_PSrc_2, eq_PSrc_3, eq_PSrc_4
    eq_QSrc_1, eq_QSrc_2, eq_QSrc_3, eq_QSrc_4
    eq_PLoad_1, eq_PLoad_2, eq_PLoad_3, eq_PLoad_4
    eq_QLoad_1, eq_QLoad_2, eq_QLoad_3, eq_QLoad_4
$IFTHEN %IncludeHarmonics% == yes
    eq_IReik_1, eq_IReik_2, eq_IReik_3, eq_IReik_4
    eq_IReki_1, eq_IReki_2, eq_IReki_3, eq_IReki_4
    eq_IImik_1, eq_IImik_2, eq_IImik_3, eq_IImik_4
    eq_IImki_1, eq_IImki_2, eq_IImki_3, eq_IImki_4
    eq_IReSrc_1, eq_IReSrc_2, eq_IReSrc_3, eq_IReSrc_4
    eq_IImSrc_1, eq_IImSrc_2, eq_IImSrc_3, eq_IImSrc_4
    eq_IReLoad_1, eq_IReLoad_2, eq_IReLoad_3, eq_IReLoad_4
    eq_IImLoad_1, eq_IImLoad_2, eq_IImLoad_3, eq_IImLoad_4
$ENDIF
    eq_PConvOutMin, eq_PConvOutMax, eq_PConvInMin, eq_PConvInMax
    eq_PConvInOut
    eq_QConvOutMin, eq_QConvOutMax, eq_QConvInMin, eq_QConvInMax
$IFTHEN %IncludeHarmonics% == yes
    eq_IReConvOutMin, eq_IReConvOutMax
    eq_IReConvInMin, eq_IReConvInMax
    eq_IImConvOutMin, eq_IImConvOutMax
    eq_IImConvInMin, eq_IImConvInMax
$ENDIF
    eq_VoltMag, eq_wC,
$IFTHEN %IncludeHarmonics% == yes
    eq_EHarmMin, eq_EHarmMax, eq_FHarmMin, eq_FHarmMax
$ENDIF
    eq_BranchAsgn, eq_SrcAsgn, eq_LoadAsgn
    eq_SrcConn, eq_LoadConn, eq_ConvIn, eq_ConvOut
$IFTHEN %ExtraCuts% == yes
    eq_PBLoss_L, eq_QBLoss_L, eq_QBLoss_U
    eq_VEF_1, eq_VEF_2, eq_VEF_3
$ENDIF
    eq_wVV, eq_wEE, eq_wEF, eq_wFF, eq_wCC
$ENDIF.mod
$ENDIF.equ
* ------------------------------------------------------------------------------
*  Algorithm-Specific Data
* ------------------------------------------------------------------------------
** Main Algorithm
$IFTHENi.algorithm %TheAlgorithm% == main
* Iteration counters
    R, RR
* Decisions
    xB_LBP, xS_LBP, xL_LBP, xC_LBP
* Objective function info.
    zLBP_LB, zLBP_UB zUBP_LB, zUBP_UB, zLB, zBest
* Algorithm convergence info.
    cvg, isFeas, rBest
* Timings
    algTiming
$ENDIF.algorithm
** Bound Tightening Algorithm
$IFTHENi.algorithm %TheAlgorithm% == bound_tightening
* Iteration counter
    R
* Bound tightening history
    boundV, boundE, boundF
    boundwVV, boundwEE, boundwEF, boundwFF, boundwC, boundwCC
    changeV, changeE, changeF
    changewVV, changewEE, changewEF, changewFF, changewC, changewCC
$ENDIF.algorithm
;
$ENDIF.coredump
