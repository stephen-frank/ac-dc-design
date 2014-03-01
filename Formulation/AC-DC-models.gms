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

$TITLE Model definitions for AC-DC model
$ONTEXT

    Filename : AC-DC-models.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    This file provides several model definitions for the mixed AC-DC
    distribution system design problem. It is designed to be included (via
    $INCLUDE) in other .gms files, and relies on the following external compile
    time flags and settings (which are documented in AC-DC-defaults.gms):
        ExtraCuts           Generate valid feasiblity cuts?
        ConvexRelaxation    Specify type of convex relaxation
        PiecewiseType       Specify piecewise relaxation subtype
        NumSeg              Specify number of segments to use in relaxations

    The models provided include:
        LBP                 Lower bounding problem
        UBP                 Upper bounding problem
        Monolith            Entire problem

$OFFTEXT

* =============================================================================
*  Models
* =============================================================================
** Lower Bounding Problem
* Uses Big-M constraints and convex relaxations of 'w' variables
* Uses extra cuts if user requested
model LBP      / eq_obj,
                 eq_PBal, eq_QBal,
$IFTHEN %IncludeHarmonics% == yes
                 eq_IReBal, eq_IImBal,
$ENDIF
                 eq_Pik_1, eq_Pik_2, eq_Pik_3, eq_Pik_4,
                 eq_Pki_1, eq_Pki_2, eq_Pki_3, eq_Pki_4,
                 eq_Qik_1, eq_Qik_2, eq_Qik_3, eq_Qik_4,
                 eq_Qki_1, eq_Qki_2, eq_Qki_3, eq_Qki_4,
                 eq_PSrc_1, eq_PSrc_2, eq_PSrc_3, eq_PSrc_4,
                 eq_QSrc_1, eq_QSrc_2, eq_QSrc_3, eq_QSrc_4,
                 eq_PLoad_1, eq_PLoad_2, eq_PLoad_3, eq_PLoad_4,
                 eq_QLoad_1, eq_QLoad_2, eq_QLoad_3, eq_QLoad_4,
$IFTHEN %IncludeHarmonics% == yes
                 eq_IReik_1, eq_IReik_2, eq_IReik_3, eq_IReik_4,
                 eq_IReki_1, eq_IReki_2, eq_IReki_3, eq_IReki_4,
                 eq_IImik_1, eq_IImik_2, eq_IImik_3, eq_IImik_4,
                 eq_IImki_1, eq_IImki_2, eq_IImki_3, eq_IImki_4,
                 eq_IReSrc_1, eq_IReSrc_2, eq_IReSrc_3, eq_IReSrc_4,
                 eq_IImSrc_1, eq_IImSrc_2, eq_IImSrc_3, eq_IImSrc_4,
                 eq_IReLoad_1, eq_IReLoad_2, eq_IReLoad_3, eq_IReLoad_4,
                 eq_IImLoad_1, eq_IImLoad_2, eq_IImLoad_3, eq_IImLoad_4,
$ENDIF
                 eq_PConvOutMin, eq_PConvOutMax, eq_PConvInMin, eq_PConvInMax,
                 eq_PConvInOut,
                 eq_QConvOutMin, eq_QConvOutMax, eq_QConvInMin, eq_QConvInMax,
$IFTHEN %IncludeHarmonics% == yes
                 eq_IReConvOutMin, eq_IReConvOutMax,
                 eq_IReConvInMin, eq_IReConvInMax,
                 eq_IImConvOutMin, eq_IImConvOutMax,
                 eq_IImConvInMin, eq_IImConvInMax,
$ENDIF
                 eq_VoltMag, eq_wC,
$IFTHEN %IncludeHarmonics% == yes
                 eq_EHarmMin, eq_EHarmMax, eq_FHarmMin, eq_FHarmMax,
$ENDIF
                 eq_BranchAsgn, eq_SrcAsgn, eq_LoadAsgn, eq_BranchLoad,
                 eq_SrcConn, eq_LoadConn, eq_ConvIn, eq_ConvOut,
$IFTHEN %ExtraCuts% == yes
                 eq_PBLoss_L, eq_QBLoss_L, eq_QBLoss_U,
                 eq_VEF_1, eq_VEF_2, eq_VEF_3,
$ENDIF
$IFTHEN %ConvexRelaxation% == linear
                 eq_wVVii_over, eq_wVVii_under,
                 eq_wEEii_over, eq_wEEii_under,
                 eq_wFFii_over, eq_wFFii_under,
                 eq_wCC_over, eq_wCC_under,
                 eq_wVVik_over_1, eq_wVVik_over_2,
                 eq_wVVik_under_1, eq_wVVik_under_2,
                 eq_wEEik_over_1, eq_wEEik_over_2,
                 eq_wEEik_under_1, eq_wEEik_under_2,
                 eq_wEFik_over_1, eq_wEFik_over_2,
                 eq_wEFik_under_1, eq_wEFik_under_2,
                 eq_wFFik_over_1, eq_wFFik_over_2,
                 eq_wFFik_under_1, eq_wFFik_under_2,
$ELSEIF %ConvexRelaxation% == quadratic
                 eq_wVVii_over, eq_wVVii_under,
                 eq_wEEii_over, eq_wEEii_under,
                 eq_wFFii_over, eq_wFFii_under,
                 eq_wCC_over, eq_wCC_under,
                 eq_wVVik_over_1, eq_wVVik_over_2, eq_wVVik_under,
                 eq_wEEik_over_1, eq_wEEik_over_2, eq_wEEik_under,
                 eq_wEFik_over_1, eq_wEFik_over_2, eq_wEFik_under,
                 eq_wFFik_over_1, eq_wFFik_over_2, eq_wFFik_under,
$ELSEIF %ConvexRelaxation% == piecewise
                 eq_lambdaV, eq_breakV, eq_wVVii_over, eq_wVVii_under,
                 eq_lambdaE, eq_breakE, eq_wEEii_over, eq_wEEii_under,
                 eq_lambdaF, eq_breakF, eq_wFFii_over, eq_wFFii_under,
                 eq_lambdaC, eq_breakC, eq_wCC_over, eq_wCC_under,
$IFTHEN.PW %PiecewiseType% == mccormick
                 eq_segsel_VV, eq_segsel_EE, eq_segsel_EF, eq_segsel_FF,
                 eq_seg_VVi, eq_seg_VVk, eq_seg_EEi, eq_seg_EEk,
                 eq_seg_EFi, eq_seg_EFk, eq_seg_FFi, eq_seg_FFk,
                 eq_seg_VVi_lo, eq_seg_VVi_up, eq_seg_VVk_lo, eq_seg_VVk_up,
                 eq_seg_EEi_lo, eq_seg_EEi_up, eq_seg_EEk_lo, eq_seg_EEk_up,
                 eq_seg_EFi_lo, eq_seg_EFi_up, eq_seg_EFk_lo, eq_seg_EFk_up,
                 eq_seg_FFi_lo, eq_seg_FFi_up, eq_seg_FFk_lo, eq_seg_FFk_up,
                 eq_wVVik_over_1, eq_wVVik_over_2,
                 eq_wVVik_under_1, eq_wVVik_under_2,
                 eq_wEEik_over_1, eq_wEEik_over_2,
                 eq_wEEik_under_1, eq_wEEik_under_2,
                 eq_wEFik_over_1, eq_wEFik_over_2,
                 eq_wEFik_under_1, eq_wEFik_under_2,
                 eq_wFFik_over_1, eq_wFFik_over_2,
                 eq_wFFik_under_1, eq_wFFik_under_2,
$ELSEIF.PW %PiecewiseType% == logtransform
                 eq_lambdaVV, eq_breakVV,
                 eq_lambdaEE, eq_breakEE,
                 eq_lambdaEF, eq_breakEF,
                 eq_lambdaFF, eq_breakFF,
                 eq_lnV_over, eq_lnV_under,
                 eq_lnE_over, eq_lnE_under,
                 eq_lnF_over, eq_lnF_under,
                 eq_lnwVV, eq_lnwVV_over, eq_lnwVV_under,
                 eq_lnwEE, eq_lnwEE_over, eq_lnwEE_under,
                 eq_lnwEF, eq_lnwEF_over, eq_lnwEF_under,
                 eq_lnwFF, eq_lnwFF_over, eq_lnwFF_under,
$ENDIF.PW
$ENDIF
               / ;

** Upper Bounding Problem
* Does not Big-M constraints (requires fixed 'x' variables)
* Uses exact reformulations of 'w' variables
* Uses extra cuts if user requested
model UBP      / eq_obj,
                 eq_PBal, eq_QBal,
$IFTHEN %IncludeHarmonics% == yes
                 eq_IReBal, eq_IImBal,
$ENDIF
                 eq_Pik, eq_Pki, eq_Qik, eq_Qki,
                 eq_PSrc, eq_QSrc,
                 eq_PLoad, eq_QLoad,
$IFTHEN %IncludeHarmonics% == yes
                 eq_IReik, eq_IReki, eq_IImik, eq_IImki,
                 eq_IReSrc, eq_IImSrc,
                 eq_IReLoad, eq_IImLoad,
$ENDIF
                 eq_PConvOutMin, eq_PConvOutMax, eq_PConvInMin, eq_PConvInMax,
                 eq_PConvInOut,
                 eq_QConvOutMin, eq_QConvOutMax, eq_QConvInMin, eq_QConvInMax,
$IFTHEN %IncludeHarmonics% == yes
                 eq_IReConvOutMin, eq_IReConvOutMax,
                 eq_IReConvInMin, eq_IReConvInMax,
                 eq_IImConvOutMin, eq_IImConvOutMax,
                 eq_IImConvInMin, eq_IImConvInMax,
$ENDIF
                 eq_VoltMag, eq_wC,
$IFTHEN %IncludeHarmonics% == yes
                 eq_EHarmMin, eq_EHarmMax, eq_FHarmMin, eq_FHarmMax,
$ENDIF
$IFTHEN %ExtraCuts% == yes
                 eq_PBLoss_L, eq_QBLoss_L, eq_QBLoss_U,
                 eq_VEF_1, eq_VEF_2, eq_VEF_3,
$ENDIF
                 eq_wVV, eq_wEE, eq_wEF, eq_wFF, eq_wCC
               / ;

** Reduced Upper Bounding Problem
* Same as UBP except excludes all harmonic equations; useful for the main
* algorithm (where the harmonic solution can be ported directly from the LBP).
model UBPr     / eq_obj,
                 eq_PBal, eq_QBal,
                 eq_Pik, eq_Pki, eq_Qik, eq_Qki,
                 eq_PSrc, eq_QSrc,
                 eq_PLoad, eq_QLoad,
                 eq_PConvOutMin, eq_PConvOutMax, eq_PConvInMin, eq_PConvInMax,
                 eq_PConvInOut,
                 eq_QConvOutMin, eq_QConvOutMax, eq_QConvInMin, eq_QConvInMax,
                 eq_VoltMag, eq_wC,
$IFTHEN %ExtraCuts% == yes
                 eq_PBLoss_L, eq_QBLoss_L, eq_QBLoss_U,
                 eq_VEF_1, eq_VEF_2, eq_VEF_3,
$ENDIF
                 eq_wVV, eq_wEE, eq_wEF, eq_wFF, eq_wCC
               / ;

** Monolith
* Uses Big-M constraints and exact reformulations of 'w' variables
* Uses extra cuts if user requested
model Monolith / eq_obj,
                 eq_PBal, eq_QBal,
$IFTHEN %IncludeHarmonics% == yes
                 eq_IReBal, eq_IImBal,
$ENDIF
                 eq_Pik_1, eq_Pik_2, eq_Pik_3, eq_Pik_4,
                 eq_Pki_1, eq_Pki_2, eq_Pki_3, eq_Pki_4,
                 eq_Qik_1, eq_Qik_2, eq_Qik_3, eq_Qik_4,
                 eq_Qki_1, eq_Qki_2, eq_Qki_3, eq_Qki_4,
                 eq_PSrc_1, eq_PSrc_2, eq_PSrc_3, eq_PSrc_4,
                 eq_QSrc_1, eq_QSrc_2, eq_QSrc_3, eq_QSrc_4,
                 eq_PLoad_1, eq_PLoad_2, eq_PLoad_3, eq_PLoad_4,
                 eq_QLoad_1, eq_QLoad_2, eq_QLoad_3, eq_QLoad_4,
$IFTHEN %IncludeHarmonics% == yes
                 eq_IReik_1, eq_IReik_2, eq_IReik_3, eq_IReik_4,
                 eq_IReki_1, eq_IReki_2, eq_IReki_3, eq_IReki_4,
                 eq_IImik_1, eq_IImik_2, eq_IImik_3, eq_IImik_4,
                 eq_IImki_1, eq_IImki_2, eq_IImki_3, eq_IImki_4,
                 eq_IReSrc_1, eq_IReSrc_2, eq_IReSrc_3, eq_IReSrc_4,
                 eq_IImSrc_1, eq_IImSrc_2, eq_IImSrc_3, eq_IImSrc_4,
                 eq_IReLoad_1, eq_IReLoad_2, eq_IReLoad_3, eq_IReLoad_4,
                 eq_IImLoad_1, eq_IImLoad_2, eq_IImLoad_3, eq_IImLoad_4,
$ENDIF
                 eq_PConvOutMin, eq_PConvOutMax, eq_PConvInMin, eq_PConvInMax,
                 eq_PConvInOut,
                 eq_QConvOutMin, eq_QConvOutMax, eq_QConvInMin, eq_QConvInMax,
$IFTHEN %IncludeHarmonics% == yes
                 eq_IReConvOutMin, eq_IReConvOutMax,
                 eq_IReConvInMin, eq_IReConvInMax,
                 eq_IImConvOutMin, eq_IImConvOutMax,
                 eq_IImConvInMin, eq_IImConvInMax,
$ENDIF
                 eq_VoltMag, eq_wC,
$IFTHEN %IncludeHarmonics% == yes
                 eq_EHarmMin, eq_EHarmMax, eq_FHarmMin, eq_FHarmMax,
$ENDIF
                 eq_BranchAsgn, eq_SrcAsgn, eq_LoadAsgn, eq_BranchLoad,
                 eq_SrcConn, eq_LoadConn, eq_ConvIn, eq_ConvOut,
$IFTHEN %ExtraCuts% == yes
                 eq_PBLoss_L, eq_QBLoss_L, eq_QBLoss_U,
                 eq_VEF_1, eq_VEF_2, eq_VEF_3,
$ENDIF
                 eq_wVV, eq_wEE, eq_wEF, eq_wFF, eq_wCC
               / ;

** Harmonic Feasibility Problem
* Includes only harmonic-related equations; designed to be used to check
* the feasibility of a particular fundamental solution
$IFTHEN %IncludeHarmonics% == yes
model HFeasbility / eq_dummy,
                    eq_IReik_1, eq_IReik_2, eq_IReik_3, eq_IReik_4,
                    eq_IReki_1, eq_IReki_2, eq_IReki_3, eq_IReki_4,
                    eq_IImik_1, eq_IImik_2, eq_IImik_3, eq_IImik_4,
                    eq_IImki_1, eq_IImki_2, eq_IImki_3, eq_IImki_4,
                    eq_IReSrc_1, eq_IReSrc_2, eq_IReSrc_3, eq_IReSrc_4,
                    eq_IImSrc_1, eq_IImSrc_2, eq_IImSrc_3, eq_IImSrc_4,
                    eq_IReLoad_1, eq_IReLoad_2, eq_IReLoad_3, eq_IReLoad_4,
                    eq_IImLoad_1, eq_IImLoad_2, eq_IImLoad_3, eq_IImLoad_4,
                    eq_IReConvOutMin, eq_IReConvOutMax,
                    eq_IReConvInMin, eq_IReConvInMax,
                    eq_IImConvOutMin, eq_IImConvOutMax,
                    eq_IImConvInMin, eq_IImConvInMax,
                    eq_EHarmMin, eq_EHarmMax, eq_FHarmMin, eq_FHarmMax,
                    eq_BranchAsgn, eq_SrcAsgn, eq_LoadAsgn,
                    eq_SrcConn, eq_LoadConn, eq_ConvIn, eq_ConvOut,
                  / ;
$ENDIF
