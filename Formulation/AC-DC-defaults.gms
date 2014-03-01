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

$TITLE Default compile-time settings for AC-DC model
$ONTEXT

    Filename : AC-DC-defaults.gms
    Author   : Stephen Frank
    Release  : March 1, 2014
    Version  : 1.1

    The various AC-DC formulation files use a number of compile-time flags,
    many of which are used in multiple files. This file consolidates the
    documentation and default values for these flags. Include this file (via
    $INCLUDE) to set default values for otherwise unspecified compile-time
    flags and settings.

    Flags and settings are listed separately in alphabetical order. Each item
    has documentation that states
        1. The purpose of the flag or setting
        2. Possible values (if applicable)
        3. The file(s) in which the flag or setting is used
    
$OFFTEXT

* ==============================================================================
*  Flags
* ==============================================================================
** BigM
* Generate big-M constraints for binary-continuous terms of the form u = y f(x),
* y binary, f(x) linear?
*   Values:     yes, no
*   Files:      AC-DC-equations-reform.gms
$IF NOT SET BigM                    $SETGLOBAL BigM                     yes

** BilinearEq
* Generate bilinear equality constraints of the form u = y f(x), y binary,
* f(x) linear?
*   Values:     yes, no
*   Files:      AC-DC-equations-reform.gms
$IF NOT SET BilinearEq              $SETGLOBAL BilinearEq               yes

** ExactReform
* Generate exact reformulation variable relationships? (These are nonconvex)
*   Values:     yes, no
*   Files:      AC-DC-equations-reform.gms
$IF NOT SET ExactReform             $SETGLOBAL ExactReform              yes

** ExportEveryIteration
* Export a GDX on every NGBD iteration?
*   Values:     yes, no
*   Files:      AC-DC-algorithm.gms
*   Note:       Caution! Generates a lot of files.
$IF NOT SET ExportEveryIteration    $SETGLOBAL ExportEveryIteration     no

** ExtraCuts
* Generate the extra valid ineqalities (feasibility cuts)?
*   Values:     yes, no
*   Files:      AC-DC-equations-reform.gms
*               AC-DC-models.gms
$IF NOT SET ExtraCuts               $SETGLOBAL ExtraCuts                yes

** FindAllFeasibleDesigns
* Find all feasible system designs (not just best ones)?
*   Values:     yes, no
*   Files:      AC-DC-algorithm.gms
$IF NOT SET FindAllFeasibleDesigns  $SETGLOBAL FindAllFeasibleDesigns   no

** LoadSets
* Load set elements on data import?
*   Values:     yes, no
*   Files:      AC-DC-data.gms
$IF NOT SET LoadSets                $SETGLOBAL LoadSets                 yes

** LoadParams
* Load parameter values on data import?
*   Values:     yes, no
*   Files:      AC-DC-data.gms
$IF NOT SET LoadParams              $SETGLOBAL LoadParams               yes

** LoadVars
* Load variable levels and values on data import?
*   Values:     yes, no
*   Files:      AC-DC-data.gms
*   Note:       If %LoadVars% == yes, then %InitVars% is ignored
$IF NOT SET LoadVars                $SETGLOBAL LoadVars                 no

** LoadBigM
* Load big M values on data import?
*   Values:     yes, no
*   Files:      AC-DC-data.gms
$IF NOT SET LoadBigM                $SETGLOBAL LoadBigM                 no

** IncludeHarmonics
* Include harmonic current equations in the formulation?
*   Values:     yes, no
*   Files:      AC-DC-big-M.gms
*               AC-DC-bound-tightening.gms
*               AC-DC-equations-reform.gms
*               AC-DC-models.gms
$IF NOT SET IncludeHarmonics        $SETGLOBAL IncludeHarmonics         yes

** InitVars
* Initialize variable levels and bounds after data import?
*   Values:     yes, no
*   Files:      AC-DC-data.gms
*   Note:       Ignored if %LoadVars% == yes
$IF NOT SET InitVars                $SETGLOBAL InitVars                 yes

* ==============================================================================
*  Settings
* ==============================================================================
** BigMSlack
* Specify an absolute extra slack for the Big M values
* (used to avoid numerical error when tight against a computed Big M value)
*   Values:     {any positive real number}
*   Files:      AC-DC-big-M.gms
*               AC-DC-bound-tightening.gms
$IF NOT SET BigMSlack               $SETGLOBAL BigMSlack                1e-3

** BoundSlack
* Specify an absolute extra slack for bound tightening
* (used to avoid numerical error when tight against a computed bound)
*   Values:     {any positive real number}
*   Files:      AC-DC-bound-tightening.gms
$IF NOT SET BoundSlack              $SETGLOBAL BoundSlack               1e-3

** ConvexRelaxation
* Specify a type of convex relaxation
*   Values:     linear = McCormick inequalities + valid linear cuts
*               quadratic = Convex quadratic under-, linear overestimator
*               piecewise = Piecewise linear + valid linear cuts
*   Files:      AC-DC-equations-relax.gms
*               AC-DC-models.gms
$IF NOT SET ConvexRelaxation        $SETGLOBAL ConvexRelaxation         linear

** CvgTol
* Specify relative convergence tolerance (for bound tightening)
*   Values:     {any positive real number}
*   Files:      AC-DC-bound-tightening.gms
$IF NOT SET CvgTol                  $SETGLOBAL CvgTol                   1e-2

** DataFile
* Specify the data file to use
*   Values:     {any valid .gdx file name}
*   Files:      AC-DC-data.gms
*   Note:       There is no default for this setting; it must be user specified
$IF NOT SET DataFile                $ABORT 'DataFile' must be specified.

** IntegerTol
* Specify tolerance for determining valid integer values during bound tightening
*   Values:     {any positive real number less than 1}
*   Files:      AC-DC-bound-tightening.gms
$IF NOT SET IntegerTol              $SETGLOBAL IntegerTol               1e-6

* ------------------------------------------------------------------------------
*  LBP Settings
* ------------------------------------------------------------------------------
** LBP_abs
* Absolute tolerance for lower bounding problem
*   Values:     {any positive real number}
*   Files:      AC-DC-algorithm.gms
$IF NOT SET LBP_abs                $SETGLOBAL LBP_abs                   0

** LBP_rel
* Relative tolerance for lower bounding problem
*   Values:     {any positive real number}
*   Files:      AC-DC-algorithm.gms
$IF NOT SET LBP_rel                $SETGLOBAL LBP_rel                   1e-3

** LBP_reslim
* Resource limit (CPU sec) for lower bounding problem
*   Values:     {any positive real number}
*   Files:      AC-DC-algorithm.gms
*   Note:       By default, this value is left unset such that the global
*               resource limit is used instead.
$IFTHEN NOT SET LBP_reslim
$LOG LBP_reslim is unset; using default resource limit for LBP.
$ENDIF

** PiecewiseType
* Specify a type of piecewise linear relaxation
*   Values:     mccormick = McCormick inequalities over segments
*               logtransform = Log transformation: bilinear terms to univariate
*   Files:      AC-DC-equations-relax.gms
*               AC-DC-models.gms
$IF NOT SET PiecewiseType           $SETGLOBAL PiecewiseType            mccormick

** NumSeg
* Set number of segments for linear cuts and piecewise linear relaxations
*   Values:     {any positive integer}
*   Files:      AC-DC-equations-relax.gms
$IF NOT SET NumSeg                  $SETGLOBAL NumSeg                   5

** MaxIter
* Maximum number of iterations for NGBD algorithm
*   Values:     {any positive integer}
*   Files:      AC-DC-algorithm.gms
*               AC-DC-defs.gms
*               AC-DC-bound-tightening.gms
*   Note:       Sets the maximum iterations for the NGBD algorithm in
*               AC-DC-algorithm.gms indirectly by limiting the size of the
*               looping set
$IF NOT SET MaxIter                 $SETGLOBAL MaxIter                  1000

** RelPath
* Specify the relative path to the AC-DC formulation files
*   Values:     {a valid system path, without terminating /}
*   Files:      AC-DC-bound-tightening.gms
*               AC-DC-data.gms,
*               AC-DC-equations-relax.gms
*   Note:       Allows GAMS to properly locate includes within other includes
$IF NOT SET RelPath                 $SETGLOBAL RelPath                  .

** TimeLimit
* Specifies a total time limit, in seconds, for bound tightening or the
* solution algorithm
*   Values:     {any positive number}
*   Files:      AC-DC-algorithm.gms
*               AC-DC-bound-tightening.gms
$IF NOT SET TimeLimit               $SETGLOBAL TimeLimit                3600

* ------------------------------------------------------------------------------
*  UBP Settings
* ------------------------------------------------------------------------------
** UBP_abs
* Absolute tolerance for upper bounding problem
*   Values:     {any positive real number}
*   Files:      AC-DC-algorithm.gms
*   Note:       Required to be nonzero for GLoMIQO
$IF NOT SET UBP_abs                $SETGLOBAL UBP_abs                   1e-6

** UBP_rel
* Relative tolerance for upper bounding problem
*   Values:     {any positive real number}
*   Files:      AC-DC-algorithm.gms
$IF NOT SET UBP_rel                $SETGLOBAL UBP_rel                   1e-2

** UBP_reslim
* Resource limit (CPU sec) for upper bounding problem
*   Values:     {any positive real number}
*   Files:      AC-DC-algorithm.gms
*   Note:       This value is specified per UBP subproblem. By default, it is
*               left unset such that the global resource limit is used instead.
$IFTHEN NOT SET UBP_reslim
$LOG UBP_reslim is unset; using default resource limit for UBP.
$ENDIF