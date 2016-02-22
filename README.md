Mahadevan-Garofalini Water Potential Implementation for LAMMPS
================================================================================

This is the code that implements the [Mahadevan & Garofalini potential][1] in
[LAMMPS][2].  This code is _not_ very thoroughly tested and is _not_ the code
used to obtain the results in any publications put forth by Garofalini or others
at the Intermolecular Science Laboratory at Rutgers.  

Although the numerics and overall dynamics of the LAMMPS code are consistent
with the original implementation by Mahadevan and Garofalini, the density of
water at 1 atm, 298 K comes out consistently higher in this LAMMPS
implementation by a small amount.  Using a tabulated 2-body potential generated
by the original Mahadevan code with LAMMPS (`pair_style table`) shows the same
behavior.  This suggests that the difference in density is caused either by the
3-body potential in LAMMPS or some difference in the thermostat or barostat,
but the exact reason has never been determined.

Nevertheless, this repository contains all of the files one would need to add
to LAMMPS.  The following files are included:

* `pair_gg_coul_wolf.cpp` and `pair_gg_coul_wolf.h` - the code to calculate the
  2-body potential with Wolf summation
* `in.ggh2o` - sample input that illustrates how to combine the custom 2-body
  with LAMMPS' own implementation of the Stillinger Weber potential to get the
  3-body interactions
* `gg.sw.kcal` - the Stillinger Weber parameters to get the correct 3-body
  interactions for the dissociative water potential.  The two-body parameters
  are all zeroed out.
* `data.knite9.0792` and `data.knite9.4000` - sample water configurations to run

Inputs that exercise the SiO2 or water-silica interactions have not been created,
but there should be a way to use LAMMPS's Stillinger-Weber potential to achieve
the same 3-body cross-terms for Si-O-H described in the Mahadevan paper.

This code was last tested with a February 2016 (`16Feb16`) release of LAMMPS and
installation should be a simple matter of

1. copying `pair_gg_coul_wolf.cpp` and `pair_gg_coul_wolf.h` to the `src`
   directory within the LAMMPS source distribution
2. adding `#include "pair_gg_coul_wolf.h"` to `src/style_pair.h`
3. `make my_system` as usual

This code has not thoroughly tested, so there may be bugs.  It is intended to be
a helpful starting point for implementing the Mahadevan-Garofalini Dissociative
Water Potential, but it comes with **no guarantees of correctness**.

[1]: http://dx.doi.org/10.1021/jp072530o
[2]: http://lammps.sandia.gov
