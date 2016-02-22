Mahadevan-Garofalini Water Potential Implementation for LAMMPS
================================================================================

This is the code that implements the Mahadevan & Garofalini potential in LAMMPS,
but it's not very thoroughly tested.  Although the numerics and overall
dynamics of the LAMMPS code are consistent with the original implementation by
Mahadevan and Garofalini, the density of water (1 atm, 298 K) came out
consistently higher in the LAMMPS implementation by a small amount.  I never
figured out why, although I did use a tabulated 2-body potential generated with
the Mahadevan code with LAMMPS and saw the same behavior.  This suggested that
the difference in density was caused either by the 3-body potential in LAMMPS or
some difference in the thermostat or barostat.

Nevertheless, this repository contains all of the files one would need to add
to LAMMPS to experiment with it yourself.  The following files are included:

* `pair_gg_coul_wolf.cpp` and `pair_gg_coul_wolf.h` - the code to calculate the
  2-body potential with Wolf summation
* `in.ggh2o` - sample input that illustrates how to combine the custom 2-body
  with LAMMPS' own implementation of the Stillinger Weber potential to get the
  3-body interactions
* `gg.sw.kcal` - the Stillinger Weber parameters to get the correct 3-body
  interactions for the dissociative water potential.  The two-body parameters
  are all zeroed out.
* `data.knite9.0792` and `data.knite9.4000` - sample water configurations to run

I never made a set of inputs to test the SiO2 or water-silica interactions but
I recall seeing a way to use LAMMPS's Stillinger Weber potential to achieve
the same 3-body cross-terms for Si-O-H described in the Mahadevan paper.

I last tested this code with a February 2014 release of LAMMPS so you should
be able to just drop the .cpp and .h files into the lammps/src directory and
do 'make' to compile them in.  I should stress again that this code has not
been used much though, so there may be bugs.  Hopefully it is a helpful
starting point though.
