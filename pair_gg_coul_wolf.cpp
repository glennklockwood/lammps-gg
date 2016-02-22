/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Glenn K. Lockwood, glock@rutgers.edu
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_gg_coul_wolf.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairGGCoulWolf::PairGGCoulWolf(LAMMPS *lmp) : Pair(lmp) {}

/* ---------------------------------------------------------------------- */

PairGGCoulWolf::~PairGGCoulWolf()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(arep);
    memory->destroy(csix);
    memory->destroy(xi);
    memory->destroy(xir);
    memory->destroy(invtwoxi);
    memory->destroy(invsqtwoxi);
    memory->destroy(invtwoxir);
    memory->destroy(eshift);
    memory->destroy(fshift);
  }
}

/* ---------------------------------------------------------------------- */

void PairGGCoulWolf::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double rsq,r2inv,r6inv,forcecoul,forcegg,factor_coul,factor_lj;
  double prefactor;
  double r;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double potgg,e_self,self_con,qisq;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;
  double picon = 2.0 / MY_PIS;

  // temporary GG-specific registers 
  double erbet2, xtrval, rdd, rdd2, erfval, rmix, rmix2;
  double uqq, fqq, udd, fdd, umix, fmix, udisp, fdisp, urep, ggz, ggz2;
  double invtwoxi_, invsqtwoxi_;

#if defined(DWP_DEBUG)
#define TERM_QQ         0
#define TERM_QDQD       1
#define TERM_QQD        2
#define TERM_QDQ        3
#define TERM_REP        4
#define TERM_DISP       5
#define TERM_SELF       6
#define TERM_SHIFT      7
#define TERM_MAX        8
  const char *labels[] = {
    "qq", "qdqd", "qqd", "qdq", "rep", "disp", "self", "shift" };
  double uterms[TERM_MAX], fterms[TERM_MAX];
  for ( ii = 0; ii < TERM_MAX; ii++ )
  {
    uterms[ii] = 0.0;
    fterms[ii] = 0.0;
  }
#endif
  // self and shifted coulombic energy
  self_con = -0.5*erfc(alf*cut_coul)/cut_coul - alf/MY_PIS;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
 
  // loop over neighbors of my atoms
#if defined(DWP_DEBUG)
  fterms[TERM_SELF] = 0.0;
#endif
  for (ii = 0; ii < inum; ii++) 
  {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    qisq = qtmp*qtmp;
    e_self = self_con * qisq * qqrd2e;
#if defined(DWP_DEBUG)
    uterms[TERM_SELF] += e_self;
#endif

    if (evflag) 
        ev_tally(i,i,nlocal,0,0.0,e_self,0.0,0.0,0.0,0.0);

    for (jj = 0; jj < jnum; jj++) 
    {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) 
      {
	r2inv = 1.0/rsq;

	if (rsq < cut_coulsq) 
        {
	  r = sqrt(rsq);
          prefactor = qqrd2e*qtmp*q[j];
          xtrval = erfc(alf*r);
          erbet2 = exp(-alf*alf*r*r);
          invtwoxi_ = invtwoxi[itype][jtype];
          invsqtwoxi_ = invsqtwoxi[itype][jtype];

          // point-point interactions
          uqq = xtrval / r;
          fqq = uqq + picon*erbet2*alf;

          // diffuse-diffuse interactions
          rdd = r * invtwoxi_;
          rdd2 = - rdd * rdd;
          erfval = 1.0 - erfc(rdd);
          udd = 0.06250 * uqq * erfval;
          fdd = udd 
              - 0.06250*(exp(rdd2)*xtrval*picon*invtwoxi_
              - erfval*picon*erbet2*alf);

          // mixed point/diffuse interactions
          rmix = r * invsqtwoxi[itype][jtype];
          rmix2 = - rmix*rmix;
          erfval = 1.0 - erfc(rmix);
          umix = -0.25 * uqq * erfval;
          fmix = umix
              + 0.25*(exp(rmix2)*xtrval*picon*invsqtwoxi_
              - picon*erfval*alf*erbet2);

          potgg = (uqq + 2.0*umix + udd - eshift[itype][jtype]) * prefactor;
#if defined(DWP_DEBUG)
          fterms[TERM_QQ]    += fqq  * prefactor / r;
          fterms[TERM_QDQ]   += fmix * prefactor / r;
          fterms[TERM_QQD]   += fmix * prefactor / r;
          fterms[TERM_QDQD]  += fdd  * prefactor / r;
          fterms[TERM_SHIFT] += fshift[itype][jtype]*prefactor;

          uterms[TERM_QQ]    += uqq  * prefactor;
          uterms[TERM_QDQ]   += umix * prefactor;
          uterms[TERM_QQD]   += umix * prefactor;
          uterms[TERM_QDQD]  += udd  * prefactor;
          uterms[TERM_SHIFT] += eshift[itype][jtype]*prefactor;
#endif
          // forcecoul is missing a factor of 1/r at this point
          forcecoul = (fqq + 2.0*fmix + fdd - fshift[itype][jtype]*r) * prefactor;

	  if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;

          // dispersion shares the same cutoff as Coulomb
          r6inv = r2inv*r2inv*r2inv;
          udisp = csix[itype][jtype] * r6inv;
          fdisp = 6.0 * udisp;          // fdisp is missing a factor of 1/r
#if defined(DWP_DEBUG)
          if ( udisp != 0.0 ) printf( "r6inv: %15.4f %15.4f = %15.4f\n",
              r6inv, csix[itype][jtype],udisp);
          uterms[TERM_DISP] += udisp;
          fterms[TERM_DISP] += fdisp / r;
#endif
	}
        else 
        {
          potgg = 0.0;
          forcecoul = 0.0;
          udisp = 0.0;
          fdisp = 0.0;
        }

// attraction should be truncated at 5.0 A
	if (rsq < cut_ljsq[itype][jtype]) 
        {
          r = sqrt(rsq);

          ggz = r * invtwoxir[itype][jtype];
          ggz2 = -1.0 * ggz * ggz;
          erfval = erfc(ggz);
          urep = arep[itype][jtype] * erfval / ggz;
          // forcegg is missing a factor of 1/r
          forcegg = urep + picon * arep[itype][jtype] * exp(ggz2);
#if defined(DWP_DEBUG)
          fterms[TERM_REP] += forcegg / r;
          uterms[TERM_REP] += urep;
#endif
	}
        else 
        {
          urep = 0.0;
          forcegg = 0.0;
        }

        forcegg += fdisp;
	fpair = (forcecoul + factor_lj*forcegg) * r2inv;

// apply forces to i and j
	f[i][0] += delx*fpair;
	f[i][1] += dely*fpair;
	f[i][2] += delz*fpair;
	if (newton_pair || j < nlocal) 
        {
	  f[j][0] -= delx*fpair;
	  f[j][1] -= dely*fpair;
	  f[j][2] -= delz*fpair;
	}
	if (eflag) 
        {
          if (rsq < cut_coulsq) 
          {
            ecoul = potgg;
	    if (factor_coul < 1.0) ecoul *= (1.0-factor_coul);
          }
          else
            ecoul = 0.0;
          if (rsq < cut_ljsq[itype][jtype])
          {
            evdwl = urep + udisp;

            evdwl *= factor_lj;
          } else evdwl = 0.0;
	}

        if ( evflag )
          ev_tally(i,j,nlocal,newton_pair,evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
  }

#if defined(DWP_DEBUG)
  if ( alf > 1.0 )      // units are in CGS
  {
    double usum = 0.0, fsum = 0.0;
    for ( ii = 0; ii < TERM_MAX; ii++ )
    {
      usum += uterms[ii];
      fsum += fterms[ii];
      printf( "%10s (%1d) %13.4E %13.4E      %13.4E %13.4E\n", 
        labels[ii], ii+1,
        uterms[ii], uterms[ii]/inum,
        fterms[ii], fterms[ii]/inum );
    }
    printf(   "%10s (%1d) %13.4E %13.4E      %13.4E %13.4E\n\n", 
      "Sum", 0, 
      usum, usum/inum, 
      fsum, fsum/inum );
  }
  else
  {
    double usum = 0.0, fsum = 0.0;
    for ( ii = 0; ii < TERM_MAX; ii++ )
    {
      usum += uterms[ii];
      fsum += fterms[ii];
      printf( "%10s (%1d) %13.4f %13.4f      %13.4f %13.4f\n", 
        labels[ii], ii+1,
        uterms[ii], uterms[ii]/inum,
        fterms[ii], fterms[ii]/inum );
    }
    printf(   "%10s (%1d) %13.4f %13.4f      %13.4f %13.4f\n\n", 
      "Sum", 0, 
      usum, usum/inum, 
      fsum, fsum/inum );
  }
#endif
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairGGCoulWolf::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,       "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,         "pair:cutsq");
  memory->create(cut_lj,n+1,n+1,        "pair:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,      "pair:cut_ljsq");
  memory->create(arep,n+1,n+1,          "pair:arep");
  memory->create(csix,n+1,n+1,          "pair:csix");
  memory->create(xi,n+1,n+1,            "pair:xi");
  memory->create(xir,n+1,n+1,           "pair:xir");
  memory->create(invtwoxi,n+1,n+1,      "pair:invtwoxi");
  memory->create(invsqtwoxi,n+1,n+1,    "pair:invsqtwoxi");
  memory->create(invtwoxir,n+1,n+1,     "pair:invtwoxir");
  memory->create(eshift,n+1,n+1,        "pair:eshift");
  memory->create(fshift,n+1,n+1,        "pair:fshift");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGGCoulWolf::settings(int narg, char **arg)
{
  if (narg < 2 || narg > 3) error->all(FLERR,"Illegal pair_style command");

  alf = force->numeric(FLERR,arg[0]);
  cut_lj_global = force->numeric(FLERR,arg[1]);
  if (narg == 2) cut_coul = cut_lj_global;
  else cut_coul = force->numeric(FLERR,arg[2]);

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j]) cut_lj[i][j] = cut_lj_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGGCoulWolf::coeff(int narg, char **arg)
{
// narg = 6 or narg = 7
  if (narg < 6 || narg > 7) 
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double arep_one = force->numeric(FLERR,arg[2]);
  double csix_one = force->numeric(FLERR,arg[3]);
  double xi_one   = force->numeric(FLERR,arg[4]);
  double xir_one  = force->numeric(FLERR,arg[5]);

  double cut_lj_one = cut_lj_global;
  if (narg == 7) cut_lj_one = force->numeric(FLERR,arg[6]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) 
  {
    for (int j = MAX(jlo,i); j <= jhi; j++) 
    {
      arep[i][j] = arep_one;
      csix[i][j] = csix_one;
      xi[i][j] = xi_one;
      xir[i][j] = xir_one;
      cut_lj[i][j] = cut_lj_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGGCoulWolf::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style gg/coul/wolf requires atom attribute q");

  int irequest = neighbor->request(this);

  cut_coulsq = cut_coul * cut_coul;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGGCoulWolf::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  double erbet2, xtrval, rdd, rdd2, erfval, rmix, rmix2;
  double uqq, fqq, udd, fdd, umix, fmix, r;
  double picon = 2.0 / MY_PIS;

  double cut = MAX(cut_lj[i][j],cut_coul);

  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];

  cut_ljsq[j][i] = cut_ljsq[i][j];
  arep[j][i] = arep[i][j];
  csix[j][i] *= -1.0;
  csix[j][i] = csix[i][j];
  xi[j][i] = xi[i][j];
  xir[j][i] = xir[i][j];

  invtwoxi[j][i] = 0.5 / xi[j][i];
  invsqtwoxi[j][i] = sqrt(0.5) / xi[j][i];
  invtwoxir[j][i] = 0.5 / xir[j][i];

  invtwoxi[i][j] = invtwoxi[j][i];
  invsqtwoxi[i][j] = invsqtwoxi[j][i];
  invtwoxir[i][j] = invtwoxir[j][i];

// calculate the Wolf shift in coulomb potential+force
  r = sqrt(cut_coulsq);
  xtrval = erfc(alf*r);
  erbet2 = exp(-alf*alf*r*r);

  // point-point interactions
  uqq = xtrval / r;
  fqq = uqq + picon*erbet2*alf;

  // diffuse-diffuse interactions
  rdd = r * invtwoxi[j][i];
  rdd2 = - rdd * rdd;
  erfval = 1.0 - erfc(rdd);
  udd = 0.06250 * uqq * erfval;
  fdd = udd 
    - 0.06250*(exp(rdd2)*xtrval*picon*invtwoxi[j][i]
      - erfval*picon*erbet2*alf);

  // mixed point/diffuse interactions
  rmix = r * invsqtwoxi[j][i];
  rmix2 = - rmix*rmix;
  erfval = 1.0 - erfc(rmix);
  umix = -0.25 * uqq * erfval;
  fmix = umix
      + 0.25*(exp(rmix2)*xtrval*picon*invsqtwoxi[j][i]
      - picon*erfval*alf*erbet2);

  eshift[j][i] = uqq + 2.0*umix + udd;
  fshift[j][i] = (fqq + 2.0*fmix + fdd)/r;
  eshift[i][j] = eshift[j][i];
  fshift[i][j] = fshift[j][i];

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGGCoulWolf::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
	fwrite(&arep[i][j],sizeof(double),1,fp);
	fwrite(&csix[i][j],sizeof(double),1,fp);
	fwrite(&xi[i][j],sizeof(double),1,fp);
	fwrite(&xir[i][j],sizeof(double),1,fp);
	fwrite(&cut_lj[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGGCoulWolf::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
	if (me == 0) {
	  fread(&arep[i][j],sizeof(double),1,fp);
	  fread(&csix[i][j],sizeof(double),1,fp);
	  fread(&xi[i][j],sizeof(double),1,fp);
	  fread(&xir[i][j],sizeof(double),1,fp);
	  fread(&cut_lj[i][j],sizeof(double),1,fp);
	}
	MPI_Bcast(&arep[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&csix[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&xi[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&xir[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGGCoulWolf::write_restart_settings(FILE *fp)
{
  fwrite(&alf,sizeof(double),1,fp);
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGGCoulWolf::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&alf,sizeof(double),1,fp);
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&alf,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   only the pair part is calculated here
------------------------------------------------------------------------- */

double PairGGCoulWolf::single(int i, int j, int itype, int jtype,
				double rsq, 
				double factor_coul, double factor_lj,
				double &fforce)
{
  double r2inv,r6inv,r,prefactor;
  double forcecoul,forcegg;
  double ecoul, evdwl;
  // temporary GG-specific registers 
  double erbet2, xtrval, rdd, rdd2, erfval, rmix, rmix2;
  double uqq, fqq, udd, fdd, umix, fmix, udisp, fdisp, urep, ggz, ggz2;
  double invtwoxi_, invsqtwoxi_;
  double picon = 2.0 / MY_PIS;

  r2inv = 1.0/rsq;

  if (rsq < cut_coulsq) 
  {
    r = sqrt(rsq);
    prefactor = force->qqrd2e * atom->q[i] * atom->q[j];
    xtrval = erfc(alf*r);
    erbet2 = exp(-alf*alf*r*r);
    invtwoxi_ = invtwoxi[itype][jtype];
    invsqtwoxi_ = invsqtwoxi[itype][jtype];

    // point-point interactions
    uqq = xtrval / r;
    fqq = uqq + picon*erbet2*alf;

    // diffuse-diffuse interactions
    rdd = r * invtwoxi_;
    rdd2 = - rdd * rdd;
    erfval = 1.0 - erfc(rdd);
    udd = 0.06250 * uqq * erfval;
    fdd = udd - 0.06250*(exp(rdd2)*xtrval*picon*invtwoxi_
        - erfval*picon*erbet2*alf);

    // mixed point/diffuse interactions
    rmix = r * invsqtwoxi_;
    rmix2 = - rmix*rmix;
    erfval = 1.0 - erfc(rmix);
    umix = -0.25 * uqq * erfval;
    fmix = umix
         + 0.25*(exp(rmix2)*xtrval*picon*invsqtwoxi_
         - picon*erfval*alf*erbet2);

    ecoul     = (uqq + 2.0*umix + udd - eshift[itype][jtype]) * prefactor;
    // forcecoul is missing a factor of 1/r at this point
    forcecoul = (fqq + 2.0*fmix + fdd - fshift[itype][jtype]*r) * prefactor;

    if (factor_coul < 1.0) 
    {
      forcecoul -= (1.0-factor_coul)*prefactor;
      ecoul -= (1.0-factor_coul)*prefactor;
    }

    // dispersion shares the same cutoff as Coulomb
    r6inv = r2inv*r2inv*r2inv;
    udisp = csix[itype][jtype] * r6inv;
    fdisp = 6.0 * udisp;          // fdisp is missing a factor of 1/r
  }
  else 
  {
    ecoul = 0.0;
    forcecoul = 0.0;
    udisp = 0.0;
    fdisp = 0.0;
  }

// attraction should be truncated at 5.0 A
  if (rsq < cut_ljsq[itype][jtype]) 
  {
    r = sqrt(rsq);

    ggz = r * invtwoxir[itype][jtype];
    ggz2 = -1.0 * ggz * ggz;
    erfval = erfc(ggz);
    urep = arep[itype][jtype] * erfval / ggz;
    // forcegg is missing a factor of 1/r
    forcegg = urep + picon * arep[itype][jtype] * exp(ggz2);
    evdwl = factor_lj*(urep + udisp);
  }
  else 
  {
    urep = 0.0;
    forcegg = 0.0;
    evdwl = 0.0;
  }

  forcegg += fdisp;
  fforce = (forcecoul + factor_lj*forcegg) * r2inv;

  return ecoul+evdwl;
}
