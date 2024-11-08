/* ----------------------------------------------------------------------
 *   This is the MLIP-LAMMPS interface
 *   MLIP is a software for Machine Learning Interatomic Potentials
 *   distributed by A. Shapeev, Skoltech (Moscow)
 *   Contributors: Evgeny Podryabinkin

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
   LAMMPS is distributed under a GNU General Public License
   and is not a part of MLIP.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Evgeny Podryabinkin (Skoltech)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <vector>
#include "pair_MLIP.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"

using namespace LAMMPS_NS;


#define MAXLINE 1024


/* ---------------------------------------------------------------------- */

PairMLIP::PairMLIP(LAMMPS *lmp) : Pair(lmp)
{
  // Edit1 begin
  centroidstressflag = CENTROID_AVAIL;
  // Edit1 end
  restartinfo = 0;
  manybody_flag = 1;

  single_enable = 0;

  inited = false;
  allocated = 0;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairMLIP::~PairMLIP()
{
  if (copymode) return;

  if (allocated) {
      memory->destroy(setflag);
      memory->destroy(cutsq);
  }
  
  if (inited) MLIP_finalize();
}

/* ---------------------------------------------------------------------- */

void PairMLIP::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

  if (mode == 0) // nbh version
  {
    double energy = 0;
    double *p_site_en = NULL;
    double **p_site_virial = NULL;
    if (eflag_atom) p_site_en = &eatom[0];
    if (vflag_atom) p_site_virial = vatom;
	// Edit2 begin
    if (cvflag_atom) p_site_virial = cvatom;
	// Edit2 end

    MLIP_calc_nbh(list->inum, 
		  list->ilist, 
		  list->numneigh, 
		  list->firstneigh,
                  atom->nlocal,
		  atom->nghost,
		  atom->x, 
		  atom->type,
		  atom->f, 
		  energy,
		  p_site_en,      // if NULL no site energy is calculated
		  p_site_virial); // if NULL no virial stress per atom is calculated

    if (eflag_global) eng_vdwl += energy;
    if (vflag_fdotr) virial_fdotr_compute();
  }
  else
  {
    double lattice[9];
    lattice[0] = domain->xprd;
    lattice[1] = 0.0;
    lattice[2] = 0.0;
    lattice[3] = domain->xy;
    lattice[4] = domain->yprd;
    lattice[5] = 0.0;
    lattice[6] = domain->xz;
    lattice[7] = domain->yz;
    lattice[8] = domain->zprd;
                          
    double en = 0.0;
    double virstr[9] = {0.0,0.0,0.0 , 0.0,0.0,0.0 , 0.0,0.0,0.0};
                            
//                                printf("PairMLIP::compute begin\n");
//                                //printf("\tPairMLIP:: CutOff = %f\n", cutoff);
//                                printf("\tPairMLIP:: inum = %d\n", list->inum);
//                                printf("\tPairMLIP:: x[0] = %e %e %e\n", x[0][0], x[0][1], x[0][2]);
//                                printf("\tPairMLIP:: lat = %e %e %e   %e %e %e   %e %e %e\n", lattice[0], lattice[1], lattice[2], lattice[3], lattice[4], lattice[5], lattice[6], lattice[7], lattice[8]);
                            
    MLIP_calc_cfg(list->inum, lattice, atom->x, atom->type, atom->tag, en, atom->f, virstr);
                             
//                                  printf("\tPairMLIP:: en = %f\n", en);
//                                  printf("\tPairMLIP:: inum = %d\n", list->inum);
//                                  printf("\tPairMLIP:: x[0] = %f %f %f\n", x[0][0], x[0][1], x[0][2]);
//                                  printf("\tPairMLIP:: f[0] = %f %f %f\n", f[0][0], f[0][1], f[0][2]);
//                                  printf("\tPairMLIP:: lat = %e %e %e   %e %e %e   %e %e %e\n", lattice[0], lattice[1], lattice[2], lattice[3], lattice[4], lattice[5], lattice[6], lattice[7], lattice[8]);
//                                  printf("PairMLIP::compute end\n");
                              
    if (eflag)
      eng_vdwl += en;
                                    
    if (vflag)
    {
      virial[0] += virstr[0];
      virial[1] += virstr[4];
      virial[2] += virstr[8];
      virial[3] += (virstr[1] + virstr[3]) / 2;
      virial[4] += (virstr[2] + virstr[6]) / 2;
      virial[5] += (virstr[5] + virstr[7]) / 2;
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMLIP::allocate()
{
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      setflag[i][j] = 1;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  allocated = 1;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMLIP::settings(int narg, char **arg)
{
  if (narg != 1 && narg != 2) 
    error->all(FLERR, "Illegal pair_style command");

  if (strlen(arg[0]) > 999)
    error->all(FLERR, "MLIP settings file name is too long");

  strcpy(MLIPsettings_filename, arg[0]);
  if (narg == 2)
    strcpy(MLIPlog_filename, arg[1]);
  else
    strcpy(MLIPlog_filename, "");
}

/* ----------------------------------------------------------------------
   set flags for type pairs
------------------------------------------------------------------------- */

void PairMLIP::coeff(int narg, char **arg)
{
  if (strcmp(arg[0],"*") || strcmp(arg[1],"*") )
    error->all(FLERR, "Incorrect args for pair coefficients");

  if (!allocated) allocate();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairMLIP::init_style()
{
  if (force->newton_pair != 1)
      error->all(FLERR, "Pair style MLIP requires Newton pair on");

  if (inited)
    MLIP_finalize();
  
  if (MLIPlog_filename[0] != '\0')
    MLIP_init(MLIPsettings_filename, MLIPlog_filename, atom->ntypes, cutoff, mode);
  else
    MLIP_init(MLIPsettings_filename, NULL, atom->ntypes, cutoff, mode);

  cutoffsq = cutoff*cutoff;
  int n = atom->ntypes;
  for (int i=1; i<=n; i++)
    for (int j=1; j<=n; j++)
      cutsq[i][j] = cutoffsq;

  if (comm->nprocs != 1 && mode == 1)
    error->all(FLERR, "MLIP settings are incompatible with parallel LAMMPS mode");
  inited = true;

  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMLIP::init_one(int i, int j)
{
  return cutoff;
}
