/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

#include "fix_LubricationFLD.h"
#include "atom.h"
#include "fix_property_atom.h"
#include "fix_scalar_transport_equation.h"
#include "force.h"
#include "group.h"
#include "math_extra.h"
#include "modify.h"
#include "pair_gran.h"
#include <stdlib.h>
#include <random>
#include "update.h"
#include "mpi_liggghts.h" 
#include "domain.h"

#include "neighbor.h"
#include "neigh_list.h"
#include "pair_gran.h"



using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

Fix_LubricationFLD::Fix_LubricationFLD(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg){

  if ((!atom->radius_flag)||(!atom->rmass_flag)) error->all(FLERR,"Fix heat/gran needs per particle radius and mass");

  if (narg < 13)
    error->fix_error(FLERR,this,"not enough arguments");

  int iarg = 3; 
  
  if(strcmp(arg[iarg++],"mu"))
    error->fix_error(FLERR,this,"expecting keyword 'mu'");
  mu = atof(arg[iarg++]);
  
  if(strcmp(arg[iarg++],"flagfld"))
    error->fix_error(FLERR,this,"expecting keyword 'flagfld'");
  flagfld = atoi(arg[iarg++]);
  
  if(strcmp(arg[iarg++],"cutinner"))
    error->fix_error(FLERR,this,"expecting keyword 'cutinner'");
  cutinner = atof(arg[iarg++]);
  
  if(strcmp(arg[iarg++],"cutoff"))
    error->fix_error(FLERR,this,"expecting keyword 'cutoff'");
  cutoff = atof(arg[iarg++]);
  
  if(strcmp(arg[iarg++],"flagVF"))
    error->fix_error(FLERR,this,"expecting keyword 'flagVF'");
  flagVF = atoi(arg[iarg++]);
  
   
}

/* ---------------------------------------------------------------------- */


/* ---------------------------------------------------------------------- */



/* ---------------------------------------------------------------------- */

void Fix_LubricationFLD::init()
{	
 
shearing = flagdeform = 0;
	// checking for deform fix. 
	for (int i = 0; i < modify->nfix; i++){
    if (strcmp(modify->fix[i]->style,"deform") == 0) {
      shearing = flagdeform = 1;
      //if (((FixDeform *) modify->fix[i])->remapflag != Domain::V_REMAP)
        //error->all(FLERR,"Using pair lubricate with inconsistent "         "fix deform remap option");
    }
    }
	
	// calucalting the vol fraction and adjusting the stokes coefficients. 
	double vol_T;
	vol_T = domain->xprd*domain->yprd*domain->zprd;
	
	double volP;
	
	 int nlocal = atom->nlocal;
	
	for (int i = 0; i < nlocal; i++){
	volP += (4.0/3.0)*M_PI*pow(atom->radius[i],3.0);
	}
	
	MPI_Allreduce(&volP,&volP,1,MPI_DOUBLE,MPI_SUM,world);
	 
	 double vol_f = volP/vol_T;
	 
    R0  = 6*M_PI*mu*(1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
    RT0 = 8*M_PI*mu*(1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
    RS0 = 20.0/3.0*M_PI*mu*(1.0 + 3.64*vol_f - 6.95*vol_f*vol_f);

	// intializing the stress tensor. 
  Ef[0][0] = Ef[0][1] = Ef[0][2] = 0.0;
  Ef[1][0] = Ef[1][1] = Ef[1][2] = 0.0;
  Ef[2][0] = Ef[2][1] = Ef[2][2] = 0.0;
	 
	 return; 
}

/* ---------------------------------------------------------------------- */

int Fix_LubricationFLD::setmask()
{
  int mask = 0;
  //mask |= INIT;
  mask |= POST_FORCE;
  //mask |= min_setup;
  return mask;
}

/* ---------------------------------------------------------------------- */



/* ---------------------------------------------------------------------- */



void Fix_LubricationFLD::post_force(int vflag) 
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz,tx,ty,tz,t2x,t2y,t2z;
  double rsq,r,h_sep,beta0,beta1,radi,radj;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3;
  double vt1,vt2,vt3,wt1,wt2,wt3,wdotn;
  double vRS0;
  double vi[3],vj[3],wi[3],wj[3],xl[3],jl[3];
  double a_sq,a_sh,a_pu;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double lamda[3],vstream[3];

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  
  
  pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
  NeighList *list = pair_gran->list;
  
  pair = static_cast<Pair*>(force->pair_match("gran", 0));
  //int newton_pair = pair->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;


// subtract streaming component of velocity, omega, angmom
  // assume fluid streaming velocity = box deformation rate
  // vstream = (ux,uy,uz)
  // ux = h_rate[0]*x + h_rate[5]*y + h_rate[4]*z
  // uy = h_rate[1]*y + h_rate[3]*z
  // uz = h_rate[2]*z
  // omega_new = omega - curl(vstream)/2
  // angmom_new = angmom - I*curl(vstream)/2
  // Ef = (grad(vstream) + (grad(vstream))^T) / 2

 /*double **v = atom->v;
 int nlocal = atom->nlocal;
 int i;

	for (i = 0; i < nlocal; i++){
		v[i][0]=0;
		v[i][1]=0;
	v[i][2]=0;
	}
*/


  if (shearing) {
    double *h_rate = domain->h_rate;
    double *h_ratelo = domain->h_ratelo;

    for (ii = 0; ii < inum; ii++) { // removing the fluid velocity from the particle velocity. 
      i = ilist[ii];
      itype = type[i];
      radi = radius[i];
      domain->x2lamda(x[i],lamda);
      vstream[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] +
        h_rate[4]*lamda[2] + h_ratelo[0];
      vstream[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
      vstream[2] = h_rate[2]*lamda[2] + h_ratelo[2];
      v[i][0] -= vstream[0];
      v[i][1] -= vstream[1];
      v[i][2] -= vstream[2];

      omega[i][0] += 0.5*h_rate[3];
      omega[i][1] -= 0.5*h_rate[4];
      omega[i][2] += 0.5*h_rate[5];
    }
	
    Ef[0][0] = h_rate[0]/domain->xprd;
    Ef[1][1] = h_rate[1]/domain->yprd;
    Ef[2][2] = h_rate[2]/domain->zprd;
    Ef[0][1] = Ef[1][0] = 0.5 * h_rate[5]/domain->yprd;
    Ef[0][2] = Ef[2][0] = 0.5 * h_rate[4]/domain->zprd;
    Ef[1][2] = Ef[2][1] = 0.5 * h_rate[3]/domain->zprd;

  }
  
	wi[0] = omega[i][0];
	wi[1] = omega[i][1];
	wi[2] = omega[i][2];
	
	for (i=0; i<nlocal; i++) // stokes forces on the particles. 
	{
		 radi = radius[i];
		
	if (flagfld) {
      f[i][0] -= R0*radi*v[i][0];
      f[i][1] -= R0*radi*v[i][1];
      f[i][2] -= R0*radi*v[i][2];
      const double radi3 = radi*radi*radi;
      torque[i][0] -= RT0*radi3*wi[0];
      torque[i][1] -= RT0*radi3*wi[1];
      torque[i][2] -= RT0*radi3*wi[2];

	

      if (shearing) {
        vRS0 = -RS0*radi3;
	
	
	Pair::v_tally_tensor(i,i,nlocal,newton_pair,
                       vRS0*Ef[0][0],vRS0*Ef[1][1],vRS0*Ef[2][2],
                      vRS0*Ef[0][1],vRS0*Ef[0][2],vRS0*Ef[1][2]);
              
      }
	
	}
	}
	
	

// calculating lubrication forces and updatin in the preforce step. 
	for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    radi = radius[i];

    // angular velocity

    wi[0] = omega[i][0];
    wi[1] = omega[i][1];
    wi[2] = omega[i][2];

    
    



for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      radj = atom->radius[j];
		  r = sqrt(rsq);
      if (r<(cutoff*(radi+radj))) {
      

        // angular momentum = I*omega = 2/5 * M*R^2 * omega

        wj[0] = omega[j][0];
        wj[1] = omega[j][1];
        wj[2] = omega[j][2];

        // xl = point of closest approach on particle i from its center

        xl[0] = -delx/r*radi;
        xl[1] = -dely/r*radi;
        xl[2] = -delz/r*radi;
        jl[0] = -delx/r*radj;
        jl[1] = -dely/r*radj;
        jl[2] = -delz/r*radj;

        // velocity at the point of closest approach on both particles
        // v = v + omega_cross_xl - Ef.xl

        // particle i

        vi[0] = v[i][0] + (wi[1]*xl[2] - wi[2]*xl[1])
                        - (Ef[0][0]*xl[0] + Ef[0][1]*xl[1] + Ef[0][2]*xl[2]);

        vi[1] = v[i][1] + (wi[2]*xl[0] - wi[0]*xl[2])
                        - (Ef[1][0]*xl[0] + Ef[1][1]*xl[1] + Ef[1][2]*xl[2]);

        vi[2] = v[i][2] + (wi[0]*xl[1] - wi[1]*xl[0])
                        - (Ef[2][0]*xl[0] + Ef[2][1]*xl[1] + Ef[2][2]*xl[2]);

        // particle j

        vj[0] = v[j][0] - (wj[1]*jl[2] - wj[2]*jl[1])
                        + (Ef[0][0]*jl[0] + Ef[0][1]*jl[1] + Ef[0][2]*jl[2]);

        vj[1] = v[j][1] - (wj[2]*jl[0] - wj[0]*jl[2])
                        + (Ef[1][0]*jl[0] + Ef[1][1]*jl[1] + Ef[1][2]*jl[2]);

        vj[2] = v[j][2] - (wj[0]*jl[1] - wj[1]*jl[0])
                        + (Ef[2][0]*jl[0] + Ef[2][1]*jl[1] + Ef[2][2]*jl[2]);

        // scalar resistances XA and YA

        h_sep = r - radi-radj;

        // if less than the minimum gap use the minimum gap instead

        if (r < cutinner*(radi+radj))
          h_sep = (cutinner-1)*(radi+radj);

        // scale h_sep by radi

        h_sep = h_sep/radi;
        beta0 = radj/radi;
        beta1 = 1.0 + beta0;

	// scalar resistances

        
         a_sq = beta0*beta0/beta1/beta1/h_sep +
            (1.0+7.0*beta0+beta0*beta0)/5.0/pow(beta1,3.0)*log(1.0/h_sep);
          a_sq += (1.0+18.0*beta0-29.0*beta0*beta0+18.0 *
                   pow(beta0,3.0)+pow(beta0,4.0))/21.0/pow(beta1,4.0) *
            h_sep*log(1.0/h_sep);
          a_sq *= 6.0*M_PI*mu*radi;
                //  a_sq=0;
          a_sh = 4.0*beta0*(2.0+beta0+2.0*beta0*beta0)/15.0/pow(beta1,3.0) *
            log(1.0/h_sep);
          a_sh += 4.0*(16.0-45.0*beta0+58.0*beta0*beta0-45.0*pow(beta0,3.0) +
                       16.0*pow(beta0,4.0))/375.0/pow(beta1,4.0) *
            h_sep*log(1.0/h_sep);
                a_sh *= 6.0*M_PI*mu*radi;
          // old invalid eq for pumping term
          // changed 29Jul16 from eq 9.25 -> 9.27 in Kim and Karilla
//          a_pu = beta0*(4.0+beta0)/10.0/beta1/beta1*log(1.0/h_sep);
//          a_pu += (32.0-33.0*beta0+83.0*beta0*beta0+43.0 *
//                   pow(beta0,3.0))/250.0/pow(beta1,3.0)*h_sep*log(1.0/h_sep);
//          a_pu *= 8.0*M_PI*mu*pow(radi,3.0);
          a_pu = 2.0*beta0/5.0/beta1*log(1.0/h_sep);
          a_pu += 2.0*(8.0+6.0*beta0+33.0*beta0*beta0)/125.0/beta1/beta1*
                   h_sep*log(1.0/h_sep);
          a_pu *= 8.0*M_PI*mu*pow(radi,3.0);
         // a_pu=0;
		 
		 

			//a_sq = 6.0*M_PI*mu*radi*(beta0*beta0/beta1/beta1/h_sep);

        // relative velocity at the point of closest approach
        // includes fluid velocity

        vr1 = vi[0] - vj[0];
        vr2 = vi[1] - vj[1];
        vr3 = vi[2] - vj[2];

        // normal component (vr.n)n

        vnnr = (vr1*delx + vr2*dely + vr3*delz)/r;
        vn1 = vnnr*delx/r;
        vn2 = vnnr*dely/r;
        vn3 = vnnr*delz/r;

        // tangential component vr - (vr.n)n

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // force due to squeeze type motion

        fx  = a_sq*vn1;
        fy  = a_sq*vn2;
        fz  = a_sq*vn3;

	
	fx = fx + a_sh*vt1;
        fy = fy + a_sh*vt2;
        fz = fz + a_sh*vt3;

	// add to total force

        f[i][0] -= fx;
        f[i][1] -= fy;
        f[i][2] -= fz;
		
		
		f[j][0] += fx;
		f[j][1] += fy;
		f[j][2] += fz;


        // torque due to this force

          tx = xl[1]*fz - xl[2]*fy;
          ty = xl[2]*fx - xl[0]*fz;
          tz = xl[0]*fy - xl[1]*fx;
		  
		  t2x = jl[1]*fz - jl[2]*fy;
          t2y = jl[2]*fx - jl[0]*fz;  			// torque due to the a_sq and a_sh force on particle j
          t2z = jl[0]*fy - jl[1]*fx;

          torque[i][0] -= tx;
          torque[i][1] -= ty;
          torque[i][2] -= tz;
		  
		  torque[j][0] -= t2x;
		  torque[j][1] -= t2y;
		  torque[j][2] -= t2z;

          // torque due to a_pu

          wdotn = ((wi[0]-wj[0])*delx + (wi[1]-wj[1])*dely +
                   (wi[2]-wj[2])*delz)/r;
          wt1 = (wi[0]-wj[0]) - wdotn*delx/r;
          wt2 = (wi[1]-wj[1]) - wdotn*dely/r;
          wt3 = (wi[2]-wj[2]) - wdotn*delz/r;

          tx = a_pu*wt1;
          ty = a_pu*wt2;
          tz = a_pu*wt3;
		  
		 t2x = -a_pu*wt1;
		 t2y = -a_pu*wt2;
		 t2z = -a_pu*wt3;

          torque[i][0] -= tx;
          torque[i][1] -= ty;
          torque[i][2] -= tz;
		  
		  torque[j][0] -= t2x;
          torque[j][1] -= t2y;
          torque[j][2] -= t2z;

	}
	}
	
	}
// restore streaming component of velocity, omega, angmom


  if (shearing) {
    double *h_rate = domain->h_rate;
    double *h_ratelo = domain->h_ratelo;

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      itype = type[i];
      radi = atom->radius[i];

      domain->x2lamda(x[i],lamda);
      vstream[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] +
        h_rate[4]*lamda[2] + h_ratelo[0];
      vstream[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
      vstream[2] = h_rate[2]*lamda[2] + h_ratelo[2];
      v[i][0] += vstream[0];
      v[i][1] += vstream[1];
      v[i][2] += vstream[2];

      omega[i][0] -= 0.5*h_rate[3];
      omega[i][1] += 0.5*h_rate[4];
      omega[i][2] -= 0.5*h_rate[5];
  	  
    }
  } 

}


/* ---------------------------------------------------------------------- */



