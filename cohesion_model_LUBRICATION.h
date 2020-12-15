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
    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2014-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_LUBRICATION,LUBRICATION,9)
#else

#ifndef COHESION_MODEL_LUBRICATION_H_
#define COHESION_MODEL_LUBRICATION_H_

#include "contact_models.h"
#include "cohesion_model_base.h"
#include <cmath>
#include <algorithm>
#include "math_extra_liggghts.h"
#include "global_properties.h"
#include "fix_property_atom.h"
#include "neighbor.h"
#include "fix.h"
#include "fix_deform.h"
#include "fix_wall.h"



namespace MODEL_PARAMS
{
    
    inline static ScalarProperty* createMinSeparationDistanceRatioLub(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* minSeparationDistanceRatioScalar = MODEL_PARAMS::createScalarProperty(registry, "minSeparationDistanceRatio", caller);
      return minSeparationDistanceRatioScalar;
    }

    inline static ScalarProperty* createMaxSeparationDistanceRatioLub(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* maxSeparationDistanceRatioScalar = MODEL_PARAMS::createScalarProperty(registry, "maxSeparationDistanceRatio", caller);
      return maxSeparationDistanceRatioScalar;
    }

    inline static ScalarProperty* createFluidViscosityLub(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* fluidViscosityScalar = MODEL_PARAMS::createScalarProperty(registry, "fluidViscosity", caller);
      return fluidViscosityScalar;
    }
	inline static ScalarProperty* createShearTermLub(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* ShearTermScalar = MODEL_PARAMS::createScalarProperty(registry, "ShearTerm", caller);
      return ShearTermScalar;
    }
	inline static ScalarProperty* createSqueezeTermLub(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* SqueezeTermScalar = MODEL_PARAMS::createScalarProperty(registry, "SqueezeTerm", caller);
      return SqueezeTermScalar;
    }
	inline static ScalarProperty* createPumpTermLub(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* PumpTermScalar = MODEL_PARAMS::createScalarProperty(registry, "PumpTerm", caller);
      return PumpTermScalar;
    }
	inline static ScalarProperty* createActivateTorqueLub(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* ActivateTorqueScalar = MODEL_PARAMS::createScalarProperty(registry, "ActivateTorque", caller);
      return ActivateTorqueScalar;
    }
	inline static ScalarProperty* createVolCorrectionLub(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* VolCorrectionScalar = MODEL_PARAMS::createScalarProperty(registry, "VolCorrection", caller);
      return VolCorrectionScalar;
    }
	inline static ScalarProperty* createFLDLub(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* FLDScalar = MODEL_PARAMS::createScalarProperty(registry, "FLD", caller);
      return FLDScalar;
    }
}

namespace LIGGGHTS {

namespace ContactModels {

  template<>
  class CohesionModel<COHESION_LUBRICATION> : public CohesionModelBase {

  public:
    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup,class ContactModelBase *cmb) :
      CohesionModelBase(lmp, hsetup, cmb),
      minSeparationDistanceRatio(0.0),
      maxSeparationDistanceRatio(0.0),
      fluidViscosity(0.),
	  ShearTerm(0.0),
	  SqueezeTerm(0.0),
	  ActivateTorque(0.0),
	  PumpTerm(0.0),
	  VolCorrection(0.0),
      history_offset(0.0)
     {
      history_offset = hsetup->add_history_value("contflag", "0");
      
	  if(cmb->is_wall())
        error->warning(FLERR,"Using cohesion model LUBRICATION for walls only supports dry walls");
    
	// Conditions for the Lubrication forces and initialization.
	
  // ensure all particles are finite-size
  // for pair hybrid, should limit test to types using the pair style

  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (radius[i] == 0.0)
      error->one(FLERR,"Pair lubricate/poly requires extended particles");

  /*int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1; */

  // set the isotropic constants that depend on the volume fraction
  // vol_T = total volume

	}

    void registerSettings(Settings& settings)
    {
        settings.registerOnOff("tangential_reduce",tangentialReduce_,false);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

    void connectToProperties(PropertyRegistry & registry)
    {
      
      registry.registerProperty("fluidViscosity", &MODEL_PARAMS::createFluidViscosityLub);
      registry.registerProperty("minSeparationDistanceRatio", &MODEL_PARAMS::createMinSeparationDistanceRatioLub);
      registry.registerProperty("maxSeparationDistanceRatio", &MODEL_PARAMS::createMaxSeparationDistanceRatioLub);
	  registry.registerProperty("ShearTerm", &MODEL_PARAMS::createShearTermLub);
	  registry.registerProperty("SqueezeTerm", &MODEL_PARAMS::createSqueezeTermLub);
	  registry.registerProperty("PumpTerm", &MODEL_PARAMS::createPumpTermLub);
	  registry.registerProperty("ActivateTorque", &MODEL_PARAMS::createActivateTorqueLub);
	  registry.registerProperty("VolCorrection", &MODEL_PARAMS::createVolCorrectionLub);
	  registry.registerProperty("FLD", &MODEL_PARAMS::createFLDLub);


      registry.connect("fluidViscosity", fluidViscosity,"cohesion_model LUBRICATION"); 
      registry.connect("minSeparationDistanceRatio", minSeparationDistanceRatio,"cohesion_model LUBRICATION");
      registry.connect("maxSeparationDistanceRatio", maxSeparationDistanceRatio,"cohesion_model LUBRICATION");
      registry.connect("ShearTerm", ShearTerm,"cohesion_model LUBRICATION");
      registry.connect("SqueezeTerm", SqueezeTerm,"cohesion_model LUBRICATION");
      registry.connect("PumpTerm", PumpTerm,"cohesion_model LUBRICATION");
      registry.connect("ActivateTorque", ActivateTorque,"cohesion_model LUBRICATION");
	  registry.connect("VolCorrection", VolCorrection,"cohesion_model LUBRICATION");
	  registry.connect("FLD", FLD,"cohesion_model LUBRICATION");
	  
	  
      ln1overMinSeparationDistanceRatio = log(1./minSeparationDistanceRatio);

      
      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"cohesion model LUBRICATION");

      neighbor->register_contact_dist_factor(maxSeparationDistanceRatio*1.1); 
      if(maxSeparationDistanceRatio < 1.0)
            error->one(FLERR,"\n\ncohesion model LUBRICATION requires maxSeparationDistanceRatio >= 1.0. Please increase this value.\n");
	  if(maxSeparationDistanceRatio < 1.0)
            error->one(FLERR,"\n\ncohesion model LUBRICATION requires maxSeparationDistanceRatio >= 1.0. Please increase this value.\n");
	  if(!(ShearTerm == 1.0 || ShearTerm ==0.0))
            error->one(FLERR,"\n\ncohesion model LUBRICATION requires ShearTerm to be Boolean. Please change.\n");
	  if(!(SqueezeTerm == 1.0 || SqueezeTerm ==0.0))
            error->one(FLERR,"\n\ncohesion model LUBRICATION requires SqueezeTerm to be Boolean. Please change.\n");
	  if(!(PumpTerm == 1.0 || PumpTerm ==0.0))
            error->one(FLERR,"\n\ncohesion model LUBRICATION requires PumpTerm to be Boolean. Please change.\n");
	  if(!(ActivateTorque == 1.0 || ActivateTorque ==0.0))
            error->one(FLERR,"\n\ncohesion model LUBRICATION requires ActivateTorque to be Boolean. Please change.\n");
	  if(!(VolCorrection == 1.0 || VolCorrection ==0.0))
            error->one(FLERR,"\n\ncohesion model LUBRICATION requires VolCorrection to be Boolean. Please change.\n");
	  if(!(FLD == 1.0 || FLD ==0.0))
            error->one(FLERR,"\n\ncohesion model LUBRICATION requires FLD to be Boolean. Please change.\n");
	
	
	
	//Lammps lubrication part
	
	
  shearing = flagdeform = flagwall = 0;
  for (int i = 0; i < modify->nfix; i++){
     if (strcmp(modify->fix[i]->style,"deform") == 0) {
      shearing = flagdeform = 1;
    }
  }

  vol_T = domain->xprd*domain->yprd*domain->zprd;
  int nlocal = atom->nlocal;

  volP = 0.0;
  for (int i = 0; i < nlocal; i++)
    volP += (4.0/3.0)*M_PI*pow(atom->radius[i],3.0);
  MPI_Allreduce(&volP,&vol_P,1,MPI_DOUBLE,MPI_SUM,world);

  vol_f = vol_P/vol_T;


			
			
			Ef[0][0] = Ef[0][1] = Ef[0][2] = 0.0;
			Ef[1][0] = Ef[1][1] = Ef[1][2] = 0.0;
			Ef[2][0] = Ef[2][1] = Ef[2][2] = 0.0;	

	

	
	}

    inline void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&) {}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {


      double * const contflag = &sidata.contact_history[history_offset];
	  double vi[3],vj[3],wi[3],wj[3],xl[3],jl[3];
      
		
		double fx,fy,fz,tx,ty,tz,t2x,t2y,t2z;
		double beta0,beta1;     // r = sqrt(sidata.rsq); dist is the h_sep in this template; r is distance between two atoms
		double wt1,wt2,wt3,wdotn;
		double vRS0;                                // for FLD and shearing is done the correction in strain terms
		double a_sq,a_sh,a_pu;
	
		double lamda[3],vstream[3];  // only requried for shearing, 
		//apparently lamda is a coordinate system and used for shearing (fix defrom systems) 
	  
	  
	  // 2 cases: (i) no bridge present, (ii) bridge active
      // for this model, bridge is active of present.




      const int i = sidata.i;
      const int j = sidata.j;
      const double radi = sidata.radi;
      const double radj = sidata.is_wall ? radi : sidata.radj;
      const double r = sqrt(sidata.rsq);
	  const double rsq = sqrt(sidata.rsq);
      const double radsum = sidata.radsum;
      // if the gap is less than minSeparationDistance set the h_sep to minSeparationDistance
	  double h_sep = (minSeparationDistanceRatio-1)*(radi+radj);

      
      const double rEff = radi*radj / (radi+radj);

      // in case of contact use h_sep = h_min
      // 
      if (VolCorrection==0) vol_f = 0;

		// set isotropic constants

		
		R0  = 6*M_PI*fluidViscosity*(1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
		RT0 = 8*M_PI*fluidViscosity*(1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
		RS0 = 20.0/3.0*M_PI*fluidViscosity*(1.0 + 3.64*vol_f - 6.95*vol_f*vol_f);
			Ef[0][0] = Ef[0][1] = Ef[0][2] = 0.0;
			Ef[1][0] = Ef[1][1] = Ef[1][2] = 0.0;
			Ef[2][0] = Ef[2][1] = Ef[2][2] = 0.0;	
		
		  if (shearing) {
			double *h_rate = domain->h_rate;
			//double *h_ratelo = domain->h_ratelo;
	
	
			Ef[0][0] = h_rate[0]/domain->xprd;
			Ef[1][1] = h_rate[1]/domain->yprd;
			Ef[2][2] = h_rate[2]/domain->zprd;
			Ef[0][1] = Ef[1][0] = 0.5 * h_rate[5]/domain->yprd;
			Ef[0][2] = Ef[2][0] = 0.5 * h_rate[4]/domain->zprd;
			Ef[1][2] = Ef[2][1] = 0.5 * h_rate[3]/domain->zprd;

			//comm->forward_comm_pair(this);
			}
	  // This section of code adjusts R0/RT0/RS0 if necessary due to changes
  // in the volume fraction as a result of fix deform or moving walls
	  
	  double dims[3], wallcoord;
  if (VolCorrection) // Flag for volume fraction corrections
    if (flagdeform){ // Possible changes in volume fraction
        for (int jk = 0; jk < 3; jk++)
          dims[jk] = domain->prd[jk];
      
      double vol_T = dims[0]*dims[1]*dims[2];
      vol_f = vol_P/vol_T;
      
        R0  = 6*M_PI*fluidViscosity*(1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
        RT0 = 8*M_PI*fluidViscosity*(1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
        RS0 = 20.0/3.0*M_PI*fluidViscosity*(1.0 + 3.64*vol_f - 6.95*vol_f*vol_f);
    }
	  
		
		
      

         double **v = atom->v;
		 
		  
		  
		  
          double  *x_i = atom->x[i];
		  double  *x_j = atom->x[j];
		  
	
          
          // calculate forces, case no collision

          // calculate the lubrication forcese based on Segante kim and Seppo J.karrila
          
          const double R2 = (radi >=radj) ? radi : radj;
          const double R2inv = 1./R2;
		  
		  
          
          

          // calculate vn and vt since not in struct
          const double rinv = 1.0 / r;
          const double dx = sidata.delta[0];
          const double dy = sidata.delta[1];
          const double dz = sidata.delta[2];

	     double const *omega_i = atom->omega[i];
         double const *omega_j = atom->omega[j];
		  
		  // angular velocity

			wi[0] = omega_i[0];
			wi[1] = omega_i[1];
			wi[2] = omega_i[2];
		  
			wj[0] = omega_j[0];
			wj[1] = omega_j[1];
			wj[2] = omega_j[2];

	
			xl[0] = -dx/r*radi;
            xl[1] = -dy/r*radi;
			xl[2] = -dz/r*radi;
			jl[0] = -dx/r*radj;
			jl[1] = -dy/r*radj;
			jl[2] = -dz/r*radj;
	
			// particle i

			vi[0] = v[i][0] + (wi[1]*xl[2] - wi[2]*xl[1])-(Ef[0][0]*xl[0] + Ef[0][1]*xl[1] + Ef[0][2]*xl[2]);
			vi[1] = v[i][1] + (wi[2]*xl[0] - wi[0]*xl[2])-(Ef[1][0]*xl[0] + Ef[1][1]*xl[1] + Ef[1][2]*xl[2]);
			vi[2] = v[i][2] + (wi[0]*xl[1] - wi[1]*xl[0])-(Ef[2][0]*xl[0] + Ef[2][1]*xl[1] + Ef[2][2]*xl[2]);

			//particle j
			vj[0] = v[j][0] -(wj[1]*jl[2] - wj[2]*jl[1])+ (Ef[0][0]*jl[0] + Ef[0][1]*jl[1] + Ef[0][2]*jl[2]);
			vj[1] = v[j][1] -(wj[2]*jl[0] - wj[0]*jl[2])+ (Ef[1][0]*jl[0] + Ef[1][1]*jl[1] + Ef[1][2]*jl[2]);
			vj[2] = v[j][2] -(wj[0]*jl[1] - wj[1]*jl[0])+ (Ef[2][0]*jl[0] + Ef[2][1]*jl[1] + Ef[2][2]*jl[2]);

	
	

	
			// rescale h_sep by (radi +radj)/2
	
			h_sep = h_sep/((radi+radj)/2);
			beta0 =	radj/radi;
			beta1 = 1.0 + beta0;

			//in this program mu =FluidViscosity;
			a_sq=0;
			a_sh=0;
			a_pu=0;	
			if(SqueezeTerm) 
			{
			a_sq = beta0*beta0/beta1/beta1/beta1/h_sep +
            beta0*(1.0+7.0*beta0+beta0*beta0)/5.0/pow(beta1,3.0)*log(1.0/h_sep);
			a_sq += (1.0+18.0*beta0-29.0*beta0*beta0+18.0 *pow(beta0,3.0)+pow(beta0,4.0))/21.0/pow(beta1,4.0) *
            h_sep*log(1.0/h_sep);
			a_sq *= 6.0*M_PI*fluidViscosity*radi;
			}
				
			if(ShearTerm)
			{
			a_sh = 4.0*beta0*(2.0+beta0+2.0*beta0*beta0)/15.0/pow(beta1,3.0) *
            log(1.0/h_sep);
			a_sh += 4.0*(16.0-45.0*beta0+58.0*beta0*beta0-45.0*pow(beta0,3.0) +
                       16.0*pow(beta0,4.0))/375.0/pow(beta1,4.0) *
            h_sep*log(1.0/h_sep);
			a_sh *= 6.0*M_PI*fluidViscosity*radi;
			}
			
			if(PumpTerm)
			{
			a_pu = 2.0*beta0/5.0/beta1*log(1.0/h_sep);
			a_pu += 2.0*(8.0+6.0*beta0+33.0*beta0*beta0)/125.0/beta1/beta1*
                   h_sep*log(1.0/h_sep);
			a_pu *= 8.0*M_PI*fluidViscosity*pow(radi,3.0);
			}

			
          
          



          const double enx = dx * rinv;
          const double eny = dy * rinv;    //i think it is not needed...
          const double enz = dz * rinv;

          // relative translational velocity
           double vr1 =  vi[0]-vj[0];
           double vr2 =  vi[1]-vj[1];
           double vr3 =  vi[2]-vj[2];

          // normal component
          double vn = vr1 * enx + vr2 * eny + vr3 * enz;
          double vn1 = vn * enx;
          double vn2 = vn * eny;
          double vn3 = vn * enz;

          // tangential component
          double vt1 = vr1 - vn1;
          double vt2 = vr2 - vn2;
          double vt3 = vr3 - vn3;

	// force due to squeeze type motion

        fx  = a_sq*vn1;
        fy  = a_sq*vn2;
        fz  = a_sq*vn3;

        // force due to all shear kind of motions

        
        fx = fx + a_sh*vt1;
        fy = fy + a_sh*vt2;
        fz = fz + a_sh*vt3;
        
			if (FLD) {
			fx -= R0*radi*v[i][0];
			fy -= R0*radi*v[i][1];
			fz -= R0*radi*v[i][2];
			double radi3 = radi*radi*radi;
			tx -=  RT0*radi3*wi[0];
			ty -=  RT0*radi3*wi[1];
			tz -=  RT0*radi3*wi[2];
			
			/*if (shearing && vflag_either) {
			vRS0 = -RS0*radi3;
			v_tally_tensor(i,i,nlocal,newton_pair,
                       vRS0*Ef[0][0],vRS0*Ef[1][1],vRS0*Ef[2][2],
                       vRS0*Ef[0][1],vRS0*Ef[0][2],vRS0*Ef[1][2]); //to add virial components
			}*/
		  }

          // relative rotational velocity
          double wr1, wr2, wr3;
          
          if(sidata.is_wall) {
            wr1 = radi * omega_i[0] * rinv;
            wr2 = radi * omega_i[1] * rinv;
            wr3 = radi * omega_i[2] * rinv;
          } else {
            wr1 = (radi * omega_i[0] + radj * omega_j[0]) * rinv;
            wr2 = (radi * omega_i[1] + radj * omega_j[1]) * rinv;
            wr3 = (radi * omega_i[2] + radj * omega_j[2]) * rinv;
          }

          // relative velocities
          const double vtr1 = vt1 - (dz * wr2 - dy * wr3);
          const double vtr2 = vt2 - (dx * wr3 - dz * wr1);
          const double vtr3 = vt3 - (dy * wr1 - dx * wr2);

          // viscous force
          

          // tangential force components
          tx = 0;
		  ty = 0;
		  tz = 0;
		  t2x=0;
		  t2y=0;	
		  t2z=0;
		  
		  
		  if (ActivateTorque) {

          // torques
		  tx = xl[1]*fz - xl[2]*fy;
          ty = xl[2]*fx - xl[0]*fz;  			// torque due to the a_sq and a_sh force
          tz = xl[0]*fy - xl[1]*fx;
			
		  t2x = jl[1]*fz - jl[2]*fy;
          t2y = jl[2]*fx - jl[0]*fz;  			// torque due to the a_sq and a_sh force
          t2z = jl[0]*fy - jl[1]*fx;
	
	
	  // torque due to a_pu
	
	      wdotn = ((wi[0]-wj[0])*dx + (wi[1]-wj[1])*dy +(wi[2]-wj[2])*dz)/r;
          wt1 = (wi[0]-wj[0]) - wdotn*dx/r;
          wt2 = (wi[1]-wj[1]) - wdotn*dy/r;
          wt3 = (wi[2]-wj[2]) - wdotn*dz/r;

			tx += a_pu*wt1;
			ty += a_pu*wt2;
			tz += a_pu*wt3;
 		 t2x += -a_pu*wt1;
		 t2y += -a_pu*wt2;
		 t2z += -a_pu*wt3;

		  }		
	

          const double tor1 = tx;
          const double tor2 = ty;
          const double tor3 = tz;
		  const double tor4 = t2x;
          const double tor5 = t2y;
          const double tor6 = t2z;
          // add to fn, Ft
          //if(tangentialReduce_) sidata.Fn += Fcapilary+FviscN;
          //sidata.Ft += ...

         
         

          sidata.has_force_update = true;

          // return resulting forces
          if(sidata.is_wall) {
            const double area_ratio = sidata.area_ratio;
            i_forces.delta_F[0] -= fx * area_ratio;
            i_forces.delta_F[1] -= fy * area_ratio;
            i_forces.delta_F[2] -= fz * area_ratio;
            i_forces.delta_torque[0] += - tor1 * area_ratio;
            i_forces.delta_torque[1] += - tor2 * area_ratio;
            i_forces.delta_torque[2] += - tor3 * area_ratio;
          } else {
            i_forces.delta_F[0] -= fx;
            i_forces.delta_F[1] -= fy;
            i_forces.delta_F[2] -= fz;
            i_forces.delta_torque[0] += - tor1; // using radius here, not contact radius
            i_forces.delta_torque[1] += - tor2;
            i_forces.delta_torque[2] += - tor3;

            j_forces.delta_F[0] += fx;
            j_forces.delta_F[1] += fy;
            j_forces.delta_F[2] += fz;
            j_forces.delta_torque[0] += - tor4; // using radius here, not contact radius
            j_forces.delta_torque[1] += - tor5;
            j_forces.delta_torque[2] += - tor6;
          }
     
    }

    void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces)
    {

      double * const contflag = &scdata.contact_history[history_offset];
	  double vi[3],vj[3],wi[3],wj[3],xl[3],jl[3];
      
		double fx,fy,fz,tx,ty,tz,t2x,t2y,t2z;
		double beta0,beta1;     // r = sqrt(scdata.rsq); dist is the h_sep in this template; r is distance between two atoms
		
		double wt1,wt2,wt3,wdotn;
		double vRS0;                                // for FLD and shearing is done the correction in strain terms
		double a_sq,a_sh,a_pu;
		
		double lamda[3],vstream[3];  // only requried for shearing, for now we shall ignore the same.
		//apparently lamda is a coordinate system and used for shearing (fix defrom systems) 
	  
	  
	  // 2 cases: (i) no bridge present, (ii) bridge active
      // for this model, bridge is active of present.

      
      // cases (i) and (ii)

      const int i = scdata.i;
      const int j = scdata.j;
      
      const double radi = scdata.radi;
      const double radj = scdata.is_wall ? radi : scdata.radj;
      const double r = sqrt(scdata.rsq);
	  const double rsq = sqrt(scdata.rsq);
      const double radsum = scdata.radsum;
      double h_sep = scdata.is_wall ? r - radi : r - (radi + radj);
      
      
      const double rEff = radi*radj / (radi+radj);
      const double distMax = (maxSeparationDistanceRatio-1)*radi+(maxSeparationDistanceRatio-1)*radj;
		//one is not sure whether we require this 
      // check if liquid bridge exists
      bool bridge_active = true; // this model has only bridge_active or not

      if (r > (maxSeparationDistanceRatio)*(radi+radj)) // case (i) when the bridge doesnot exists.
      {
        bridge_active = false;
      }
      
	  
		if (VolCorrection==0) vol_f = 0;

		// set isotropic constants

		
		R0  = 6*M_PI*fluidViscosity*(1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
		RT0 = 8*M_PI*fluidViscosity*(1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
		RS0 = 20.0/3.0*M_PI*fluidViscosity*(1.0 + 3.64*vol_f - 6.95*vol_f*vol_f);
			Ef[0][0] = Ef[0][1] = Ef[0][2] = 0.0;
			Ef[1][0] = Ef[1][1] = Ef[1][2] = 0.0;
			Ef[2][0] = Ef[2][1] = Ef[2][2] = 0.0;	
		
      if (shearing) {
    double *h_rate = domain->h_rate;
    //double *h_ratelo = domain->h_ratelo;
	
	
	Ef[0][0] = h_rate[0]/domain->xprd;
    Ef[1][1] = h_rate[1]/domain->yprd;
    Ef[2][2] = h_rate[2]/domain->zprd;
    Ef[0][1] = Ef[1][0] = 0.5 * h_rate[5]/domain->yprd;
    Ef[0][2] = Ef[2][0] = 0.5 * h_rate[4]/domain->zprd;
    Ef[1][2] = Ef[2][1] = 0.5 * h_rate[3]/domain->zprd;

	//comm->forward_comm_pair(this);
	  }
	  
	// This section of code adjusts R0/RT0/RS0 if necessary due to changes
  // in the volume fraction as a result of fix deform or moving walls
	  
	  double dims[3], wallcoord;
  if (VolCorrection) // Flag for volume fraction corrections
    if (flagdeform){ // Possible changes in volume fraction
        for (int jk = 0; jk < 3; jk++)
          dims[jk] = domain->prd[jk];
      
      double vol_T = dims[0]*dims[1]*dims[2];
      vol_f = vol_P/vol_T;
      
        R0  = 6*M_PI*fluidViscosity*(1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
        RT0 = 8*M_PI*fluidViscosity*(1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
        RS0 = 20.0/3.0*M_PI*fluidViscosity*(1.0 + 3.64*vol_f - 6.95*vol_f*vol_f);
    }
	  
	  
      // case (ii)   // when the bridge_is_active
      if(bridge_active)
      {
             double **v = atom->v;
		 
		  
		  
		  
          double  *x_i = atom->x[i];
		  double  *x_j = atom->x[j];
		  
	
          
          // calculate forces, case no collision

          // calculate the lubrication forcese based on Segante kim and Seppo J.karrila
          
           double R2 = (radi >=radj) ? radi : radj;
           double R2inv = 1./R2;
		  
		  
          
          

          // calculate vn and vt since not in struct
           double rinv = 1.0 / r;
           double dx = scdata.delta[0];
           double dy = scdata.delta[1];
           double dz = scdata.delta[2];

	    double  *omega_i = atom->omega[i];
        double  *omega_j = atom->omega[j];
		  
		  // angular velocity

			wi[0] = omega_i[0];
			wi[1] = omega_i[1];
			wi[2] = omega_i[2];
		  
			wj[0] = omega_j[0];
			wj[1] = omega_j[1];
			wj[2] = omega_j[2];

	
		xl[0] = -dx/r*radi;
        xl[1] = -dy/r*radi;
        xl[2] = -dz/r*radi;
        jl[0] = -dx/r*radj;
        jl[1] = -dy/r*radj;
        jl[2] = -dz/r*radj;

	
	
	// particle i

    vi[0] = v[i][0] + (wi[1]*xl[2] - wi[2]*xl[1])-(Ef[0][0]*xl[0] + Ef[0][1]*xl[1] + Ef[0][2]*xl[2]);
	vi[1] = v[i][1] + (wi[2]*xl[0] - wi[0]*xl[2])-(Ef[1][0]*xl[0] + Ef[1][1]*xl[1] + Ef[1][2]*xl[2]);
	vi[2] = v[i][2] + (wi[0]*xl[1] - wi[1]*xl[0])-(Ef[2][0]*xl[0] + Ef[2][1]*xl[1] + Ef[2][2]*xl[2]);

	//particle j
	vj[0] = v[j][0] -(wj[1]*jl[2] - wj[2]*jl[1])+ (Ef[0][0]*jl[0] + Ef[0][1]*jl[1] + Ef[0][2]*jl[2]);
	vj[1] = v[j][1] -(wj[2]*jl[0] - wj[0]*jl[2])+ (Ef[1][0]*jl[0] + Ef[1][1]*jl[1] + Ef[1][2]*jl[2]);
	vj[2] = v[j][2] -(wj[0]*jl[1] - wj[1]*jl[0])+ (Ef[2][0]*jl[0] + Ef[2][1]*jl[1] + Ef[2][2]*jl[2]);

	

	// if the gap is less than minSeparationDistance set the h_sep to minSeparationDistance
	if (r < minSeparationDistanceRatio*(radi+radj))
	 h_sep = (minSeparationDistanceRatio-1)*(radi+radj);

	
	// rescale h_sep by (radi +radj)/2
	
	h_sep = h_sep/((radi+radj)/2);
	beta0 =	radj/radi;
	beta1 = 1.0 + beta0;

	//in this program mu =FluidViscosity;
			a_sq=0;
			a_sh=0;
			a_pu=0;	
			if(SqueezeTerm) 
			{
			a_sq = beta0*beta0/beta1/beta1/beta1/h_sep+
            beta0*(1.0+7.0*beta0+beta0*beta0)/5.0/pow(beta1,3.0)*log(1.0/h_sep);
			a_sq += (1.0+18.0*beta0-29.0*beta0*beta0+18.0 *pow(beta0,3.0)+pow(beta0,4.0))/21.0/pow(beta1,4.0) *
            h_sep*log(1.0/h_sep);
			a_sq *= 6.0*M_PI*fluidViscosity*radi;
			}
				
			if(ShearTerm)
			{
			a_sh = 4.0*beta0*(2.0+beta0+2.0*beta0*beta0)/15.0/pow(beta1,3.0)*log(1.0/h_sep);
			a_sh += 4.0*(16.0-45.0*beta0+58.0*beta0*beta0-45.0*pow(beta0,3.0)+16.0*pow(beta0,4.0))/375.0/pow(beta1,4.0) *
            h_sep*log(1.0/h_sep);
			a_sh *= 6.0*M_PI*fluidViscosity*radi;
			}
			
			if(PumpTerm)
			{
			a_pu = 2.0*beta0/5.0/beta1*log(1.0/h_sep);
			a_pu += 2.0*(8.0+6.0*beta0+33.0*beta0*beta0)/125.0/beta1/beta1*
                   h_sep*log(1.0/h_sep);
			a_pu *= 8.0*M_PI*fluidViscosity*pow(radi,3.0);
			}
			
          const double enx = dx * rinv;
          const double eny = dy * rinv;    //i think it is not needed...
          const double enz = dz * rinv;

          // relative translational velocity
         double vr1 = vi[0]-vj[0];//scdata.v_i[0] - scdata.v_j[0];
         double vr2 = vi[1]-vj[1];//scdata.v_i[1] - scdata.v_j[1];
         double vr3 = vi[2]-vj[2];//scdata.v_i[2] - scdata.v_j[2];

          // normal component
          const double vn = vr1 * enx + vr2 * eny + vr3 * enz;
           double vn1 = vn * enx;
           double vn2 = vn * eny;
           double vn3 = vn * enz;

          // tangential component
           double vt1 = vr1 - vn1;
           double vt2 = vr2 - vn2;
           double vt3 = vr3 - vn3;

	// force due to squeeze type motion

        fx  = a_sq*vn1;
        fy  = a_sq*vn2;
        fz  = a_sq*vn3;

        // force due to all shear kind of motions

        
        fx = fx + a_sh*vt1;
        fy = fy + a_sh*vt2;
        fz = fz + a_sh*vt3;
        
		if (FLD) {
			fx -= R0*radi*v[i][0];
			fy -= R0*radi*v[i][1];
			fz -= R0*radi*v[i][2];
			double radi3 = radi*radi*radi;
			tx -=  RT0*radi3*wi[0];
			ty -=  RT0*radi3*wi[1];
			tz -=  RT0*radi3*wi[2];
			
		}
			
          // relative rotational velocity
          double wr1, wr2, wr3;
          
          if(scdata.is_wall) {
            wr1 = radi * omega_i[0] * rinv;
            wr2 = radi * omega_i[1] * rinv;
            wr3 = radi * omega_i[2] * rinv;
          } else {
            wr1 = (radi * omega_i[0] + radj * omega_j[0]) * rinv;
            wr2 = (radi * omega_i[1] + radj * omega_j[1]) * rinv;
            wr3 = (radi * omega_i[2] + radj * omega_j[2]) * rinv;
          }

          // relative velocities
          const double vtr1 = vt1 - (dz * wr2 - dy * wr3);
          const double vtr2 = vt2 - (dx * wr3 - dz * wr1);
          const double vtr3 = vt3 - (dy * wr1 - dx * wr2);

          // viscous force
          
		tx=0;
		ty=0;	
		tz=0;
		t2x=0;
		t2y=0;	
		t2z=0;
		
		
		
		if(ActivateTorque){
          // torques
		  tx = xl[1]*fz - xl[2]*fy;
          ty = xl[2]*fx - xl[0]*fz;  			// torque due to the a_sq and a_sh force
          tz = xl[0]*fy - xl[1]*fx;
		  
		  t2x = jl[1]*fz - jl[2]*fy;
          t2y = jl[2]*fx - jl[0]*fz;  			// torque due to the a_sq and a_sh force
          t2z = jl[0]*fy - jl[1]*fx;
	
	
	  // torque due to a_pu
	
		  wdotn = ((wi[0]-wj[0])*dx + (wi[1]-wj[1])*dy +(wi[2]-wj[2])*dz)/r;
          wt1 = (wi[0]-wj[0]) - wdotn*dx/r;
          wt2 = (wi[1]-wj[1]) - wdotn*dy/r;
          wt3 = (wi[2]-wj[2]) - wdotn*dz/r;

		 tx += a_pu*wt1;
		 ty += a_pu*wt2;
		 tz += a_pu*wt3;
		 
		 t2x += -a_pu*wt1;
		 t2y += -a_pu*wt2;
		 t2z += -a_pu*wt3;
		 
		}	
	

          const double tor1 = tx;
          const double tor2 = ty;
          const double tor3 = tz;
		  const double tor4 = t2x;
          const double tor5 = t2y;
          const double tor6 = t2z;

          // add to fn, Ft
          //if(tangentialReduce_) scdata.Fn += Fcapilary+FviscN;
          //scdata.Ft += ...

         
         

          scdata.has_force_update = true;

          // return resulting forces
          if(scdata.is_wall) {
            const double area_ratio = scdata.area_ratio;
            i_forces.delta_F[0] -= fx * area_ratio;
            i_forces.delta_F[1] -= fy * area_ratio;
            i_forces.delta_F[2] -= fz * area_ratio;
            i_forces.delta_torque[0] += - tor1 * area_ratio;
            i_forces.delta_torque[1] += - tor2 * area_ratio;
            i_forces.delta_torque[2] += - tor3 * area_ratio;
          } else {
            i_forces.delta_F[0] -= fx;
            i_forces.delta_F[1] -= fy;
            i_forces.delta_F[2] -= fz;
            i_forces.delta_torque[0] += - tor1; // using radius here, not contact radius
            i_forces.delta_torque[1] += - tor2;
            i_forces.delta_torque[2] += - tor3;

            j_forces.delta_F[0] += fx;
            j_forces.delta_F[1] += fy;
            j_forces.delta_F[2] += fz;
            j_forces.delta_torque[0] += - tor4; // using radius here, not contact radius
            j_forces.delta_torque[1] += - tor5;
            j_forces.delta_torque[2] += - tor6;
          }
    
	  }
      // 
     // end of case (iii)
    }

  private:
    double minSeparationDistanceRatio, maxSeparationDistanceRatio, fluidViscosity,ShearTerm,SqueezeTerm,ActivateTorque,PumpTerm,VolCorrection,FLD; 
    double ln1overMinSeparationDistanceRatio;
	double shearing,flagdeform,flagwall;
  double wallfix,vol_T;
  double wallcoord;
  double Ef[3][3];
  double R0,RT0,RS0;
  double vol_P,vol_f;
  double volP;
    int history_offset;
    bool tangentialReduce_;
  };
}
}
#endif // COHESION_MODEL_Lubrication_poly
#endif
