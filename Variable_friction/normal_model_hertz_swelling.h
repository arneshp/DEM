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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef NORMAL_MODEL
NORMAL_MODEL(HERTZ_SWELLING,hertz_swelling,13)
#else
#ifndef NORMAL_MODEL_HERTZ_SWELLING_H_
#define NORMAL_MODEL_HERTZ_SWELLING_H_
#include "global_properties.h"
#include "fix_property_atom.h"
#include <cmath>
#include "normal_model_base.h"
#include "fix_mesh_surface.h"

namespace MODEL_PARAMS
{
    
    inline static ScalarProperty* createRoughnessRatioswell(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      ScalarProperty* RoughnessRatioScalar = MODEL_PARAMS::createScalarProperty(registry, "roughnessRatio", caller);
      return RoughnessRatioScalar;
    }
	
}






namespace LIGGGHTS {

namespace ContactModels
{
  class ContactModelBase;

  template<>
  class NormalModel<HERTZ_SWELLING> : public NormalModelBase
  {
  public:
    NormalModel(LAMMPS * lmp, IContactHistorySetup* hsetup, class ContactModelBase *c) :
      NormalModelBase(lmp, hsetup, c),
      Yh(NULL),
      Yl(NULL),
      nu(NULL),
      betaeff(NULL),
      limitForce(false),
      displayedSettings(false),
      heating(false),
      heating_track(false),
      elastic_potential_offset_(-1),
      elasticpotflag_(false),
      fix_dissipated_(NULL),
      dissipatedflag_(false),
      overlap_offset_(0.0),
      disable_when_bonded_(false),
      bond_history_offset_(-1),
      dissipation_history_offset_(-1),
      cmb(c)
    {
      
    }

    void registerSettings(Settings & settings)
    {
      
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("limitForce", limitForce);
      settings.registerOnOff("heating_normal_hertz",heating,false);
      settings.registerOnOff("heating_tracking",heating_track,false);
      settings.registerOnOff("computeElasticPotential", elasticpotflag_, false);
      settings.registerOnOff("computeDissipatedEnergy", dissipatedflag_, false);
      settings.registerOnOff("disableNormalWhenBonded", disable_when_bonded_, false);
      //TODO error->one(FLERR,"TODO here also check if right surface model used");
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb)
    {
        if (elasticpotflag_)
        {
            elastic_potential_offset_ = cmb->get_history_offset("elastic_potential_normal");
            if (elastic_potential_offset_ == -1)
            {
                elastic_potential_offset_ = hsetup->add_history_value("elastic_potential_normal", "0");
                hsetup->add_history_value("elastic_force_normal_0", "1");
                hsetup->add_history_value("elastic_force_normal_1", "1");
                hsetup->add_history_value("elastic_force_normal_2", "1");
                hsetup->add_history_value("elastic_torque_normal_i_0", "0");
                hsetup->add_history_value("elastic_torque_normal_i_1", "0");
                hsetup->add_history_value("elastic_torque_normal_i_2", "0");
                hsetup->add_history_value("elastic_torque_normal_j_0", "0");
                hsetup->add_history_value("elastic_torque_normal_j_1", "0");
                hsetup->add_history_value("elastic_torque_normal_j_2", "0");
                if (cmb->is_wall())
                    hsetup->add_history_value("elastic_potential_wall", "0");
                cmb->add_history_offset("elastic_potential_normal", elastic_potential_offset_);
            }
        }
        if (dissipatedflag_)
        {
            if (cmb->is_wall())
            {
                fix_dissipated_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dissipated_energy_wall", "property/atom", "vector", 0, 0, "dissipated energy"));
                dissipation_history_offset_ = cmb->get_history_offset("dissipation_force");
                if (!dissipation_history_offset_)
                    error->one(FLERR, "Internal error: Could not find dissipation history offset");
            }
            else
                fix_dissipated_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dissipated_energy", "property/atom", "vector", 0, 0, "dissipated energy"));
            if (!fix_dissipated_)
                error->one(FLERR, "Surface model has not registered dissipated_energy fix");
        }
        if (disable_when_bonded_)
        {
            bond_history_offset_ = cmb->get_history_offset("bond_contactflag");
            if (bond_history_offset_ < 0)
                error->one(FLERR, "Could not find bond history offset");
            overlap_offset_ = hsetup->add_history_value("overlap_offset", "0");
        }
    }

    void connectToProperties(PropertyRegistry & registry)
    {
      
      registry.registerProperty("youngsModulushigh", &MODEL_PARAMS::createYoungsModulushigh,"model hertz_swelling");
      registry.registerProperty("youngsModuluslow", &MODEL_PARAMS::createYoungsModuluslow,"model hertz_swelling");
      registry.registerProperty("poissonsRatio", &MODEL_PARAMS::createPoissonsRatio,"model hertz_swelling");
      registry.registerProperty("betaeff", &MODEL_PARAMS::createBetaEff,"model hertz_swelling");

	  registry.registerProperty("roughnessRatio", &MODEL_PARAMS::createRoughnessRatioswell);

      //registry.connect("Yeff", Yeff,"model hertz");
      //registry.connect("Geff", Geff,"model hertz");
      
      
      registry.connect("youngsModulushigh", Yh, "model hertz_swelling");
      registry.connect("youngsModuluslow", Yl, "model hertz_swelling");
      registry.connect("poissonsRatio", nu, "model hertz_swelling");
      registry.connect("betaeff", betaeff,"model hertz_swelling");
      registry.connect("roughnessRatio", roughnessRatio,"model hertz_swelling");

      // enlarge contact distance flag in case of elastic energy computation
      // to ensure that surfaceClose is called after a contact
      if (elasticpotflag_)
          //set neighbor contact_distance_factor here
          neighbor->register_contact_dist_factor(1.1*roughnessRatio);
    }

    // effective exponent for stress-strain relationship
    
    inline double stressStrainExponent()
    {
      return 1.5;
    }

    void dissipateElasticPotential(SurfacesCloseData &scdata)
    {
        if (elasticpotflag_)
        {
            double * const elastic_energy = &scdata.contact_history[elastic_potential_offset_];
            if (scdata.is_wall)
            {
                // we need to calculate half an integration step which was left over to ensure no energy loss, but only for the elastic energy. The dissipation part is handled in fix_wall_gran_base.h.
                double delta[3];
                scdata.fix_mesh->triMesh()->get_global_vel(delta);
                vectorScalarMult3D(delta, update->dt);
                // -= because force is in opposite direction
                // no *dt as delta is v*dt of the contact position
                elastic_energy[0] -= (delta[0]*(elastic_energy[1]) +
                                      delta[1]*(elastic_energy[2]) +
                                      delta[2]*(elastic_energy[3]))*0.5
                                     // from previous half step
                                     + elastic_energy[10];
                elastic_energy[10] = 0.0;
            }
            elastic_energy[1] = 0.0;
            elastic_energy[2] = 0.0;
            elastic_energy[3] = 0.0;
            elastic_energy[4] = 0.0;
            elastic_energy[5] = 0.0;
            elastic_energy[6] = 0.0;
            elastic_energy[7] = 0.0;
            elastic_energy[8] = 0.0;
            elastic_energy[9] = 0.0;
        }
    }

    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      if (sidata.contact_flags) *sidata.contact_flags |= CONTACT_NORMAL_MODEL;
      const bool update_history = sidata.computeflag && sidata.shearupdate;
      
	
      const int i = sidata.i;
      const int j = sidata.j;
      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      const double radi = sidata.radi;
      const double radj = sidata.radj;
      double reff = sidata.is_wall ? radi : (radi*radj/(radi+radj));

	
	
	
	double Yi,Yj;
	

      #ifdef SUPERQUADRIC_ACTIVE_FLAG
      if(sidata.is_non_spherical && atom->superquadric_flag)
        reff = sidata.reff;
      #endif
      const double meff=sidata.meff;

      if(sidata.deltan < 0)
        error->one(FLERR, "sidata.deltan < 0!");
      const double sqrtval = sqrt(reff*sidata.deltan);
      #ifdef LIGGGHTS_DEBUG
        if(std::isnan(sqrtval))
          error->one(FLERR, "sqrtval is NaN!");
      #endif

      if (disable_when_bonded_ && update_history && sidata.deltan < sidata.contact_history[overlap_offset_])
        sidata.contact_history[overlap_offset_] = sidata.deltan;
      const double deltan = disable_when_bonded_ ? fmax(sidata.deltan-sidata.contact_history[overlap_offset_], 0.0) : sidata.deltan;
      
       
      
      fix_Sdim = static_cast<FixPropertyAtom*>(modify->find_fix_property("Sdim","property/atom","scalar",0,0,"Sdim"));
      Sdim = fix_Sdim->vector_atom; 

	
	
      Yi = Yh[itype]*(1-Sdim[i]) + (Sdim[i])*Yl[itype];
      Yj = Yh[jtype]*(1-Sdim[j]) + (Sdim[j])*Yl[jtype];
      
      
      
      Yeff =1/( (1-nu[itype]*nu[itype])/Yi + (1-nu[jtype]*nu[jtype])/Yj); 
        
      Geff =1/( 2*(2-nu[itype])*(1+nu[itype])/Yi + 2*(2-nu[jtype])*(1+nu[jtype])/Yj);      

      const double Sn=2.*Yeff*sqrtval;
      const double St=8.*Geff*sqrtval;

      double kn=4./3.*Yeff*sqrtval;
      double kt=St;
      const double sqrtFiveOverSix = 0.91287092917527685576161630466800355658790782499663875;
      const double gamman=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(Sn*meff);
      const double gammat= tangential_damping ? -2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(St*meff) : 0.0;
      
      if(!displayedSettings)
      {
        displayedSettings = true;

        /*
        if(limitForce)
            if(0 == comm->me) fprintf(screen," NormalModel<HERTZ_STIFFNESS>: will limit normal force.\n");
        */
      }
      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

      const double Fn_damping = -gamman*sidata.vn;
      const double Fn_contact = kn*deltan;
      double Fn = Fn_damping + Fn_contact;

      //limit force to avoid the artefact of negative repulsion force
      if(limitForce && (Fn<0.0) )
      {
          Fn = 0.0;
      }

      sidata.Fn = Fn;
      sidata.kn = kn;
      sidata.kt = kt;
      sidata.gamman = gamman;
      sidata.gammat = gammat;

      #ifdef NONSPHERICAL_ACTIVE_FLAG
          double torque_i[3] = {0., 0., 0.};
          double Fn_i[3] = { Fn * sidata.en[0], Fn * sidata.en[1], Fn * sidata.en[2]};
          if(sidata.is_non_spherical) {
            double xci[3];
            vectorSubtract3D(sidata.contact_point, atom->x[sidata.i], xci);
            vectorCross3D(xci, Fn_i, torque_i);
          }
          
      #endif

      // apply normal force
      if (!disable_when_bonded_ || sidata.contact_history[bond_history_offset_] < 0.5)
      {

        if(heating)
        {
            const double mj = sidata.is_wall ? sidata.mi : sidata.mj;
            const double E_therm = fabs((-sidata.vn - update->dt*Fn*0.5*(1.0/sidata.mi + 1.0/mj))*Fn_damping);
            sidata.P_diss += E_therm; 
            if(heating_track && sidata.is_wall)
                cmb->tally_pw(E_therm ,sidata.i,jtype,0);
            if(heating_track && !sidata.is_wall)
                cmb->tally_pp(E_therm ,sidata.i,sidata.j,0);
        }

        // energy balance terms
        if (update_history)
        {
            // compute increment in elastic potential
            if (elasticpotflag_)
            {
                double * const elastic_energy = &sidata.contact_history[elastic_potential_offset_];
                // correct for wall influence
                if (sidata.is_wall)
                {
                    double delta[3];
                    sidata.fix_mesh->triMesh()->get_global_vel(delta);
                    vectorScalarMult3D(delta, update->dt);
                    // -= because force is in opposite direction
                    // no *dt as delta is v*dt of the contact position
                    //printf("pela %e %e %e %e\n",  update->get_cur_time()-update->dt, deb, -sidata.radj, deb-sidata.radj);
                    elastic_energy[0] -= (delta[0]*elastic_energy[1] +
                                          delta[1]*elastic_energy[2] +
                                          delta[2]*elastic_energy[3])*0.5
                                         // from previous half step
                                         + elastic_energy[10];
                    elastic_energy[10] = -(delta[0]*Fn_contact*sidata.en[0] +
                                           delta[1]*Fn_contact*sidata.en[1] +
                                           delta[2]*Fn_contact*sidata.en[2])*0.5;
                }
                elastic_energy[1] = -Fn_contact*sidata.en[0];
                elastic_energy[2] = -Fn_contact*sidata.en[1];
                elastic_energy[3] = -Fn_contact*sidata.en[2];
                elastic_energy[4] = 0.0;
                elastic_energy[5] = 0.0;
                elastic_energy[6] = 0.0;
                elastic_energy[7] = 0.0;
                elastic_energy[8] = 0.0;
                elastic_energy[9] = 0.0;
            }
            // compute increment in dissipated energy
            if (dissipatedflag_)
            {
                double * const * const dissipated = fix_dissipated_->array_atom;
                double * const dissipated_i = dissipated[sidata.i];
                double * const dissipated_j = dissipated[sidata.j];
                const double F_diss = -Fn_damping;
                dissipated_i[1] += sidata.en[0]*F_diss;
                dissipated_i[2] += sidata.en[1]*F_diss;
                dissipated_i[3] += sidata.en[2]*F_diss;
                if (sidata.j < atom->nlocal && !sidata.is_wall)
                {
                    dissipated_j[1] -= sidata.en[0]*F_diss;
                    dissipated_j[2] -= sidata.en[1]*F_diss;
                    dissipated_j[3] -= sidata.en[2]*F_diss;
                }
                else if (sidata.is_wall)
                {
                    double * const diss_force = &sidata.contact_history[dissipation_history_offset_];
                    diss_force[0] -= sidata.en[0]*F_diss;
                    diss_force[1] -= sidata.en[1]*F_diss;
                    diss_force[2] -= sidata.en[2]*F_diss;
                }
            }
            #ifdef NONSPHERICAL_ACTIVE_FLAG
            if ((dissipatedflag_ || elasticpotflag_) && sidata.is_non_spherical)
                error->one(FLERR,"Dissipation and elastic potential do not compute torque influence for nonspherical particles");
            #endif
        }

        if(sidata.is_wall) {
          const double Fn_ = Fn * sidata.area_ratio;
          i_forces.delta_F[0] += Fn_ * sidata.en[0];
          i_forces.delta_F[1] += Fn_ * sidata.en[1];
          i_forces.delta_F[2] += Fn_ * sidata.en[2];
          #ifdef NONSPHERICAL_ACTIVE_FLAG
                  if(sidata.is_non_spherical) {
                    //for non-spherical particles normal force can produce torque!
                    i_forces.delta_torque[0] += torque_i[0];
                    i_forces.delta_torque[1] += torque_i[1];
                    i_forces.delta_torque[2] += torque_i[2];
                  }
          #endif
        } else {
          i_forces.delta_F[0] += sidata.Fn * sidata.en[0];
          i_forces.delta_F[1] += sidata.Fn * sidata.en[1];
          i_forces.delta_F[2] += sidata.Fn * sidata.en[2];

          j_forces.delta_F[0] += -i_forces.delta_F[0];
          j_forces.delta_F[1] += -i_forces.delta_F[1];
          j_forces.delta_F[2] += -i_forces.delta_F[2];
          #ifdef NONSPHERICAL_ACTIVE_FLAG
                  if(sidata.is_non_spherical) {
                    //for non-spherical particles normal force can produce torque!
                    double xcj[3], torque_j[3];
                    double Fn_j[3] = { -Fn_i[0], -Fn_i[1], -Fn_i[2]};
                    vectorSubtract3D(sidata.contact_point, atom->x[sidata.j], xcj);
                    vectorCross3D(xcj, Fn_j, torque_j);

                    i_forces.delta_torque[0] += torque_i[0];
                    i_forces.delta_torque[1] += torque_i[1];
                    i_forces.delta_torque[2] += torque_i[2];

                    j_forces.delta_torque[0] += torque_j[0];
                    j_forces.delta_torque[1] += torque_j[1];
                    j_forces.delta_torque[2] += torque_j[2];
                  }
          #endif
        }
      }
      else if (update_history)
      {
        sidata.contact_history[overlap_offset_] = sidata.deltan;
        dissipateElasticPotential(sidata);
      }
    }

    void surfacesClose(SurfacesCloseData &scdata, ForceData& i_forces, ForceData&  j_forces)
    {
		const int i = scdata.i;
		const int j = scdata.j;
		const int itype = scdata.itype;
        const int jtype = scdata.jtype;
		
		const double radi = scdata.radi;
		const double radj = scdata.radj;
		const double r = sqrt(scdata.rsq);
		
		double h = r -(radi+radj);
		
		
		if (r< (roughnessRatio*(radi+radj)))
		{
			 fix_Sdim = static_cast<FixPropertyAtom*>(modify->find_fix_property("Sdim","property/atom","scalar",0,0,"Sdim"));
			 Sdim = fix_Sdim->vector_atom; 

				double Yi,Yj;
	
			  Yi = Yh[itype]*(1-Sdim[i]) + (Sdim[i])*Yl[itype];
			  Yj = Yh[jtype]*(1-Sdim[j]) + (Sdim[j])*Yl[jtype];
			  
			  
			  Yeff =1/( (1-nu[itype]*nu[itype])/Yi + (1-nu[jtype]*nu[jtype])/Yj); 				
			  Geff =1/( 2*(2-nu[itype])*(1+nu[itype])/Yi + 2*(2-nu[jtype])*(1+nu[jtype])/Yj);      
		
			   double Cvi,Cvj,Lcbari, Lcbarj;
			   double deltacriticalbari, deltacriticalbarj;
			   
			   Cvi =1.234 + 1.256*nu[itype];
			   Cvj =1.234 + 1.256*nu[jtype];
				
				Lcbari = 8.88*nu[itype] - 10.13*(nu[itype]*nu[itype] + 0.089);
				Lcbarj = 8.88*nu[jtype] - 10.13*(nu[jtype]*nu[jtype] + 0.089);
				
				deltacriticalbari = 6.82 - 7.83*(nu[itype]*nu[itype] + 0.0586); 
				deltacriticalbarj = 6.82 - 7.83*(nu[jtype]*nu[jtype] + 0.0586); 
		
				
				double deltani,deltanj; 
				double deltacriticali, deltacriticalj;
				double Lci, Lcj;
				double betai, betaj;
				
				
				
				
				
				
				deltani = r- (1+roughnessRatio)*(radi) -radj; 
				deltanj = r- (1+roughnessRatio)*(radj) -radi; 
					 
					 
					 deltacriticali = deltacriticalbari*(roughnessRatio*radi)*(M_PI*Cvi*(1-nu[itype]*nu[itype])*(1/1000));
					 deltacriticalj = deltacriticalbarj*(roughnessRatio*radj)*(M_PI*Cvj*(1-nu[jtype]*nu[jtype])*(1/1000));
					 
					 Lci = Lcbari*pow(M_PI,3)*(Yi/6000)*pow(Cvi,3)*pow(((roughnessRatio*radi)*(2*(1-nu[itype]*nu[itype])/1000)),2);
					 Lcj = Lcbarj*pow(M_PI,3)*(Yj/6000)*pow(Cvj,3)*pow(((roughnessRatio*radj)*(2*(1-nu[jtype]*nu[jtype])/1000)),2);
					 
					 betai = 0.174 + 0.08*nu[itype];
					 betaj = 0.174 + 0.08*nu[jtype];
					 
				double Fni[3],Fnj[3];
				
				if(deltani < deltacriticali)
				{
						Fni[0] = -Lci*pow((fabs(deltani)/deltacriticali),1.5)*scdata.delta[0]/r;
						Fni[1] = -Lci*pow((fabs(deltani)/deltacriticali),1.5)*scdata.delta[1]/r;
						Fni[2] = -Lci*pow((fabs(deltani)/deltacriticali),1.5)*scdata.delta[2]/r;
				} else if ((r/2)< roughnessRatio*(radi)) {
					Fni[0] = -Lci*(pow((fabs(deltani)/deltacriticali),1.5))*(1- exp( 1/(1- pow( (fabs(deltani)/deltacriticali),betai))))*scdata.delta[0]/r;
					Fni[1] = -Lci*(pow((fabs(deltani)/deltacriticali),1.5))*(1- exp( 1/(1- pow( (fabs(deltani)/deltacriticali),betai))))*scdata.delta[1]/r;
					Fni[2] = -Lci*(pow((fabs(deltani)/deltacriticali),1.5))*(1- exp( 1/(1- pow( (fabs(deltani)/deltacriticali),betai))))*scdata.delta[2]/r;
				} else 
					{
						Fni[0] = 0;
						Fni[1] = 1;
						Fni[2] = 2;
					}
				
					 
				if(deltanj < deltacriticalj)
				{
						Fnj[0] = -Lcj*pow((fabs(deltanj)/deltacriticalj),1.5)*scdata.delta[0]/r;
						Fnj[1] = -Lcj*pow((fabs(deltanj)/deltacriticalj),1.5)*scdata.delta[1]/r;
						Fnj[2] = -Lcj*pow((fabs(deltanj)/deltacriticalj),1.5)*scdata.delta[2]/r;
				} else if ((r/2)< roughnessRatio*(radj)) 
				{
					Fnj[0] = -Lcj*(pow((fabs(deltanj)/deltacriticalj),1.5))*(1- exp( 1/(1- pow( (fabs(deltanj)/deltacriticalj),betaj))))*scdata.delta[0]/r;
					Fnj[1] = -Lcj*(pow((fabs(deltanj)/deltacriticalj),1.5))*(1- exp( 1/(1- pow( (fabs(deltanj)/deltacriticalj),betaj))))*scdata.delta[1]/r;
					Fnj[2] = -Lcj*(pow((fabs(deltanj)/deltacriticalj),1.5))*(1- exp( 1/(1- pow( (fabs(deltanj)/deltacriticalj),betaj))))*scdata.delta[2]/r;
				}	 else 
				{
					Fnj[0]=0;
					Fnj[1]=1;
					Fnj[2]=2;
				}
				
				double Fn[3];
				
				Fn[0] = Fni[0] + Fnj[0];
				Fn[1] = Fni[1] + Fnj[1];
				Fn[2] = Fni[2] + Fnj[2];
				
				
				double kt; 
				kt = (2/7)*fabs(pow((Fn[0]*Fn[0] + Fn[1]*Fn[1]+Fn[2]*Fn[2]),0.5)/(deltani+deltanj));
				
				kt /= force->nktv2p;
				//scdata.kt = kt; 
	
          i_forces.delta_F[0] += Fn[0];
          i_forces.delta_F[1] += Fn[1];
          i_forces.delta_F[2] += Fn[2];

          j_forces.delta_F[0] += -i_forces.delta_F[0];
          j_forces.delta_F[1] += -i_forces.delta_F[1];
          j_forces.delta_F[2] += -i_forces.delta_F[2];
          
        
		
		}
		
        if (scdata.contact_flags)
            *scdata.contact_flags |= CONTACT_NORMAL_MODEL;
        dissipateElasticPotential(scdata);
		
    }

    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

  protected:
    
    double ** betaeff;
    double *Yh;
    double *Yl;
    double *nu;
    
    double *Sdim;
    double Yeff;
    double Geff;
    double roughnessRatio;


    bool tangential_damping;
    bool limitForce;
    bool displayedSettings;
    bool heating;
    bool heating_track;
    int elastic_potential_offset_;
    bool elasticpotflag_;
    FixPropertyAtom *fix_dissipated_;
    bool dissipatedflag_;
    int overlap_offset_;
    bool disable_when_bonded_;
    int bond_history_offset_;
    int dissipation_history_offset_;
    class ContactModelBase *cmb;
    
    class FixPropertyAtom* fix_Sdim; 
    
  };

}

}
#endif
#endif
