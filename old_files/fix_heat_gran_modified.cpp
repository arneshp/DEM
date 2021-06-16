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

#include "fix_heat_gran_modified.h"
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

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixHeatGran_modified::FixHeatGran_modified(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg){

  if ((!atom->radius_flag)||(!atom->rmass_flag)) error->all(FLERR,"Fix heat/gran needs per particle radius and mass");

  if (narg < 13)
    error->fix_error(FLERR,this,"not enough arguments");

  int iarg = 3;

  if(strcmp(arg[iarg++],"initial_temperature"))
    error->fix_error(FLERR,this,"expecting keyword 'initial_temperature'");
  T0 = atof(arg[iarg++]);
  
  if(strcmp(arg[iarg++],"K0"))
    error->fix_error(FLERR,this,"expecting keyword 'K0'");
  K0 = atof(arg[iarg++]);
  
  if(strcmp(arg[iarg++],"tau_ref"))
    error->fix_error(FLERR,this,"expecting keyword 'tau_ref'");
  tau_ref = atof(arg[iarg++]);

  if(strcmp(arg[iarg++],"Min_temperature"))
    error->fix_error(FLERR,this,"expecting keyword 'Min_temperature'");
  Tmin = atof(arg[iarg++]);

  if(strcmp(arg[iarg++],"Reference_temperature"))
    error->fix_error(FLERR,this,"expecting keyword 'Reference_temperature'");
  Tref = atof(arg[iarg++]);

  fix_tdelayelapsed = NULL;
  fix_alpha = NULL;
  fix_temp = fix_heatFlux = fix_heatSource = NULL;
  fix_ste = NULL;
  fix_directionalHeatFlux = NULL;
  peratom_flag = 1;      
  size_peratom_cols = 0; 
  peratom_freq = 1;

  scalar_flag = 1; 
  global_freq = 1; 

  cpl = NULL;
}

/* ---------------------------------------------------------------------- */

void FixHeatGran_modified::post_create()
{
  // register directional flux
  fix_directionalHeatFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("directionalHeatFlux","property/atom","vector",3,0,this->style,false));
  if(!fix_directionalHeatFlux)
  {
    const char* fixarg[11];
    fixarg[0]="directionalHeatFlux";   // This part is similar to defining a fix statement from the input script. so this fix defines the variable directional heat flux by default. 
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="directionalHeatFlux";
    fixarg[4]="vector";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    fixarg[9]="0.";
    fixarg[10]="0.";
    fix_directionalHeatFlux = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);  // This is the line that is adding the new atom/property.  
  }

 if(!fix_tdelayelapsed)
  {
    const char* fixarg[9];
    fixarg[0]="tdelayelapsed";   // This part is similar to defining a fix statement from the input script. so this fix defines the variable tdelayelapsed by default. 
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="tdelayelapsed";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    
    fix_tdelayelapsed = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);  // This is the line that is adding the new atom/property.  
  }

 if(!fix_alpha)
  {
    const char* fixarg[9];
    fixarg[0]="alpha";   // This part is similar to defining a fix statement from the input script. so this fix defines the variable alpha-difficulty of swelling parameter by default. 
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="alpha";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    
    fix_alpha = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);  // This is the line that is adding the new atom/property.  
  }

 if(!fix_Sdim)
  {
    const char* fixarg[9];
    fixarg[0]="Sdim";   // This part is similar to defining a fix statement from the input script. so this fix defines the variable Sdim-Nondimensional diameter  by default. 
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="Sdim";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
   
    fix_Sdim = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);  // This is the line that is adding the new atom/property.  
  }


 if(!fix_Swratio)
  {
    const char* fixarg[9];
    fixarg[0]="Swratio";   // This part is similar to defining a fix statement from the input script. so this fix defines the variable Swratio-Swelling ratio by default. 
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="Swratio";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    
    fix_Swratio= modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);  // This is the line that is adding the new atom/property.  
  }

 if(!fix_iradi)
  {
    const char* fixarg[9];
    fixarg[0]="iradi";   // This part is similar to defining a fix statement from the input script. so this fix defines the variable iradi-initial radius by default. 
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="iradi";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0001";
    
    fix_iradi = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);  // This is the line that is adding the new atom/property.  
  }

  fix_ste = modify->find_fix_scalar_transport_equation("heattransfer");  // This line is checking if there already exists a heattransfer property, this apparently is defined in the scalar_transport_equation "sub-"fix. 
  if(!fix_ste)
  {
    const char * newarg[15];
    newarg[0] = "ste_heattransfer";
    newarg[1] = group->names[igroup];
    newarg[2] = "transportequation/scalar";
    newarg[3] = "equation_id";
    newarg[4] = "heattransfer";
    newarg[5] = "quantity";
    newarg[6] = "Temp";
    newarg[7] = "default_value";
    char arg8[30];
    sprintf(arg8,"%f",T0);     // sprintf basically pushes the value of T0 into arg8. T0 is the inital temperature of the particles, (all same temperature). 
    newarg[8] = arg8;
    newarg[9] = "flux_quantity";
    newarg[10] = "heatFlux";
    newarg[11] = "source_quantity";
    newarg[12] = "heatSource";
    newarg[13] = "capacity_quantity";
    newarg[14] = "thermalCapacity";
    modify->add_fix(15,(char**)newarg);   // adding a new fix variable to the modify class. in the add_fix function has provisions to add all the fixes defined in the input script.  
  }

  fix_temp = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",0,0,style));   // static_cast is to convert the output into a <FixPropertyAtom*>, all these properties are defined in the fix by default. 
  fix_heatFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFlux","property/atom","scalar",0,0,style));
  fix_heatSource = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatSource","property/atom","scalar",0,0,style));
  fix_directionalHeatFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("directionalHeatFlux","property/atom","vector",0,0,style));

  fix_tdelayelapsed = static_cast<FixPropertyAtom*>(modify->find_fix_property("tdelayelapsed","property/atom", "scalar",0,0,style));	
  fix_alpha = static_cast<FixPropertyAtom*>(modify->find_fix_property("alpha", "property/atom", "scalar",0,0,style));	
  fix_Sdim = static_cast<FixPropertyAtom*>(modify->find_fix_property("Sdim", "property/atom", "scalar",0,0,style));
  fix_Swratio = static_cast<FixPropertyAtom*>(modify->find_fix_property("Swratio", "property/atom", "scalar",0,0,style));
  fix_iradi = static_cast<FixPropertyAtom*>(modify->find_fix_property("iradi", "property/atom", "scalar",0,0,style));

  if(!fix_temp || !fix_heatFlux || !fix_heatSource || !fix_directionalHeatFlux || !fix_tdelayelapsed || !fix_alpha || !fix_Sdim || !fix_Swratio || !fix_iradi)
    error->one(FLERR,"internal error");
}

/* ---------------------------------------------------------------------- */

void FixHeatGran_modified::updatePtrs()
{
  Temp = fix_temp->vector_atom;
  vector_atom = Temp; 

  heatFlux = fix_heatFlux->vector_atom;              // each of these newly created fixes are classes and have peratom variables as vector if the quantity is scalar or array if they are vector quantities, (x,y,z)
  heatSource = fix_heatSource->vector_atom;		// 	Looks like fix_heatFlux , fix_heatSource, fix_directionalHeatFlux are all "sub-" fixes of fix_temp. 
  directionalHeatFlux = fix_directionalHeatFlux->array_atom;
  alpha = fix_alpha->vector_atom;
  tdelayelapsed = fix_tdelayelapsed->vector_atom;
  Sdim = fix_Sdim->vector_atom;
  Swratio = fix_Swratio->vector_atom;
  iradi = fix_iradi->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixHeatGran_modified::init()
{
  
  if (!atom->radius_flag || !atom->rmass_flag)
    error->fix_error(FLERR,this,"must use a granular atom style ");

    // check if a fix of this style already exists
  if(modify->n_fixes_style(style) > 1) 			// style is nothing but heat/gran for this fix. other fixes will have 
    error->fix_error(FLERR,this,"cannot have more than one fix of this style"); // similarly different style. 

  if(!force->pair_match("gran", 0))
    error->fix_error(FLERR,this,"needs a granular pair style to be used");

  pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
  history_flag = pair_gran->is_history();

  fix_ste = modify->find_fix_scalar_transport_equation("heattransfer");    // checks if the fix_ste with "heattransfer" is defined or not.  The governing equations for the heat transfer is defined in that sub-fix. 
  if(!fix_ste) error->fix_error(FLERR,this,"needs a fix transportequation/scalar to work with");

  fix_temp = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",0,0,style));
  fix_heatFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFlux","property/atom","scalar",0,0,style));
  fix_heatSource = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatSource","property/atom","scalar",0,0,style));
  fix_directionalHeatFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("directionalHeatFlux","property/atom","vector",0,0,style));
  
  fix_tdelayelapsed = static_cast<FixPropertyAtom*>(modify->find_fix_property("tdelayelapsed","property/atom", "scalar",0,0,style));	
  fix_alpha = static_cast<FixPropertyAtom*>(modify->find_fix_property("alpha", "property/atom", "scalar",0,0,style));
  fix_Sdim = static_cast<FixPropertyAtom*>(modify->find_fix_property("Sdim", "property/atom", "scalar",0,0,style));
  fix_Swratio = static_cast<FixPropertyAtom*>(modify->find_fix_property("Swratio", "property/atom", "scalar",0,0,style));	
  fix_iradi = static_cast<FixPropertyAtom*>(modify->find_fix_property("iradi", "property/atom", "scalar",0,0,style));

  if(!fix_temp || !fix_heatFlux || !fix_heatSource || !fix_directionalHeatFlux || !fix_tdelayelapsed || !fix_alpha)
    error->one(FLERR,"internal error");

  updatePtrs();
  
   int nlocal = atom->nlocal;   // This is all the atoms on the current processor. 
   double *radius = atom->radius;
   std::default_random_engine generator;
  std::normal_distribution<double> distribution(2.34,0.33);
  std::uniform_real_distribution<> dis(0.00001, 1.0);
   
   for (int i = 0; i < nlocal; i++) 
   {
     tdelayelapsed[i] = 0.; // setting the elapsed delaytime in the reference time frame for all the granules to zero.
     Sdim[i] = 0.; // setting intial non-dimensional diameters to 0.  
     *iradi++ = *radius++;
     
     Swratio[i] = distribution(generator);
     alpha[i] = dis(generator);   
   }
}

/* ---------------------------------------------------------------------- */

int FixHeatGran_modified::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHeatGran_modified::initial_integrate(int vflag)
{
  updatePtrs();
   //reset heat flux
  //sources are not reset
  //int *mask = atom->mask;
  int nlocal = atom->nlocal;   // This is all the atoms on the current processor. 
  double *radius = atom->radius;
  

  for (int i = 0; i < nlocal; i++) // This loop is setting the initial fluxes as zeros. to start integrating in the further timesteps. 
  {
     //if (mask[i] & groupbit)
     {
        directionalHeatFlux[i][0] = 0.;
        directionalHeatFlux[i][1] = 0.;
        directionalHeatFlux[i][2] = 0.;
     }
  }
 
  //update ghosts
  fix_directionalHeatFlux->do_forward_comm();
 }

/* ---------------------------------------------------------------------- */

double FixHeatGran_modified::compute_scalar()
{
    return fix_ste->compute_scalar(); // The fix_ste sub-fix is returning scalar. i dont know what it is. This is the last impelemented method in this fix, so the action goes to the sub-fix. 
}


void FixHeatGran_modified::pre_force(int vflag) 
{
  updatePtrs();
 double *radius = atom->radius;
 double tdelay_ref;
 double dt = update->dt;
 int nlocal = atom->nlocal;
   for (int i=0; i<nlocal; i++)
   {
       tdelay_ref = -log(1-alpha[i])*tau_ref;
   	if (tdelay_ref<tdelayelapsed[i])
   	{
   	  Sdim[i] = Sdim[i]+ K0*pow((1-alpha[i]),0.5)*pow((Temp[i]-Tmin)/(Tref-Tmin),2)*(1-Sdim[i])*dt ;
   	  radius[i] = iradi[i]*(Swratio[i]-1)*Sdim[i] + iradi[i];
   	    	}
   	else 
   	{ if (Temp[i]>Tmin) 
   	   tdelayelapsed[i] = tdelayelapsed[i] + dt*(Temp[i]-Tmin)/(Tref-Tmin);
   	}
   }
   updatePtrs();
}


/* ---------------------------------------------------------------------- */

void FixHeatGran_modified::cpl_evaluate(class ComputePairGranLocal * cpl){

  char *mystyle = style;
  char *emsg = new char[100];
  sprintf(emsg, "Fix %s does not implement cpl_evaluate().\n", mystyle);  // I dont understand this method. may this is not relevant for this fix and the error message says it is not implemented for this fix.
  error->all(FLERR, emsg);

}

/* ---------------------------------------------------------------------- */

void FixHeatGran_modified::register_compute_pair_local(class ComputePairGranLocal *ptr){

  char *mystyle = style;
  char *emsg = new char[100];
  sprintf(emsg, "Fix %s does not implement register_compute_pair_local().\n", mystyle);
  error->all(FLERR, emsg);

}

/* ---------------------------------------------------------------------- */

void FixHeatGran_modified::unregister_compute_pair_local(class ComputePairGranLocal *ptr){

  char *mystyle = style;
  char *emsg = new char[100];
  sprintf(emsg, "Fix %s does not implement unregister_compute_pair_local().\n", mystyle);
  error->all(FLERR, emsg);

}
