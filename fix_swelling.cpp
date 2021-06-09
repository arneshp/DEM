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

#include "fix_swelling.h"
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

// caution: This fix should only be defined after fix definition of gran conduction. 

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

Fix_swelling::Fix_swelling(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg){

  if ((!atom->radius_flag)||(!atom->rmass_flag)) error->all(FLERR,"Fix heat/gran needs per particle radius and mass");

  if (narg < 11)
    error->fix_error(FLERR,this,"not enough arguments");

  int iarg = 3; 
  
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
  fix_temp = NULL;
  fix_Sdim = NULL;
  fix_Swratio = NULL;
  fix_iradi = NULL;
  fix_densityp = NULL;
  fix_densityf = NULL;
  
  peratom_flag = 1;      
  size_peratom_cols = 0; 
  peratom_freq = 1;

  scalar_flag = 1; 
  global_freq = 1; 

  rad_mass_vary_flag = 1;
}

/* ---------------------------------------------------------------------- */

void Fix_swelling::post_create()
{
  // register directional flux
 

 if(!fix_tdelayelapsed)
  {
    const char* fixarg[9];
    fixarg[0]="tdelayelapsed";   // This part is similar to defining a fix statement from the input script. so this fix defines the variable tdelayelapsed by default. 
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="tdelayelapsed";
    fixarg[4]="scalar";
    fixarg[5]="yes";
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
    fixarg[5]="yes";
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
    fixarg[5]="yes";
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
    fixarg[5]="yes";
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
    fixarg[5]="yes";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0001";
    
    fix_iradi = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);  // This is the line that is adding the new atom/property.  
  }

  
  fix_temp = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",0,0,style));   // static_cast is to convert the output into a <FixPropertyAtom*>, all these properties are defined in the fix by default. 

  fix_tdelayelapsed = static_cast<FixPropertyAtom*>(modify->find_fix_property("tdelayelapsed","property/atom", "scalar",0,0,style));	
  fix_alpha = static_cast<FixPropertyAtom*>(modify->find_fix_property("alpha", "property/atom", "scalar",0,0,style));	
  fix_Sdim = static_cast<FixPropertyAtom*>(modify->find_fix_property("Sdim", "property/atom", "scalar",0,0,style));
  fix_Swratio = static_cast<FixPropertyAtom*>(modify->find_fix_property("Swratio", "property/atom", "scalar",0,0,style));
  fix_iradi = static_cast<FixPropertyAtom*>(modify->find_fix_property("iradi", "property/atom", "scalar",0,0,style));
  fix_densityp = static_cast<FixPropertyGlobal*>(modify->find_fix_property("densityp", "property/global", "scalar",0,0,style));
  fix_densityf = static_cast<FixPropertyGlobal*>(modify->find_fix_property("densityf", "property/global", "scalar",0,0,style));

  if(!fix_densityp || !fix_densityf)
  error->one(FLERR, "fix swelling requires global properties density of particle 'densityp' and density of fluid 'densityf' to be defined");  

  if(!fix_temp || !fix_tdelayelapsed || !fix_alpha || !fix_Sdim || !fix_Swratio || !fix_iradi)
    error->one(FLERR,"internal error");
}

/* ---------------------------------------------------------------------- */

void Fix_swelling::updatePtrs()
{
  Temp = fix_temp->vector_atom; // Temp is a temperature of each atom, it is a vector. 
  vector_atom = Temp; 

  alpha = fix_alpha->vector_atom;
  tdelayelapsed = fix_tdelayelapsed->vector_atom;
  Sdim = fix_Sdim->vector_atom;
  Swratio = fix_Swratio->vector_atom;
  iradi = fix_iradi->vector_atom;
  
}

/* ---------------------------------------------------------------------- */

void Fix_swelling::init()
{
  
  if (!atom->radius_flag || !atom->rmass_flag)
    error->fix_error(FLERR,this,"must use a granular atom style ");

    // check if a fix of this style already exists
  //if(modify->n_fixes_style(style) > 1) 			// style is nothing but heat/gran for this fix. other fixes will have 
  // error->fix_error(FLERR,this,"cannot have more than one fix of this style"); // similarly different style. 

  if(!force->pair_match("gran", 0))
    error->fix_error(FLERR,this,"needs a granular pair style to be used");

  pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
  history_flag = pair_gran->is_history();

 
  fix_temp = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",0,0,style));
  
  
  fix_tdelayelapsed = static_cast<FixPropertyAtom*>(modify->find_fix_property("tdelayelapsed","property/atom", "scalar",0,0,style));	
  fix_alpha = static_cast<FixPropertyAtom*>(modify->find_fix_property("alpha", "property/atom", "scalar",0,0,style));
  fix_Sdim = static_cast<FixPropertyAtom*>(modify->find_fix_property("Sdim", "property/atom", "scalar",0,0,style));
  fix_Swratio = static_cast<FixPropertyAtom*>(modify->find_fix_property("Swratio", "property/atom", "scalar",0,0,style));	
  fix_iradi = static_cast<FixPropertyAtom*>(modify->find_fix_property("iradi", "property/atom", "scalar",0,0,style));

  if(!fix_temp || !fix_tdelayelapsed || !fix_alpha)
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

int Fix_swelling::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */



/* ---------------------------------------------------------------------- */



void Fix_swelling::post_integrate() 
{
  updatePtrs();
 double *radius = atom->radius;
 double tdelay_ref;
 double dt = update->dt;
 int nlocal = atom->nlocal;
 
 double excess_vol, mass_old;
 
 double densityp, densityf;

 densityp = fix_densityp->get_values()[0];
 densityf = fix_densityf->get_values()[0]; 
 
 
   for (int i=0; i<nlocal; i++)
   {
       tdelay_ref = -log(1-alpha[i])*tau_ref;
   	if (tdelay_ref<tdelayelapsed[i])
   	{
   	  Sdim[i] = Sdim[i]+ K0*pow((1-alpha[i]),0.5)*pow((Temp[i]-Tmin)/(Tref-Tmin),2)*(1-Sdim[i])*dt ;
   	  radius[i] = iradi[i]*(Swratio[i]-1)*Sdim[i] + iradi[i];
   	  excess_vol = 4.0*3.14159/3.0 * atom->radius[i]*atom->radius[i]*atom->radius[i] - 4.0*3.14159/3.0*iradi[i]*iradi[i]*iradi[i];
   	  
   	  mass_old = atom->rmass[i];
   	  
   	  atom->rmass[i] = 4.0*3.14159/3.0 * iradi[i]*iradi[i]*iradi[i]*densityp + excess_vol*densityf;
   	  atom->density[i] = atom->rmass[i]/(4.0*3.14159/3.0 * atom->radius[i]*atom->radius[i]*atom->radius[i]);   	   
   	    	
   	  atom->v[i][0] = atom->v[i][0]*sqrt(mass_old/(atom->rmass[i]));
   	  atom->v[i][1] = atom->v[i][1]*sqrt(mass_old/(atom->rmass[i]));
   	  atom->v[i][2] = atom->v[i][2]*sqrt(mass_old/(atom->rmass[i])); 
   	    	
   	    	}
   	else 
   	{ if (Temp[i]>Tmin) 
   	   tdelayelapsed[i] = tdelayelapsed[i] + dt*(Temp[i]-Tmin)/(Tref-Tmin);
   	}
   }
   updatePtrs();
}


/* ---------------------------------------------------------------------- */



