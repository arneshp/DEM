#### DEM
LIGGGHTS files used for DEM simulation.

Lubrication implementation in LIGGGHTS Based on Ball and melrose 1997. 
implementation similar to https://lammps.sandia.gov/doc/pair_lubricate.html 
For NEMD simulations use fix_LubricateFLD.cpp and for simulation without moving walls and non-deforming boundaries such as for CFDEM coupling one can use both fix_LubricationFLD.cpp as well as cohension_Lubrication.h 


fix swelling 
is the swelling model of starch granules based on the paper 
https://www.sciencedirect.com/science/article/pii/S2213329120300150 

This fix (swelling) needs fix heat/gran to be defined ealier, as temperature variable is calculated and derived from the same.   

The parameter difficulty of swelling "alpha", Non-dimensional diameter "Sdim" are the parameters that can be output using the
dump command in the input script. 

The input parameters of the model are particle "initial_temperautre", "K0" (rate constant of the kinetic model), "tau_ref" (time constant for 
the delay in swelling of granules), "Min_temperature" (minimum temperature for swelling), "Reference_temperature" (time constant reference temperature). 

The particle radius is increased based on a parameter, ratio of swelling which is from a normal distribution from N(2.34,0.33) 
This is internally hard coded in the fix. 

Please note that the stochastic parameters and Sdim are reset and reinitialised everytime the fix swelling is redefined irrespective of the values in the restart file. Thus for validation runs for such as parallel vs serial. Please set them to a fixed value (say mean value, ~0.5 for alpha, etc) and validate.


fix_LubricationFLD.cpp can be used to implement the lubrication contact model popularly known as Ball and melrose. 
This file is similar to pair_lubricate_poly.cpp in lammps. 

The following functions v_tally_tensor() and v_tally_xyz() are needed to compute the virial stress tensors for NEMD simulations. 

added v_tally_tensor() and v_tally_xyz() 
These are needed in the fix_LubricationFLD.cpp 

The stress contribution due to lubrication forces is calculated here, which can independently be accessed in the compute pressure command with keyword fix.
compute pressure with keyword pair gives the contributions due to the contact models (say Hertz).


######################################################################


work-in progress to implement variable friction contact model proposed by lorby et al 2019. 

