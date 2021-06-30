Lubrication contact model is implemented as described in https://docs.lammps.org/pair_lubricate.html
Ball and melrose approximation is used to approximate stokes forces to the near-field lubrication and Fast Lubrication dynamics.

fix implementation (i.e. fix_lubricationFLD.cpp) is accurate for less-edward boundary conditions simulations as virial tensors are not 
updated in the cohesion_model_LUBRICATION.h implementation for the fix deform case. 

However, the forces and torques are nearly identitical for both the implementation. The added advantage of fix impelmentation is the contribution of lubrication forces to the viscosity of the suspension can be extrated using the fix keyword in the compute pressure command. 

Please note the Fix::v_tally_tensor(..) is added in the fix.h and fix.cpp files. 
This is necessary for updating the stresslets due to the fast lubrication forces. 

The test cases are a pair of particles or 3 particles interacting with each other due to the lubrication forces when the box is deformed. 

The equivalence with lammps is shown. 


