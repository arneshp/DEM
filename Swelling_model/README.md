Swelling model as described in the paper Palanisamy et al. 2020 is implemented here, https://doi.org/10.1016/j.foostr.2020.100150 .

It is impelmented as fix. 
Swelling fix requires one to define the heat/gran fix to be defined earlier. 
Refer the documentation here to read more about the same https://www.cfdem.com/media/DEM/docu/fix_heat_gran_conduction.html#fix-heat-gran-conduction-command. 

This is because the temperauture peratom property is defined and derived from this fix (heat/gran).

The typical implemenation of this fix is as follows. 

fix         swelling all swelling K0 0.0023 tau_ref 190. Min_temperature 333. Reference_temperature 335.  

The test cases include the stationary particles with constant temperature of particles case. 
The impelmetation of the particle growth is compared with the matlab script for this simple case. 

Second test case is the particles placed linearly and have motion in a single direction. 

Third test case is particles with linear velocity with a small offset in the second dimension. temperature zoning is also implemented  using "set             region halfbed property/atom Temp 353 ". 

final example is for test the parallel and serial runs.
This example has init and Run files. 
 
The stochastic paramters of the model are reset when the restart file is re-read. 
Thus to test this example manually, set the values of all the atoms to alpha = 0.5 and Swratio = 2.34 in the fix_swelling.cpp 
(lines 272 and 273) and recompile the code. 

 Now run in.init file to generate restart file. (This contains locations of particles). 

Then run the in.runparallel and in.runserial with the required number of processors. 
compare stats such as kinetic energy vs time. 



