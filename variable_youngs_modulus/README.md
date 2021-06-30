noramal_model_hertz_swelling impelments the variation of young's modulus of the particles based on the swelling degree. 

Two property creators are added to the global_properties.cpp 
VectorProperty* createYoungsModuluslow and VectorProperty* createYoungsModulushigh 

This the lower and upper bound of the young's modulus. 

The young's modulus is determined based on the swelling degree of the particle.

When the "swelling degree" Sdim is zero. (no swelling), the young's modulus is YoungsModulushigh 
and when the "swelling degree" Sdim is 1 (fully swollen), the young's modulus is YoungsModuluslow.

For the intermediate stages, 
Y = Yh*(1-Sdim) + Yl*(Sdim). 


This model requries the defintion of these global pertype properties along with other pre-requisties of hertz contact model. 

YoungsModulushigh
YoungsModuluslow 

Also the peratom quantity degree of swelling Sdim needs to be defined. 
This can be done iin two ways either by defining the swelling model(ideal way) with fix swelling which in-turn requires heat/gran fix defintion or 
by simply setting the peratom Sdim to a fixed quantity like in the test cases. 

 fix Sdim all property/atom Sdim scalar yes yes no 0.
 
 
The test case shows a pair of particles coliding, based the value of Sdim the Young's modulus is changed and the forces change leading to different trajectories. 
This is tested by taking the edge case scenarios and setting the young's modulus manually in the LIGGGTHS Hertz impelmentation at the hard and soft limits. 
The responses are tested with the setting the Sdim to 0 (hard limit) and Sdim to 1 (soft limit). 

checkout the input scripts.

 
 


