# DEM
LIGGGHTS files used for DEM simulation.

Lubrication implementation in LIGGGHTS Based on Ball and melrose 1997. 
implementation similar to https://lammps.sandia.gov/doc/pair_lubricate.html 

fix heat/gran/modified 
is the swelling model of starch granules based on the paper 
https://www.sciencedirect.com/science/article/pii/S2213329120300150 

The parameter difficulty of swelling "alpha", Non-dimensional diameter "Sdim" are the parameters that can be output using the
dump command in the input script. 

The input parameters of the model are particle "initial_temperautre", "K0" (rate constant of the kinetic model), "tau_ref" (time constant for 
the delay in swelling of granules), "Min_temperature" (minimum temperature for swelling), "Reference_temperature" (time constant reference temperature). 

The particle radius is increased based on a parameter, ratio of swelling which is from a normal distribution from N(2.34,0.33) 
This is internally hard coded in the fix. 
