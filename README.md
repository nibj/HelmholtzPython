# HelmholtzPython
## Netgen code for 2D Helmholtz equation
Note that this code is supplied without any guarantees.  You cannot assume it is correct, reliable or fit for purpose.  You
are responsible!

To use this code
1. Install netgen from https://ngsolve.org (I prefer to use pip)
2. Run the code using
  *python3 helmholtz_dirichlet.py* (thats what works on my mac!)
   or
   *jupyter notebook helmholtz_penetrable.ipynb*

Note that 
1. geometries_dirichlet.py contains several geometries for simple impenentrable scatterers
2. helmholtz_dirichlet.py contains the ngspy code for the Dirichlet problem
3. geometries_penetrable.py contains several geometries for simple penetrable scatterers
4. helmholtz_penetrable.ipynb contains the ngspy code for the penetrable medium problem with some comments
