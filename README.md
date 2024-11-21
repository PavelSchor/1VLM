Introduction

1VLM is an implementation of Vortex Lattice Method as described in [1] 
It is fully implemented in Python, without any external FORTRAN or C library.  It relies on NUMBA to provide reasonable performance.
![Alt text](doc/intro.png?raw=true "Example aircraft")

Features
- Steady and unsteady aerodynamic analyses of complex configurations
- Groups of panels (for loads analyses)
- Deflections of control surfaces
- Basic trim algorithm, see example4 (to be coded in FDM later)
- Simple flight dynamic solver with euler integration method [2]
- Exports to VTK and JSON

TODO
- implement proper wake relaxation and variable stregth of wake vortex rings
- investigate / implement Runge Kutta integration scheme in FDM

References

  [1]: J. Katz and A. Plotkin, “Low Speed Aerodynamics,” 2nd Edition, Cambridge University Press, Cambridge, 2001.
  
  [2]: STEVENS, Brian L. and Frank L. LEWIS. c2003. Aircraft control and simulation. 2nd ed. Hoboken: Wiley. ISBN 0-471-37145-9.
