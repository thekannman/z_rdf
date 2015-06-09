# z_rdf
This program calculates radial distribution functions (including various weighted variations) from molecular dynamics simulations.

THIS PROGRAM HAS NOT BEEN THOROUGHLY TESTED!!!

Current limitations include:
* Only trajectories from the GROMACS simulation package can be used as input.

The following libraries are required:
* GROMACS XTC library for reading position/velocity files. 
  http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library
* The Armadillo library for matrix calculations. http://arma.sourceforge.net/
* The Boost libraries for input option parsing, 4D tensor calculations, lexical casting, and output formatting. 
  http://www.boost.org/
