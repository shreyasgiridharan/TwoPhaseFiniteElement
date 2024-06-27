# Two Phase Finite element code
A simple FE code following the theory of porous media. 

The main subroutines are written in fem_2D.f90. Use that as a starting point. 

An example input file is provided in input.dat. 

The original program accepts input file only with the name 'input.dat'. You can however change that to whatever you want in fem_2D.f90.

Don't use the Neutral loading surface model for the moment, as it has not been thoroughly tested. The model is based off of Prof. Dieter Stolle's idea (stolle@mcmaster.ca). 

You will need GiD to pre- and post process the simulation.

If you encounted any errors, bugs or improvements, please contact shreyas.giridharan@gmail.com.
