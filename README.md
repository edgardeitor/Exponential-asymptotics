# Exponential-asymptotics
General code to run the computation to study homoclinic snaking in general reaction-diffusion equations.

Two steps are necessary to obtain results from this code. First, one must run the Python code 'Beyond-all-order asymptotics.py' inside the folder 'demos'. To run it, one must create a folder 'foo' with a Python script 'foo.py' inside following the same format as in the other demos in that folder. That code will compute all the coefficients of the parametric expansion, determine whether the expansion leads to the existence of a Maxwell point, and let the user know whether the expansion can be carried out with the parameters given by default.

After this, one can go to the folder called 'Iteration' to set the Mathematica code up to compute an estimate of $c_r^{[0]}$ (see Equation (46) in arXiv preprint arXiv:2501.02698). The only things necessary for a user to set that code up are on the second and third cells. They are related to the system one is analyzing. Some Mathematica codes are ready to be run in that folder, and I recommend the user run them (with a lower value of 'nmax', to ensure the code does not take long) in order to understand how the iteration works.

The code works in a simple way, but it is long because of the number of steps that need to be followed. See details in Appendix D of arXiv preprint arXiv:2501.02698, or the documentation within the file functions.py inside the folder demos.

Feel free to contact me if there are any problems with setting anything up: edgardo.villar-sepulveda@bristol.ac.uk
