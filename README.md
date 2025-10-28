# Exponential-asymptotics
General code to run the computation to study homoclinic snaking in general reaction-diffusion equations.

Two steps are necessary to obtain results from this code. First, one must run the Python code 'Beyond-all-order asymptotics.py' inside the folder 'demos'. To run it, one must create a folder 'foo' with a Python script 'foo.py' inside following the same format as in the other demos in that folder. That code will compute all the coefficients of the parametric expansion, determine whether the expansion leads to the existence of a Maxwell point, and let the user know whether the expansion can be carried out with the parameters given by default.

After this, one can go to the folder called 'Iteration' to set the Mathematica code up to compute an estimate of $c_r^{[0]}$. The only things necessary for a user to set said code up are at the beginning, in the first four blocks of code. They are related to the system one is analyzing. Some Mathematica codes are ready to be run in that folder, and I recommend the user run them (with a lower value of 'nmax', to ensure the code does not take long) in order to understand how the iteration works.

The code works in a simple way but it is long because of the number of steps that need to be followed. Specifically, as solving linear systems of equations might be cumbersome if the coefficients are algebraically challenging, some extra classes have been created to generate matrices and vectors with corresponding dummy variables that let Python solve all the linear systems of equations easily. Check the documentation within the file 'functions.py'.

Feel free to contact me if there are any problems with setting anything up: edgardo.villar-sepulveda@bristol.ac.uk
