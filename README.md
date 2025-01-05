# Exponential-asymptotics
General code to run the computation to study homoclinic snaking in general reaction-diffusion equations.

There are two parts to run this code. First, one must run the Python code 'Beyond-all-order asymptotics.py' that is inside the folder 'demos'. To run it, one must create a folder 'foo' with a Python script 'foo.py' inside following the same format as in the other demos. That code will compute all the coefficients needed for the expansion and determine whether the expansion or some parameter values provided are right or wrong.

After this, one can go to the folder called 'Iteration' to set the Mathematica code up to compute an estimate of $c_2^{[0]}$. The only things that are necessary for a user to change are at the beginning. They are related to the system one is analyzing. The only extra thing that is needed from the user is the estimation of $A_3$, which is located in the cell called 'Iterative solver'. Some codes ready to be run can be found in the corresponding folder and I recommend the user to run them if any questions. Also, feel free to get in touch with me if there are any problems with setting everything up.
