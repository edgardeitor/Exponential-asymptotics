This folder contains the Python code to compute the amplitude equations at orders 5, 6, and 7, for a general n-component system centred at a Maxwell point, and finds the asymptotic expansion of this point up to order 7. The main Python code is called 'Beyond-all-order asymptotics.py', and it requires to be in the same folder with the file 'functions.py' and, to run the code to compute the amplitude equations for an example foo, you need a file foo.py within the folder foo. This file is needed to define the system to be studied. The structure of the code is always the same, and each of the demos have the same variables defined within each file. Specifically, the file needs to contain the definitions of

-a dictionary, parameters, that specifies all the parameters of the system that are not being expanded, together with their values;
-a list, unevaluatedparameters, that specifies all the parameters that are required not to be evaluated by Python (the recommendation is to make this list empty. Otherwise, the code will take longer to run);
-the strings, par1 and par2, which are the names of the parameters $a$ and $b$ to be expanded;
-a string, muval, which is the value of the square of the wavenumber, k^2;
-a list of strings, equilibrium, which contains the algebraic coordinates of the steady state of the system;
-a list of strings, var, which is the vector of coordinates of the system;
-a list of lists, diffmatrix, which encodes the diffusion matrix of the system;
-a vector of strings, kinetics, which corresponds to the reaction of the system;
-a real number, tol, which corresponds to the threshold that is accepted to check for Turing bifurcation and codimension-two conditions;
-a string, phiunit, which is a boolean variable to state whether the vector $\phi_1^{[1]}$ is required to be a unit vector;
-a boolean variable, simp, which states whether the algebraic expressions are to be simplified or not;
-a boolean variable, simple, which states whether the equations for $\alpha_{i, 2}$ are solved explicitly or not (this variable is relevant only when the coefficients $a_1, b_1, a_3, b_3, a_5, b_5$ have not been set to be zero from the start);
-and a function, extrapars(extraparvals), that edits the dictionary that evaluates the parameters of the expansion following, for example, the structure extraparvals[aNF[1]] = 0, to say that a_1 = 0, for your parameter a.

To run the code, you need to have the following packages installed

IPython
os
shutil
sys
sympy
mpmath
math
