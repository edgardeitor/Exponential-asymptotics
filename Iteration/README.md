This folder contains the Mathematica files to perform the late-order expansion to obtain the variable K_2. Specifically, each folder has just one file that needs to be run cell by cell in order to obtain an approximation of the coefficient $c_r^{[0]}$ explained in the article arXiv preprint arXiv:2501.02698 (Equation (46)).

To run the code, you have to set up the system in Mathematica on the second and third cells by following any of the examples. You need to run the Python code before any of these scripts because the values of $\xi$, $\eta$, $A_1$, and $A_3$ are needed to run these files.

Specifically, the code needs the following variables to run correctly:

-nmax, which represents the maximum order to be computed by the code;

-var, which represents an array with the variables of the system;

-par, which represents an array with the parameters of the system;

-nvar, which represents the number of variables of the system and is computed automatically by using the length of the array \texttt{var};

-fixedpar, which represents an array to evaluate the parameters of the system that are fixed;

-parval, which represents an array to evaluate the values of the parameters of the expansion at the codimension-two bifurcation point;

-$\eta$val, which represents the value of $\eta$ obtained from the Python code;

-$\xi$val, which represents the value of $\xi$ obtained from the Python code;

-$\kappa$val, which represents the value of $\kappa$ to be considered in the computation of the limit defining $c_r^{[0]}$;

-C11eval, which represents an array to evaluate the first-order expansion of $A_1$ at $X_0$;

-A3eval, which represents an array to evaluate the first-order expansion of $A_3$ at $X_0$.

-kinetics, which represents an array with the vector field, $\mbf f$, defined by the system \eqref{geneq};

-diffmatrix, which represents a two-dimensional array for the diffusion matrix, $\hat D$, defined by the system \eqref{geneq};

After setting the system up and choosing the maximum order you would like to expand your equation up to, you can just run the code and then plot the result.

The folders are named according to the models shown in the paper:

-'Swift-Hohenberg': Swift-Hohenberg 2-3.

-'SHDM': Swift-Hohenberg 3-5

-'Schnakenberg': Modified Schnakenberg system

-'Brusselator': Brusselator system

-'Bru 4': 4-component Brusselator
