class Vector:
    '''This class defines a vector and a dummy version of it to fasten code.'''
    def __init__(self, name):
        '''
        This function generates a vector as a 1D matrix, and defines its coordinates as symbolic dummy variables.

        :param name: This variable represents the name that will be used to call the coordinates of the vector,
        which follow the structure (name1, name2, ...).
        '''
        self.dummy = []
        self.actualcoord = []
        for varnum in range(nvar):
            self.dummy.append(symbols(name + '_' + str(varnum)))
            self.actualcoord.append(symbols(name + '_' + str(varnum)))
        self.dummy = Matrix(self.dummy)
            
class matrix:
    '''This class defines a sympy Matrix and a dummy version of it to fasten the code.'''
    def __init__(self, name, ncol = nvar, refmatrix = ones(nvar)):
        '''
        This function generates a matrix, and defines its coordinates as symbolic dummy variables.
        An extra argument can be passed to store the actual coordinates of the matrix that can be used for
        later substitution

        :param name: This variable represents the name that will be used to call the coordinates of the matrix,
        which follow the structure ((name11, name12, ...), (name21, name22, ...), ...).
        :param refmatrix: This variable represents the reference matrix to be defined as the actual coordinates
        of the matrix for later substitution if provided.
        '''
        self.dummy = zeros(ncol)
        if refmatrix!=ones(nvar):
            self.actualcoord = refmatrix
        for row in range(ncol):
            for col in range(ncol):
                if refmatrix[row,col]!=0:
                    self.dummy[row, col] = symbols(name + '_' + str(row+1) + '^' + str(col + 1))
                    
def first_order(ind_a, ind_b, u):
    '''
    This function computes the first order Taylor expansion of f_{ind_a, ind_b} at 0 in the direction of u.
    
    :param ind_a: This variable represents the index for the expansion of the vector field in terms of the parameter a
    :param ind_b: This variable represents the index for the expansion of the vector field in terms of the parameter b
    :param u: This variable is a vector that represents the direction in which the first order expansion is being applied

    :return: This function returns a vector that represents the first order Taylor expansion of f_{ind_a, ind_b} at 0 in the direction of u
    '''
    SS = zeros(nvar, 1)
    firstorderderivatives = globals()[f'firstorderderivatives{ind_a}{ind_b}']
    for i in range(nvar):
        SS = Add(SS, Mul(u.dummy[i], firstorderderivatives[i]))
    return SS
    
def second_order(ind_a, ind_b, u, v):
    '''
    This function computes the second order Taylor expansion of f_{ind_a, ind_b} at 0 in the directions of u, v
    
    :param ind_a: This variable represents the index for the expansion of the vector field in terms of the parameter a
    :param ind_b: This variable represents the index for the expansion of the vector field in terms of the parameter b
    :param u: This variable is a vector that represents the first direction in which the second order expansion is being applied
    :param v: This variable is a vector that represents the second direction in which the second order expansion is being applied

    :return: This function returns a vector that represents the second order Taylor expansion of f_{ind_a, ind_b} at 0 in the directions of u, v
    '''
    DS = zeros(nvar, 1)
    secondorderderivatives = globals()[f'secondorderderivatives{ind_a}{ind_b}']
    for i in range(nvar):
        for ii in range(nvar):
            DS = Add(DS, Mul(u.dummy[i], v.dummy[ii], secondorderderivatives[i][ii]))
    return Mul(DS, Pow(math.factorial(2), - 1))
    
def third_order(ind_a, ind_b, u, v, w):
    '''
    This function computes the third order Taylor expansion of f_{ind_a, ind_b} at 0 in the directions of u, v, w
    
    :param ind_a: This variable represents the index for the expansion of the vector field in terms of the parameter a
    :param ind_b: This variable represents the index for the expansion of the vector field in terms of the parameter b
    :param u: This variable is a vector that represents the first direction in which the third order expansion is being applied
    :param v: This variable is a vector that represents the second direction in which the third order expansion is being applied
    :param w: This variable is a vector that represents the third direction in which the third order expansion is being applied

    :return: This function returns a vector that represents the third order Taylor expansion of f_{ind_a, ind_b} at 0 in the directions of u, v, w
    '''
    TS = zeros(nvar, 1)
    thirdorderderivatives = globals()[f'thirdorderderivatives{ind_a}{ind_b}']
    for i in range(nvar):
        for ii in range(nvar):
            for iii in range(nvar):
                TS = Add(TS, Mul(u.dummy[i], v.dummy[ii], w.dummy[iii], thirdorderderivatives[i][ii][iii]))
    return Mul(TS, Pow(math.factorial(3), - 1))

def fourth_order(ind_a, ind_b, u, v, w, r):
    '''
    This function computes the fourth order Taylor expansion of f_{ind_a, ind_b} at 0 in the directions of u, v, w, r
    
    :param ind_a: This variable represents the index for the expansion of the vector field in terms of the parameter a
    :param ind_b: This variable represents the index for the expansion of the vector field in terms of the parameter b
    :param u: This variable is a vector that represents the first direction in which the fourth order expansion is being applied
    :param v: This variable is a vector that represents the second direction in which the fourth order expansion is being applied
    :param w: This variable is a vector that represents the third direction in which the fourth order expansion is being applied
    :param r: This variable is a vector that represents the fourth direction in which the fourth order expansion is being applied

    :return: This function returns a vector that represents the fourth order Taylor expansion of f_{ind_a, ind_b} at 0 in the directions of u, v, w, r
    '''
    Q4S = zeros(nvar, 1)
    fourthorderderivatives = globals()[f'fourthorderderivatives{ind_a}{ind_b}']
    for i in range(nvar):
        for ii in range(nvar):
            for iii in range(nvar):
                for iv in range(nvar):
                    Q4S = Add(Q4S, Mul(u.dummy[i], v.dummy[ii], w.dummy[iii], r.dummy[iv],
                                       fourthorderderivatives[i][ii][iii][iv]))
    return Mul(Q4S, Pow(math.factorial(4), - 1))

def fifth_order(ind_a, ind_b, u, v, w, r, s):
    '''
    This function computes the fifth order Taylor expansion of f_{ind_a, ind_b} at 0 in the directions of u, v, w, r, s
    
    :param ind_a: This variable represents the index for the expansion of the vector field in terms of the parameter a
    :param ind_b: This variable represents the index for the expansion of the vector field in terms of the parameter b
    :param u: This variable is a vector that represents the first direction in which the fifth order expansion is being applied
    :param v: This variable is a vector that represents the second direction in which the fifth order expansion is being applied
    :param w: This variable is a vector that represents the third direction in which the fifth order expansion is being applied
    :param r: This variable is a vector that represents the fourth direction in which the fifth order expansion is being applied
    :param s: This variable is a vector that represents the fifth direction in which the fifth order expansion is being applied

    :return: This function returns a vector that represents the fifth order Taylor expansion of f_{ind_a, ind_b} at 0 in the directions of u, v, w, r, s
    '''
    Q5S = zeros(nvar,1)
    fifthorderderivatives = globals()[f'fifthorderderivatives{ind_a}{ind_b}']
    for i in range(nvar):
        for ii in range(nvar):
            for iii in range(nvar):
                for iv in range(nvar):
                    for v5 in range(nvar):
                        Q5S = Add(Q5S, Mul(u.dummy[i], v.dummy[ii], w.dummy[iii], r.dummy[iv], s.dummy[v5],
                                           fifthorderderivatives[i][ii][iii][iv][v5]))
    return Mul(Q5S, Pow(math.factorial(5), - 1))

def sixth_order(ind_a, ind_b, u, v, w, r, s, t):
    '''
    This function computes the sixth order Taylor expansion of f_{ind_a, ind_b} at 0 in the directions of u, v, w, r, s, t
    
    :param ind_a: This variable represents the index for the expansion of the vector field in terms of the parameter a
    :param ind_b: This variable represents the index for the expansion of the vector field in terms of the parameter b
    :param u: This variable is a vector that represents the first direction in which the sixth order expansion is being applied
    :param v: This variable is a vector that represents the second direction in which the sixth order expansion is being applied
    :param w: This variable is a vector that represents the third direction in which the sixth order expansion is being applied
    :param r: This variable is a vector that represents the fourth direction in which the sixth order expansion is being applied
    :param s: This variable is a vector that represents the fifth direction in which the sixth order expansion is being applied
    :param t: This variable is a vector that represents the sixth direction in which the sixth order expansion is being applied

    :return: This function returns a vector that represents the sixth order Taylor expansion of f_{ind_a, ind_b} at 0 in the directions of u, v, w, r, s, t
    '''
    S6S = zeros(nvar, 1)
    sixthorderderivatives = globals()[f'sixthorderderivatives{ind_a}{ind_b}']
    for i in range(nvar):
        for ii in range(nvar):
            for iii in range(nvar):
                for iv in range(nvar):
                    for v5 in range(nvar):
                        for vi in range(nvar):
                            S6S = Add(S6S, Mul(u.dummy[i], v.dummy[ii], w.dummy[iii], r.dummy[iv], s.dummy[v5],
                                               t.dummy[vi], sixthorderderivatives[i][ii][iii][iv][v5][vi]))
    return Mul(S6S, Pow(math.factorial(6), - 1))

def seventh_order(ind_a, ind_b, u, v, w, r, s, t, x):
    '''
    This function computes the seventh order Taylor expansion of f_{ind_a, ind_b} at 0 in the directions of u, v, w, r, s, t, x
    
    :param ind_a: This variable represents the index for the expansion of the vector field in terms of the parameter a
    :param ind_b: This variable represents the index for the expansion of the vector field in terms of the parameter b
    :param u: This variable is a vector that represents the first direction in which the seventh order expansion is being applied
    :param v: This variable is a vector that represents the second direction in which the seventh order expansion is being applied
    :param w: This variable is a vector that represents the third direction in which the seventh order expansion is being applied
    :param r: This variable is a vector that represents the fourth direction in which the seventh order expansion is being applied
    :param s: This variable is a vector that represents the fifth direction in which the seventh order expansion is being applied
    :param t: This variable is a vector that represents the sixth direction in which the seventh order expansion is being applied
    :param x: This variable is a vector that represents the seventh direction in which the seventh order expansion is being applied

    :return: This function returns a vector that represents the seventh order Taylor expansion of f_{ind_a, ind_b} at 0 in the directions of u, v, w, r, s, t, x
    '''
    S7S = zeros(nvar, 1)
    seventhorderderivatives = globals()[f'seventhorderderivatives{ind_a}{ind_b}']
    for i in range(nvar):
        for ii in range(nvar):
            for iii in range(nvar):
                for iv in range(nvar):
                    for v5 in range(nvar):
                        for vi in range(nvar):
                            for vii in range(nvar):
                                S7S = Add(S7S, Mul(u.dummy[i], v.dummy[ii], w.dummy[iii], r.dummy[iv],
                                                   s.dummy[v5], t.dummy[vi], x.dummy[vii],
                                                   seventhorderderivatives[i][ii][iii][iv][v5][vi][vii]))
    return Mul(S7S, Pow(math.factorial(7), - 1))

def dummyvareval(vector, negativeRHS, coefmat):
    '''
    This function evaluates the dummy variables used to solve linear systems algebraically present on the variable vector
    
    :param vector: This variable represents the vector for which the dummy variables are being replaced
    :param negativeRHS: This variable represents a vector with the actual coordinates of the right hand side of the equation solved by the vector
    :param coefmat: This variable represents the matrix of coefficients of the system solved by the vector

    :return: This function returns the variable vector after evaluating all the dummy variables present on it
    '''
    for row in range(nvar):
        vector = vector.subs(negativeRHS.dummy[row],negativeRHS.actualcoord[row])
        for col in range(nvar):
            vector = vector.subs(coefmat.dummy[row,col],coefmat.actualcoord[row,col])
    return vector

def linearsolver(vector, negativeRHS, coefmat):
    '''
    This function solves a linear system of equations with an invertible matrix of coefficients

    :param vector: This variable represents the vector for which the dummy variables are being replaced
    :param negativeRHS: This variable represents a vector with the actual coordinates of the right hand side of the equation solved by the vector
    :param coefmat: This variable represents the matrix of coefficients of the system solved by the vector

    :return: This function returns the solution to a regular linear system of equations with an invertible matrix of coefficients
    '''
    vector.actualcoord = linsolve(Add(Mul(coefmat.dummy, vector.dummy), negativeRHS.dummy), list(vector.dummy))
    vector.actualcoord = transpose(Matrix(list(vector.actualcoord)))
    vector.actualcoord = dummyvareval(vector.actualcoord, negativeRHS, coefmat)
    return vector

def kernel_determination(vector, coefmat, criticalcol, coefsubmatrix, submatrixrows, submatrixcols):
    '''
    This function finds a vector that spans the kernel of a non-invertible matrix with a 1-dimensional kernel

    :param vector: This variable represents the vector that solves the system
    :param coefmat: This variable represents the noninvertible matrix of coefficients of the degenerate system of equations
    :param criticalcol: This variable represents the index of the column that can be removed from coefmat to obtain an invertible submatrix
    :param coefsubmatrix: This variable represents the invertible submatrix of coefficients obtained after removing a specific row and column from coefmat
    :param submatrixrows: This variable represents the indices of the rows of coefsubmatrix
    :param submatrixcols: This variable represents the indices of the columns of coefsubmatrix

    :return: This function returns a vector that spans the kernel of a non-invertible matrix with a 1-dimensional kernel
    '''   
    vector.actualcoord[criticalcol] = 1
    
    auxiliaryterm, = linsolve(Add(Mul(coefsubmatrix.dummy, Matrix(vector.actualcoord).extract(submatrixcols, [0])),
                                  coefmat.dummy.extract(submatrixrows, [criticalcol])),
                              list(Matrix(vector.actualcoord).extract(submatrixcols, [0])))
    
    vector.actualcoord[0:criticalcol] = auxiliaryterm[0:criticalcol]
    vector.actualcoord[criticalcol+1:nvar] = auxiliaryterm[criticalcol:nvar-1]
    
    vector.actualcoord = Matrix(vector.actualcoord)
        
    for row in range(nvar):
        for col in range(nvar):
            vector.actualcoord = vector.actualcoord.subs(coefmat.dummy[row, col], coefmat.actualcoord[row, col])
            if row<nvar-1 and col<nvar-1:
                vector.actualcoord = vector.actualcoord.subs(coefsubmatrix.dummy[row, col],
                                                             coefsubmatrix.actualcoord[row, col])
    
    return vector

def critical_linearsolver(vector, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols):
    '''
    This function solves a degenerate system of equations with a non-invertible matrix of coefficients with a 1-dimensional kernel

    :param vector: This variable represents the vector that solves the system
    :param negativeRHS: This variable represents the right hand side of the system of equations
    :param criticalcol: This variable represents the number of the column that is being removed to obtain an invertible submatrix from the noninvertible matrix of coefficients
    :param submatrixrows: This variable represents the indices of the rows of coefsubmatrix
    :param submatrixcols: This variable represents the indices of the columns of coefsubmatrix

    :return: This function returns a solution to a degneerate system of linear equations with a non-invertible matrix of coefficients with a 1-dimensional kernel
    '''
    vector.actualcoord[criticalcol] = 0
    
    auxiliaryterm, = linsolve(Add(Mul(coefsubmatrix.dummy, Matrix(vector.actualcoord).extract(submatrixcols, [0])),
                                  negativeRHS.dummy.extract(submatrixrows, [0])),
                              list(Matrix(vector.actualcoord).extract(submatrixcols, [0])))
        
    vector.actualcoord[0:criticalcol] = auxiliaryterm[0:criticalcol]
    vector.actualcoord[criticalcol+1:nvar] = auxiliaryterm[criticalcol:nvar-1]
    
    vector.actualcoord = Matrix(vector.actualcoord)
    
    for row in range(nvar):
        vector.actualcoord = vector.actualcoord.subs(negativeRHS.dummy[row], negativeRHS.actualcoord[row])
        for col in range(nvar - 1):
            if row<nvar - 1:
                vector.actualcoord = vector.actualcoord.subs(coefsubmatrix.dummy[row, col],
                                                             coefsubmatrix.actualcoord[row, col])
                
    return vector

def evaluation_dict(vector):
    '''
    This function generates a dictionary used to evaluate the dummy coordinates of a vector

    :param vector: This variable represents the vector that is to be evaluated

    :return: This function returns a dictionary that can be used to evaluate the dummy variables associated with the variable vector
    '''
    actual_dict = dict(zip(vector.dummy,
                           simplify(vector.actualcoord.subs(muNF, muval).subs(parameters).subs(extraparvals))))
    return actual_dict

def origin_translation(equilibrium, kinetics):
    '''
    This function shifts the kinetics of the system to make 0 the steady state under analysis

    :param equilibrium: This variable represents the coordinates of the steady state to be studied in the original system
    :param kinetics: This variable represents the reaction part of the reaction-diffusion equation that is being studied

    :return: This function returns the coordinates of the new steady state and kinetics after making 0 the steady state under study
    '''
    for varnum in range(nvar):
        kinetics = kinetics.subs(var[varnum], var[varnum] + equilibrium[varnum])
    equilibrium = [0]*nvar
    return equilibrium, kinetics
    
def check_equilibrium(equilibium, kinetics):
    '''
    This function checks that the coordinates provided for the equilibrium of the system correspond to a homogeneous steady state of the system

    :param equilibrium: This variable represents the coordinates of the steady state of the original system provided by the user
    :param kinetics: This variable represents the reaction part of the reaction-diffusion equation that is being studied

    :return: This function returns a boolean variable that is True if the steady state provided solve the equation kinetics = 0 and False, otherwise.
    '''
    for varnum in range(nvar):
        kinetics = kinetics.subs(var[varnum], equilibrium[varnum])
    kinetics = simplify(kinetics)
    if kinetics==Matrix([0]*nvar):
        return True
    else:
        return False
    
def parameter_translation(kinetics, par1, extrapara0, par2, extraparb0):
    '''
    This function shifts the parameters of the kinetics of the system to make (par1, par2) = (0, 0) the codimension-two bifurcation point

    :param kinetics: This variable represents the reaction part of the reaction-diffusion equation that is being studied
    :param par1: This variable represents the parameter a that is being shifted
    :param extrapara0: This variable represents the value of par1 at the codimension-two bifurcation point
    :param par2: This variable represents the parameter b that is being shifted
    :param extraparb0: This variable represents the value of par2 at the codimension-two bifurcation point

    :return: This function returns the kinetics of the system after shifting them to make (par1, par2) = (0, 0) the codimension-two bifurcation point
    '''
    kinetics = kinetics.subs(par1, Add(par1, extrapara0))
    kinetics = kinetics.subs(par2, Add(par2, extraparb0))
    return kinetics

def firstordereval(expr):
    '''
    This function evaluates the wavenumber, parameters, and first order eigenvectors present in the variable expr

    :param expr: This variable represents the expression that is going to be evaluated

    :return: This function returns the evaluated expression in terms of the wavenumber, parameters, and first order eigenvectors
    '''
    expr = expr.subs(phiNF_eval).subs(psiNF_eval).subs(muNF, muval).subs(parameters).subs(expandingparvals)
    return expr

def evaluation(vector):
    '''
    This function evaluates all the symbolic variables present in the variable vector

    :param vector: This variable represents the vector that is going to be evaluated

    :return: This function returns the variable vector after all the symbolic variables present on it have been evaluated
    '''
    vector.actualcoord = firstordereval(vector.actualcoord)
    try:
        vector.actualcoord = vector.actualcoord.subs(W02NF_eval).subs(W12NF_eval).subs(W22NF_eval)
    except:
        return vector
    try:
        vector.actualcoord = vector.actualcoord.subs(W03NF_eval).subs(W13NF_eval).subs(W123NF_eval)\
            .subs(W133NF_eval).subs(W23NF_eval).subs(W33NF_eval)
    except:
        return vector
    try:
        vector.actualcoord = vector.actualcoord.subs(W04NF_eval).subs(W024NF_eval).subs(W034NF_eval)\
            .subs(W14NF_eval).subs(W124NF_eval).subs(W134NF_eval).subs(W24NF_eval).subs(W224NF_eval)\
                .subs(W234NF_eval).subs(W34NF_eval).subs(W44NF_eval)
    except:
        return vector
    try:
        vector.actualcoord = vector.actualcoord.subs(W05NF_eval).subs(W025NF_eval).subs(W035NF_eval)\
            .subs(W15NF_eval).subs(W125NF_eval).subs(W135NF_eval).subs(W145NF_eval).subs(W155NF_eval)\
                .subs(W165NF_eval).subs(W175NF_eval).subs(W25NF_eval).subs(W225NF_eval).subs(W235NF_eval)\
                    .subs(W35NF_eval).subs(W325NF_eval).subs(W335NF_eval)
    except:
        return vector
    try:
        vector.actualcoord = vector.actualcoord.subs(W06NF_eval).subs(W026NF_eval).subs(W036NF_eval)\
            .subs(W046NF_eval).subs(W056NF_eval).subs(W066NF_eval).subs(W076NF_eval).subs(W16NF_eval)\
                .subs(W126NF_eval).subs(W136NF_eval).subs(W146NF_eval).subs(W156NF_eval).subs(W166NF_eval)\
                    .subs(W176NF_eval).subs(W26NF_eval).subs(W226NF_eval).subs(W236NF_eval).subs(W246NF_eval)\
                        .subs(W256NF_eval).subs(W266NF_eval).subs(W276NF_eval).subs(W286NF_eval)
    except:
        return vector
    return vector

def evaluation_alpha(alpha):
    '''
    This function evaluates all the symbolic variables present in the alpha variable

    :param alpha: This variable represents the alpha coefficient that is going to be evaluated

    :return: This function returns the variable alpha after all the symbolic variables present on it have been evaluated
    '''
    alpha = firstordereval(alpha)
    
    alpha = alpha.subs(W02NF_eval).subs(W12NF_eval).subs(W22NF_eval)
    
    alpha = alpha.subs(W03NF_eval).subs(W13NF_eval).subs(W123NF_eval).subs(W133NF_eval).subs(W23NF_eval)\
        .subs(W33NF_eval)
        
    alpha = alpha.subs(W04NF_eval).subs(W024NF_eval).subs(W034NF_eval).subs(W14NF_eval).subs(W124NF_eval)\
        .subs(W134NF_eval).subs(W24NF_eval).subs(W224NF_eval).subs(W234NF_eval).subs(W34NF_eval)\
            .subs(W44NF_eval)
            
    alpha = alpha.subs(W05NF_eval).subs(W025NF_eval).subs(W035NF_eval).subs(W15NF_eval).subs(W125NF_eval)\
        .subs(W135NF_eval).subs(W145NF_eval).subs(W155NF_eval).subs(W165NF_eval).subs(W175NF_eval)\
            .subs(W25NF_eval).subs(W225NF_eval).subs(W235NF_eval).subs(W35NF_eval).subs(W325NF_eval)\
                .subs(W335NF_eval)
                
    alpha = alpha.subs(W06NF_eval).subs(W026NF_eval).subs(W036NF_eval).subs(W046NF_eval).subs(W056NF_eval)\
        .subs(W066NF_eval).subs(W076NF_eval).subs(W16NF_eval).subs(W126NF_eval).subs(W136NF_eval)\
            .subs(W146NF_eval).subs(W156NF_eval).subs(W166NF_eval).subs(W176NF_eval).subs(W26NF_eval)\
                .subs(W226NF_eval).subs(W236NF_eval).subs(W246NF_eval).subs(W256NF_eval).subs(W266NF_eval)\
                    .subs(W276NF_eval).subs(W286NF_eval).subs(extraparvals)
                    
    return simplify(alpha)