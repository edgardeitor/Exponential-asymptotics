class Vector:
    '''This class defines a vector with two components. A dummy version of it and
    its actual coordinates'''
    def __init__(self,name):
        self.dummy = []
        self.actualcoord = []
        for varnum in range(nvar):
            self.dummy.append(symbols(name + '_' + str(varnum)))
            self.actualcoord.append(symbols(name + '_' + str(varnum)))
        self.dummy = Matrix(self.dummy)
            
class matrix:
    '''This class defines a sympy Matrix with two components. A dummy version of it and
    its actual coordinates. Here you can provide the actual matrix to be saves into the
    latter one'''
    def __init__(self, name, ncol=nvar, refmatrix=ones(nvar)):
        self.dummy = zeros(ncol)
        if refmatrix!=ones(nvar):
            self.actualcoord = refmatrix
        for row in range(ncol):
            for col in range(ncol):
                if refmatrix[row,col]!=0:
                    self.dummy[row,col] = symbols(name + '_' + str(row+1) + '^' + str(col+1))
                    
def first_order(ind_a, ind_b, u):
    '''This function is a coded version of F_1'''
    SS = zeros(nvar, 1)
    firstorderderivatives = globals()[f'firstorderderivatives{ind_a}{ind_b}']
    for i in range(nvar):
        SS = Add(SS, Mul(u.dummy[i], firstorderderivatives[i]))
    return SS
    
def second_order(ind_a, ind_b, u, v):
    '''This function is a coded version of F_2'''
    DS = zeros(nvar, 1)
    secondorderderivatives = globals()[f'secondorderderivatives{ind_a}{ind_b}']
    for i in range(nvar):
        for ii in range(nvar):
            DS = Add(DS, Mul(u.dummy[i], v.dummy[ii], secondorderderivatives[i][ii]))
    return Mul(DS, Pow(math.factorial(2), - 1))
    
def third_order(ind_a, ind_b, u, v, w):
    '''This function is a coded version of F_3'''
    TS = zeros(nvar, 1)
    thirdorderderivatives = globals()[f'thirdorderderivatives{ind_a}{ind_b}']
    for i in range(nvar):
        for ii in range(nvar):
            for iii in range(nvar):
                TS = Add(TS, Mul(u.dummy[i], v.dummy[ii], w.dummy[iii], thirdorderderivatives[i][ii][iii]))
    return Mul(TS, Pow(math.factorial(3), - 1))

def fourth_order(ind_a, ind_b, u, v, w, r):
    '''This function is a coded version of F_4'''
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
    '''This function is a coded version of F_5'''
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
    '''This function is a coded version of F_6'''
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
    '''This function is a coded version of F_7'''
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
    '''This function evaluates all the dummy variables used to solve a linear
    system in a simple way'''
    for row in range(nvar):
        vector = vector.subs(negativeRHS.dummy[row],negativeRHS.actualcoord[row])
        for col in range(nvar):
            vector = vector.subs(coefmat.dummy[row,col],coefmat.actualcoord[row,col])
    return vector

def linearsolver(vector, negativeRHS, coefmat):
    '''This function solves a linear system with dummy variables and then uses
    the previous function to evaluate the dummy variables'''
    vector.actualcoord = linsolve(Add(Mul(coefmat.dummy, vector.dummy), negativeRHS.dummy), list(vector.dummy))
    vector.actualcoord = transpose(Matrix(list(vector.actualcoord)))
    vector.actualcoord = dummyvareval(vector.actualcoord, negativeRHS, coefmat)
    return vector

def kernel_determination(vector, coefmat, criticalcol, coefsubmatrix, submatrixrows, submatrixcols):    
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
    actual_dict = dict(zip(vector.dummy,
                           simplify(vector.actualcoord.subs(muNF, muval).subs(parameters).subs(extraparvals))))
    return actual_dict

def origin_translation(equilibrium, kinetics):
    for varnum in range(nvar):
        kinetics = kinetics.subs(var[varnum], var[varnum] + equilibrium[varnum])
    equilibrium = [0]*nvar
    return equilibrium, kinetics
    
def check_equilibrium(equilibium, kinetics):
    for varnum in range(nvar):
        kinetics = kinetics.subs(var[varnum], equilibrium[varnum])
    kinetics = simplify(kinetics)
    if kinetics==Matrix([0]*nvar):
        return True
    else:
        return False
    
def parameter_translation(kinetics, par1, extrapara0, par2, extraparb0):
    kinetics = kinetics.subs(par1, Add(par1, extrapara0))
    kinetics = kinetics.subs(par2, Add(par2, extraparb0))
    return kinetics

def firstordereval(expr):
    expr = expr.subs(phiNF_eval).subs(psiNF_eval).subs(muNF, muval).subs(parameters).subs(expandingparvals)
    return expr

def evaluation(vector):
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