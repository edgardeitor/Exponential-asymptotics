from IPython import get_ipython
get_ipython().run_line_magic('reset', '-sf')

import os
import shutil
from sys import exit
from sympy import *
init_printing()
from sympy.solvers import solve
from mpmath import findroot
import math

modelname = input('Enter the name of the system you would like to analyze: ')

if not os.path.isdir(modelname):
    print('The directory of that system was not found. Create it first and place ' + 
          'the data file inside named in the same way.')
    exit()
else:
    os.chdir(os.getcwd() + '\\' + modelname)

try:
    exec(open(modelname + '.py').read())
except:
    print('The file ' + modelname + '.py could not be run')
    exit()
    
nvar = len(var)

try:
    exec(open(os.path.dirname(os.path.realpath(__file__)) + '\\functions.py').read())
except:
    print('File functions.py is not in the same folder as the script you are running')
    exit()

file = open('List of Variables.txt', 'w')

try:
    var
except:
    print('Variables were not provided')
    exit()

for varnum in range(nvar):
    try:
        exec(var[varnum] + ' = symbols(var[varnum], real=True)')
        var[varnum] = eval(var[varnum])
    except:
        print('The script could not define ' + (var[varnum]) + ' as a variable')
        exit()
    file.write(latex(var[varnum]) + '\n')

try:
    parameters
except:
    print('Parameters were not provided. The script will assume that there are no parameters.')
    parameters=[]

npar = len(parameters) + len(unevaluatedparameters)

if npar>0:
    newparameters = dict()
    for key in parameters.keys():
        try:
            exec(key + ' = symbols(key, real = True)')
            newparameters[eval(key)] = parameters[key]
        except:
            print('The script could not define your variable ' + key + ' as a variable')
            exit()
    parameters = newparameters

for parnum in range(len(unevaluatedparameters)):
    exec(unevaluatedparameters[parnum] + ' = symbols(unevaluatedparameters[parnum], real=True)')
    unevaluatedparameters[parnum] = eval(unevaluatedparameters[parnum])
    
try:
    diffmatrix
except:
    print('The diffusion matrix was not provided correctly.')
    exit()
    
exec(par1 + ' = symbols(par1, real = True)')
exec(par2 + ' = symbols(par2, real = True)')

par1 = eval(par1)
par2 = eval(par2)

expandingparvals = {par1: 0, par2: 0}

for varnum in range(nvar):
    try:
        equilibrium[varnum] = eval(equilibrium[varnum])
    except:
        pass

try:
    for row in range(nvar):
        for col in range(nvar):
            diffmatrix[row][col] = eval(str(diffmatrix[row][col]).replace('^', '**'))
except:
    print('The diffusion matrix is not a function of the parameters of the system')
    exit()
        
diffmatrix = Matrix(diffmatrix)

# Kinetics

try:
    kinetics
except:
    print('Kinetics were not provided.')
    exit()
    
for functionnumber in range(nvar):
    try:
        kinetics[functionnumber] = eval(str(kinetics[functionnumber]).replace('^', '**'))
    except:
        print('The expression ' + kinetics[functionnumber] + ' is not a function of the parameters of your system')
        exit()

kinetics = Matrix(kinetics)

if not check_equilibrium(equilibrium, kinetics):
    print('The coordinates of the equilibrium do not solve the trivial system.')
    exit()

if equilibrium!=[0]*nvar:
    equilibrium, kinetics = origin_translation(equilibrium, kinetics)
    
aNF = []
bNF = []

for counter in range(7):
    aNF.append(symbols(str(par1) + 'NF' + str(counter), real = True))
    bNF.append(symbols(str(par2) + 'NF' + str(counter), real = True))

extraparvals = {}

extraparvals = extrapars(extraparvals)

kinetics = parameter_translation(kinetics, par1, aNF[0], par2, bNF[0])

file = open('Kinetics.txt', 'w')
for functionnumber in range(nvar):
    file.write(latex(kinetics[functionnumber]) + '\n')
file.close()

jacobianmat = kinetics.jacobian(var).subs(par1, 0).subs(par2, 0)
jacobianmat = eqsubstitution(jacobianmat)

# Beyond-all-order asymptotics

muNF = symbols(r'\mu_{NF}', real=true)

try:
    muval = eval(muval)
except:
    pass

for counter in range(5):
    exec(f'coefmat{counter} = matrix("coef{counter}mat", nvar, Add(jacobianmat, Mul(- Pow(counter, 2), muNF, diffmatrix)))')

TC1 = coefmat1.dummy.det()

for row in range(nvar):
    for col in range(nvar):
        TC1 = TC1.subs(coefmat1.dummy[row, col], coefmat1.actualcoord[row, col])

for parnum1 in range(7):
    for parnum2 in range(7):
        exec('f' + str(parnum1) + str(parnum2) + '= simplify(Mul(Pow(Mul(math.factorial(parnum1), ' + 
             'math.factorial(parnum2)), -1), diff(kinetics, par1, parnum1, par2, parnum2))' +
             '.subs(par1, 0).subs(par2, 0))')

negativeRHS = Vector('negativeRHS')

for parnum1 in range(7):
    for parnum2 in range(7):
        exec(f'firstorderderivatives{parnum1}{parnum2} = list()')
        if parnum1 + parnum2<=5:
            exec(f'secondorderderivatives{parnum1}{parnum2} = list()')
            if parnum1 + parnum2<=4:
                exec(f'thirdorderderivatives{parnum1}{parnum2} = list()')
                if parnum1 + parnum2<=3:
                    exec(f'fourthorderderivatives{parnum1}{parnum2} = list()')
                    if parnum1 + parnum2<=2:
                        exec(f'fifthorderderivatives{parnum1}{parnum2} = list()')
                        if parnum1 + parnum2<=1:
                            exec(f'sixthorderderivatives{parnum1}{parnum2} = list()')
                            if parnum1 + parnum2==0:
                                exec(f'seventhorderderivatives{parnum1}{parnum2} = list()')
        for counter1 in range(nvar):
            exec(f'firstorderderivatives{parnum1}{parnum2}.append(diff(f{parnum1}{parnum2}, var[{counter1}]))')
            if parnum1 + parnum2<=5:
                exec(f'secondorderderivatives{parnum1}{parnum2}.append(list())')
                if parnum1 + parnum2<=4:
                    exec(f'thirdorderderivatives{parnum1}{parnum2}.append(list())')
                    if parnum1 + parnum2<=3:
                        exec(f'fourthorderderivatives{parnum1}{parnum2}.append(list())')
                        if parnum1 + parnum2<=2:
                            exec(f'fifthorderderivatives{parnum1}{parnum2}.append(list())')
                            if parnum1 + parnum2<=1:
                                exec(f'sixthorderderivatives{parnum1}{parnum2}.append(list())')
                                if parnum1 + parnum2==0:
                                    exec(f'seventhorderderivatives{parnum1}{parnum2}.append(list())')
            for counter2 in range(nvar):
                if parnum1 + parnum2<=5:
                    exec(f'secondorderderivatives{parnum1}{parnum2}[{counter1}].append(diff(firstorderderivatives{parnum1}{parnum2}[{counter1}], var[{counter2}]))')
                    if parnum1 + parnum2<=4:
                        exec(f'thirdorderderivatives{parnum1}{parnum2}[{counter1}].append(list())')
                        if parnum1 + parnum2<=3:
                            exec(f'fourthorderderivatives{parnum1}{parnum2}[{counter1}].append(list())')
                            if parnum1 + parnum2<=2:
                                exec(f'fifthorderderivatives{parnum1}{parnum2}[{counter1}].append(list())')
                                if parnum1 + parnum2<=1:
                                    exec(f'sixthorderderivatives{parnum1}{parnum2}[{counter1}].append(list())')
                                    if parnum1 + parnum2==0:
                                        exec(f'seventhorderderivatives{parnum1}{parnum2}[{counter1}].append(list())')
                for counter3 in range(nvar):
                    if parnum1 + parnum2<=4:
                        exec(f'thirdorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}].append(diff(secondorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}], var[{counter3}]))')
                        if parnum1 + parnum2<=3:
                            exec(f'fourthorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}].append(list())')
                            if parnum1 + parnum2<=2:
                                exec(f'fifthorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}].append(list())')
                                if parnum1 + parnum2<=1:
                                    exec(f'sixthorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}].append(list())')
                                    if parnum1 + parnum2==0:
                                        exec(f'seventhorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}].append(list())')
                    for counter4 in range(nvar):
                        if parnum1 + parnum2<=3:
                            exec(f'fourthorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}].append(diff(thirdorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}], var[{counter4}]))')
                            if parnum1 + parnum2<=2:
                                exec(f'fifthorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}].append(list())')
                                if parnum1 + parnum2<=1:
                                    exec(f'sixthorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}].append(list())')
                                    if parnum1 + parnum2==0:
                                        exec(f'seventhorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}].append(list())')
                        for counter5 in range(nvar):
                            if parnum1 + parnum2<=2:
                                exec(f'fifthorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}][{counter4}].append(diff(fourthorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}][{counter4}], var[{counter5}]))')
                                if parnum1 + parnum2<=1:
                                    exec(f'sixthorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}][{counter4}].append(list())')
                                    if parnum1 + parnum2==0:
                                        exec(f'seventhorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}][{counter4}].append(list())')
                            for counter6 in range(nvar):
                                if parnum1 + parnum2<=1:
                                    exec(f'sixthorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}][{counter4}][{counter5}].append(diff(fifthorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}][{counter4}][{counter5}], var[{counter6}]))')
                                    if parnum1 + parnum2==0:
                                        exec(f'seventhorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}][{counter4}][{counter5}].append(list())')
                                        for counter7 in range(nvar):
                                            exec(f'seventhorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}][{counter4}][{counter5}][{counter6}].append(diff(sixthorderderivatives{parnum1}{parnum2}[{counter1}][{counter2}][{counter3}][{counter4}][{counter5}][{counter6}], var[{counter7}]))')
                            
phiNF = Vector('phi^NF')
psiNF = Vector('psi^NF')

# Order 2 vectors

for counter1 in range(3):
    exec(f'W{counter1}2NF = Vector("W{counter1}2^NF")')
    
# Order 3 vectors
    
for counter1 in range(4):
    if counter1!=1:
        exec(f'W{counter1}3NF = Vector("W{counter1}3^NF")')
    else:
        for counter2 in range(1, 4):
            if counter2==1:
                exec(f'W{counter1}3NF = Vector("W{counter1}3^NF")')
            else:
                exec(f'W{counter1}{counter2}3NF = Vector("W{counter1}{counter2}3^NF")')

# Order 4 vectors
        
for counter1 in range(5):
    if counter1<3:
        for counter2 in range(1, 4):
            if counter2==1:
                exec(f'W{counter1}4NF = Vector("W{counter1}4^NF")')
            else:
                exec(f'W{counter1}{counter2}4NF = Vector("W{counter1}{counter2}4^NF")')
    else:
        exec(f'W{counter1}4NF = Vector("W{counter1}4^NF")')

# Order 5 vectors

for counter1 in range(4):
    if counter1!=1:
        for counter2 in range(1, 4):
            if counter2==1:
                exec(f'W{counter1}5NF = Vector("W{counter1}5^NF")')
            else:
                exec(f'W{counter1}{counter2}5NF = Vector("W{counter1}{counter2}5^NF")')
    else:
        for counter2 in range(1, 8):
            if counter2==1:
                exec(f'W{counter1}5NF = Vector("W{counter1}5^NF")')
            else:
                exec(f'W{counter1}{counter2}5NF = Vector("W{counter1}{counter2}5^NF")')

# Order 6 vectors
    
for counter1 in range(3):
    if counter1!=2:
        for counter2 in range(1, 8):
            if counter2==1:
                exec(f'W{counter1}6NF = Vector("W{counter1}6^NF")')
            else:
                exec(f'W{counter1}{counter2}6NF = Vector("W{counter1}{counter2}6^NF")')
    else:
        for counter2 in range(1, 9):
            if counter2==1:
                exec(f'W{counter1}6NF = Vector("W{counter1}6^NF")')
            else:
                exec(f'W{counter1}{counter2}6NF = Vector("W{counter1}{counter2}6^NF")')

getout = 0

for row in range(nvar):
    for col in range(nvar):
        submatrixrows = list(range(nvar))
        submatrixcols = list(range(nvar))
        submatrixrows.remove(row)
        submatrixcols.remove(col)
        invertiblesubmatrix = coefmat1.actualcoord.extract(submatrixrows, submatrixcols)
        submatrixeval = invertiblesubmatrix
        submatrixeval = submatrixeval.subs(muNF, muval).subs(parameters)
        for varnum in range(nvar):
            submatrixeval = submatrixeval.subs(var[varnum], equilibrium[varnum])
        try:
            if abs(N(submatrixeval.det())) > tol:
                phiNF.actualcoord[col] = 1
                criticalrow = row
                criticalcol = col
                getout = 1
                break
        except:
            phiNF.actualcoord[col] = 1
            criticalrow = row
            criticalcol = col
            getout = 1
            break
    if getout==1:
        break
    
coefsubmatrix = matrix('dummysubmatrix', nvar - 1, invertiblesubmatrix)

phiNF = kernel_determination(phiNF, coefmat1, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

if phiunit=='y':
    phiNF.actualcoord = Mul(Pow(sqrt(phiNF.actualcoord.dot(phiNF.actualcoord)), - 1), phiNF.actualcoord)

Tcoefsubmatrix = matrix('Tdummysubmatrix', nvar - 1, transpose(invertiblesubmatrix))
Tcoefmat1 = matrix('Tcoefmat', nvar, transpose(coefmat1.actualcoord))

psiNF = kernel_determination(psiNF, Tcoefmat1, criticalrow, Tcoefsubmatrix, submatrixcols, submatrixrows)

for row in range(nvar):
    for col in range(nvar):
        psiNF.actualcoord = psiNF.actualcoord.subs(coefmat1.dummy[row, col], coefmat1.actualcoord[row, col])
        if row<nvar-1 and col<nvar-1:
            psiNF.actualcoord = psiNF.actualcoord.subs(coefsubmatrix.dummy[row, col],
                                                       coefsubmatrix.actualcoord[row, col])

phiNF_eval = evaluation_dict(phiNF)
psiNF_eval = evaluation_dict(psiNF)

print('First order ready')

DS_phiphi_00 = second_order(0, 0, phiNF, phiNF)

negativeRHS.actualcoord = DS_phiphi_00

W02NF = linearsolver(W02NF, negativeRHS, coefmat0)

SS_phi_10 = first_order(1, 0, phiNF)
SS_phi_01 = first_order(0, 1, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_phi_10), Mul(bNF[1], SS_phi_01))

equation1 = psiNF.dummy.dot(negativeRHS.actualcoord)

equation1 = simplify(firstordereval(equation1))

if N(equation1.subs(extraparvals))!=0:
    if simp=='y':
        tempvar = solve(equation1.subs(extraparvals), dict = True)[0]
        for par in tempvar.keys():
            tempvar[par] = simplify(tempvar[par])
        extraparvals = extraparvals | tempvar
    else:
        extraparvals = extraparvals | solve(equation1.subs(extraparvals), dict = True)[0]
    equation1 = equation1.subs(extraparvals)

negativeRHS.actualcoord = negativeRHS.actualcoord.subs(extraparvals)

W12NF = critical_linearsolver(W12NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

negativeRHS.actualcoord = DS_phiphi_00

W22NF = linearsolver(W22NF, negativeRHS, coefmat2)

W02NF = evaluation(W02NF)
W12NF = evaluation(W12NF)
W22NF = evaluation(W22NF)

W02NF_eval = evaluation_dict(W02NF)
W12NF_eval = evaluation_dict(W12NF)
W22NF_eval = evaluation_dict(W22NF)

print('Second order ready')

SS_W02_10 = first_order(1, 0, W02NF)
SS_W02_01 = first_order(0, 1, W02NF)
DS_phiW12_00 = second_order(0, 0, phiNF, W12NF)
DS_phiphi_10 = second_order(1, 0, phiNF, phiNF)
DS_phiphi_01 = second_order(0, 1, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W02_10), Mul(bNF[1], SS_W02_01), Mul(2, DS_phiW12_00),
                              Mul(aNF[1], DS_phiphi_10), Mul(bNF[1], DS_phiphi_01)).subs(extraparvals)

W03NF = linearsolver(W03NF, negativeRHS, coefmat0)

SS_W12_10 = first_order(1, 0, W12NF)
SS_W12_01 = first_order(0, 1, W12NF)
SS_phi_20 = first_order(2, 0, phiNF)
SS_phi_11 = first_order(1, 1, phiNF)
SS_phi_02 = first_order(0, 2, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W12_10), Mul(bNF[1], SS_W12_01), Mul(aNF[2], SS_phi_10),
                              Mul(bNF[2], SS_phi_01), Mul(Pow(aNF[1], 2), SS_phi_20),
                              Mul(aNF[1], bNF[1], SS_phi_11),
                              Mul(Pow(bNF[1], 2), SS_phi_02)).subs(extraparvals)

equation2 = psiNF.dummy.dot(negativeRHS.actualcoord)

equation2 = firstordereval(equation2)
equation2 = simplify(equation2.subs(W12NF_eval))

W13NF = critical_linearsolver(W13NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

# if equation2!=0:
#     if len(solve(equation2, dict = True))==0:
#         print('There are no parameters to find a codimension-two point.')
#         exit()
#     if len(solve(equation2, dict = True))==1:
#         extraparvals = extraparvals | solve(equation2, dict = True)[0]
#     else:
#         print("You have to choose a solution:")
#         print(solve(equation2.subs(extraparvals), dict = True))
#         solchoice = input("Which solution would you like to use? ")
#         while True:
#             try:
#                 int(solchoice)
#                 extraparvals = extraparvals | solve(equation2.subs(extraparvals), dict = True)[int(solchoice) - 1]
#                 break
#             except:
#                 solchoice = input("You did not provide a valid integer. Which solution would you like to use? ")
    
# if N(equation2)!=0 and equation2.subs(aNF[1], 0)==0:
#     extraparvals[aNF[1]] = 0
# if N(equation2)!=0 and equation2.subs(bNF[1], 0)==0:
#     extraparvals[bNF[1]] = 0
# if N(equation2)!=0 and equation2.subs(aNF[2], 0)==0:
#     extraparvals[aNF[2]] = 0
# if N(equation2)!=0 and equation2.subs(bNF[2], 0)==0:
#     extraparvals[bNF[2]] = 0

DS_phiW02_00 = second_order(0, 0, phiNF, W02NF)
DS_phiW22_00 = second_order(0, 0, phiNF, W22NF)
TS_phiphiphi_00 = third_order(0, 0, phiNF, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(4, DS_phiW02_00), Mul(2, DS_phiW22_00), Mul(3, TS_phiphiphi_00))

Cod2 = psiNF.dummy.dot(negativeRHS.actualcoord)
    
W123NF = critical_linearsolver(W123NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

negativeRHS.actualcoord = Mul(2, sqrt(muNF), diffmatrix, phiNF.dummy)

TC2 = psiNF.dummy.dot(negativeRHS.actualcoord)

W133NF = critical_linearsolver(W133NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

SS_W22_10 = first_order(1, 0, W22NF)
SS_W22_01 = first_order(0, 1, W22NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W22_10), Mul(bNF[1], SS_W22_01), Mul(2, DS_phiW12_00),
                              Mul(aNF[1], DS_phiphi_10), Mul(bNF[1], DS_phiphi_01)).subs(extraparvals)

W23NF = linearsolver(W23NF, negativeRHS, coefmat2)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW22_00), TS_phiphiphi_00)

W33NF = linearsolver(W33NF, negativeRHS, coefmat3)

W03NF = evaluation(W03NF)
W13NF = evaluation(W13NF)
W123NF = evaluation(W123NF)
W133NF = evaluation(W133NF)
W23NF = evaluation(W23NF)
W33NF = evaluation(W33NF)

W03NF_eval = evaluation_dict(W03NF)
W13NF_eval = evaluation_dict(W13NF)
W123NF_eval = evaluation_dict(W123NF)
W133NF_eval = evaluation_dict(W133NF)
W23NF_eval = evaluation_dict(W23NF)
W33NF_eval = evaluation_dict(W33NF)

print('Third order ready')

SS_W03_10 = first_order(1, 0, W03NF)
SS_W03_01 = first_order(0, 1, W03NF)
SS_W02_20 = first_order(2, 0, W02NF)
SS_W02_11 = first_order(1, 1, W02NF)
SS_W02_02 = first_order(0, 2, W02NF)
DS_phiW13_00 = second_order(0, 0, phiNF, W13NF)
DS_W12W12_00 = second_order(0, 0, W12NF, W12NF)
DS_phiW12_10 = second_order(1, 0, phiNF, W12NF)
DS_phiW12_01 = second_order(0, 1, phiNF, W12NF)
DS_phiphi_20 = second_order(2, 0, phiNF, phiNF)
DS_phiphi_11 = second_order(1, 1, phiNF, phiNF)
DS_phiphi_02 = second_order(0, 2, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W03_10), Mul(bNF[1], SS_W03_01), Mul(aNF[2], SS_W02_10),
                              Mul(bNF[2], SS_W02_01), Mul(Pow(aNF[1], 2), SS_W02_20),
                              Mul(aNF[1], bNF[1], SS_W02_11), Mul(Pow(bNF[1], 2), SS_W02_02),
                              Mul(2, DS_phiW13_00), DS_W12W12_00,
                              Mul(2, aNF[1], DS_phiW12_10), Mul(2, bNF[1], DS_phiW12_01),
                              Mul(aNF[2], DS_phiphi_10), Mul(bNF[2], DS_phiphi_01),
                              Mul(Pow(aNF[1], 2), DS_phiphi_20), Mul(aNF[1], bNF[1], DS_phiphi_11),
                              Mul(Pow(bNF[1], 2), DS_phiphi_02)).subs(extraparvals)

W04NF = linearsolver(W04NF, negativeRHS, coefmat0)

DS_phiW123_00 = second_order(0, 0, phiNF, W123NF)
DS_W02W02_00 = second_order(0, 0, W02NF, W02NF)
DS_W22W22_00 = second_order(0, 0, W22NF, W22NF)
TS_phiphiW02_00 = third_order(0, 0, phiNF, phiNF, W02NF)
TS_phiphiW22_00 = third_order(0, 0, phiNF, phiNF, W22NF)
Q4S_phiphiphiphi_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW123_00), Mul(2, DS_W02W02_00), DS_W22W22_00,
                              Mul(6, TS_phiphiW02_00), Mul(3, TS_phiphiW22_00),
                              Mul(3, Q4S_phiphiphiphi_00))

W024NF = linearsolver(W024NF, negativeRHS, coefmat0)

DS_phiW133_00 = second_order(0, 0, phiNF, W133NF)

negativeRHS.actualcoord = Mul(2, DS_phiW133_00)

W034NF = linearsolver(W034NF, negativeRHS, coefmat0)

SS_W13_10 = first_order(1, 0, W13NF)
SS_W13_01 = first_order(0, 1, W13NF)
SS_W12_20 = first_order(2, 0, W12NF)
SS_W12_11 = first_order(1, 1, W12NF)
SS_W12_02 = first_order(0, 2, W12NF)
SS_phi_30 = first_order(3, 0, phiNF)
SS_phi_21 = first_order(2, 1, phiNF)
SS_phi_12 = first_order(1, 2, phiNF)
SS_phi_03 = first_order(0, 3, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W13_10), Mul(bNF[1], SS_W13_01), Mul(aNF[2], SS_W12_10),
                              Mul(bNF[2], SS_W12_01), Mul(aNF[3], SS_phi_10), Mul(bNF[3], SS_phi_01),
                              Mul(Pow(aNF[1], 2), SS_W12_20), Mul(aNF[1], bNF[1], SS_W12_11),
                              Mul(Pow(bNF[1], 2), SS_W12_02), Mul(2, aNF[1], aNF[2], SS_phi_20),
                              Mul(Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), SS_phi_11),
                              Mul(2, bNF[1], bNF[2], SS_phi_02), Mul(Pow(aNF[1], 3), SS_phi_30),
                              Mul(Pow(aNF[1], 2), bNF[1], SS_phi_21),
                              Mul(aNF[1], Pow(bNF[1], 2), SS_phi_12),
                              Mul(Pow(bNF[1], 3), SS_phi_03)).subs(extraparvals)
            
equation3 = psiNF.dummy.dot(negativeRHS.actualcoord)

equation3 = firstordereval(equation3)
equation3 = equation3.subs(W12NF_eval)
equation3 = simplify(equation3.subs(W13NF_eval))

negativeRHS.actualcoord = negativeRHS.actualcoord.subs(extraparvals)

W14NF = critical_linearsolver(W14NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

SS_W123_10 = first_order(1, 0, W123NF)
SS_W123_01 = first_order(0, 1, W123NF)
DS_phiW03_00 = second_order(0, 0, phiNF, W03NF)
DS_phiW23_00 = second_order(0, 0, phiNF, W23NF)
DS_W12W02_00 = second_order(0, 0, W12NF, W02NF)
DS_W12W22_00 = second_order(0, 0, W12NF, W22NF)
DS_phiW02_10 = second_order(1, 0, phiNF, W02NF)
DS_phiW22_10 = second_order(1, 0, phiNF, W22NF)
DS_phiW02_01 = second_order(0, 1, phiNF, W02NF)
DS_phiW22_01 = second_order(0, 1, phiNF, W22NF)
TS_phiphiW12_00 = third_order(0, 0, phiNF, phiNF, W12NF)
TS_phiphiphi_10 = third_order(1, 0, phiNF, phiNF, phiNF)
TS_phiphiphi_01 = third_order(0, 1, phiNF, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W123_10), Mul(bNF[1], SS_W123_01),
                              Mul(4, DS_phiW03_00), Mul(2, DS_phiW23_00),
                              Mul(4, DS_W12W02_00), Mul(2, DS_W12W22_00),
                              Mul(4, aNF[1], DS_phiW02_10), Mul(2, aNF[1], DS_phiW22_10),
                              Mul(4, bNF[1], DS_phiW02_01), Mul(2, bNF[1], DS_phiW22_01),
                              Mul(9, TS_phiphiW12_00), Mul(3, aNF[1], TS_phiphiphi_10),
                              Mul(3, bNF[1], TS_phiphiphi_01)).subs(extraparvals)
            
equation4 = psiNF.dummy.dot(negativeRHS.actualcoord)

equation4 = firstordereval(equation4)
equation4 = equation4.subs(W02NF_eval).subs(W12NF_eval).subs(W22NF_eval)
equation4 = simplify(equation4.subs(W03NF_eval).subs(W123NF_eval).subs(W23NF_eval))

if N(equation4.subs(extraparvals))!=0:
    if simp=='y':
        tempvar = solve(equation4.subs(extraparvals), dict = True)[0]
        for par in tempvar.keys():
            tempvar[par] = simplify(tempvar[par])
        extraparvals = extraparvals | tempvar
    else:
        extraparvals = extraparvals | solve(equation4.subs(extraparvals), dict = True)[0]
    equation4 = equation4.subs(extraparvals)
    
try:
    extraparvals[aNF[1]] = extraparvals[aNF[1]].subs(extraparvals)
except:
    pass

try:
    extraparvals[bNF[1]] = extraparvals[bNF[1]].subs(extraparvals)
except:
    pass

if equation3.subs(extraparvals)!=0:
    if len(solve(equation3, dict = True))==1:
        if simp=='y':
            tempvar = solve(equation3.subs(extraparvals), dict = True)[0]
            for par in tempvar.keys():
                tempvar[par] = simplify(tempvar[par])
            extraparvals = extraparvals | tempvar
        else:
            extraparvals = extraparvals | solve(equation3.subs(extraparvals), dict = True)[0]
        equation3 = 0
    else:
        print("You have to choose a solution to solve " + latex(equation3) + "= 0:")
        print(solve(equation3.subs(extraparvals), dict = True))
        solchoice = input("Which solution would you like to use? ")
        while True:
            try:
                if simp=='y':
                    tempvar = solve(equation3.subs(extraparvals), dict = True)[int(solchoice) - 1]
                    for par in tempvar.keys():
                        tempvar[par] = simplify(tempvar[par])
                    extraparvals = extraparvals | tempvar
                else:
                    extraparvals = extraparvals | solve(equation3.subs(extraparvals), dict = True)[int(solchoice) - 1]
                break
            except:
                solchoice = input("You did not provide a valid integer. Which solution would you like to use? ")
    
try:
    extraparvals[aNF[1]] = extraparvals[aNF[1]].subs(extraparvals)
except:
    pass

try:
    extraparvals[bNF[1]] = extraparvals[bNF[1]].subs(extraparvals)
except:
    pass

negativeRHS.actualcoord = negativeRHS.actualcoord.subs(extraparvals)

W124NF = critical_linearsolver(W124NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

SS_W133_10 = first_order(1, 0, W133NF)
SS_W133_01 = first_order(0, 1, W133NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W133_10), Mul(bNF[1], SS_W133_01),
                              Mul(2, sqrt(muNF), diffmatrix, W12NF.dummy)).subs(extraparvals)
            
equation5 = psiNF.dummy.dot(negativeRHS.actualcoord)

equation5 = firstordereval(equation5)
equation5 = equation5.subs(W12NF_eval)
equation5 = simplify(equation5.subs(W133NF_eval))

if N(equation5.subs(extraparvals))!=0:
    if simp=='y':
        tempvar = solve(equation5.subs(extraparvals), dict = True)[0]
        for par in tempvar.keys():
            tempvar[par] = simplify(tempvar[par])
        extraparvals = extraparvals | tempvar
    else:
        extraparvals = extraparvals | solve(equation5.subs(extraparvals), dict = True)[0]
    equation5 = equation5.subs(extraparvals)
    
try:
    extraparvals[aNF[1]] = extraparvals[aNF[1]].subs(extraparvals)
except:
    pass

try:
    extraparvals[bNF[1]] = extraparvals[bNF[1]].subs(extraparvals)
except:
    pass

W134NF = critical_linearsolver(W134NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

SS_W23_10 = first_order(1, 0, W23NF)
SS_W23_01 = first_order(0, 1, W23NF)
SS_W22_20 = first_order(2, 0, W22NF)
SS_W22_11 = first_order(1, 1, W22NF)
SS_W22_02 = first_order(0, 2, W22NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W23_10), Mul(bNF[1], SS_W23_01), Mul(aNF[2], SS_W22_10),
                              Mul(bNF[2], SS_W22_01), Mul(Pow(aNF[1], 2), SS_W22_20),
                              Mul(aNF[1], bNF[1], SS_W22_11), Mul(Pow(bNF[1], 2), SS_W22_02),
                              Mul(2, DS_phiW13_00), DS_W12W12_00,
                              Mul(2, aNF[1], DS_phiW12_10), Mul(2, bNF[1], DS_phiW12_01),
                              Mul(aNF[2], DS_phiphi_10), Mul(bNF[2], DS_phiphi_01),                        
                              Mul(Pow(aNF[1], 2), DS_phiphi_20), Mul(aNF[1], bNF[1], DS_phiphi_11),
                              Mul(Pow(bNF[1], 2), DS_phiphi_02)).subs(extraparvals)

W24NF = linearsolver(W24NF, negativeRHS, coefmat2)

DS_phiW33_00 = second_order(0, 0, phiNF, W33NF)
DS_W02W22_00 = second_order(0, 0, W02NF, W22NF)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW123_00), Mul(2, DS_phiW33_00), Mul(4, DS_W02W22_00),
                              Mul(6, TS_phiphiW02_00), Mul(6, TS_phiphiW22_00), Mul(4, Q4S_phiphiphiphi_00))

W224NF = linearsolver(W224NF, negativeRHS, coefmat2)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW133_00), Mul(8, sqrt(muNF), diffmatrix, W22NF.dummy))

W234NF = linearsolver(W234NF, negativeRHS, coefmat2)

SS_W33_10 = first_order(1, 0, W33NF)
SS_W33_01 = first_order(0, 1, W33NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W33_10), Mul(bNF[1], SS_W33_01),
                              Mul(2, DS_phiW23_00), Mul(2, DS_W12W22_00),
                              Mul(2, aNF[1], DS_phiW22_10), Mul(2, bNF[1], DS_phiW22_01),
                              Mul(3, TS_phiphiW12_00), Mul(aNF[1], TS_phiphiphi_10),
                              Mul(bNF[1], TS_phiphiphi_01)).subs(extraparvals)

W34NF = linearsolver(W34NF, negativeRHS, coefmat3)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW33_00), DS_W22W22_00, Mul(3, TS_phiphiW22_00),
                              Q4S_phiphiphiphi_00)

W44NF = linearsolver(W44NF, negativeRHS, coefmat4)

W04NF = evaluation(W04NF)
W024NF = evaluation(W024NF)
W034NF = evaluation(W034NF)
W14NF = evaluation(W14NF)
W124NF = evaluation(W124NF)
W134NF = evaluation(W134NF)
W24NF = evaluation(W24NF)
W224NF = evaluation(W224NF)
W234NF = evaluation(W234NF)
W34NF = evaluation(W34NF)
W44NF = evaluation(W44NF)

W04NF_eval = evaluation_dict(W04NF)
W024NF_eval = evaluation_dict(W024NF)
W034NF_eval = evaluation_dict(W034NF)
W14NF_eval = evaluation_dict(W14NF)
W124NF_eval = evaluation_dict(W124NF)
W134NF_eval = evaluation_dict(W134NF)
W24NF_eval = evaluation_dict(W24NF)
W224NF_eval = evaluation_dict(W224NF)
W234NF_eval = evaluation_dict(W234NF)
W34NF_eval = evaluation_dict(W34NF)
W44NF_eval = evaluation_dict(W44NF)

print('Fourth order ready')

# Final evaluation

if equation2.subs(extraparvals)!=0:
    if len(solve(equation2.subs(extraparvals), dict = True))==0:
        print('There are no parameters to produce a codimension-two point')
        exit()
    elif len(solve(equation2.subs(extraparvals), dict = True))==1:
        if simp=='y':
            tempvar = solve(equation2.subs(extraparvals), dict = True)[0]
            for par in tempvar.keys():
                tempvar[par] = simplify(tempvar[par])
            extraparvals = extraparvals | tempvar
        else:
            extraparvals = extraparvals | solve(equation2.subs(extraparvals), dict = True)[0]
    elif len(solve(equation2.subs(extraparvals), dict = True))>1:
        print("You have to choose a solution:")
        print(solve(equation2.subs(extraparvals), dict = True))
        solchoice = input("Which solution would you like to use? ")
        while True:
            try:
                if simp=='y':
                    tempvar = solve(equation2.subs(extraparvals), dict = True)[int(solchoice) - 1]
                    for par in tempvar.keys():
                        tempvar[par] = simplify(tempvar[par])
                    extraparvals = extraparvals | tempvar
                else:
                    extraparvals = extraparvals | solve(equation2.subs(extraparvals), dict = True)[int(solchoice) - 1]
                break
            except:
                solchoice = input("You did not provide a valid integer. Which solution would you like to use? ")

W03NF_eval = evaluation_dict(W03NF)
W13NF_eval = evaluation_dict(W13NF)
W23NF_eval = evaluation_dict(W23NF)

W04NF_eval = evaluation_dict(W04NF)
W24NF_eval = evaluation_dict(W24NF)

TC1 = simplify(firstordereval(TC1).subs(extraparvals))
try:
    if abs(TC1)>tol:
        print('You did not provide a Turing bifurcation point')
except:
    print('You did not provide a bifurcation point. You have to solve\n')
    print(str(TC1) + ' = 0\n')
    
TC2 = simplify(firstordereval(TC2).subs(extraparvals))
try:
    if abs(TC2)>tol:
        print('You did not provide a Turing bifurcation point')
except:
    print('You did not provide a bifurcation point. You have to solve\n')
    print(str(TC2) + ' = 0\n')
    
Cod2 = simplify(firstordereval(Cod2))
Cod2 = Cod2.subs(W02NF_eval).subs(W22NF_eval)
Cod2 = simplify(firstordereval(Cod2).subs(extraparvals))
try:
    if abs(Cod2)>tol:
        print('You did not provide a Codimension-two point')
except:
    print('You did not provide a codimension-two bifurcation point. You have to solve\n')
    print(str(Cod2) + ' = 0\n')
    
equation1 = simplify(equation1.subs(extraparvals))
equation2 = simplify(equation2.subs(extraparvals))
try:
    equation3 = simplify(equation3.subs(extraparvals))
except:
    pass
equation4 = simplify(equation4.subs(extraparvals))
equation5 = simplify(equation5.subs(extraparvals))
    
try:
    if abs(equation1)>tol or abs(equation2)>tol or abs(equation3)>tol or abs(equation4)>tol or abs(equation5)>tol:
        print('Something went wrong. Review extra parameters.')
    else:
        print('The expansion to meet solvability conditions has been found.')
except:
    pass
    
SS_W04_10 = first_order(1, 0, W04NF)
SS_W04_01 = first_order(0, 1, W04NF)
SS_W03_20 = first_order(2, 0, W03NF)
SS_W03_11 = first_order(1, 1, W03NF)
SS_W03_02 = first_order(0, 2, W03NF)
SS_W02_30 = first_order(3, 0, W02NF)
SS_W02_21 = first_order(2, 1, W02NF)
SS_W02_12 = first_order(1, 2, W02NF)
SS_W02_03 = first_order(0, 3, W02NF)
DS_phiW14_00 = second_order(0, 0, phiNF, W14NF)
DS_W12W13_00 = second_order(0, 0, W12NF, W13NF)
DS_phiW13_10 = second_order(1, 0, phiNF, W13NF)
DS_phiW13_01 = second_order(0, 1, phiNF, W13NF)
DS_W12W12_10 = second_order(1, 0, W12NF, W12NF)
DS_W12W12_01 = second_order(0, 1, W12NF, W12NF)
DS_phiW12_20 = second_order(2, 0, phiNF, W12NF)
DS_phiW12_11 = second_order(1, 1, phiNF, W12NF)
DS_phiW12_02 = second_order(0, 2, phiNF, W12NF)
DS_phiphi_30 = second_order(3, 0, phiNF, phiNF)
DS_phiphi_21 = second_order(2, 1, phiNF, phiNF)
DS_phiphi_12 = second_order(1, 2, phiNF, phiNF)
DS_phiphi_03 = second_order(0, 3, phiNF, phiNF)
    
negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W04_10), Mul(bNF[1], SS_W04_01),
                              Mul(aNF[2], SS_W03_10), Mul(bNF[2], SS_W03_01),
                              Mul(aNF[3], SS_W02_10), Mul(bNF[3], SS_W02_01),
                              Mul(Pow(aNF[1], 2), SS_W03_20), Mul(aNF[1], bNF[1], SS_W03_11),
                              Mul(Pow(bNF[1], 2), SS_W03_02), Mul(2, aNF[1], aNF[2], SS_W02_20),
                              Mul(Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), SS_W02_11),
                              Mul(2, bNF[1], bNF[2], SS_W02_02), Mul(Pow(aNF[1], 3), SS_W02_30),
                              Mul(Pow(aNF[1], 2), bNF[1], SS_W02_21),
                              Mul(aNF[1], Pow(bNF[1], 2), SS_W02_12),
                              Mul(Pow(bNF[1], 3), SS_W02_03),
                              Mul(2, DS_phiW14_00), Mul(2, DS_W12W13_00),
                              Mul(2, aNF[1], DS_phiW13_10), Mul(2, bNF[1], DS_phiW13_01),
                              Mul(aNF[1], DS_W12W12_10), Mul(bNF[1], DS_W12W12_01),
                              Mul(2, aNF[2], DS_phiW12_10), Mul(2, bNF[2], DS_phiW12_01),
                              Mul(aNF[3], DS_phiphi_10), Mul(bNF[3], DS_phiphi_01),
                              Mul(2, Pow(aNF[1], 2), DS_phiW12_20),
                              Mul(2, aNF[1], bNF[1], DS_phiW12_11),
                              Mul(2, Pow(bNF[1], 2), DS_phiW12_02),
                              Mul(2, aNF[1], aNF[2], DS_phiphi_20),
                              Mul(Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), DS_phiphi_11),
                              Mul(2, bNF[1], bNF[2], DS_phiphi_02), Mul(Pow(aNF[1], 3), DS_phiphi_30),
                              Mul(Pow(aNF[1], 2), bNF[1], DS_phiphi_21),
                              Mul(aNF[1], Pow(bNF[1], 2), DS_phiphi_12),
                              Mul(Pow(bNF[1], 3), DS_phiphi_03)).subs(extraparvals)

W05NF = linearsolver(W05NF, negativeRHS, coefmat0)

SS_W024_10 = first_order(1, 0, W024NF)
SS_W024_01 = first_order(0, 1, W024NF)
DS_phiW124_00 = second_order(0, 0, phiNF, W124NF)
DS_W02W03_00 = second_order(0, 0, W02NF, W03NF)
DS_W12W123_00 = second_order(0, 0, W12NF, W123NF)
DS_W22W23_00 = second_order(0, 0, W22NF, W23NF)
DS_phiW123_10 = second_order(1, 0, phiNF, W123NF)
DS_phiW123_01 = second_order(0, 1, phiNF, W123NF)
DS_W02W02_10 = second_order(1, 0, W02NF, W02NF)
DS_W02W02_01 = second_order(0, 1, W02NF, W02NF)
DS_W22W22_10 = second_order(1, 0, W22NF, W22NF)
DS_W22W22_01 = second_order(0, 1, W22NF, W22NF)
TS_phiphiW03_00 = third_order(0, 0, phiNF, phiNF, W03NF)
TS_phiphiW23_00 = third_order(0, 0, phiNF, phiNF, W23NF)
TS_phiW12W02_00 = third_order(0, 0, phiNF, W12NF, W02NF)
TS_phiW12W22_00 = third_order(0, 0, phiNF, W12NF, W22NF)
TS_phiphiW02_10 = third_order(1, 0, phiNF, phiNF, W02NF)
TS_phiphiW22_10 = third_order(1, 0, phiNF, phiNF, W22NF)
TS_phiphiW02_01 = third_order(0, 1, phiNF, phiNF, W02NF)
TS_phiphiW22_01 = third_order(0, 1, phiNF, phiNF, W22NF)
Q4S_phiphiphiW12_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W12NF)
Q4S_phiphiphiphi_10 = fourth_order(1, 0, phiNF, phiNF, phiNF, phiNF)
Q4S_phiphiphiphi_01 = fourth_order(0, 1, phiNF, phiNF, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W024_10), Mul(bNF[1], SS_W024_01),
                              Mul(2, DS_phiW124_00), Mul(4, DS_W02W03_00),
                              Mul(2, DS_W12W123_00), Mul(2, DS_W22W23_00),
                              Mul(2, aNF[1], DS_phiW123_10), Mul(2, bNF[1], DS_phiW123_01),
                              Mul(2, aNF[1], DS_W02W02_10), Mul(2, bNF[1], DS_W02W02_01),
                              Mul(aNF[1], DS_W22W22_10), Mul(bNF[1], DS_W22W22_01),
                              Mul(6, TS_phiphiW03_00), Mul(3, TS_phiphiW23_00),
                              Mul(12, TS_phiW12W02_00), Mul(6, TS_phiW12W22_00),
                              Mul(6, aNF[1], TS_phiphiW02_10), Mul(3, aNF[1], TS_phiphiW22_10),
                              Mul(6, bNF[1], TS_phiphiW02_01), Mul(3, bNF[1], TS_phiphiW22_01),
                              Mul(12, Q4S_phiphiphiW12_00), Mul(3, aNF[1], Q4S_phiphiphiphi_10),
                              Mul(3, bNF[1], Q4S_phiphiphiphi_01)).subs(extraparvals)

W025NF = linearsolver(W025NF, negativeRHS, coefmat0)

SS_W034_10 = first_order(1, 0, W034NF)
SS_W034_01 = first_order(0, 1, W034NF)
DS_phiW134_00 = second_order(0, 0, phiNF, W134NF)
DS_W12W133_00 = second_order(0, 0, W12NF, W133NF)
DS_phiW133_10 = second_order(1, 0, phiNF, W133NF)
DS_phiW133_01 = second_order(0, 1, phiNF, W133NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W034_10), Mul(bNF[1], SS_W034_01),
                              Mul(2, DS_phiW134_00), Mul(2, DS_W12W133_00),
                              Mul(2, aNF[1], DS_phiW133_10),
                              Mul(2, bNF[1], DS_phiW133_01)).subs(extraparvals)

W035NF = linearsolver(W035NF, negativeRHS, coefmat0)
    
negativeRHS.actualcoord = Add(Mul(- 2, sqrt(muNF), diffmatrix, W133NF.dummy),
                              Mul(diffmatrix, phiNF.dummy))

alpha1 = psiNF.dummy.dot(negativeRHS.actualcoord)

negativeRHS.actualcoord = Add(negativeRHS.actualcoord, - Mul(alpha1, psiNF.dummy,
                                                             Pow(psiNF.dummy.dot(psiNF.dummy), - 1)))

W15NF = critical_linearsolver(W15NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

SS_W134_10 = first_order(1, 0, W134NF)
SS_W134_01 = first_order(0, 1, W134NF)
SS_W133_20 = first_order(2, 0, W133NF)
SS_W133_11 = first_order(1, 1, W133NF)
SS_W133_02 = first_order(0, 2, W133NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W134_10), Mul(bNF[1], SS_W134_01),
                              Mul(aNF[2], SS_W133_10), Mul(bNF[2], SS_W133_01),
                              Mul(Pow(aNF[1], 2), SS_W133_20), Mul(aNF[1], bNF[1], SS_W133_11),
                              Mul(Pow(bNF[1], 2), SS_W133_02),
                              Mul(2, sqrt(muNF), diffmatrix, W13NF.dummy)).subs(extraparvals)

alpha2 = psiNF.dummy.dot(negativeRHS.actualcoord)

negativeRHS.actualcoord = Add(negativeRHS.actualcoord, - Mul(alpha2, psiNF.dummy,
                                                             Pow(psiNF.dummy.dot(psiNF.dummy), - 1)))

W125NF = critical_linearsolver(W125NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

DS_phiW034_00 = second_order(0, 0, phiNF, W034NF)
DS_phiW234_00 = second_order(0, 0, phiNF, W234NF)
DS_W02W133_00 = second_order(0, 0, W02NF, W133NF)
TS_phiphiW133_00 = third_order(0, 0, phiNF, phiNF, W133NF)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW034_00), Mul(2, DS_phiW234_00),
                              Mul(4, DS_W02W133_00), Mul(6, TS_phiphiW133_00),
                              Mul(4, sqrt(muNF), diffmatrix, W123NF.dummy)).subs(extraparvals)

alpha3 = psiNF.dummy.dot(negativeRHS.actualcoord)

negativeRHS.actualcoord = Add(negativeRHS.actualcoord, - Mul(alpha3, psiNF.dummy,
                                                             Pow(psiNF.dummy.dot(psiNF.dummy), - 1)))

W135NF = critical_linearsolver(W135NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

DS_W22W133_00 = second_order(0, 0, W22NF, W133NF)

negativeRHS.actualcoord =  Add(Mul(- 2, DS_phiW034_00), Mul(- 2, DS_W22W133_00),
                               Mul(- 3, TS_phiphiW133_00),
                               Mul(2, sqrt(muNF), diffmatrix, W123NF.dummy)).subs(extraparvals)

alpha4 = psiNF.dummy.dot(negativeRHS.actualcoord)

negativeRHS.actualcoord = Add(negativeRHS.actualcoord, - Mul(alpha4, psiNF.dummy,
                                                             Pow(psiNF.dummy.dot(psiNF.dummy), - 1)))

W145NF = critical_linearsolver(W145NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

SS_W14_10 = first_order(1, 0, W14NF)
SS_W14_01 = first_order(0, 1, W14NF)
SS_W13_20 = first_order(2, 0, W13NF)
SS_W13_11 = first_order(1, 1, W13NF)
SS_W13_02 = first_order(0, 2, W13NF)
SS_W12_30 = first_order(3, 0, W12NF)
SS_W12_21 = first_order(2, 1, W12NF)
SS_W12_12 = first_order(1, 2, W12NF)
SS_W12_03 = first_order(0, 3, W12NF)
SS_phi_40 = first_order(4, 0, phiNF)
SS_phi_31 = first_order(3, 1, phiNF)
SS_phi_22 = first_order(2, 2, phiNF)
SS_phi_13 = first_order(1, 3, phiNF)
SS_phi_04 = first_order(0, 4, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W14_10), Mul(bNF[1], SS_W14_01),
                              Mul(aNF[2], SS_W13_10), Mul(bNF[2], SS_W13_01),
                              Mul(aNF[3], SS_W12_10), Mul(bNF[3], SS_W12_01),
                              Mul(aNF[4], SS_phi_10), Mul(bNF[4], SS_phi_01),                              
                              Mul(Pow(aNF[1], 2), SS_W13_20), Mul(aNF[1], bNF[1], SS_W13_11),
                              Mul(Pow(bNF[1], 2), SS_W13_02), Mul(2, aNF[1], aNF[2], SS_W12_20),
                              Mul(Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), SS_W12_11),
                              Mul(2, bNF[1], bNF[2], SS_W12_02), Mul(2, aNF[1], aNF[3], SS_phi_20),
                              Mul(Add(Mul(aNF[1], bNF[3]), Mul(aNF[3], bNF[1])), SS_phi_11),
                              Mul(2, bNF[1], bNF[3], SS_phi_02), Mul(Pow(aNF[2], 2), SS_phi_20),
                              Mul(aNF[2], bNF[2], SS_phi_11), Mul(Pow(bNF[2], 2), SS_phi_02),
                              Mul(Pow(aNF[1], 3), SS_W12_30), Mul(Pow(aNF[1], 2), bNF[1], SS_W12_21),
                              Mul(aNF[1], Pow(bNF[1], 2), SS_W12_12), Mul(Pow(bNF[1], 3), SS_W12_03),
                              Mul(3, Pow(aNF[1], 2), aNF[2], SS_phi_30),
                              Mul(2, aNF[1], aNF[2], bNF[1], SS_phi_21),
                              Mul(aNF[2], Pow(bNF[1], 2), SS_phi_12),
                              Mul(Pow(aNF[1], 2), bNF[2], SS_phi_21),
                              Mul(2, aNF[1], bNF[1], bNF[2], SS_phi_12),
                              Mul(3, Pow(bNF[1], 2), bNF[2], SS_phi_03),
                              Mul(Pow(aNF[1], 4), SS_phi_40), Mul(Pow(aNF[1], 3), bNF[1], SS_phi_31),
                              Mul(Pow(aNF[1], 2), Pow(bNF[1], 2), SS_phi_22),
                              Mul(aNF[1], Pow(bNF[1], 3), SS_phi_13),
                              Mul(Pow(bNF[1], 4), SS_phi_04)).subs(extraparvals)

alpha5 = psiNF.dummy.dot(negativeRHS.actualcoord)

negativeRHS.actualcoord = Add(negativeRHS.actualcoord, - Mul(alpha5, psiNF.dummy,
                                                             Pow(psiNF.dummy.dot(psiNF.dummy), - 1)))

W155NF = critical_linearsolver(W155NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

SS_W124_10 = first_order(1, 0, W124NF)
SS_W124_01 = first_order(0, 1, W124NF)
SS_W123_20 = first_order(2, 0, W123NF)
SS_W123_11 = first_order(1, 1, W123NF)
SS_W123_02 = first_order(0, 2, W123NF)
DS_phiW04_00 = second_order(0, 0, phiNF, W04NF)
DS_phiW24_00 = second_order(0, 0, phiNF, W24NF)
DS_W12W03_00 = second_order(0, 0, W12NF, W03NF)
DS_W12W23_00 = second_order(0, 0, W12NF, W23NF)
DS_W13W02_00 = second_order(0, 0, W13NF, W02NF)
DS_W13W22_00 = second_order(0, 0, W13NF, W22NF)
DS_phiW03_10 = second_order(1, 0, phiNF, W03NF)
DS_phiW23_10 = second_order(1, 0, phiNF, W23NF)
DS_phiW03_01 = second_order(0, 1, phiNF, W03NF)
DS_phiW23_01 = second_order(0, 1, phiNF, W23NF)
DS_W12W02_10 = second_order(1, 0, W12NF, W02NF)
DS_W12W22_10 = second_order(1, 0, W12NF, W22NF)
DS_W12W02_01 = second_order(0, 1, W12NF, W02NF)
DS_W12W22_01 = second_order(0, 1, W12NF, W22NF)
DS_phiW02_20 = second_order(2, 0, phiNF, W02NF)
DS_phiW22_20 = second_order(2, 0, phiNF, W22NF)
DS_phiW02_11 = second_order(1, 1, phiNF, W02NF)
DS_phiW22_11 = second_order(1, 1, phiNF, W22NF)
DS_phiW02_02 = second_order(0, 2, phiNF, W02NF)
DS_phiW22_02 = second_order(0, 2, phiNF, W22NF)
TS_phiphiW13_00 = third_order(0, 0, phiNF, phiNF, W13NF)
TS_phiW12W12_00 = third_order(0, 0, phiNF, W12NF, W12NF)
TS_phiphiW12_10 = third_order(1, 0, phiNF, phiNF, W12NF)
TS_phiphiW12_01 = third_order(0, 1, phiNF, phiNF, W12NF)
TS_phiphiphi_20 = third_order(2, 0, phiNF, phiNF, phiNF)
TS_phiphiphi_11 = third_order(1, 1, phiNF, phiNF, phiNF)
TS_phiphiphi_02 = third_order(0, 2, phiNF, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W124_10), Mul(bNF[1], SS_W124_01),
                              Mul(aNF[2], SS_W123_10), Mul(bNF[2], SS_W123_01),
                              Mul(Pow(aNF[1], 2), SS_W123_20), Mul(aNF[1], bNF[1], SS_W123_11),
                              Mul(Pow(bNF[1], 2), SS_W123_02), Mul(4, DS_phiW04_00),
                              Mul(2, DS_phiW24_00), Mul(4, DS_W12W03_00),
                              Mul(2, DS_W12W23_00), Mul(4, DS_W13W02_00),
                              Mul(2, DS_W13W22_00), Mul(4, aNF[1], DS_phiW03_10),
                              Mul(2, aNF[1], DS_phiW23_10), Mul(4, bNF[1], DS_phiW03_01),
                              Mul(2, bNF[1], DS_phiW23_01), Mul(4, aNF[1], DS_W12W02_10),
                              Mul(2, aNF[1], DS_W12W22_10), Mul(4, bNF[1], DS_W12W02_01),
                              Mul(2, bNF[1], DS_W12W22_01), Mul(4, aNF[2], DS_phiW02_10),
                              Mul(2, aNF[2], DS_phiW22_10), Mul(4, bNF[2], DS_phiW02_01),
                              Mul(2, bNF[2], DS_phiW22_01), Mul(4, Pow(aNF[1], 2), DS_phiW02_20),
                              Mul(2, Pow(aNF[1], 2), DS_phiW22_20), Mul(4, aNF[1], bNF[1], DS_phiW02_11),
                              Mul(2, aNF[1], bNF[1], DS_phiW22_11), Mul(4, Pow(bNF[1], 2), DS_phiW02_02),
                              Mul(2, Pow(bNF[1], 2), DS_phiW22_02), Mul(9, TS_phiphiW13_00),
                              Mul(9, TS_phiW12W12_00), Mul(9, aNF[1], TS_phiphiW12_10),
                              Mul(9, bNF[1], TS_phiphiW12_01), Mul(3, aNF[2], TS_phiphiphi_10),
                              Mul(3, bNF[2], TS_phiphiphi_01), Mul(3, Pow(aNF[1], 2), TS_phiphiphi_20),
                              Mul(3, aNF[1], bNF[1], TS_phiphiphi_11),
                              Mul(3, Pow(bNF[1], 2), TS_phiphiphi_02)).subs(extraparvals)

alpha6 = psiNF.dummy.dot(negativeRHS.actualcoord)

negativeRHS.actualcoord = Add(negativeRHS.actualcoord, - Mul(alpha6, psiNF.dummy,
                                                             Pow(psiNF.dummy.dot(psiNF.dummy), - 1)))

W165NF = critical_linearsolver(W165NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

DS_phiW024_00 = second_order(0, 0, phiNF, W024NF)
DS_phiW224_00 = second_order(0, 0, phiNF, W224NF)
DS_W22W33_00 = second_order(0, 0, W22NF, W33NF)
DS_W123W02_00 = second_order(0, 0, W123NF, W02NF)
DS_W123W22_00 = second_order(0, 0, W123NF, W22NF)
TS_phiphiW123_00 = third_order(0, 0, phiNF, phiNF, W123NF)
TS_phiphiW33_00 = third_order(0, 0, phiNF, phiNF, W33NF)
TS_phiW02W02_00 = third_order(0, 0, phiNF, W02NF, W02NF)
TS_phiW22W02_00 = third_order(0, 0, phiNF, W22NF, W02NF)
TS_phiW22W22_00 = third_order(0, 0, phiNF, W22NF, W22NF)
Q4S_phiphiphiW02_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W02NF)
Q4S_phiphiphiW22_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W22NF)
Q5S_phiphiphiphiphi_00 = fifth_order(0, 0, phiNF, phiNF, phiNF, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(4, DS_phiW024_00), Mul(2, DS_phiW224_00),
                              Mul(2, DS_W22W33_00), Mul(4, DS_W123W02_00),
                              Mul(2, DS_W123W22_00), Mul(9, TS_phiphiW123_00),
                              Mul(3, TS_phiphiW33_00), Mul(12, TS_phiW02W02_00),
                              Mul(12, TS_phiW22W02_00), Mul(6, TS_phiW22W22_00),
                              Mul(24, Q4S_phiphiphiW02_00), Mul(16, Q4S_phiphiphiW22_00),
                              Mul(10, Q5S_phiphiphiphiphi_00)).subs(extraparvals)

alpha7 = psiNF.dummy.dot(negativeRHS.actualcoord)

negativeRHS.actualcoord = Add(negativeRHS.actualcoord, - Mul(alpha7, psiNF.dummy,
                                                             Pow(psiNF.dummy.dot(psiNF.dummy), - 1)))

W175NF = critical_linearsolver(W175NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

SS_W24_10 = first_order(1, 0, W24NF)
SS_W24_01 = first_order(0, 1, W24NF)
SS_W23_20 = first_order(2, 0, W23NF)
SS_W23_11 = first_order(1, 1, W23NF)
SS_W23_02 = first_order(0, 2, W23NF)
SS_W22_30 = first_order(3, 0, W22NF)
SS_W22_21 = first_order(2, 1, W22NF)
SS_W22_12 = first_order(1, 2, W22NF)
SS_W22_03 = first_order(0, 3, W22NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W24_10), Mul(bNF[1], SS_W24_01),
                              Mul(aNF[2], SS_W23_10), Mul(bNF[2], SS_W23_01),
                              Mul(aNF[3], SS_W22_10), Mul(bNF[3], SS_W22_01),
                              Mul(Pow(aNF[1], 2), SS_W23_20), Mul(aNF[1], bNF[1], SS_W23_11),
                              Mul(Pow(bNF[1], 2), SS_W23_02), Mul(2, aNF[1], aNF[2], SS_W22_20),
                              Mul(Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), SS_W22_11),
                              Mul(2, bNF[1], bNF[2], SS_W22_02), Mul(Pow(aNF[1], 3), SS_W22_30),
                              Mul(Pow(aNF[1], 2), bNF[1], SS_W22_21),
                              Mul(aNF[1], Pow(bNF[1], 2), SS_W22_12),
                              Mul(Pow(bNF[1], 3), SS_W22_03), Mul(2, DS_phiW14_00),
                              Mul(2, DS_W12W13_00), Mul(2, aNF[1], DS_phiW13_10),
                              Mul(2, bNF[1], DS_phiW13_01), Mul(aNF[1], DS_W12W12_10),
                              Mul(bNF[1], DS_W12W12_01), Mul(2, aNF[2], DS_phiW12_10),
                              Mul(2, bNF[2], DS_phiW12_01), Mul(aNF[3], DS_phiphi_10),
                              Mul(bNF[3], DS_phiphi_01), Mul(2, Pow(aNF[1], 2), DS_phiW12_20),
                              Mul(2, aNF[1], bNF[1], DS_phiW12_11),
                              Mul(2, Pow(bNF[1], 2), DS_phiW12_02),
                              Mul(2, aNF[1], aNF[2], DS_phiphi_20),
                              Mul(Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), DS_phiphi_11),
                              Mul(2, bNF[1], bNF[2], DS_phiphi_02),
                              Mul(Pow(aNF[1], 3), DS_phiphi_30), Mul(Pow(aNF[1], 2), bNF[1], DS_phiphi_21),
                              Mul(aNF[1], Pow(bNF[1], 2), DS_phiphi_12),
                              Mul(Pow(bNF[1], 3), DS_phiphi_03)).subs(extraparvals)

W25NF = linearsolver(W25NF, negativeRHS, coefmat2)

SS_W224_10 = first_order(1, 0, W224NF)
SS_W224_01 = first_order(0, 1, W224NF)
DS_phiW34_00 = second_order(0, 0, phiNF, W34NF)
DS_W02W23_00 = second_order(0, 0, W02NF, W23NF)
DS_W12W33_00 = second_order(0, 0, W12NF, W33NF)
DS_W22W03_00 = second_order(0, 0, W22NF, W03NF)
DS_phiW33_10 = second_order(1, 0, phiNF, W33NF)
DS_phiW33_01 = second_order(0, 1, phiNF, W33NF)
DS_W02W22_10 = second_order(1, 0, W02NF, W22NF)
DS_W02W22_01 = second_order(0, 1, W02NF, W22NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W224_10), Mul(bNF[1], SS_W224_01),
                              Mul(2, DS_phiW124_00), Mul(2, DS_phiW34_00),
                              Mul(4, DS_W02W23_00), Mul(2, DS_W12W123_00),
                              Mul(2, DS_W12W33_00), Mul(4, DS_W22W03_00),
                              Mul(2, aNF[1], DS_phiW123_10), Mul(2, aNF[1], DS_phiW33_10),
                              Mul(2, bNF[1], DS_phiW123_01), Mul(2, bNF[1], DS_phiW33_01),
                              Mul(4, aNF[1], DS_W02W22_10), Mul(4, bNF[1], DS_W02W22_01),
                              Mul(6, TS_phiphiW03_00), Mul(6, TS_phiphiW23_00),
                              Mul(12, TS_phiW12W02_00), Mul(12, TS_phiW12W22_00),
                              Mul(6, aNF[1], TS_phiphiW02_10), Mul(6, aNF[1], TS_phiphiW22_10),
                              Mul(6, bNF[1], TS_phiphiW02_01), Mul(6, bNF[1], TS_phiphiW22_01),
                              Mul(16, Q4S_phiphiphiW12_00), Mul(4, aNF[1], Q4S_phiphiphiphi_10),
                              Mul(4, bNF[1], Q4S_phiphiphiphi_01)).subs(extraparvals)

W225NF = linearsolver(W225NF, negativeRHS, coefmat2)

SS_W234_10 = first_order(1, 0, W234NF)
SS_W234_01 = first_order(0, 1, W234NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W234_10), Mul(bNF[1], SS_W234_01),
                              Mul(2, DS_phiW134_00), Mul(2, DS_W12W133_00),
                              Mul(2, aNF[1], DS_phiW133_10), Mul(2, bNF[1], DS_phiW133_01),
                              Mul(8, sqrt(muNF), diffmatrix, W23NF.dummy)).subs(extraparvals)

W235NF = linearsolver(W235NF, negativeRHS, coefmat2)

SS_W34_10 = first_order(1, 0, W34NF)
SS_W34_01 = first_order(0, 1, W34NF)
SS_W33_20 = first_order(2, 0, W33NF)
SS_W33_11 = first_order(1, 1, W33NF)
SS_W33_02 = first_order(0, 2, W33NF)
DS_W22W13_00 = second_order(0, 0, W22NF, W13NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W34_10), Mul(bNF[1], SS_W34_01),
                              Mul(aNF[2], SS_W33_10), Mul(bNF[2], SS_W33_01),
                              Mul(Pow(aNF[1], 2), SS_W33_20), Mul(aNF[1], bNF[1], SS_W33_11),
                              Mul(Pow(bNF[1], 2), SS_W33_02), Mul(2, DS_phiW24_00),
                              Mul(2, DS_W12W23_00), Mul(2, DS_W22W13_00),
                              Mul(2, aNF[1], DS_phiW23_10), Mul(2, bNF[1], DS_phiW23_01),
                              Mul(2, aNF[1], DS_W12W22_10), Mul(2, bNF[1], DS_W12W22_01),
                              Mul(2, aNF[2], DS_phiW22_10), Mul(2, bNF[2], DS_phiW22_01),
                              Mul(2, Pow(aNF[1], 2), DS_phiW22_20),
                              Mul(2, aNF[1], bNF[1], DS_phiW22_11),
                              Mul(2, Pow(bNF[1], 2), DS_phiW22_02), Mul(3, TS_phiphiW13_00),
                              Mul(3, TS_phiW12W12_00), Mul(3, aNF[1], TS_phiphiW12_10),
                              Mul(3, bNF[1], TS_phiphiW12_01), Mul(aNF[2], TS_phiphiphi_10),
                              Mul(bNF[2], TS_phiphiphi_01), Mul(Pow(aNF[1], 2), TS_phiphiphi_20),
                              Mul(aNF[1], bNF[1], TS_phiphiphi_11),
                              Mul(Pow(bNF[1], 2), TS_phiphiphi_02)).subs(extraparvals)

W35NF = linearsolver(W35NF, negativeRHS, coefmat3)

DS_phiW44_00 = second_order(0, 0, phiNF, W44NF)
DS_W02W33_00 = second_order(0, 0, W02NF, W33NF)
DS_W22W123_00 = second_order(0, 0, W22NF, W123NF)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW224_00), Mul(2, DS_phiW44_00),
                              Mul(4, DS_W02W33_00), Mul(2, DS_W22W123_00),
                              Mul(3, TS_phiphiW123_00), Mul(6, TS_phiphiW33_00),
                              Mul(12, TS_phiW22W02_00), Mul(3, TS_phiW22W22_00),
                              Mul(8, Q4S_phiphiphiW02_00), Mul(12, Q4S_phiphiphiW22_00),
                              Mul(5, Q5S_phiphiphiphiphi_00)).subs(extraparvals)

W325NF = linearsolver(W325NF, negativeRHS, coefmat3)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW234_00), Mul(2, DS_W22W133_00),
                              Mul(3, TS_phiphiW133_00),
                              Mul(18, sqrt(muNF), diffmatrix, W33NF.dummy)).subs(extraparvals)

W335NF = linearsolver(W335NF, negativeRHS, coefmat3)

W05NF = evaluation(W05NF)
W025NF = evaluation(W025NF)
W035NF = evaluation(W035NF)
W15NF = evaluation(W15NF)
W125NF = evaluation(W125NF)
W135NF = evaluation(W135NF)
W145NF = evaluation(W145NF)
W155NF = evaluation(W155NF)
W165NF = evaluation(W165NF)
W175NF = evaluation(W175NF)
W25NF = evaluation(W25NF)
W225NF = evaluation(W225NF)
W235NF = evaluation(W235NF)
W35NF = evaluation(W35NF)
W325NF = evaluation(W325NF)
W335NF = evaluation(W335NF)

W05NF_eval = evaluation_dict(W05NF)
W025NF_eval = evaluation_dict(W025NF)
W035NF_eval = evaluation_dict(W035NF)
W15NF_eval = evaluation_dict(W15NF)
W125NF_eval = evaluation_dict(W125NF)
W135NF_eval = evaluation_dict(W135NF)
W145NF_eval = evaluation_dict(W145NF)
W155NF_eval = evaluation_dict(W155NF)
W165NF_eval = evaluation_dict(W165NF)
W175NF_eval = evaluation_dict(W175NF)
W25NF_eval = evaluation_dict(W25NF)
W225NF_eval = evaluation_dict(W225NF)
W235NF_eval = evaluation_dict(W235NF)
W35NF_eval = evaluation_dict(W35NF)
W325NF_eval = evaluation_dict(W325NF)
W335NF_eval = evaluation_dict(W335NF)

print('Fifth order ready')

SS_W05_10 = first_order(1, 0, W05NF)
SS_W05_01 = first_order(0, 1, W05NF)
SS_W04_20 = first_order(2, 0, W04NF)
SS_W04_11 = first_order(1, 1, W04NF)
SS_W04_02 = first_order(0, 2, W04NF)
SS_W03_30 = first_order(3, 0, W03NF)
SS_W03_21 = first_order(2, 1, W03NF)
SS_W03_12 = first_order(1, 2, W03NF)
SS_W03_03 = first_order(0, 3, W03NF)
SS_W02_40 = first_order(4, 0, W02NF)
SS_W02_31 = first_order(3, 1, W02NF)
SS_W02_22 = first_order(2, 2, W02NF)
SS_W02_13 = first_order(1, 3, W02NF)
SS_W02_04 = first_order(0, 4, W02NF)
DS_phiW155_00 = second_order(0, 0, phiNF, W155NF)
DS_W12W14_00 = second_order(0, 0, W12NF, W14NF)
DS_W13W13_00 = second_order(0, 0, W13NF, W13NF)
DS_phiW14_10 = second_order(1, 0, phiNF, W14NF)
DS_phiW14_01 = second_order(0, 1, phiNF, W14NF)
DS_W12W13_10 = second_order(1, 0, W12NF, W13NF)
DS_W12W13_01 = second_order(0, 1, W12NF, W13NF)
DS_phiW13_20 = second_order(2, 0, phiNF, W13NF)
DS_phiW13_11 = second_order(1, 1, phiNF, W13NF)
DS_phiW13_02 = second_order(0, 2, phiNF, W13NF)
DS_W12W12_20 = second_order(2, 0, W12NF, W12NF)
DS_W12W12_11 = second_order(1, 1, W12NF, W12NF)
DS_W12W12_02 = second_order(0, 2, W12NF, W12NF)
DS_phiW12_30 = second_order(3, 0, phiNF, W12NF)
DS_phiW12_21 = second_order(2, 1, phiNF, W12NF)
DS_phiW12_12 = second_order(1, 2, phiNF, W12NF)
DS_phiW12_03 = second_order(0, 3, phiNF, W12NF)
DS_phiphi_40 = second_order(4, 0, phiNF, phiNF)
DS_phiphi_31 = second_order(3, 1, phiNF, phiNF)
DS_phiphi_22 = second_order(2, 2, phiNF, phiNF)
DS_phiphi_13 = second_order(1, 3, phiNF, phiNF)
DS_phiphi_04 = second_order(0, 4, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W05_10), Mul(bNF[1], SS_W05_01),
                              Mul(aNF[2], SS_W04_10), Mul(bNF[2], SS_W04_01),
                              Mul(aNF[3], SS_W03_10), Mul(bNF[3], SS_W03_01),
                              Mul(aNF[4], SS_W02_10), Mul(bNF[4], SS_W02_01),
                              Mul(Pow(aNF[1], 2), SS_W04_20), Mul(aNF[1], bNF[1], SS_W04_11),
                              Mul(Pow(bNF[1], 2), SS_W04_02), Mul(2, aNF[1], aNF[2], SS_W03_20),
                              Mul(Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), SS_W03_11),
                              Mul(2, bNF[1], bNF[2], SS_W03_02), Mul(2, aNF[1], aNF[3], SS_W02_20),
                              Mul(Add(Mul(aNF[1], bNF[3]), Mul(aNF[3], bNF[1])), SS_W02_11),
                              Mul(2, bNF[1], bNF[3], SS_W02_02), Mul(Pow(aNF[2], 2), SS_W02_20),
                              Mul(aNF[2], bNF[2], SS_W02_11), Mul(Pow(bNF[2], 2), SS_W02_02),
                              Mul(Pow(aNF[1], 3), SS_W03_30), Mul(Pow(aNF[1], 2), bNF[1], SS_W03_21),
                              Mul(aNF[1], Pow(bNF[1], 2), SS_W03_12), Mul(Pow(bNF[1], 3), SS_W03_03),
                              Mul(3, Pow(aNF[1], 2), aNF[2], SS_W02_30),
                              Mul(2, aNF[1], aNF[2], bNF[1], SS_W02_21),
                              Mul(aNF[2], Pow(bNF[1], 2), SS_W02_12),
                              Mul(Pow(aNF[1], 2), bNF[2], SS_W02_21),
                              Mul(2, aNF[1], bNF[1], bNF[2], SS_W02_12),
                              Mul(3, Pow(bNF[1], 2), bNF[2], SS_W02_03),
                              Mul(Pow(aNF[1], 4), SS_W02_40),
                              Mul(Pow(aNF[1], 3), bNF[1], SS_W02_31),
                              Mul(Pow(aNF[1], 2), Pow(bNF[1], 2), SS_W02_22),
                              Mul(aNF[1], Pow(bNF[1], 3), SS_W02_13), Mul(Pow(bNF[1], 4), SS_W02_04),
                              Mul(2, DS_phiW155_00), Mul(2, DS_W12W14_00),
                              DS_W13W13_00, Mul(2, aNF[1], DS_phiW14_10),
                              Mul(2, bNF[1], DS_phiW14_01), Mul(2, aNF[1], DS_W12W13_10),
                              Mul(2, bNF[1], DS_W12W13_01), Mul(2, aNF[2], DS_phiW13_10),
                              Mul(2, bNF[2], DS_phiW13_01), Mul(aNF[2], DS_W12W12_10),
                              Mul(bNF[2], DS_W12W12_01), Mul(2, aNF[3], DS_phiW12_10),
                              Mul(2, bNF[3], DS_phiW12_01), Mul(aNF[4], DS_phiphi_10),
                              Mul(bNF[4], DS_phiphi_01), Mul(2, Pow(aNF[1], 2), DS_phiW13_20),
                              Mul(2, aNF[1], bNF[1], DS_phiW13_11),
                              Mul(2, Pow(bNF[1], 2), DS_phiW13_02),
                              Mul(Pow(aNF[1], 2), DS_W12W12_20), Mul(aNF[1], bNF[1], DS_W12W12_11),
                              Mul(Pow(bNF[1], 2), DS_W12W12_02), Mul(4, aNF[1], aNF[2], DS_phiW12_20),
                              Mul(2, Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), DS_phiW12_11),
                              Mul(4, bNF[1], bNF[2], DS_phiW12_02),
                              Mul(2, aNF[1], aNF[3], DS_phiphi_20),
                              Mul(Add(Mul(aNF[1], bNF[3]), Mul(aNF[3], bNF[1])), DS_phiphi_11),
                              Mul(2, bNF[1], bNF[3], DS_phiphi_02),
                              Mul(Pow(aNF[2], 2), DS_phiphi_20),
                              Mul(aNF[2], bNF[2], DS_phiphi_11), Mul(Pow(bNF[2], 2), DS_phiphi_02),
                              Mul(2, Pow(aNF[1], 3), DS_phiW12_30),
                              Mul(2, Pow(aNF[1], 2), bNF[1], DS_phiW12_21),
                              Mul(2, aNF[1], Pow(bNF[1], 2), DS_phiW12_12),
                              Mul(2, Pow(bNF[1], 3), DS_phiW12_03),
                              Mul(3, Pow(aNF[1], 2), aNF[2], DS_phiphi_30),
                              Mul(2, aNF[1], aNF[2], bNF[1], DS_phiphi_21),
                              Mul(aNF[2], Pow(bNF[1], 2), DS_phiphi_12),
                              Mul(Pow(aNF[1], 2), bNF[2], DS_phiphi_21),
                              Mul(2, aNF[1], bNF[1], bNF[2], DS_phiphi_12),
                              Mul(3, Pow(bNF[1], 2), bNF[2], DS_phiphi_03),
                              Mul(Pow(aNF[1], 4), DS_phiphi_40),
                              Mul(Pow(aNF[1], 3), bNF[1], DS_phiphi_31),
                              Mul(Pow(aNF[1], 2), Pow(bNF[1], 2), DS_phiphi_22),
                              Mul(aNF[1], Pow(bNF[1], 3), DS_phiphi_13),
                              Mul(Pow(bNF[1], 4), DS_phiphi_04)).subs(extraparvals)

W06NF = linearsolver(W06NF, negativeRHS, coefmat0)

SS_W025_10 = first_order(1, 0, W025NF)
SS_W025_01 = first_order(0, 1, W025NF)
SS_W024_20 = first_order(2, 0, W024NF)
SS_W024_11 = first_order(1, 1, W024NF)
SS_W024_02 = first_order(0, 2, W024NF)
DS_phiW165_00 = second_order(0, 0, phiNF, W165NF)
DS_W02W04_00 = second_order(0, 0, W02NF, W04NF)
DS_W12W124_00 = second_order(0, 0, W12NF, W124NF)
DS_W22W24_00 = second_order(0, 0, W22NF, W24NF)
DS_W03W03_00 = second_order(0, 0, W03NF, W03NF)
DS_W13W123_00 = second_order(0, 0, W13NF, W123NF)
DS_W23W23_00 = second_order(0, 0, W23NF, W23NF)
DS_phiW124_10 = second_order(1, 0, phiNF, W124NF)
DS_phiW124_01 = second_order(0, 1, phiNF, W124NF)
DS_W02W03_10 = second_order(1, 0, W02NF, W03NF)
DS_W02W03_01 = second_order(0, 1, W02NF, W03NF)
DS_W12W123_10 = second_order(1, 0, W12NF, W123NF)
DS_W12W123_01 = second_order(0, 1, W12NF, W123NF)
DS_W22W23_10 = second_order(1, 0, W22NF, W23NF)
DS_W22W23_01 = second_order(0, 1, W22NF, W23NF)
DS_phiW123_20 = second_order(2, 0, phiNF, W123NF)
DS_phiW123_11 = second_order(1, 1, phiNF, W123NF)
DS_phiW123_02 = second_order(0, 2, phiNF, W123NF)
DS_W02W02_20 = second_order(2, 0, W02NF, W02NF)
DS_W02W02_11 = second_order(1, 1, W02NF, W02NF)
DS_W02W02_02 = second_order(0, 2, W02NF, W02NF)
DS_W22W22_20 = second_order(2, 0, W22NF, W22NF)
DS_W22W22_11 = second_order(1, 1, W22NF, W22NF)
DS_W22W22_02 = second_order(0, 2, W22NF, W22NF)
TS_phiphiW04_00 = third_order(0, 0, phiNF, phiNF, W04NF)
TS_phiphiW24_00 = third_order(0, 0, phiNF, phiNF, W24NF)
TS_phiW12W03_00 = third_order(0, 0, phiNF, W12NF, W03NF)
TS_phiW12W23_00 = third_order(0, 0, phiNF, W12NF, W23NF)
TS_phiW13W02_00 = third_order(0, 0, phiNF, W13NF, W02NF)
TS_phiW13W22_00 = third_order(0, 0, phiNF, W13NF, W22NF)
TS_W12W12W02_00 = third_order(0, 0, W12NF, W12NF, W02NF)
TS_W12W12W22_00 = third_order(0, 0, W12NF, W12NF, W22NF)
TS_phiphiW03_10 = third_order(1, 0, phiNF, phiNF, W03NF)
TS_phiphiW23_10 = third_order(1, 0, phiNF, phiNF, W23NF)
TS_phiphiW03_01 = third_order(0, 1, phiNF, phiNF, W03NF)
TS_phiphiW23_01 = third_order(0, 1, phiNF, phiNF, W23NF)
TS_phiW12W02_10 = third_order(1, 0, phiNF, W12NF, W02NF)
TS_phiW12W22_10 = third_order(1, 0, phiNF, W12NF, W22NF)
TS_phiW12W02_01 = third_order(0, 1, phiNF, W12NF, W02NF)
TS_phiW12W22_01 = third_order(0, 1, phiNF, W12NF, W22NF)
TS_phiphiW02_20 = third_order(2, 0, phiNF, phiNF, W02NF)
TS_phiphiW22_20 = third_order(2, 0, phiNF, phiNF, W22NF)
TS_phiphiW02_11 = third_order(1, 1, phiNF, phiNF, W02NF)
TS_phiphiW22_11 = third_order(1, 1, phiNF, phiNF, W22NF)
TS_phiphiW02_02 = third_order(0, 2, phiNF, phiNF, W02NF)
TS_phiphiW22_02 = third_order(0, 2, phiNF, phiNF, W22NF)
Q4S_phiphiphiW13_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W13NF)
Q4S_phiphiW12W12_00 = fourth_order(0, 0, phiNF, phiNF, W12NF, W12NF)
Q4S_phiphiphiW12_10 = fourth_order(1, 0, phiNF, phiNF, phiNF, W12NF)
Q4S_phiphiphiW12_01 = fourth_order(0, 1, phiNF, phiNF, phiNF, W12NF)
Q4S_phiphiphiphi_20 = fourth_order(2, 0, phiNF, phiNF, phiNF, phiNF)
Q4S_phiphiphiphi_11 = fourth_order(1, 1, phiNF, phiNF, phiNF, phiNF)
Q4S_phiphiphiphi_02 = fourth_order(0, 2, phiNF, phiNF, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W025_10), Mul(bNF[1], SS_W025_01),
                              Mul(aNF[2], SS_W024_10), Mul(bNF[2], SS_W024_01),
                              Mul(Pow(aNF[1], 2), SS_W024_20), Mul(aNF[1], bNF[1], SS_W024_11),
                              Mul(Pow(bNF[1], 2), SS_W024_02), Mul(2, DS_phiW165_00),
                              Mul(4, DS_W02W04_00), Mul(2, DS_W12W124_00), Mul(2, DS_W22W24_00),
                              Mul(2, DS_W03W03_00), Mul(2, DS_W13W123_00),
                              DS_W23W23_00, Mul(2, aNF[1], DS_phiW124_10),
                              Mul(2, bNF[1], DS_phiW124_01), Mul(4, aNF[1], DS_W02W03_10),
                              Mul(4, bNF[1], DS_W02W03_01), Mul(2, aNF[1], DS_W12W123_10),
                              Mul(2, bNF[1], DS_W12W123_01), Mul(2, aNF[1], DS_W22W23_10),
                              Mul(2, bNF[1], DS_W22W23_01), Mul(2, aNF[2], DS_phiW123_10),
                              Mul(2, bNF[2], DS_phiW123_01), Mul(2, aNF[2], DS_W02W02_10),
                              Mul(2, bNF[2], DS_W02W02_01), Mul(aNF[2], DS_W22W22_10),
                              Mul(bNF[2], DS_W22W22_01), Mul(2, Pow(aNF[1], 2), DS_phiW123_20),
                              Mul(2, aNF[1], bNF[1], DS_phiW123_11), Mul(2, Pow(bNF[1], 2), DS_phiW123_02),
                              Mul(2, Pow(aNF[1], 2), DS_W02W02_20), Mul(2, aNF[1], bNF[1], DS_W02W02_11),
                              Mul(2, Pow(bNF[1], 2), DS_W02W02_02), Mul(Pow(aNF[1], 2), DS_W22W22_20),
                              Mul(aNF[1], bNF[1], DS_W22W22_11), Mul(Pow(bNF[1], 2), DS_W22W22_02),
                              Mul(6, TS_phiphiW04_00), Mul(3, TS_phiphiW24_00),
                              Mul(12, TS_phiW12W03_00), Mul(6, TS_phiW12W23_00),
                              Mul(12, TS_phiW13W02_00), Mul(6, TS_phiW13W22_00),
                              Mul(6, TS_W12W12W02_00), Mul(3, TS_W12W12W22_00),
                              Mul(6, aNF[1], TS_phiphiW03_10), Mul(3, aNF[1], TS_phiphiW23_10),
                              Mul(6, bNF[1], TS_phiphiW03_01), Mul(3, bNF[1], TS_phiphiW23_01),
                              Mul(12, aNF[1], TS_phiW12W02_10), Mul(6, aNF[1], TS_phiW12W22_10),
                              Mul(12, bNF[1], TS_phiW12W02_01), Mul(6, bNF[1], TS_phiW12W22_01),
                              Mul(6, aNF[2], TS_phiphiW02_10), Mul(3, aNF[2], TS_phiphiW22_10),
                              Mul(6, bNF[2], TS_phiphiW02_01), Mul(3, bNF[2], TS_phiphiW22_01),
                              Mul(6, Pow(aNF[1], 2), TS_phiphiW02_20), Mul(3, Pow(aNF[1], 2), TS_phiphiW22_20),
                              Mul(6, aNF[1], bNF[1], TS_phiphiW02_11), Mul(3, aNF[1], bNF[1], TS_phiphiW22_11),
                              Mul(6, Pow(bNF[1], 2), TS_phiphiW02_02), Mul(3, Pow(bNF[1], 2), TS_phiphiW22_02),
                              Mul(12, Q4S_phiphiphiW13_00), Mul(18, Q4S_phiphiW12W12_00),
                              Mul(12, aNF[1], Q4S_phiphiphiW12_10), Mul(12, bNF[1], Q4S_phiphiphiW12_01),
                              Mul(3, aNF[2], Q4S_phiphiphiphi_10), Mul(3, bNF[2], Q4S_phiphiphiphi_01),
                              Mul(3, Pow(aNF[1], 2), Q4S_phiphiphiphi_20),
                              Mul(3, aNF[1], bNF[1], Q4S_phiphiphiphi_11),
                              Mul(3, Pow(bNF[1], 2), Q4S_phiphiphiphi_02)).subs(extraparvals)

W026NF = linearsolver(W026NF, negativeRHS, coefmat0)

DS_phiW175_00 = second_order(0, 0, phiNF, W175NF)
DS_W02W024_00 = second_order(0, 0, W02NF, W024NF)
DS_W22W224_00 = second_order(0, 0, W22NF, W224NF)
DS_W123W123_00 = second_order(0, 0, W123NF, W123NF)
DS_W33W33_00 = second_order(0, 0, W33NF, W33NF)
TS_phiphiW024_00 = third_order(0, 0, phiNF, phiNF, W024NF)
TS_phiphiW224_00 = third_order(0, 0, phiNF, phiNF, W224NF)
TS_phiW22W33_00 = third_order(0, 0, phiNF, W22NF, W33NF)
TS_phiW123W02_00 = third_order(0, 0, phiNF, W123NF, W02NF)
TS_phiW123W22_00 = third_order(0, 0, phiNF, W123NF, W22NF)
TS_W02W02W02_00 = third_order(0, 0, W02NF, W02NF, W02NF)
TS_W02W22W22_00 = third_order(0, 0, W02NF, W22NF, W22NF)
Q4S_phiphiphiW123_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W123NF)
Q4S_phiphiphiW33_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W33NF)
Q4S_phiphiW02W02_00 = fourth_order(0, 0, phiNF, phiNF, W02NF, W02NF)
Q4S_phiphiW22W02_00 = fourth_order(0, 0, phiNF, phiNF, W22NF, W02NF)
Q4S_phiphiW22W22_00 = fourth_order(0, 0, phiNF, phiNF, W22NF, W22NF)
Q5S_phiphiphiphiW02_00 = fifth_order(0, 0, phiNF, phiNF, phiNF, phiNF, W02NF)
Q5S_phiphiphiphiW22_00 = fifth_order(0, 0, phiNF, phiNF, phiNF, phiNF, W22NF)
S6S_phiphiphiphiphiphi_00 = sixth_order(0, 0, phiNF, phiNF, phiNF, phiNF, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW175_00), Mul(4, DS_W02W024_00),
                              Mul(2, DS_W22W224_00), DS_W123W123_00,
                              DS_W33W33_00, Mul(6, TS_phiphiW024_00),
                              Mul(3, TS_phiphiW224_00), Mul(6, TS_phiW22W33_00),
                              Mul(12, TS_phiW123W02_00), Mul(6, TS_phiW123W22_00),
                              Mul(4, TS_W02W02W02_00), Mul(6, TS_W02W22W22_00),
                              Mul(12, Q4S_phiphiphiW123_00), Mul(4, Q4S_phiphiphiW33_00),
                              Mul(24, Q4S_phiphiW02W02_00), Mul(24, Q4S_phiphiW22W02_00),
                              Mul(12, Q4S_phiphiW22W22_00), Mul(30, Q5S_phiphiphiphiW02_00),
                              Mul(20, Q5S_phiphiphiphiW22_00),
                              Mul(10, S6S_phiphiphiphiphiphi_00)).subs(extraparvals)

W036NF = linearsolver(W036NF, negativeRHS, coefmat0)

DS_W133W133_00 = second_order(0, 0, W133NF, W133NF)

negativeRHS.actualcoord = Add(DS_W133W133_00, Mul(2, diffmatrix, W02NF.dummy))

W046NF = linearsolver(W046NF, negativeRHS, coefmat0)

SS_W035_10 = first_order(1, 0, W035NF)
SS_W035_01 = first_order(0, 1, W035NF)
SS_W034_20 = first_order(2, 0, W034NF)
SS_W034_11 = first_order(1, 1, W034NF)
SS_W034_02 = first_order(0, 2, W034NF)
DS_phiW125_00 = second_order(0, 0, phiNF, W125NF)
DS_W12W134_00 = second_order(0, 0, W12NF, W134NF)
DS_W13W133_00 = second_order(0, 0, W13NF, W133NF)
DS_phiW134_10 = second_order(1, 0, phiNF, W134NF)
DS_phiW134_01 = second_order(0, 1, phiNF, W134NF)
DS_W12W133_10 = second_order(1, 0, W12NF, W133NF)
DS_W12W133_01 = second_order(0, 1, W12NF, W133NF)
DS_phiW133_20 = second_order(2, 0, phiNF, W133NF)
DS_phiW133_11 = second_order(1, 1, phiNF, W133NF)
DS_phiW133_02 = second_order(0, 2, phiNF, W133NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W035_10), Mul(bNF[1], SS_W035_01),
                              Mul(aNF[2], SS_W034_10), Mul(bNF[2], SS_W034_01),
                              Mul(Pow(aNF[1], 2), SS_W034_20), Mul(aNF[1], bNF[1], SS_W034_11),
                              Mul(Pow(bNF[1], 2), SS_W034_02), Mul(2, DS_phiW125_00),
                              Mul(2, DS_W12W134_00), Mul(2, DS_W13W133_00),
                              Mul(2, aNF[1], DS_phiW134_10), Mul(2, bNF[1], DS_phiW134_01),
                              Mul(2, aNF[1], DS_W12W133_10), Mul(2, bNF[1], DS_W12W133_01),
                              Mul(2, aNF[2], DS_phiW133_10), Mul(2, bNF[2], DS_phiW133_01),
                              Mul(2, Pow(aNF[1], 2), DS_phiW133_20), Mul(2, aNF[1], bNF[1], DS_phiW133_11),
                              Mul(2, Pow(bNF[1], 2), DS_phiW133_02)).subs(extraparvals)

W056NF = linearsolver(W056NF, negativeRHS, coefmat0)

DS_phiW135_00 = second_order(0, 0, phiNF, W135NF)
DS_phiW145_00 = second_order(0, 0, phiNF, W145NF)
DS_W02W034_00 = second_order(0, 0, W02NF, W034NF)
DS_W22W234_00 = second_order(0, 0, W22NF, W234NF)
DS_W123W133_00 = second_order(0, 0, W123NF, W133NF)
TS_phiphiW034_00 = third_order(0, 0, phiNF, phiNF, W034NF)
TS_phiphiW234_00 = third_order(0, 0, phiNF, phiNF, W234NF)
TS_phiW133W02_00 = third_order(0, 0, phiNF, W133NF, W02NF)
TS_phiW133W22_00 = third_order(0, 0, phiNF, W133NF, W22NF)
Q4S_phiphiphiW133_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W133NF)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW135_00), Mul(- 2, DS_phiW145_00),
                              Mul(4, DS_W02W034_00), Mul(2, DS_W22W234_00),
                              Mul(2, DS_W123W133_00), Mul(6, TS_phiphiW034_00),
                              Mul(3, TS_phiphiW234_00), Mul(12, TS_phiW133W02_00),
                              Mul(6, TS_phiW133W22_00), Mul(12, Q4S_phiphiphiW133_00)).subs(extraparvals)

W066NF = linearsolver(W066NF, negativeRHS, coefmat0)

DS_phiW15_00 = second_order(0, 0, phiNF, W15NF)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW15_00), Mul(2, diffmatrix, W02NF.dummy))

W076NF = linearsolver(W076NF, negativeRHS, coefmat0)

SS_W15_10 = first_order(1, 0, W15NF)
SS_W15_01 = first_order(0, 1, W15NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W15_10), Mul(bNF[1], SS_W15_01),
                              Mul(- 2, sqrt(muNF), diffmatrix, W134NF.dummy),
                              Mul(diffmatrix, W12NF.dummy)).subs(extraparvals)

alpha12 = psiNF.dummy.dot(negativeRHS.actualcoord)

negativeRHS.actualcoord = Add(negativeRHS.actualcoord, - Mul(alpha12, psiNF.dummy,
                                                             Pow(psiNF.dummy.dot(psiNF.dummy), - 1)))

W16NF = critical_linearsolver(W16NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

SS_W125_10 = first_order(1, 0, W125NF)
SS_W125_01 = first_order(0, 1, W125NF)
SS_W134_20 = first_order(2, 0, W134NF)
SS_W134_11 = first_order(1, 1, W134NF)
SS_W134_02 = first_order(0, 2, W134NF)
SS_W133_30 = first_order(3, 0, W133NF)
SS_W133_21 = first_order(2, 1, W133NF)
SS_W133_12 = first_order(1, 2, W133NF)
SS_W133_03 = first_order(0, 3, W133NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W125_10), Mul(bNF[1], SS_W125_01),
                              Mul(aNF[2], SS_W134_10), Mul(bNF[2], SS_W134_01),
                              Mul(aNF[3], SS_W133_10), Mul(bNF[3], SS_W133_01),
                              Mul(Pow(aNF[1], 2), SS_W134_20), Mul(aNF[1], bNF[1], SS_W134_11),
                              Mul(Pow(bNF[1], 2), SS_W134_02), Mul(2, aNF[1], aNF[2], SS_W133_20),
                              Mul(Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), SS_W133_11),
                              Mul(2, bNF[1], bNF[2], SS_W133_02), Mul(Pow(aNF[1], 3), SS_W133_30),
                              Mul(Pow(aNF[1], 2), bNF[1], SS_W133_21),
                              Mul(aNF[1], Pow(bNF[1], 2), SS_W133_12),
                              Mul(Pow(bNF[1], 3), SS_W133_03),
                              Mul(2, sqrt(muNF), diffmatrix, W14NF.dummy)).subs(extraparvals)

alpha22 = psiNF.dummy.dot(negativeRHS.actualcoord)

negativeRHS.actualcoord = Add(negativeRHS.actualcoord, - Mul(alpha22, psiNF.dummy,
                                                             Pow(psiNF.dummy.dot(psiNF.dummy), - 1)))

W126NF = critical_linearsolver(W126NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

SS_W135_10 = first_order(1, 0, W135NF)
SS_W135_01 = first_order(0, 1, W135NF)
DS_phiW035_00 = second_order(0, 0, phiNF, W035NF)
DS_phiW235_00 = second_order(0, 0, phiNF, W235NF)
DS_W02W134_00 = second_order(0, 0, W02NF, W134NF)
DS_W12W034_00 = second_order(0, 0, W12NF, W034NF)
DS_W12W234_00 = second_order(0, 0, W12NF, W234NF)
DS_W03W133_00 = second_order(0, 0, W03NF, W133NF)
DS_phiW034_10 = second_order(1, 0, phiNF, W034NF)
DS_phiW234_10 = second_order(1, 0, phiNF, W234NF)
DS_phiW034_01 = second_order(0, 1, phiNF, W034NF)
DS_phiW234_01 = second_order(0, 1, phiNF, W234NF)
DS_W02W133_10 = second_order(1, 0, W02NF, W133NF)
DS_W02W133_01 = second_order(0, 1, W02NF, W133NF)
TS_phiphiW134_00 = third_order(0, 0, phiNF, phiNF, W134NF)
TS_phiW12W133_00 = third_order(0, 0, phiNF, W12NF, W133NF)
TS_phiphiW133_10 = third_order(1, 0, phiNF, phiNF, W133NF)
TS_phiphiW133_01 = third_order(0, 1, phiNF, phiNF, W133NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W135_10), Mul(bNF[1], SS_W135_01),
                              Mul(2, DS_phiW035_00), Mul(2, DS_phiW235_00),
                              Mul(4, DS_W02W134_00), Mul(2, DS_W12W034_00),
                              Mul(2, DS_W12W234_00), Mul(4, DS_W03W133_00),
                              Mul(2, aNF[1], DS_phiW034_10), Mul(2, aNF[1], DS_phiW234_10),
                              Mul(2, bNF[1], DS_phiW034_01), Mul(2, bNF[1], DS_phiW234_01),
                              Mul(4, aNF[1], DS_W02W133_10), Mul(4, bNF[1], DS_W02W133_01),
                              Mul(6, TS_phiphiW134_00), Mul(12, TS_phiW12W133_00),
                              Mul(6, aNF[1], TS_phiphiW133_10), Mul(6, bNF[1], TS_phiphiW133_01),
                              Mul(4, sqrt(muNF), diffmatrix, W124NF.dummy)).subs(extraparvals)

alpha32 = psiNF.dummy.dot(negativeRHS.actualcoord)

negativeRHS.actualcoord = Add(negativeRHS.actualcoord, - Mul(alpha32, psiNF.dummy,
                                                             Pow(psiNF.dummy.dot(psiNF.dummy), - 1)))

W136NF = critical_linearsolver(W136NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

SS_W145_10 = first_order(1, 0, W145NF)
SS_W145_01 = first_order(0, 1, W145NF)
DS_W22W134_00 = second_order(0, 0, W22NF, W134NF)
DS_W23W133_00 = second_order(0, 0, W23NF, W133NF)
DS_W22W133_10 = second_order(1, 0, W22NF, W133NF)
DS_W22W133_01 = second_order(0, 1, W22NF, W133NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W145_10), Mul(bNF[1], SS_W145_01),
                              Mul(- 2, DS_phiW035_00), Mul(- 2, DS_W12W034_00),
                              Mul(- 2, DS_W22W134_00), Mul(- 2, DS_W23W133_00),
                              Mul(- 2, aNF[1], DS_phiW034_10), Mul(- 2, bNF[1], DS_phiW034_01),
                              Mul(- 2, aNF[1], DS_W22W133_10), Mul(- 2, bNF[1], DS_W22W133_01),
                              Mul(- 3, TS_phiphiW134_00), Mul(- 6, TS_phiW12W133_00),
                              Mul(- 3, aNF[1], TS_phiphiW133_10), Mul(- 3, bNF[1], TS_phiphiW133_01),
                              Mul(2, sqrt(muNF), diffmatrix, W124NF.dummy)).subs(extraparvals)

alpha42 = psiNF.dummy.dot(negativeRHS.actualcoord)

negativeRHS.actualcoord = Add(negativeRHS.actualcoord, - Mul(alpha42, psiNF.dummy,
                                                             Pow(psiNF.dummy.dot(psiNF.dummy), - 1)))

W146NF = critical_linearsolver(W146NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

SS_W155_10 = first_order(1, 0, W155NF)
SS_W155_01 = first_order(0, 1, W155NF)
SS_W14_20 = first_order(2, 0, W14NF)
SS_W14_11 = first_order(1, 1, W14NF)
SS_W14_02 = first_order(0, 2, W14NF)
SS_W13_30 = first_order(3, 0, W13NF)
SS_W13_21 = first_order(2, 1, W13NF)
SS_W13_12 = first_order(1, 2, W13NF)
SS_W13_03 = first_order(0, 3, W13NF)
SS_W12_40 = first_order(4, 0, W12NF)
SS_W12_31 = first_order(3, 1, W12NF)
SS_W12_22 = first_order(2, 2, W12NF)
SS_W12_13 = first_order(1, 3, W12NF)
SS_W12_04 = first_order(0, 4, W12NF)
SS_phi_50 = first_order(5, 0, phiNF)
SS_phi_41 = first_order(4, 1, phiNF)
SS_phi_32 = first_order(3, 2, phiNF)
SS_phi_23 = first_order(2, 3, phiNF)
SS_phi_14 = first_order(1, 4, phiNF)
SS_phi_05 = first_order(0, 5, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W155_10), Mul(bNF[1], SS_W155_01),
                              Mul(aNF[2], SS_W14_10), Mul(bNF[2], SS_W14_01),
                              Mul(aNF[3], SS_W13_10), Mul(bNF[3], SS_W13_01),
                              Mul(aNF[4], SS_W12_10), Mul(bNF[4], SS_W12_01),
                              Mul(aNF[5], SS_phi_10), Mul(bNF[5], SS_phi_01),
                              Mul(Pow(aNF[1], 2), SS_W14_20), Mul(aNF[1], bNF[1], SS_W14_11),
                              Mul(Pow(bNF[1], 2), SS_W14_02), Mul(2, aNF[1], aNF[2], SS_W13_20),
                              Mul(Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), SS_W13_11),
                              Mul(2, bNF[1], bNF[2], SS_W13_02), Mul(2, aNF[1], aNF[3], SS_W12_20),
                              Mul(Add(Mul(aNF[1], bNF[3]), Mul(aNF[3], bNF[1])), SS_W12_11),
                              Mul(2, bNF[1], bNF[3], SS_W12_02), Mul(2, aNF[1], aNF[4], SS_phi_20),
                              Mul(Add(Mul(aNF[1], bNF[4]), Mul(aNF[4], bNF[1])), SS_phi_11),
                              Mul(2, bNF[1], bNF[4], SS_phi_02), Mul(Pow(aNF[2], 2), SS_W12_20),
                              Mul(aNF[2], bNF[2], SS_W12_11), Mul(Pow(bNF[2], 2), SS_W12_02),
                              Mul(2, aNF[2], aNF[3], SS_phi_20),
                              Mul(Add(Mul(aNF[2], bNF[3]), Mul(aNF[3], bNF[2])), SS_phi_11),
                              Mul(2, bNF[2], bNF[3], SS_phi_02), Mul(Pow(aNF[1], 3), SS_W13_30),
                              Mul(Pow(aNF[1], 2), bNF[1], SS_W13_21),
                              Mul(aNF[1], Pow(bNF[1], 2), SS_W13_12), Mul(Pow(bNF[1], 3), SS_W13_03),
                              Mul(3, Pow(aNF[1], 2), aNF[2], SS_W12_30),
                              Mul(2, aNF[1], aNF[2], bNF[1], SS_W12_21),
                              Mul(aNF[2], Pow(bNF[1], 2), SS_W12_12),
                              Mul(Pow(aNF[1], 2), bNF[2], SS_W12_21),
                              Mul(2, aNF[1], bNF[1], bNF[2], SS_W12_12),
                              Mul(3, Pow(bNF[1], 2), bNF[2], SS_W12_03),
                              Mul(3, Pow(aNF[1], 2), aNF[3], SS_phi_30),
                              Mul(2, aNF[1], aNF[3], bNF[1], SS_phi_21),
                              Mul(aNF[3], Pow(bNF[1], 2), SS_phi_12),
                              Mul(Pow(aNF[1], 2), bNF[3], SS_phi_21),
                              Mul(2, aNF[1], bNF[1], bNF[3], SS_phi_12),
                              Mul(3, Pow(bNF[1], 2), bNF[3], SS_phi_03),
                              Mul(3, aNF[1], Pow(aNF[2], 2), SS_phi_30),
                              Mul(2, aNF[1], aNF[2], bNF[2], SS_phi_21),
                              Mul(aNF[1], Pow(bNF[2], 2), SS_phi_12),
                              Mul(Pow(aNF[2], 2), bNF[1], SS_phi_21),
                              Mul(2, aNF[2], bNF[1], bNF[2], SS_phi_12),
                              Mul(3, bNF[1], Pow(bNF[2], 2), SS_phi_03),
                              Mul(Pow(aNF[1], 4), SS_W12_40), Mul(Pow(aNF[1], 3), bNF[1], SS_W12_31),
                              Mul(Pow(aNF[1], 2), Pow(bNF[1], 2), SS_W12_22),
                              Mul(aNF[1], Pow(bNF[1], 3), SS_W12_13),
                              Mul(Pow(bNF[1], 4), SS_W12_04),
                              Mul(4, Pow(aNF[1], 3), aNF[2], SS_phi_40),
                              Mul(3, Pow(aNF[1], 2), aNF[2], bNF[1], SS_phi_31),
                              Mul(2, aNF[1], aNF[2], Pow(bNF[1], 2), SS_phi_22),
                              Mul(aNF[2], Pow(bNF[1], 3), SS_phi_13),
                              Mul(Pow(aNF[1], 3), bNF[2], SS_phi_31),
                              Mul(2, Pow(aNF[1], 2), bNF[1], bNF[2], SS_phi_22),
                              Mul(3, aNF[1], Pow(bNF[1], 2), bNF[2], SS_phi_13),
                              Mul(4, Pow(bNF[1], 3), bNF[2], SS_phi_04),
                              Mul(Pow(aNF[1], 5), SS_phi_50),
                              Mul(Pow(aNF[1], 4), bNF[1], SS_phi_41),
                              Mul(Pow(aNF[1], 3), Pow(bNF[1], 2), SS_phi_32),
                              Mul(Pow(aNF[1], 2), Pow(bNF[1], 3), SS_phi_23),
                              Mul(aNF[1], Pow(bNF[1], 4), SS_phi_14),
                              Mul(Pow(bNF[1], 5), SS_phi_05)).subs(extraparvals)

alpha52 = psiNF.dummy.dot(negativeRHS.actualcoord)

negativeRHS.actualcoord = Add(negativeRHS.actualcoord, - Mul(alpha52, psiNF.dummy,
                                                             Pow(psiNF.dummy.dot(psiNF.dummy), - 1)))

W156NF = critical_linearsolver(W156NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

SS_W165_10 = first_order(1, 0, W165NF)
SS_W165_01 = first_order(0, 1, W165NF)
SS_W124_20 = first_order(2, 0, W124NF)
SS_W124_11 = first_order(1, 1, W124NF)
SS_W124_02 = first_order(0, 2, W124NF)
SS_W123_30 = first_order(3, 0, W123NF)
SS_W123_21 = first_order(2, 1, W123NF)
SS_W123_12 = first_order(1, 2, W123NF)
SS_W123_03 = first_order(0, 3, W123NF)
DS_phiW05_00 = second_order(0, 0, phiNF, W05NF)
DS_phiW25_00 = second_order(0, 0, phiNF, W25NF)
DS_W12W04_00 = second_order(0, 0, W12NF, W04NF)
DS_W12W24_00 = second_order(0, 0, W12NF, W24NF)
DS_W13W03_00 = second_order(0, 0, W13NF, W03NF)
DS_W13W23_00 = second_order(0, 0, W13NF, W23NF)
DS_W14W02_00 = second_order(0, 0, W14NF, W02NF)
DS_W14W22_00 = second_order(0, 0, W14NF, W22NF)
DS_phiW04_10 = second_order(1, 0, phiNF, W04NF)
DS_phiW24_10 = second_order(1, 0, phiNF, W24NF)
DS_phiW04_01 = second_order(0, 1, phiNF, W04NF)
DS_phiW24_01 = second_order(0, 1, phiNF, W24NF)
DS_W12W03_10 = second_order(1, 0, W12NF, W03NF)
DS_W12W23_10 = second_order(1, 0, W12NF, W23NF)
DS_W12W03_01 = second_order(0, 1, W12NF, W03NF)
DS_W12W23_01 = second_order(0, 1, W12NF, W23NF)
DS_W13W02_10 = second_order(1, 0, W13NF, W02NF)
DS_W13W22_10 = second_order(1, 0, W13NF, W22NF)
DS_W13W02_01 = second_order(0, 1, W13NF, W02NF)
DS_W13W22_01 = second_order(0, 1, W13NF, W22NF)
DS_phiW03_20 = second_order(2, 0, phiNF, W03NF)
DS_phiW23_20 = second_order(2, 0, phiNF, W23NF)
DS_phiW03_11 = second_order(1, 1, phiNF, W03NF)
DS_phiW23_11 = second_order(1, 1, phiNF, W23NF)
DS_phiW03_02 = second_order(0, 2, phiNF, W03NF)
DS_phiW23_02 = second_order(0, 2, phiNF, W23NF)
DS_W12W02_20 = second_order(2, 0, W12NF, W02NF)
DS_W12W22_20 = second_order(2, 0, W12NF, W22NF)
DS_W12W02_11 = second_order(1, 1, W12NF, W02NF)
DS_W12W22_11 = second_order(1, 1, W12NF, W22NF)
DS_W12W02_02 = second_order(0, 2, W12NF, W02NF)
DS_W12W22_02 = second_order(0, 2, W12NF, W22NF)
DS_phiW02_30 = second_order(3, 0, phiNF, W02NF)
DS_phiW22_30 = second_order(3, 0, phiNF, W22NF)
DS_phiW02_21 = second_order(2, 1, phiNF, W02NF)
DS_phiW22_21 = second_order(2, 1, phiNF, W22NF)
DS_phiW02_12 = second_order(1, 2, phiNF, W02NF)
DS_phiW22_12 = second_order(1, 2, phiNF, W22NF)
DS_phiW02_03 = second_order(0, 3, phiNF, W02NF)
DS_phiW22_03 = second_order(0, 3, phiNF, W22NF)
TS_phiphiW14_00 = third_order(0, 0, phiNF, phiNF, W14NF)
TS_phiW12W13_00 = third_order(0, 0, phiNF, W12NF, W13NF)
TS_W12W12W12_00 = third_order(0, 0, W12NF, W12NF, W12NF)
TS_phiphiW13_10 = third_order(1, 0, phiNF, phiNF, W13NF)
TS_phiphiW13_01 = third_order(0, 1, phiNF, phiNF, W13NF)
TS_phiW12W12_10 = third_order(1, 0, phiNF, W12NF, W12NF)
TS_phiW12W12_01 = third_order(0, 1, phiNF, W12NF, W12NF)
TS_phiphiW12_20 = third_order(2, 0, phiNF, phiNF, W12NF)
TS_phiphiW12_11 = third_order(1, 1, phiNF, phiNF, W12NF)
TS_phiphiW12_02 = third_order(0, 2, phiNF, phiNF, W12NF)
TS_phiphiphi_30 = third_order(3, 0, phiNF, phiNF, phiNF)
TS_phiphiphi_21 = third_order(2, 1, phiNF, phiNF, phiNF)
TS_phiphiphi_12 = third_order(1, 2, phiNF, phiNF, phiNF)
TS_phiphiphi_03 = third_order(0, 3, phiNF, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W165_10), Mul(bNF[1], SS_W165_01),
                              Mul(aNF[2], SS_W124_10), Mul(bNF[2], SS_W124_01),
                              Mul(aNF[3], SS_W123_10), Mul(bNF[3], SS_W123_01),
                              Mul(Pow(aNF[1], 2), SS_W124_20), Mul(aNF[1], bNF[1], SS_W124_11),
                              Mul(Pow(bNF[1], 2), SS_W124_02), Mul(2, aNF[1], aNF[2], SS_W123_20),
                              Mul(Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), SS_W123_11),
                              Mul(2, bNF[1], bNF[2], SS_W123_02), Mul(Pow(aNF[1], 3), SS_W123_30),
                              Mul(Pow(aNF[1], 2), bNF[1], SS_W123_21),
                              Mul(aNF[1], Pow(bNF[1], 2), SS_W123_12),
                              Mul(Pow(bNF[1], 3), SS_W123_03), Mul(4, DS_phiW05_00),
                              Mul(2, DS_phiW25_00), Mul(4, DS_W12W04_00),
                              Mul(2, DS_W12W24_00), Mul(4, DS_W13W03_00),
                              Mul(2, DS_W13W23_00), Mul(4, DS_W14W02_00),
                              Mul(2, DS_W14W22_00), Mul(4, aNF[1], DS_phiW04_10),
                              Mul(2, aNF[1], DS_phiW24_10), Mul(4, bNF[1], DS_phiW04_01),
                              Mul(2, bNF[1], DS_phiW24_01), Mul(4, aNF[1], DS_W12W03_10),
                              Mul(2, aNF[1], DS_W12W23_10), Mul(4, bNF[1], DS_W12W03_01),
                              Mul(2, bNF[1], DS_W12W23_01), Mul(4, aNF[1], DS_W13W02_10),
                              Mul(2, aNF[1], DS_W13W22_10), Mul(4, bNF[1], DS_W13W02_01),
                              Mul(2, bNF[1], DS_W13W22_01), Mul(4, aNF[2], DS_phiW03_10),
                              Mul(2, aNF[2], DS_phiW23_10), Mul(4, bNF[2], DS_phiW03_01),
                              Mul(2, bNF[2], DS_phiW23_01), Mul(4, aNF[2], DS_W12W02_10),
                              Mul(2, aNF[2], DS_W12W22_10), Mul(4, bNF[2], DS_W12W02_01),
                              Mul(2, bNF[2], DS_W12W22_01), Mul(4, aNF[3], DS_phiW02_10),
                              Mul(2, aNF[3], DS_phiW22_10), Mul(4, bNF[3], DS_phiW02_01),
                              Mul(2, bNF[3], DS_phiW22_01), Mul(4, Pow(aNF[1], 2), DS_phiW03_20),
                              Mul(2, Pow(aNF[1], 2), DS_phiW23_20),
                              Mul(4, aNF[1], bNF[1], DS_phiW03_11),
                              Mul(2, aNF[1], bNF[1], DS_phiW23_11),
                              Mul(4, Pow(bNF[1], 2), DS_phiW03_02),
                              Mul(2, Pow(bNF[1], 2), DS_phiW23_02),
                              Mul(4, Pow(aNF[1], 2), DS_W12W02_20),
                              Mul(2, Pow(aNF[1], 2), DS_W12W22_20),
                              Mul(4, aNF[1], bNF[1], DS_W12W02_11),
                              Mul(2, aNF[1], bNF[1], DS_W12W22_11),
                              Mul(4, Pow(bNF[1], 2), DS_W12W02_02),
                              Mul(2, Pow(bNF[1], 2), DS_W12W22_02),
                              Mul(8, aNF[1], aNF[2], DS_phiW02_20),
                              Mul(4, aNF[1], aNF[2], DS_phiW22_20),
                              Mul(4, Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), DS_phiW02_11),
                              Mul(2, Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), DS_phiW22_11),
                              Mul(8, bNF[1], bNF[2], DS_phiW02_02),
                              Mul(4, bNF[1], bNF[2], DS_phiW22_02),
                              Mul(4, Pow(aNF[1], 3), DS_phiW02_30),
                              Mul(2, Pow(aNF[1], 3), DS_phiW22_30),
                              Mul(4, Pow(aNF[1], 2), bNF[1], DS_phiW02_21),
                              Mul(2, Pow(aNF[1], 2), bNF[1], DS_phiW22_21),
                              Mul(4, aNF[1], Pow(bNF[1], 2), DS_phiW02_12),
                              Mul(2, aNF[1], Pow(bNF[1], 2), DS_phiW22_12),
                              Mul(4, Pow(bNF[1], 3), DS_phiW02_03),
                              Mul(2, Pow(bNF[1], 3), DS_phiW22_03), Mul(9, TS_phiphiW14_00),
                              Mul(18, TS_phiW12W13_00), Mul(3, TS_W12W12W12_00),
                              Mul(9, aNF[1], TS_phiphiW13_10), Mul(9, bNF[1], TS_phiphiW13_01),
                              Mul(9, aNF[1], TS_phiW12W12_10), Mul(9, bNF[1], TS_phiW12W12_01),
                              Mul(9, aNF[2], TS_phiphiW12_10), Mul(9, bNF[2], TS_phiphiW12_01),
                              Mul(3, aNF[3], TS_phiphiphi_10), Mul(3, bNF[3], TS_phiphiphi_01),
                              Mul(9, Pow(aNF[1], 2), TS_phiphiW12_20),
                              Mul(9, aNF[1], bNF[1], TS_phiphiW12_11),
                              Mul(9, Pow(bNF[1], 2), TS_phiphiW12_02),
                              Mul(6, aNF[1], aNF[2], TS_phiphiphi_20),
                              Mul(3, Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), TS_phiphiphi_11),
                              Mul(6, bNF[1], bNF[2], TS_phiphiphi_02),
                              Mul(3, Pow(aNF[1], 3), TS_phiphiphi_30),
                              Mul(3, Pow(aNF[1], 2), bNF[1], TS_phiphiphi_21),
                              Mul(3, aNF[1], Pow(bNF[1], 2), TS_phiphiphi_12),
                              Mul(3, Pow(bNF[1], 3), TS_phiphiphi_03)).subs(extraparvals)

alpha62 = psiNF.dummy.dot(negativeRHS.actualcoord)

negativeRHS.actualcoord = Add(negativeRHS.actualcoord, - Mul(alpha62, psiNF.dummy,
                                                             Pow(psiNF.dummy.dot(psiNF.dummy), - 1)))

W166NF = critical_linearsolver(W166NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

SS_W175_10 = first_order(1, 0, W175NF)
SS_W175_01 = first_order(0, 1, W175NF)
DS_phiW025_00 = second_order(0, 0, phiNF, W025NF)
DS_phiW225_00 = second_order(0, 0, phiNF, W225NF)
DS_W12W024_00 = second_order(0, 0, W12NF, W024NF)
DS_W12W224_00 = second_order(0, 0, W12NF, W224NF)
DS_W22W34_00 = second_order(0, 0, W22NF, W34NF)
DS_W123W03_00 = second_order(0, 0, W123NF, W03NF)
DS_W123W23_00 = second_order(0, 0, W123NF, W23NF)
DS_W23W33_00 = second_order(0, 0, W23NF, W33NF)
DS_W124W02_00 = second_order(0, 0, W124NF, W02NF)
DS_W124W22_00 = second_order(0, 0, W124NF, W22NF)
DS_phiW024_10 = second_order(1, 0, phiNF, W024NF)
DS_phiW224_10 = second_order(1, 0, phiNF, W224NF)
DS_phiW024_01 = second_order(0, 1, phiNF, W024NF)
DS_phiW224_01 = second_order(0, 1, phiNF, W224NF)
DS_W22W33_10 = second_order(1, 0, W22NF, W33NF)
DS_W22W33_01 = second_order(0, 1, W22NF, W33NF)
DS_W123W02_10 = second_order(1, 0, W123NF, W02NF)
DS_W123W22_10 = second_order(1, 0, W123NF, W22NF)
DS_W123W02_01 = second_order(0, 1, W123NF, W02NF)
DS_W123W22_01 = second_order(0, 1, W123NF, W22NF)
TS_phiphiW124_00 = third_order(0, 0, phiNF, phiNF, W124NF)
TS_phiphiW34_00 = third_order(0, 0, phiNF, phiNF, W34NF)
TS_phiW12W123_00 = third_order(0, 0, phiNF, W12NF, W123NF)
TS_phiW12W33_00 = third_order(0, 0, phiNF, W12NF, W33NF)
TS_phiW03W02_00 = third_order(0, 0, phiNF, W03NF, W02NF)
TS_phiW03W22_00 = third_order(0, 0, phiNF, W03NF, W22NF)
TS_phiW23W02_00 = third_order(0, 0, phiNF, W23NF, W02NF)
TS_phiW23W22_00 = third_order(0, 0, phiNF, W23NF, W22NF)
TS_W12W02W02_00 = third_order(0, 0, W12NF, W02NF, W02NF)
TS_W12W02W22_00 = third_order(0, 0, W12NF, W02NF, W22NF)
TS_W12W22W22_00 = third_order(0, 0, W12NF, W22NF, W22NF)
TS_phiphiW123_10 = third_order(1, 0, phiNF, phiNF, W123NF)
TS_phiphiW33_10 = third_order(1, 0, phiNF, phiNF, W33NF)
TS_phiphiW123_01 = third_order(0, 1, phiNF, phiNF, W123NF)
TS_phiphiW33_01 = third_order(0, 1, phiNF, phiNF, W33NF)
TS_phiW02W02_10 = third_order(1, 0, phiNF, W02NF, W02NF)
TS_phiW02W02_01 = third_order(0, 1, phiNF, W02NF, W02NF)
TS_phiW22W02_10 = third_order(1, 0, phiNF, W22NF, W02NF)
TS_phiW22W22_10 = third_order(1, 0, phiNF, W22NF, W22NF)
TS_phiW22W02_01 = third_order(0, 1, phiNF, W22NF, W02NF)
TS_phiW22W22_01 = third_order(0, 1, phiNF, W22NF, W22NF)
Q4S_phiphiphiW03_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W03NF)
Q4S_phiphiphiW23_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W23NF)
Q4S_phiphiW12W02_00 = fourth_order(0, 0, phiNF, phiNF, W12NF, W02NF)
Q4S_phiphiW12W22_00 = fourth_order(0, 0, phiNF, phiNF, W12NF, W22NF)
Q4S_phiphiphiW02_10 = fourth_order(1, 0, phiNF, phiNF, phiNF, W02NF)
Q4S_phiphiphiW22_10 = fourth_order(1, 0, phiNF, phiNF, phiNF, W22NF)
Q4S_phiphiphiW02_01 = fourth_order(0, 1, phiNF, phiNF, phiNF, W02NF)
Q4S_phiphiphiW22_01 = fourth_order(0, 1, phiNF, phiNF, phiNF, W22NF)
Q5S_phiphiphiphiW12_00 = fifth_order(0, 0, phiNF, phiNF, phiNF, phiNF, W12NF)
Q5S_phiphiphiphiphi_10 = fifth_order(1, 0, phiNF, phiNF, phiNF, phiNF, phiNF)
Q5S_phiphiphiphiphi_01 = fifth_order(0, 1, phiNF, phiNF, phiNF, phiNF, phiNF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W175_10), Mul(bNF[1], SS_W175_01),
                              Mul(4, DS_phiW025_00), Mul(2, DS_phiW225_00),
                              Mul(4, DS_W12W024_00), Mul(2, DS_W12W224_00),
                              Mul(2, DS_W22W34_00), Mul(4, DS_W123W03_00),
                              Mul(2, DS_W123W23_00), Mul(2, DS_W23W33_00),
                              Mul(4, DS_W124W02_00), Mul(2, DS_W124W22_00),
                              Mul(4, aNF[1], DS_phiW024_10), Mul(2, aNF[1], DS_phiW224_10),
                              Mul(4, bNF[1], DS_phiW024_01), Mul(2, bNF[1], DS_phiW224_01),
                              Mul(2, aNF[1], DS_W22W33_10), Mul(2, bNF[1], DS_W22W33_01),
                              Mul(4, aNF[1], DS_W123W02_10), Mul(2, aNF[1], DS_W123W22_10),
                              Mul(4, bNF[1], DS_W123W02_01), Mul(2, bNF[1], DS_W123W22_01),
                              Mul(9, TS_phiphiW124_00), Mul(3, TS_phiphiW34_00),
                              Mul(18, TS_phiW12W123_00), Mul(6, TS_phiW12W33_00),
                              Mul(24, TS_phiW03W02_00), Mul(12, TS_phiW03W22_00),
                              Mul(12, TS_phiW23W02_00), Mul(12, TS_phiW23W22_00),
                              Mul(12, TS_W12W02W02_00), Mul(12, TS_W12W02W22_00),
                              Mul(6, TS_W12W22W22_00), Mul(9, aNF[1], TS_phiphiW123_10),
                              Mul(3, aNF[1], TS_phiphiW33_10), Mul(9, bNF[1], TS_phiphiW123_01),
                              Mul(3, bNF[1], TS_phiphiW33_01), Mul(12, aNF[1], TS_phiW02W02_10),
                              Mul(12, bNF[1], TS_phiW02W02_01), Mul(12, aNF[1], TS_phiW22W02_10),
                              Mul(6, aNF[1], TS_phiW22W22_10), Mul(12, bNF[1], TS_phiW22W02_01),
                              Mul(6, bNF[1], TS_phiW22W22_01), Mul(24, Q4S_phiphiphiW03_00),
                              Mul(16, Q4S_phiphiphiW23_00), Mul(72, Q4S_phiphiW12W02_00),
                              Mul(48, Q4S_phiphiW12W22_00), Mul(24, aNF[1], Q4S_phiphiphiW02_10),
                              Mul(16, aNF[1], Q4S_phiphiphiW22_10),
                              Mul(24, bNF[1], Q4S_phiphiphiW02_01),
                              Mul(16, bNF[1], Q4S_phiphiphiW22_01), Mul(50, Q5S_phiphiphiphiW12_00),
                              Mul(10, aNF[1], Q5S_phiphiphiphiphi_10),
                              Mul(10, bNF[1], Q5S_phiphiphiphiphi_01)).subs(extraparvals)

alpha72 = psiNF.dummy.dot(negativeRHS.actualcoord)

negativeRHS.actualcoord = Add(negativeRHS.actualcoord, - Mul(alpha72, psiNF.dummy,
                                                             Pow(psiNF.dummy.dot(psiNF.dummy), - 1)))

W176NF = critical_linearsolver(W176NF, negativeRHS, criticalcol, coefsubmatrix, submatrixrows, submatrixcols)

SS_W25_10 = first_order(1, 0, W25NF)
SS_W25_01 = first_order(0, 1, W25NF)
SS_W24_20 = first_order(2, 0, W24NF)
SS_W24_11 = first_order(1, 1, W24NF)
SS_W24_02 = first_order(0, 2, W24NF)
SS_W23_30 = first_order(3, 0, W23NF)
SS_W23_21 = first_order(2, 1, W23NF)
SS_W23_12 = first_order(1, 2, W23NF)
SS_W23_03 = first_order(0, 3, W23NF)
SS_W22_40 = first_order(4, 0, W22NF)
SS_W22_31 = first_order(3, 1, W22NF)
SS_W22_22 = first_order(2, 2, W22NF)
SS_W22_13 = first_order(1, 3, W22NF)
SS_W22_04 = first_order(0, 4, W22NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W25_10), Mul(bNF[1], SS_W25_01),
                              Mul(aNF[2], SS_W24_10), Mul(bNF[2], SS_W24_01),
                              Mul(aNF[3], SS_W23_10), Mul(bNF[3], SS_W23_01),
                              Mul(aNF[4], SS_W22_10), Mul(bNF[4], SS_W22_01),
                              Mul(Pow(aNF[1], 2), SS_W24_20), Mul(aNF[1], bNF[1], SS_W24_11),
                              Mul(Pow(bNF[1], 2), SS_W24_02), Mul(2, aNF[1], aNF[2], SS_W23_20),
                              Mul(Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), SS_W23_11),
                              Mul(2, bNF[1], bNF[2], SS_W23_02), Mul(2, aNF[1], aNF[3], SS_W22_20),
                              Mul(Add(Mul(aNF[1], bNF[3]), Mul(aNF[3], bNF[1])), SS_W22_11),
                              Mul(2, bNF[1], bNF[3], SS_W22_02), Mul(Pow(aNF[2], 2), SS_W22_20),
                              Mul(aNF[2], bNF[2], SS_W22_11), Mul(Pow(bNF[2], 2), SS_W22_02),
                              Mul(Pow(aNF[1], 3), SS_W23_30), Mul(Pow(aNF[1], 2), bNF[1], SS_W23_21),
                              Mul(aNF[1], Pow(bNF[1], 2), SS_W23_12), Mul(Pow(bNF[1], 3), SS_W23_03),
                              Mul(3, Pow(aNF[1], 2), aNF[2], SS_W22_30),
                              Mul(2, aNF[1], aNF[2], bNF[1], SS_W22_21),
                              Mul(aNF[2], Pow(bNF[1], 2), SS_W22_12),
                              Mul(Pow(aNF[1], 2), bNF[2], SS_W22_21),
                              Mul(2, aNF[1], bNF[1], bNF[2], SS_W22_12),
                              Mul(3, Pow(bNF[1], 2), bNF[2], SS_W22_03),
                              Mul(Pow(aNF[1], 4), SS_W22_40), Mul(Pow(aNF[1], 3), bNF[1], SS_W22_31),
                              Mul(Pow(aNF[1], 2), Pow(bNF[1], 2), SS_W22_22),
                              Mul(aNF[1], Pow(bNF[1], 3), SS_W22_13), Mul(Pow(bNF[1], 4), SS_W22_04),
                              Mul(2, DS_phiW155_00), Mul(2, DS_W12W14_00),
                              DS_W13W13_00, Mul(2, aNF[1], DS_phiW14_10),
                              Mul(2, bNF[1], DS_phiW14_01), Mul(2, aNF[1], DS_W12W13_10),
                              Mul(2, bNF[1], DS_W12W13_01), Mul(2, aNF[2], DS_phiW13_10),
                              Mul(2, bNF[2], DS_phiW13_01), Mul(aNF[2], DS_W12W12_10),
                              Mul(bNF[2], DS_W12W12_01), Mul(2, aNF[3], DS_phiW12_10),
                              Mul(2, bNF[3], DS_phiW12_01), Mul(aNF[4], DS_phiphi_10),
                              Mul(bNF[4], DS_phiphi_01), Mul(2, Pow(aNF[1], 2), DS_phiW13_20),
                              Mul(2, aNF[1], bNF[1], DS_phiW13_11),
                              Mul(2, Pow(bNF[1], 2), DS_phiW13_02),
                              Mul(Pow(aNF[1], 2), DS_W12W12_20), Mul(aNF[1], bNF[1], DS_W12W12_11),
                              Mul(Pow(bNF[1], 2), DS_W12W12_02), Mul(4, aNF[1], aNF[2], DS_phiW12_20),
                              Mul(2, Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), DS_phiW12_11),
                              Mul(4, bNF[1], bNF[2], DS_phiW12_02),
                              Mul(2, aNF[1], aNF[3], DS_phiphi_20),
                              Mul(Add(Mul(aNF[1], bNF[3]), Mul(aNF[3], bNF[1])), DS_phiphi_11),
                              Mul(2, bNF[1], bNF[3], DS_phiphi_02), Mul(Pow(aNF[2], 2), DS_phiphi_20),
                              Mul(aNF[2], bNF[2], DS_phiphi_11), Mul(Pow(bNF[2], 2), DS_phiphi_02),
                              Mul(2, Pow(aNF[1], 3), DS_phiW12_30),
                              Mul(2, Pow(aNF[1], 2), bNF[1], DS_phiW12_21),
                              Mul(2, aNF[1], Pow(bNF[1], 2), DS_phiW12_12),
                              Mul(2, Pow(bNF[1], 3), DS_phiW12_03),
                              Mul(3, Pow(aNF[1], 2), aNF[2], DS_phiphi_30),
                              Mul(2, aNF[1], aNF[2], bNF[1], DS_phiphi_21),
                              Mul(aNF[2], Pow(bNF[1], 2), DS_phiphi_12),
                              Mul(Pow(aNF[1], 2), bNF[2], DS_phiphi_21),
                              Mul(2, aNF[1], bNF[1], bNF[2], DS_phiphi_12),
                              Mul(3, Pow(bNF[1], 2), bNF[2], DS_phiphi_03),
                              Mul(Pow(aNF[1], 4), DS_phiphi_40),
                              Mul(Pow(aNF[1], 3), bNF[1], DS_phiphi_31),
                              Mul(Pow(aNF[1], 2), Pow(bNF[1], 2), DS_phiphi_22),
                              Mul(aNF[1], Pow(bNF[1], 3), DS_phiphi_13),
                              Mul(Pow(bNF[1], 4), DS_phiphi_04)).subs(extraparvals)

W26NF = linearsolver(W26NF, negativeRHS, coefmat2)

SS_W225_10 = first_order(1, 0, W225NF)
SS_W225_01 = first_order(0, 1, W225NF)
SS_W224_20 = first_order(2, 0, W224NF)
SS_W224_11 = first_order(1, 1, W224NF)
SS_W224_02 = first_order(0, 2, W224NF)
DS_phiW35_00 = second_order(0, 0, phiNF, W35NF)
DS_W02W24_00 = second_order(0, 0, W02NF, W24NF)
DS_W12W34_00 = second_order(0, 0, W12NF, W34NF)
DS_W22W04_00 = second_order(0, 0, W22NF, W04NF)
DS_W03W23_00 = second_order(0, 0, W03NF, W23NF)
DS_W13W33_00 = second_order(0, 0, W13NF, W33NF)
DS_phiW34_10 = second_order(1, 0, phiNF, W34NF)
DS_phiW34_01 = second_order(0, 1, phiNF, W34NF)
DS_W02W23_10 = second_order(1, 0, W02NF, W23NF)
DS_W02W23_01 = second_order(0, 1, W02NF, W23NF)
DS_W12W33_10 = second_order(1, 0, W12NF, W33NF)
DS_W12W33_01 = second_order(0, 1, W12NF, W33NF)
DS_W22W03_10 = second_order(1, 0, W22NF, W03NF)
DS_W22W03_01 = second_order(0, 1, W22NF, W03NF)
DS_phiW33_20 = second_order(2, 0, phiNF, W33NF)
DS_phiW33_11 = second_order(1, 1, phiNF, W33NF)
DS_phiW33_02 = second_order(0, 2, phiNF, W33NF)
DS_W02W22_20 = second_order(2, 0, W02NF, W22NF)
DS_W02W22_11 = second_order(1, 1, W02NF, W22NF)
DS_W02W22_02 = second_order(0, 2, W02NF, W22NF)
TS_phiW13W02_00 = third_order(0, 0, phiNF, W13NF, W02NF)
TS_phiW13W22_00 = third_order(0, 0, phiNF, W13NF, W22NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W225_10), Mul(bNF[1], SS_W225_01),
                              Mul(aNF[2], SS_W224_10), Mul(bNF[2], SS_W224_01),
                              Mul(Pow(aNF[1], 2), SS_W224_20), Mul(aNF[1], bNF[1], SS_W224_11),
                              Mul(Pow(bNF[1], 2), SS_W224_02), Mul(2, DS_phiW165_00),
                              Mul(2, DS_phiW35_00), Mul(4, DS_W02W24_00),
                              Mul(2, DS_W12W124_00), Mul(2, DS_W12W34_00),
                              Mul(4, DS_W22W04_00), Mul(4, DS_W03W23_00),
                              Mul(2, DS_W13W123_00), Mul(2, DS_W13W33_00),
                              Mul(2, aNF[1], DS_phiW124_10), Mul(2, aNF[1], DS_phiW34_10),
                              Mul(2, bNF[1], DS_phiW124_01), Mul(2, bNF[1], DS_phiW34_01),
                              Mul(4, aNF[1], DS_W02W23_10), Mul(4, bNF[1], DS_W02W23_01),
                              Mul(2, aNF[1], DS_W12W123_10), Mul(2, aNF[1], DS_W12W33_10),
                              Mul(2, bNF[1], DS_W12W123_01), Mul(2, bNF[1], DS_W12W33_01),
                              Mul(4, aNF[1], DS_W22W03_10), Mul(4, bNF[1], DS_W22W03_01),
                              Mul(2, aNF[2], DS_phiW123_10), Mul(2, aNF[2], DS_phiW33_10),
                              Mul(2, bNF[2], DS_phiW123_01), Mul(2, bNF[2], DS_phiW33_01),
                              Mul(4, aNF[2], DS_W02W22_10), Mul(4, bNF[2], DS_W02W22_01),
                              Mul(2, Pow(aNF[1], 2), DS_phiW123_20), Mul(2, Pow(aNF[1], 2), DS_phiW33_20),
                              Mul(2, aNF[1], bNF[1], DS_phiW123_11), Mul(2, aNF[1], bNF[1], DS_phiW33_11),
                              Mul(2, Pow(bNF[1], 2), DS_phiW123_02), Mul(2, Pow(bNF[1], 2), DS_phiW33_02),
                              Mul(4, Pow(aNF[1], 2), DS_W02W22_20), Mul(4, aNF[1], bNF[1], DS_W02W22_11),
                              Mul(4, Pow(bNF[1], 2), DS_W02W22_02), Mul(6, TS_phiphiW04_00),
                              Mul(6, TS_phiphiW24_00), Mul(12, TS_phiW12W03_00),
                              Mul(12, TS_phiW12W23_00), Mul(12, TS_phiW13W02_00),
                              Mul(12, TS_phiW13W22_00), Mul(6, TS_W12W12W02_00),
                              Mul(6, TS_W12W12W22_00), Mul(6, aNF[1], TS_phiphiW03_10),
                              Mul(6, aNF[1], TS_phiphiW23_10), Mul(6, bNF[1], TS_phiphiW03_01),
                              Mul(6, bNF[1], TS_phiphiW23_01), Mul(12, aNF[1], TS_phiW12W02_10),
                              Mul(12, aNF[1], TS_phiW12W22_10), Mul(12, bNF[1], TS_phiW12W02_01),
                              Mul(12, bNF[1], TS_phiW12W22_01), Mul(6, aNF[2], TS_phiphiW02_10),
                              Mul(6, aNF[2], TS_phiphiW22_10), Mul(6, bNF[2], TS_phiphiW02_01),
                              Mul(6, bNF[2], TS_phiphiW22_01), Mul(6, Pow(aNF[1], 2), TS_phiphiW02_20),
                              Mul(6, Pow(aNF[1], 2), TS_phiphiW22_20), Mul(6, aNF[1], bNF[1], TS_phiphiW02_11),
                              Mul(6, aNF[1], bNF[1], TS_phiphiW22_11), Mul(6, Pow(bNF[1], 2), TS_phiphiW02_02),
                              Mul(6, Pow(bNF[1], 2), TS_phiphiW22_02), Mul(16, Q4S_phiphiphiW13_00),
                              Mul(24, Q4S_phiphiW12W12_00), Mul(16, aNF[1], Q4S_phiphiphiW12_10),
                              Mul(16, bNF[1], Q4S_phiphiphiW12_01), Mul(4, aNF[2], Q4S_phiphiphiphi_10),
                              Mul(4, bNF[2], Q4S_phiphiphiphi_01), Mul(4, Pow(aNF[1], 2), Q4S_phiphiphiphi_20),
                              Mul(4, aNF[1], bNF[1], Q4S_phiphiphiphi_11),
                              Mul(4, Pow(bNF[1], 2), Q4S_phiphiphiphi_02)).subs(extraparvals)

W226NF = linearsolver(W226NF, negativeRHS, coefmat2)

DS_phiW325_00 = second_order(0, 0, phiNF, W325NF)
DS_W02W224_00 = second_order(0, 0, W02NF, W224NF)
DS_W22W024_00 = second_order(0, 0, W22NF, W024NF)
DS_W22W44_00 = second_order(0, 0, W22NF, W44NF)
DS_W133W33_00 = second_order(0, 0, W133NF, W33NF)
TS_phiphiW44_00 = third_order(0, 0, phiNF, phiNF, W44NF)
TS_phiW02W123_00 = third_order(0, 0, phiNF, W02NF, W123NF)
TS_phiW02W33_00 = third_order(0, 0, phiNF, W02NF, W33NF)
TS_phiW22W123_00 = third_order(0, 0, phiNF, W22NF, W123NF)
TS_W02W02W22_00 = third_order(0, 0, W02NF, W02NF, W22NF)
TS_W22W22W22_00 = third_order(0, 0, W22NF, W22NF, W22NF)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW175_00), Mul(2, DS_phiW325_00),
                              Mul(4, DS_W02W224_00), Mul(4, DS_W22W024_00),
                              Mul(2, DS_W22W44_00), DS_W123W123_00,
                              Mul(2, DS_W133W33_00), Mul(6, TS_phiphiW024_00),
                              Mul(6, TS_phiphiW224_00), Mul(3, TS_phiphiW44_00),
                              Mul(12, TS_phiW02W123_00), Mul(12, TS_phiW02W33_00),
                              Mul(12, TS_phiW22W123_00), Mul(6, TS_phiW22W33_00),
                              Mul(12, TS_W02W02W22_00), Mul(3, TS_W22W22W22_00),
                              Mul(16, Q4S_phiphiphiW123_00), Mul(12, Q4S_phiphiphiW33_00),
                              Mul(24, Q4S_phiphiW02W02_00), Mul(48, Q4S_phiphiW22W02_00),
                              Mul(18, Q4S_phiphiW22W22_00), Mul(40, Q5S_phiphiphiphiW02_00),
                              Mul(35, Q5S_phiphiphiphiW22_00),
                              Mul(15, S6S_phiphiphiphiphiphi_00)).subs(extraparvals)

W236NF = linearsolver(W236NF, negativeRHS, coefmat2)

negativeRHS.actualcoord = Add(- DS_W133W133_00, Mul(- 4, sqrt(muNF), diffmatrix, W234NF.dummy),
                              Mul(2, diffmatrix, W22NF.dummy)).subs(extraparvals)

W246NF = linearsolver(W246NF, negativeRHS, coefmat2)

SS_W235_10 = first_order(1, 0, W235NF)
SS_W235_01 = first_order(0, 1, W235NF)
SS_W234_20 = first_order(2, 0, W234NF)
SS_W234_11 = first_order(1, 1, W234NF)
SS_W234_02 = first_order(0, 2, W234NF)

negativeRHS.actualcoord = Add(Mul(aNF[1], SS_W235_10), Mul(bNF[1], SS_W235_01),
                              Mul(aNF[2], SS_W234_10), Mul(bNF[2], SS_W234_01),
                              Mul(Pow(aNF[1], 2), SS_W234_20), Mul(aNF[1], bNF[1], SS_W234_11),
                              Mul(Pow(bNF[1], 2), SS_W234_02), Mul(2, DS_phiW125_00),
                              Mul(2, DS_W12W134_00), Mul(2, DS_W13W133_00),
                              Mul(2, aNF[1], DS_phiW134_10), Mul(2, bNF[1], DS_phiW134_01),
                              Mul(2, aNF[1], DS_W12W133_10), Mul(2, bNF[1], DS_W12W133_01),
                              Mul(2, aNF[2], DS_phiW133_10), Mul(2, bNF[2], DS_phiW133_01),
                              Mul(2, Pow(aNF[1], 2), DS_phiW133_20), Mul(2, aNF[1], bNF[1], DS_phiW133_11),
                              Mul(2, Pow(bNF[1], 2), DS_phiW133_02),
                              Mul(8, sqrt(muNF), diffmatrix, W24NF.dummy)).subs(extraparvals)

W256NF = linearsolver(W256NF, negativeRHS, coefmat2)

DS_phiW335_00 = second_order(0, 0, phiNF, W335NF)
DS_W02W234_00 = second_order(0, 0, W02NF, W234NF)
DS_W22W034_00 = second_order(0, 0, W22NF, W034NF)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW135_00), Mul(2, DS_phiW335_00),
                              Mul(4, DS_W02W234_00), Mul(2, DS_W22W034_00),
                              Mul(2, DS_W123W133_00), Mul(3, TS_phiphiW034_00),
                              Mul(6, TS_phiphiW234_00), Mul(12, TS_phiW133W02_00),
                              Mul(6, TS_phiW133W22_00), Mul(12, Q4S_phiphiphiW133_00),
                              Mul(12, sqrt(muNF), diffmatrix, W224NF.dummy)).subs(extraparvals)

W266NF = linearsolver(W266NF, negativeRHS, coefmat2)

TS_phiW22W133_00 = third_order(0, 0, phiNF, W22NF, W133NF)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW145_00), Mul(- 2, DS_W22W034_00),
                              Mul(- 2, DS_W133W33_00), Mul(- 3, TS_phiphiW034_00),
                              Mul(- 6, TS_phiW22W133_00), Mul(- 4, Q4S_phiphiphiW133_00),
                              Mul(4, sqrt(muNF), diffmatrix, W224NF.dummy)).subs(extraparvals)

W276NF = linearsolver(W276NF, negativeRHS, coefmat2)

negativeRHS.actualcoord = Add(Mul(2, DS_phiW15_00), Mul(- 4, sqrt(muNF), diffmatrix, W234NF.dummy),
                              Mul(2, diffmatrix, W22NF.dummy)).subs(extraparvals)

W286NF = linearsolver(W286NF, negativeRHS, coefmat2)

W06NF = evaluation(W06NF)
W026NF = evaluation(W026NF)
W036NF = evaluation(W036NF)
W046NF = evaluation(W046NF)
W056NF = evaluation(W056NF)
W066NF = evaluation(W066NF)
W076NF = evaluation(W076NF)
W16NF = evaluation(W16NF)
W126NF = evaluation(W126NF)
W136NF = evaluation(W136NF)
W146NF = evaluation(W146NF)
W156NF = evaluation(W156NF)
W166NF = evaluation(W166NF)
W176NF = evaluation(W176NF)
W26NF = evaluation(W26NF)
W226NF = evaluation(W226NF)
W236NF = evaluation(W236NF)
W246NF = evaluation(W246NF)
W256NF = evaluation(W256NF)
W266NF = evaluation(W266NF)
W276NF = evaluation(W276NF)
W286NF = evaluation(W286NF)

W06NF_eval = evaluation_dict(W06NF)
W026NF_eval = evaluation_dict(W026NF)
W036NF_eval = evaluation_dict(W036NF)
W046NF_eval = evaluation_dict(W046NF)
W056NF_eval = evaluation_dict(W056NF)
W066NF_eval = evaluation_dict(W066NF)
W076NF_eval = evaluation_dict(W076NF)
W16NF_eval = evaluation_dict(W16NF)
W126NF_eval = evaluation_dict(W126NF)
W136NF_eval = evaluation_dict(W136NF)
W146NF_eval = evaluation_dict(W146NF)
W156NF_eval = evaluation_dict(W156NF)
W166NF_eval = evaluation_dict(W166NF)
W176NF_eval = evaluation_dict(W176NF)
W26NF_eval = evaluation_dict(W26NF)
W226NF_eval = evaluation_dict(W226NF)
W236NF_eval = evaluation_dict(W236NF)
W246NF_eval = evaluation_dict(W246NF)
W256NF_eval = evaluation_dict(W256NF)
W266NF_eval = evaluation_dict(W266NF)
W276NF_eval = evaluation_dict(W276NF)
W286NF_eval = evaluation_dict(W286NF)

print('Sixth order ready')

alpha13 = psiNF.dummy.dot(Add(Mul(2, sqrt(muNF), diffmatrix, W15NF.dummy),
                              Mul(diffmatrix, W133NF.dummy))).subs(extraparvals)

SS_W16_10 = first_order(1, 0, W16NF)
SS_W16_01 = first_order(0, 1, W16NF)
SS_W15_20 = first_order(2, 0, W15NF)
SS_W15_11 = first_order(1, 1, W15NF)
SS_W15_02 = first_order(0, 2, W15NF)

alpha23 = psiNF.dummy.dot(Add(Mul(aNF[1], SS_W16_10), Mul(bNF[1], SS_W16_01),
                              Mul(aNF[2], SS_W15_10), Mul(bNF[2], SS_W15_01),
                              Mul(Pow(aNF[1], 2), SS_W15_20), Mul(aNF[1], bNF[1], SS_W15_11),
                              Mul(Pow(bNF[1], 2), SS_W15_02), Mul(- 2, sqrt(muNF), diffmatrix, W125NF.dummy),
                              Mul(diffmatrix, W13NF.dummy))).subs(extraparvals)

DS_phiW076_00 = second_order(0, 0, phiNF, W076NF)
DS_phiW286_00 = second_order(0, 0, phiNF, W286NF)
DS_W02W15_00 = second_order(0, 0, W02NF, W15NF)
TS_phiphiW15_00 = third_order(0, 0, phiNF, phiNF, W15NF)

alpha33 = psiNF.dummy.dot(Add(Mul(2, DS_phiW076_00), Mul(2, DS_phiW286_00),
                              Mul(4, DS_W02W15_00), Mul(6, TS_phiphiW15_00),
                              Mul(- 2, sqrt(muNF), diffmatrix, W135NF.dummy),
                              Mul(2, diffmatrix, W123NF.dummy))).subs(extraparvals)

DS_W22W15_00 = second_order(0, 0, W22NF, W15NF)
TS_phiphiW15_00 = third_order(0, 0, phiNF, phiNF, W15NF)

alpha43 = psiNF.dummy.dot(Add(Mul(2, DS_phiW076_00), Mul(2, DS_W22W15_00),
                              Mul(3, TS_phiphiW15_00), Mul(- 2, sqrt(muNF), diffmatrix, W145NF.dummy),
                              Mul(diffmatrix, W123NF.dummy))).subs(extraparvals)

SS_W126_10 = first_order(1, 0, W126NF)
SS_W126_01 = first_order(0, 1, W126NF)
SS_W125_20 = first_order(2, 0, W125NF)
SS_W125_11 = first_order(1, 1, W125NF)
SS_W125_02 = first_order(0, 2, W125NF)
SS_W134_30 = first_order(3, 0, W134NF)
SS_W134_21 = first_order(2, 1, W134NF)
SS_W134_12 = first_order(1, 2, W134NF)
SS_W134_03 = first_order(0, 3, W134NF)
SS_W133_40 = first_order(4, 0, W133NF)
SS_W133_31 = first_order(3, 1, W133NF)
SS_W133_22 = first_order(2, 2, W133NF)
SS_W133_13 = first_order(1, 3, W133NF)
SS_W133_04 = first_order(0, 4, W133NF)

alpha53 = psiNF.dummy.dot(Add(Mul(aNF[1], SS_W126_10), Mul(bNF[1], SS_W126_01),
                              Mul(aNF[2], SS_W125_10), Mul(bNF[2], SS_W125_01),
                              Mul(aNF[3], SS_W134_10), Mul(bNF[3], SS_W134_01),
                              Mul(aNF[4], SS_W133_10), Mul(bNF[4], SS_W133_01),
                              Mul(Pow(aNF[1], 2), SS_W125_20), Mul(aNF[1], bNF[1], SS_W125_11),
                              Mul(Pow(bNF[1], 2), SS_W125_02), Mul(2, aNF[1], aNF[2], SS_W134_20),
                              Mul(Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), SS_W134_11),
                              Mul(2, bNF[1], bNF[2], SS_W134_02), Mul(2, aNF[1], aNF[3], SS_W133_20),
                              Mul(Add(Mul(aNF[1], bNF[3]), Mul(aNF[3], bNF[1])), SS_W133_11),
                              Mul(2, bNF[1], bNF[3], SS_W133_02), Mul(Pow(aNF[2], 2), SS_W133_20),
                              Mul(aNF[2], bNF[2], SS_W133_11), Mul(Pow(bNF[2], 2), SS_W133_02),
                              Mul(Pow(aNF[1], 3), SS_W134_30), Mul(Pow(aNF[1], 2), bNF[1], SS_W134_21),
                              Mul(aNF[1], Pow(bNF[1], 2), SS_W134_12), Mul(Pow(bNF[1], 3), SS_W134_03),
                              Mul(3, Pow(aNF[1], 2), aNF[2], SS_W133_30),
                              Mul(2, aNF[1], aNF[2], bNF[1], SS_W133_21),
                              Mul(aNF[2], Pow(bNF[1], 2), SS_W133_12),
                              Mul(Pow(aNF[1], 2), bNF[2], SS_W133_21),
                              Mul(2, aNF[1], bNF[1], bNF[2], SS_W133_12),
                              Mul(3, Pow(bNF[1], 2), bNF[2], SS_W133_03),
                              Mul(Pow(aNF[1], 4), SS_W133_40),
                              Mul(Pow(aNF[1], 3), bNF[1], SS_W133_31),
                              Mul(Pow(aNF[1], 2), Pow(bNF[1], 2), SS_W133_22),
                              Mul(aNF[1], Pow(bNF[1], 3), SS_W133_13),
                              Mul(Pow(bNF[1], 4), SS_W133_04),
                              Mul(2, sqrt(muNF), diffmatrix, W155NF.dummy))).subs(extraparvals)

SS_W136_10 = first_order(1, 0, W136NF)
SS_W136_01 = first_order(0, 1, W136NF)
SS_W135_20 = first_order(2, 0, W135NF)
SS_W135_11 = first_order(1, 1, W135NF)
SS_W135_02 = first_order(0, 2, W135NF)
DS_phiW056_00 = second_order(0, 0, phiNF, W056NF)
DS_phiW256_00 = second_order(0, 0, phiNF, W256NF)
DS_W02W125_00 = second_order(0, 0, W02NF, W125NF)
DS_W12W035_00 = second_order(0, 0, W12NF, W035NF)
DS_W12W235_00 = second_order(0, 0, W12NF, W235NF)
DS_W03W134_00 = second_order(0, 0, W03NF, W134NF)
DS_W13W034_00 = second_order(0, 0, W13NF, W034NF)
DS_W13W234_00 = second_order(0, 0, W13NF, W234NF)
DS_W133W04_00 = second_order(0, 0, W133NF, W04NF)
DS_phiW035_10 = second_order(1, 0, phiNF, W035NF)
DS_phiW235_10 = second_order(1, 0, phiNF, W235NF)
DS_phiW035_01 = second_order(0, 1, phiNF, W035NF)
DS_phiW235_01 = second_order(0, 1, phiNF, W235NF)
DS_W02W134_10 = second_order(1, 0, W02NF, W134NF)
DS_W02W134_01 = second_order(0, 1, W02NF, W134NF)
DS_W12W034_10 = second_order(1, 0, W12NF, W034NF)
DS_W12W234_10 = second_order(1, 0, W12NF, W234NF)
DS_W12W034_01 = second_order(0, 1, W12NF, W034NF)
DS_W12W234_01 = second_order(0, 1, W12NF, W234NF)
DS_W03W133_10 = second_order(1, 0, W03NF, W133NF)
DS_W03W133_01 = second_order(0, 1, W03NF, W133NF)
DS_phiW034_20 = second_order(2, 0, phiNF, W034NF)
DS_phiW234_20 = second_order(2, 0, phiNF, W234NF)
DS_phiW034_11 = second_order(1, 1, phiNF, W034NF)
DS_phiW234_11 = second_order(1, 1, phiNF, W234NF)
DS_phiW034_02 = second_order(0, 2, phiNF, W034NF)
DS_phiW234_02 = second_order(0, 2, phiNF, W234NF)
DS_W02W133_20 = second_order(2, 0, W02NF, W133NF)
DS_W02W133_11 = second_order(1, 1, W02NF, W133NF)
DS_W02W133_02 = second_order(0, 2, W02NF, W133NF)
TS_phiphiW125_00 = third_order(0, 0, phiNF, phiNF, W125NF)
TS_phiW12W134_00 = third_order(0, 0, phiNF, W12NF, W134NF)
TS_phiW13W133_00 = third_order(0, 0, phiNF, W13NF, W133NF)
TS_W12W12W133_00 = third_order(0, 0, W12NF, W12NF, W133NF)
TS_phiphiW134_10 = third_order(1, 0, phiNF, phiNF, W134NF)
TS_phiphiW134_01 = third_order(0, 1, phiNF, phiNF, W134NF)
TS_phiW12W133_10 = third_order(1, 0, phiNF, W12NF, W133NF)
TS_phiW12W133_01 = third_order(0, 1, phiNF, W12NF, W133NF)
TS_phiphiW133_20 = third_order(2, 0, phiNF, phiNF, W133NF)
TS_phiphiW133_11 = third_order(1, 1, phiNF, phiNF, W133NF)
TS_phiphiW133_02 = third_order(0, 2, phiNF, phiNF, W133NF)

alpha63 = psiNF.dummy.dot(Add(Mul(aNF[1], SS_W136_10), Mul(bNF[1], SS_W136_01),
                          Mul(aNF[2], SS_W135_10), Mul(bNF[2], SS_W135_01),
                          Mul(Pow(aNF[1], 2), SS_W135_20), Mul(aNF[1], bNF[1], SS_W135_11),
                          Mul(Pow(bNF[1], 2), SS_W135_02), Mul(2, DS_phiW056_00),
                          Mul(2, DS_phiW256_00), Mul(4, DS_W02W125_00),
                          Mul(2, DS_W12W035_00), Mul(2, DS_W12W235_00),
                          Mul(4, DS_W03W134_00), Mul(2, DS_W13W034_00),
                          Mul(2, DS_W13W234_00), Mul(4, DS_W133W04_00),
                          Mul(2, aNF[1], DS_phiW035_10), Mul(2, aNF[1], DS_phiW235_10),
                          Mul(2, bNF[1], DS_phiW035_01), Mul(2, bNF[1], DS_phiW235_01),
                          Mul(4, aNF[1], DS_W02W134_10), Mul(4, bNF[1], DS_W02W134_01),
                          Mul(2, aNF[1], DS_W12W034_10), Mul(2, aNF[1], DS_W12W234_10),
                          Mul(2, bNF[1], DS_W12W034_01), Mul(2, bNF[1], DS_W12W234_01),
                          Mul(4, aNF[1], DS_W03W133_10), Mul(4, bNF[1], DS_W03W133_01),
                          Mul(2, aNF[2], DS_phiW034_10), Mul(2, aNF[2], DS_phiW234_10),
                          Mul(2, bNF[2], DS_phiW034_10), Mul(2, bNF[2], DS_phiW234_01),
                          Mul(4, aNF[2], DS_W02W133_10), Mul(4, bNF[2], DS_W02W133_01),
                          Mul(2, Pow(aNF[1], 2), DS_phiW034_20),
                          Mul(2, Pow(aNF[1], 2), DS_phiW234_20),
                          Mul(2, aNF[1], bNF[1], DS_phiW034_11),
                          Mul(2, aNF[1], bNF[1], DS_phiW234_11),
                          Mul(2, Pow(bNF[1], 2), DS_phiW034_02),
                          Mul(2, Pow(bNF[1], 2), DS_phiW234_02),
                          Mul(4, Pow(aNF[1], 2), DS_W02W133_20),
                          Mul(4, aNF[1], bNF[1], DS_W02W133_11),
                          Mul(4, Pow(bNF[1], 2), DS_W02W133_02), Mul(6, TS_phiphiW125_00),
                          Mul(12, TS_phiW12W134_00), Mul(12, TS_phiW13W133_00),
                          Mul(6, TS_W12W12W133_00), Mul(6, aNF[1], TS_phiphiW134_10),
                          Mul(6, bNF[1], TS_phiphiW134_01), Mul(12, aNF[1], TS_phiW12W133_10),
                          Mul(12, bNF[1], TS_phiW12W133_01), Mul(6, aNF[2], TS_phiphiW133_10),
                          Mul(6, bNF[2], TS_phiphiW133_01), Mul(6, Pow(aNF[1], 2), TS_phiphiW133_20),
                          Mul(6, aNF[1], bNF[1], TS_phiphiW133_11), Mul(6, Pow(bNF[1], 2), TS_phiphiW133_02),
                          Mul(4, sqrt(muNF), diffmatrix, W165NF.dummy))).subs(extraparvals)

DS_phiW066_00 = second_order(0, 0, phiNF, W066NF)
DS_phiW266_00 = second_order(0, 0, phiNF, W266NF)
DS_W02W135_00 = second_order(0, 0, W02NF, W135NF)
DS_W22W145_00 = second_order(0, 0, W22NF, W145NF)
DS_W22W335_00 = second_order(0, 0, W22NF, W335NF)
DS_W123W034_00 = second_order(0, 0, W123NF, W034NF)
DS_W123W234_00 = second_order(0, 0, W123NF, W234NF)
DS_W133W024_00 = second_order(0, 0, W133NF, W024NF)
TS_phiphiW135_00 = third_order(0, 0, phiNF, phiNF, W135NF)
TS_phiphiW145_00 = third_order(0, 0, phiNF, phiNF, W145NF)
TS_phiphiW335_00 = third_order(0, 0, phiNF, phiNF, W335NF)
TS_phiW02W034_00 = third_order(0, 0, phiNF, W02NF, W034NF)
TS_phiW02W234_00 = third_order(0, 0, phiNF, W02NF, W234NF)
TS_phiW22W034_00 = third_order(0, 0, phiNF, W22NF, W034NF)
TS_phiW22W234_00 = third_order(0, 0, phiNF, W22NF, W234NF)
TS_phiW123W133_00 = third_order(0, 0, phiNF, W123NF, W133NF)
TS_W02W02W133_00 = third_order(0, 0, W02NF, W02NF, W133NF)
TS_W22W22W133_00 = third_order(0, 0, W22NF, W22NF, W133NF)
Q4S_phiphiphiW034_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W034NF)
Q4S_phiphiphiW234_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W234NF)
Q4S_phiphiW133W02_00 = fourth_order(0, 0, phiNF, phiNF, W133NF, W02NF)
Q4S_phiphiW133W22_00 = fourth_order(0, 0, phiNF, phiNF, W133NF, W22NF)
Q5S_phiphiphiphiW133_00 = fifth_order(0, 0, phiNF, phiNF, phiNF, phiNF, W133NF)

alpha73 = psiNF.dummy.dot(Add(Mul(2, DS_phiW066_00), Mul(2, DS_phiW266_00),
                              Mul(4, DS_W02W135_00), Mul(- 2, DS_W22W145_00),
                              Mul(2, DS_W22W335_00), Mul(2, DS_W123W034_00),
                              Mul(2, DS_W123W234_00), Mul(4, DS_W133W024_00),
                              Mul(6, TS_phiphiW135_00), Mul(- 3, TS_phiphiW145_00),
                              Mul(3, TS_phiphiW335_00), Mul(12, TS_phiW02W034_00),
                              Mul(12, TS_phiW02W234_00), Mul(6, TS_phiW22W034_00),
                              Mul(6, TS_phiW22W234_00), Mul(12, TS_phiW123W133_00),
                              Mul(12, TS_W02W02W133_00), Mul(6, TS_W22W22W133_00),
                              Mul(12, Q4S_phiphiphiW034_00), Mul(12, Q4S_phiphiphiW234_00),
                              Mul(48, Q4S_phiphiW133W02_00), Mul(24, Q4S_phiphiW133W22_00),
                              Mul(30, Q5S_phiphiphiphiW133_00),
                              Mul(6, sqrt(muNF), diffmatrix, W175NF.dummy))).subs(extraparvals)

SS_W146_10 = first_order(1, 0, W146NF)
SS_W146_01 = first_order(0, 1, W146NF)
SS_W145_20 = first_order(2, 0, W145NF)
SS_W145_11 = first_order(1, 1, W145NF)
SS_W145_02 = first_order(0, 2, W145NF)
DS_W12W135_00 = second_order(0, 0, W12NF, W135NF)
DS_W22W125_00 = second_order(0, 0, W22NF, W125NF)
DS_W23W134_00 = second_order(0, 0, W23NF, W134NF)
DS_W24W133_00 = second_order(0, 0, W24NF, W133NF)
DS_W22W134_10 = second_order(1, 0, W22NF, W134NF)
DS_W22W134_01 = second_order(0, 1, W22NF, W134NF)
DS_W133W23_10 = second_order(1, 0, W133NF, W23NF)
DS_W133W23_01 = second_order(0, 1, W133NF, W23NF)
DS_W22W133_20 = second_order(2, 0, W22NF, W133NF)
DS_W22W133_11 = second_order(1, 1, W22NF, W133NF)
DS_W22W133_02 = second_order(0, 2, W22NF, W133NF)

alpha83 = psiNF.dummy.dot(Add(Mul(aNF[1], SS_W146_10), Mul(bNF[1], SS_W146_01),
                              Mul(aNF[2], SS_W145_10), Mul(bNF[2], SS_W145_01),
                              Mul(Pow(aNF[1], 2), SS_W145_20), Mul(aNF[1], bNF[1], SS_W145_11),
                              Mul(Pow(bNF[1], 2), SS_W145_02), Mul(- 2, DS_phiW056_00),
                              Mul(- 2, DS_W12W135_00), Mul(- 2, DS_W22W125_00),
                              Mul(- 2, DS_W13W034_00), Mul(- 2, DS_W23W134_00),
                              Mul(- 2, DS_W24W133_00), Mul(- 2, aNF[1], DS_phiW035_10),
                              Mul(- 2, bNF[1], DS_phiW035_01), Mul(- 2, aNF[1], DS_W12W034_10),
                              Mul(- 2, bNF[1], DS_W12W034_01), Mul(- 2, aNF[1], DS_W22W134_10),
                              Mul(- 2, bNF[1], DS_W22W134_01), Mul(- 2, aNF[1], DS_W133W23_10),
                              Mul(- 2, bNF[1], DS_W133W23_01), Mul(- 2, aNF[2], DS_phiW034_10),
                              Mul(- 2, bNF[2], DS_phiW034_01), Mul(- 2, aNF[2], DS_W22W133_10),
                              Mul(- 2, bNF[2], DS_W22W133_01),
                              Mul(- 2, Pow(aNF[1], 2), DS_phiW034_20),
                              Mul(- 2, aNF[1], bNF[1], DS_phiW034_11),
                              Mul(- 2, Pow(bNF[1], 2), DS_phiW034_02),
                              Mul(- 2, Pow(aNF[1], 2), DS_W22W133_20),
                              Mul(- 2, aNF[1], bNF[1], DS_W22W133_11),
                              Mul(- 2, Pow(bNF[1], 2), DS_W22W133_02), Mul(- 3, TS_phiphiW125_00),
                              Mul(- 6, TS_phiW12W134_00), Mul(- 6, TS_phiW13W133_00),
                              Mul(- 3, TS_W12W12W133_00), Mul(- 3, aNF[1], TS_phiphiW134_10),
                              Mul(- 3, bNF[1], TS_phiphiW134_01), Mul(- 6, aNF[1], TS_phiW12W133_10),
                              Mul(- 6, bNF[1], TS_phiW12W133_01), Mul(- 3, aNF[2], TS_phiphiW133_10),
                              Mul(- 3, bNF[2], TS_phiphiW133_01),
                              Mul(- 3, Pow(aNF[1], 2), TS_phiphiW133_20),
                              Mul(- 3, aNF[1], bNF[1], TS_phiphiW133_11),
                              Mul(- 3, Pow(bNF[1], 2), TS_phiphiW133_02),
                              Mul(2, sqrt(muNF), diffmatrix, W165NF.dummy))).subs(extraparvals)

DS_phiW276_00 = second_order(0, 0, phiNF, W276NF)
DS_W02W145_00 = second_order(0, 0, W02NF, W145NF)
DS_W22W135_00 = second_order(0, 0, W22NF, W135NF)
DS_W133W224_00 = second_order(0, 0, W133NF, W224NF)
DS_W33W234_00 = second_order(0, 0, W33NF, W234NF)
TS_phiW133W123_00 = third_order(0, 0, phiNF, W133NF, W123NF)
TS_phiW133W33_00 = third_order(0, 0, phiNF, W133NF, W33NF)
TS_W02W22W133_00 = third_order(0, 0, W02NF, W22NF, W133NF)

alpha93 = psiNF.dummy.dot(Add(Mul(- 2, DS_phiW066_00), Mul(2, DS_phiW276_00), Mul(4, DS_W02W145_00),
                              Mul(- 2, DS_W22W135_00), Mul(- 2, DS_W123W034_00),
                              Mul(- 2, DS_W133W224_00), Mul(- 2, DS_W33W234_00),
                              Mul(- 3, TS_phiphiW135_00), Mul(6, TS_phiphiW145_00),
                              Mul(- 12, TS_phiW02W034_00), Mul(- 6, TS_phiW22W034_00),
                              Mul(- 6, TS_phiW22W234_00), Mul(- 6, TS_phiW133W123_00),
                              Mul(- 6, TS_phiW133W33_00), Mul(- 12, TS_W02W22W133_00),
                              Mul(- 12, Q4S_phiphiphiW034_00), Mul(- 4, Q4S_phiphiphiW234_00),
                              Mul(- 24, Q4S_phiphiW133W02_00), Mul(- 24, Q4S_phiphiW133W22_00),
                              Mul(- 20, Q5S_phiphiphiphiW133_00),
                              Mul(4, sqrt(muNF), diffmatrix, W175NF.dummy))).subs(extraparvals)

DS_phiW246_00 = second_order(0, 0, phiNF, W246NF)
DS_W133W034_00 = second_order(0, 0, W133NF, W034NF)
TS_phiW133W133_00 = third_order(0, 0, phiNF, W133NF, W133NF)

alpha103 = psiNF.dummy.dot(Add(Mul(2, DS_phiW246_00), Mul(- 2, DS_W133W034_00),
                           Mul(- 3, TS_phiW133W133_00),
                           Mul(- 2, sqrt(muNF), diffmatrix, W135NF.dummy),
                           Mul(2, diffmatrix, W123NF.dummy))).subs(extraparvals)

DS_phiW046_00 = second_order(0, 0, phiNF, W046NF)
DS_W133W234_00 = second_order(0, 0, W133NF, W234NF)

alpha113 = psiNF.dummy.dot(Add(Mul(4, DS_phiW046_00), Mul(2, DS_W133W034_00),
                               Mul(2, DS_W133W234_00), Mul(6, TS_phiW133W133_00),
                               Mul(- 2, sqrt(muNF), diffmatrix, W135NF.dummy),
                               Mul(- 4, sqrt(muNF), diffmatrix, W145NF.dummy),
                               Mul(4, diffmatrix, W123NF.dummy))).subs(extraparvals)

SS_W156_10 = first_order(1, 0, W156NF)
SS_W156_01 = first_order(0, 1, W156NF)
SS_W155_20 = first_order(2, 0, W155NF)
SS_W155_11 = first_order(1, 1, W155NF)
SS_W155_02 = first_order(0, 2, W155NF)
SS_W14_30 = first_order(3, 0, W14NF)
SS_W14_21 = first_order(2, 1, W14NF)
SS_W14_12 = first_order(1, 2, W14NF)
SS_W14_03 = first_order(0, 3, W14NF)
SS_W13_40 = first_order(4, 0, W13NF)
SS_W13_31 = first_order(3, 1, W13NF)
SS_W13_22 = first_order(2, 2, W13NF)
SS_W13_13 = first_order(1, 3, W13NF)
SS_W13_04 = first_order(0, 4, W13NF)
SS_W12_50 = first_order(5, 0, W12NF)
SS_W12_41 = first_order(4, 1, W12NF)
SS_W12_32 = first_order(3, 2, W12NF)
SS_W12_23 = first_order(2, 3, W12NF)
SS_W12_14 = first_order(1, 4, W12NF)
SS_W12_05 = first_order(0, 5, W12NF)
SS_phi_60 = first_order(6, 0, phiNF)
SS_phi_51 = first_order(5, 1, phiNF)
SS_phi_42 = first_order(4, 2, phiNF)
SS_phi_33 = first_order(3, 3, phiNF)
SS_phi_24 = first_order(2, 4, phiNF)
SS_phi_15 = first_order(1, 5, phiNF)
SS_phi_06 = first_order(0, 6, phiNF)

alpha123 = psiNF.dummy.dot(Add(Mul(aNF[1], SS_W156_10), Mul(bNF[1], SS_W156_01),
                               Mul(aNF[2], SS_W155_10), Mul(bNF[2], SS_W155_01),
                               Mul(aNF[3], SS_W14_10), Mul(bNF[3], SS_W14_01),
                               Mul(aNF[4], SS_W13_10), Mul(bNF[4], SS_W13_01),
                               Mul(aNF[5], SS_W12_10), Mul(bNF[5], SS_W12_01),
                               Mul(aNF[6], SS_phi_10), Mul(bNF[6], SS_phi_01),
                               Mul(Pow(aNF[1], 2), SS_W155_20), Mul(aNF[1], bNF[1], SS_W155_11),
                               Mul(Pow(bNF[1], 2), SS_W155_02), Mul(2, aNF[1], aNF[2], SS_W14_20),
                               Mul(Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), SS_W14_11),
                               Mul(2, bNF[1], bNF[2], SS_W14_02), Mul(2, aNF[1], aNF[3], SS_W13_20),
                               Mul(Add(Mul(aNF[1], bNF[3]), Mul(aNF[3], bNF[1])), SS_W13_11),
                               Mul(2, bNF[1], bNF[3], SS_W13_02), Mul(2, aNF[1], aNF[4], SS_W12_20),
                               Mul(Add(Mul(aNF[1], bNF[4]), Mul(aNF[4], bNF[1])), SS_W12_11),
                               Mul(2, bNF[1], bNF[4], SS_W12_02), Mul(2, aNF[1], aNF[5], SS_phi_20),
                               Mul(Add(Mul(aNF[1], bNF[5]), Mul(aNF[5], bNF[1])), SS_phi_11),
                               Mul(2, bNF[1], bNF[5], SS_phi_02), Mul(Pow(aNF[2], 2), SS_W13_20),
                               Mul(aNF[2], bNF[2], SS_W13_11), Mul(Pow(bNF[2], 2), SS_W13_02),
                               Mul(2, aNF[2], aNF[3], SS_W12_20),
                               Mul(Add(Mul(aNF[2], bNF[3]), Mul(aNF[3], bNF[2])), SS_W12_11),
                               Mul(2, bNF[2], bNF[3], SS_W12_02), Mul(2, aNF[2], aNF[4], SS_phi_20),
                               Mul(Add(Mul(aNF[2], bNF[4]), Mul(aNF[4], bNF[2])), SS_phi_11),
                               Mul(2, bNF[2], bNF[4], SS_phi_02), Mul(Pow(aNF[3], 2), SS_phi_20),
                               Mul(aNF[3], bNF[3], SS_phi_11), Mul(Pow(bNF[3], 2), SS_phi_02),
                               Mul(Pow(aNF[1], 3), SS_W14_30), Mul(Pow(aNF[1], 2), bNF[1], SS_W14_21),
                               Mul(aNF[1], Pow(bNF[1], 2), SS_W14_12), Mul(Pow(bNF[1], 3), SS_W14_03),
                               Mul(3, Pow(aNF[1], 2), aNF[2], SS_W13_30),
                               Mul(2, aNF[1], aNF[2], bNF[1], SS_W13_21),
                               Mul(aNF[2], Pow(bNF[1], 2), SS_W13_12),
                               Mul(Pow(aNF[1], 2), bNF[2], SS_W13_21),
                               Mul(2, aNF[1], bNF[1], bNF[2], SS_W13_12),
                               Mul(3, Pow(bNF[1], 2), bNF[2], SS_W13_03),
                               Mul(3, Pow(aNF[1], 2), aNF[3], SS_W12_30),
                               Mul(2, aNF[1], aNF[3], bNF[1], SS_W12_21),
                               Mul(aNF[3], Pow(bNF[1], 2), SS_W12_12),
                               Mul(Pow(aNF[1], 2), bNF[3], SS_W12_21),
                               Mul(2, aNF[1], bNF[1], bNF[3], SS_W12_12),
                               Mul(3, Pow(bNF[1], 2), bNF[3], SS_W12_03),
                               Mul(3, Pow(aNF[1], 2), aNF[4], SS_phi_30),
                               Mul(2, aNF[1], aNF[4], bNF[1], SS_phi_21),
                               Mul(aNF[4], Pow(bNF[1], 2), SS_phi_12),
                               Mul(Pow(aNF[1], 2), bNF[4], SS_phi_21),
                               Mul(2, aNF[1], bNF[1], bNF[4], SS_phi_12),
                               Mul(3, Pow(bNF[1], 2), bNF[4], SS_phi_03),
                               Mul(3, aNF[1], Pow(aNF[2], 2), SS_W12_30),
                               Mul(2, aNF[1], aNF[2], bNF[2], SS_W12_21),
                               Mul(aNF[1], Pow(bNF[2], 2), SS_W12_12),
                               Mul(Pow(aNF[2], 2), bNF[1], SS_W12_21),
                               Mul(2, aNF[2], bNF[1], bNF[2], SS_W12_12),
                               Mul(3, bNF[1], Pow(bNF[2], 2), SS_W12_03),
                               Mul(Pow(aNF[2], 3), SS_phi_30), Mul(Pow(aNF[2], 2), bNF[2], SS_phi_21),
                               Mul(aNF[2], Pow(bNF[2], 2), SS_phi_12), Mul(Pow(bNF[2], 3), SS_phi_03),
                               Mul(6, aNF[1], aNF[2], aNF[3], SS_phi_30),
                               Mul(2, Add(Mul(aNF[2], aNF[3], bNF[1]), Mul(aNF[1], aNF[3], bNF[2]),
                                          Mul(aNF[1], aNF[2], bNF[3])), SS_phi_21),
                               Mul(2, Add(Mul(aNF[1], bNF[2], bNF[3]), Mul(aNF[2], bNF[1], bNF[3]),
                                          Mul(aNF[3], bNF[1], bNF[2])), SS_phi_12),
                               Mul(6, bNF[1], bNF[2], bNF[3], SS_phi_03),
                               Mul(Pow(aNF[1], 4), SS_W13_40),
                               Mul(Pow(aNF[1], 3), bNF[1], SS_W13_31),
                               Mul(Pow(aNF[1], 2), Pow(bNF[1], 2), SS_W13_22),
                               Mul(aNF[1], Pow(bNF[1], 3), SS_W13_13),
                               Mul(Pow(bNF[1], 4), SS_W13_04),
                               Mul(4, Pow(aNF[1], 3), aNF[2], SS_W12_40),
                               Mul(3, Pow(aNF[1], 2), aNF[2], bNF[1], SS_W12_31),
                               Mul(2, aNF[1], aNF[2], Pow(bNF[1], 2), SS_W12_22),
                               Mul(aNF[2], Pow(bNF[1], 3), SS_W12_13),
                               Mul(Pow(aNF[1], 3), bNF[2], SS_W12_31),
                               Mul(2, Pow(aNF[1], 2), bNF[1], bNF[2], SS_W12_22),
                               Mul(3, aNF[1], Pow(bNF[1], 2), bNF[2], SS_W12_13),
                               Mul(4, Pow(bNF[1], 3), bNF[2], SS_W12_04),
                               Mul(4, Pow(aNF[1], 3), aNF[3], SS_phi_40),
                               Mul(3, Pow(aNF[1], 2), aNF[3], bNF[1], SS_phi_31),
                               Mul(2, aNF[1], aNF[3], Pow(bNF[1], 2), SS_phi_22),
                               Mul(aNF[3], Pow(bNF[1], 3), SS_phi_13),
                               Mul(Pow(aNF[1], 3), bNF[3], SS_phi_31),
                               Mul(2, Pow(aNF[1], 2), bNF[1], bNF[3], SS_phi_22),
                               Mul(3, aNF[1], Pow(bNF[1], 2), bNF[3], SS_phi_13),
                               Mul(4, Pow(bNF[1], 3), bNF[3], SS_phi_04),
                               Mul(6, Pow(aNF[1], 2), Pow(aNF[2], 2), SS_phi_40),
                               Mul(3, aNF[1], Pow(aNF[2], 2), bNF[1], SS_phi_31),
                               Mul(Pow(aNF[2], 2), Pow(bNF[1], 2), SS_phi_22),
                               Mul(3, Pow(aNF[1], 2), aNF[2], bNF[2], SS_phi_31),
                               Mul(4, aNF[1], aNF[2], bNF[1], bNF[2], SS_phi_22),
                               Mul(3, aNF[2], Pow(bNF[1], 2), bNF[2], SS_phi_13),
                               Mul(Pow(aNF[1], 2), Pow(bNF[2], 2), SS_phi_22),
                               Mul(3, aNF[1], bNF[1], Pow(bNF[2], 2), SS_phi_13),
                               Mul(6, Pow(bNF[1], 2), Pow(bNF[2], 2), SS_phi_04),
                               Mul(Pow(aNF[1], 5), SS_W12_50), Mul(Pow(aNF[1], 4), bNF[1], SS_W12_41),
                               Mul(Pow(aNF[1], 3), Pow(bNF[1], 2), SS_W12_32),
                               Mul(Pow(aNF[1], 2), Pow(bNF[1], 3), SS_W12_23),
                               Mul(aNF[1], Pow(bNF[1], 4), SS_W12_14), Mul(Pow(bNF[1], 5), SS_W12_05),
                               Mul(5, Pow(aNF[1], 4), aNF[2], SS_phi_50),
                               Mul(4, Pow(aNF[1], 3), aNF[2], bNF[1], SS_phi_41),
                               Mul(3, Pow(aNF[1], 2), aNF[2], Pow(bNF[1], 2), SS_phi_32),
                               Mul(2, aNF[1], aNF[2], Pow(bNF[1], 3), SS_phi_23),
                               Mul(aNF[2], Pow(bNF[1], 4), SS_phi_14),
                               Mul(Pow(aNF[1], 4), bNF[2], SS_phi_41),
                               Mul(2, Pow(aNF[1], 3), bNF[1], bNF[2], SS_phi_32),
                               Mul(3, Pow(aNF[1], 2), Pow(bNF[1], 2), bNF[2], SS_phi_23),
                               Mul(4, aNF[1], Pow(bNF[1], 3), bNF[2], SS_phi_14),
                               Mul(5, Pow(bNF[1], 4), bNF[2], SS_phi_05),
                               Mul(Pow(aNF[1], 6), SS_phi_60), Mul(Pow(aNF[1], 5), bNF[1], SS_phi_51),
                               Mul(Pow(aNF[1], 4), Pow(bNF[1], 2), SS_phi_42),
                               Mul(Pow(aNF[1], 3), Pow(bNF[1], 3), SS_phi_33),
                               Mul(Pow(aNF[1], 2), Pow(bNF[1], 4), SS_phi_24),
                               Mul(aNF[1], Pow(bNF[1], 5), SS_phi_15),
                               Mul(Pow(bNF[1], 6), SS_phi_06))).subs(extraparvals)

SS_W166_10 = first_order(1, 0, W166NF)
SS_W166_01 = first_order(0, 1, W166NF)
SS_W165_20 = first_order(2, 0, W165NF)
SS_W165_11 = first_order(1, 1, W165NF)
SS_W165_02 = first_order(0, 2, W165NF)
SS_W124_30 = first_order(3, 0, W124NF)
SS_W124_21 = first_order(2, 1, W124NF)
SS_W124_12 = first_order(1, 2, W124NF)
SS_W124_03 = first_order(0, 3, W124NF)
SS_W123_40 = first_order(4, 0, W123NF)
SS_W123_31 = first_order(3, 1, W123NF)
SS_W123_22 = first_order(2, 2, W123NF)
SS_W123_13 = first_order(1, 3, W123NF)
SS_W123_04 = first_order(0, 4, W123NF)
DS_phiW06_00 = second_order(0, 0, phiNF, W06NF)
DS_phiW26_00 = second_order(0, 0, phiNF, W26NF)
DS_W12W05_00 = second_order(0, 0, W12NF, W05NF)
DS_W12W25_00 = second_order(0, 0, W12NF, W25NF)
DS_W13W04_00 = second_order(0, 0, W13NF, W04NF)
DS_W13W24_00 = second_order(0, 0, W13NF, W24NF)
DS_W14W03_00 = second_order(0, 0, W14NF, W03NF)
DS_W14W23_00 = second_order(0, 0, W14NF, W23NF)
DS_W155W02_00 = second_order(0, 0, W155NF, W02NF)
DS_W155W22_00 = second_order(0, 0, W155NF, W22NF)
DS_phiW05_10 = second_order(1, 0, phiNF, W05NF)
DS_phiW25_10 = second_order(1, 0, phiNF, W25NF)
DS_phiW05_01 = second_order(0, 1, phiNF, W05NF)
DS_phiW25_01 = second_order(0, 1, phiNF, W25NF)
DS_W12W04_10 = second_order(1, 0, W12NF, W04NF)
DS_W12W24_10 = second_order(1, 0, W12NF, W24NF)
DS_W12W04_01 = second_order(0, 1, W12NF, W04NF)
DS_W12W24_01 = second_order(0, 1, W12NF, W24NF)
DS_W13W03_10 = second_order(1, 0, W13NF, W03NF)
DS_W13W23_10 = second_order(1, 0, W13NF, W23NF)
DS_W13W03_01 = second_order(0, 1, W13NF, W03NF)
DS_W13W23_01 = second_order(0, 1, W13NF, W23NF)
DS_W14W02_10 = second_order(1, 0, W14NF, W02NF)
DS_W14W22_10 = second_order(1, 0, W14NF, W22NF)
DS_W14W02_01 = second_order(0, 1, W14NF, W02NF)
DS_W14W22_01 = second_order(0, 1, W14NF, W22NF)
DS_phiW04_20 = second_order(2, 0, phiNF, W04NF)
DS_phiW24_20 = second_order(2, 0, phiNF, W24NF)
DS_phiW04_11 = second_order(1, 1, phiNF, W04NF)
DS_phiW24_11 = second_order(1, 1, phiNF, W24NF)
DS_phiW04_02 = second_order(0, 2, phiNF, W04NF)
DS_phiW24_02 = second_order(0, 2, phiNF, W24NF)
DS_W12W03_20 = second_order(2, 0, W12NF, W03NF)
DS_W12W23_20 = second_order(2, 0, W12NF, W23NF)
DS_W12W03_11 = second_order(1, 1, W12NF, W03NF)
DS_W12W23_11 = second_order(1, 1, W12NF, W23NF)
DS_W12W03_02 = second_order(0, 2, W12NF, W03NF)
DS_W12W23_02 = second_order(0, 2, W12NF, W23NF)
DS_W13W02_20 = second_order(2, 0, W13NF, W02NF)
DS_W13W22_20 = second_order(2, 0, W13NF, W22NF)
DS_W13W02_11 = second_order(1, 1, W13NF, W02NF)
DS_W13W22_11 = second_order(1, 1, W13NF, W22NF)
DS_W13W02_02 = second_order(0, 2, W13NF, W02NF)
DS_W13W22_02 = second_order(0, 2, W13NF, W22NF)
DS_phiW03_30 = second_order(3, 0, phiNF, W03NF)
DS_phiW23_30 = second_order(3, 0, phiNF, W23NF)
DS_phiW03_21 = second_order(2, 1, phiNF, W03NF)
DS_phiW23_21 = second_order(2, 1, phiNF, W23NF)
DS_phiW03_12 = second_order(1, 2, phiNF, W03NF)
DS_phiW23_12 = second_order(1, 2, phiNF, W23NF)
DS_phiW03_03 = second_order(0, 3, phiNF, W03NF)
DS_phiW23_03 = second_order(0, 3, phiNF, W23NF)
DS_W12W02_30 = second_order(3, 0, W12NF, W02NF)
DS_W12W22_30 = second_order(3, 0, W12NF, W22NF)
DS_W12W02_21 = second_order(2, 1, W12NF, W02NF)
DS_W12W22_21 = second_order(2, 1, W12NF, W22NF)
DS_W12W02_12 = second_order(1, 2, W12NF, W02NF)
DS_W12W22_12 = second_order(1, 2, W12NF, W22NF)
DS_W12W02_03 = second_order(0, 3, W12NF, W02NF)
DS_W12W22_03 = second_order(0, 3, W12NF, W22NF)
DS_phiW02_40 = second_order(4, 0, phiNF, W02NF)
DS_phiW22_40 = second_order(4, 0, phiNF, W22NF)
DS_phiW02_31 = second_order(3, 1, phiNF, W02NF)
DS_phiW22_31 = second_order(3, 1, phiNF, W22NF)
DS_phiW02_22 = second_order(2, 2, phiNF, W02NF)
DS_phiW22_22 = second_order(2, 2, phiNF, W22NF)
DS_phiW02_13 = second_order(1, 3, phiNF, W02NF)
DS_phiW22_13 = second_order(1, 3, phiNF, W22NF)
DS_phiW02_04 = second_order(0, 4, phiNF, W02NF)
DS_phiW22_04 = second_order(0, 4, phiNF, W22NF)
TS_phiphiW155_00 = third_order(0, 0, phiNF, phiNF, W155NF)
TS_phiW12W14_00 = third_order(0, 0, phiNF, W12NF, W14NF)
TS_phiW13W13_00 = third_order(0, 0, phiNF, W13NF, W13NF)
TS_W12W12W13_00 = third_order(0, 0, W12NF, W12NF, W13NF)
TS_phiphiW14_10 = third_order(1, 0, phiNF, phiNF, W14NF)
TS_phiphiW14_01 = third_order(0, 1, phiNF, phiNF, W14NF)
TS_phiW12W13_10 = third_order(1, 0, phiNF, W12NF, W13NF)
TS_phiW12W13_01 = third_order(0, 1, phiNF, W12NF, W13NF)
TS_W12W12W12_10 = third_order(1, 0, W12NF, W12NF, W12NF)
TS_W12W12W12_01 = third_order(0, 1, W12NF, W12NF, W12NF)
TS_phiphiW13_20 = third_order(2, 0, phiNF, phiNF, W13NF)
TS_phiphiW13_11 = third_order(1, 1, phiNF, phiNF, W13NF)
TS_phiphiW13_02 = third_order(0, 2, phiNF, phiNF, W13NF)
TS_phiW12W12_20 = third_order(2, 0, phiNF, W12NF, W12NF)
TS_phiW12W12_11 = third_order(1, 1, phiNF, W12NF, W12NF)
TS_phiW12W12_02 = third_order(0, 2, phiNF, W12NF, W12NF)
TS_phiphiW12_30 = third_order(3, 0, phiNF, phiNF, W12NF)
TS_phiphiW12_21 = third_order(2, 1, phiNF, phiNF, W12NF)
TS_phiphiW12_12 = third_order(1, 2, phiNF, phiNF, W12NF)
TS_phiphiW12_03 = third_order(0, 3, phiNF, phiNF, W12NF)
TS_phiphiphi_40 = third_order(4, 0, phiNF, phiNF, phiNF)
TS_phiphiphi_31 = third_order(3, 1, phiNF, phiNF, phiNF)
TS_phiphiphi_22 = third_order(2, 2, phiNF, phiNF, phiNF)
TS_phiphiphi_13 = third_order(1, 3, phiNF, phiNF, phiNF)
TS_phiphiphi_04 = third_order(0, 4, phiNF, phiNF, phiNF)

alpha133 = psiNF.dummy.dot(Add(Mul(aNF[1], SS_W166_10), Mul(bNF[1], SS_W166_01),
                               Mul(aNF[2], SS_W165_10), Mul(bNF[2], SS_W165_01),
                               Mul(aNF[3], SS_W124_10), Mul(bNF[3], SS_W124_01),
                               Mul(aNF[4], SS_W123_10), Mul(bNF[4], SS_W123_01),
                               Mul(Pow(aNF[1], 2), SS_W165_20), Mul(aNF[1], bNF[1], SS_W165_11),
                               Mul(Pow(bNF[1], 2), SS_W165_02), Mul(2, aNF[1], aNF[2], SS_W124_20),
                               Mul(Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), SS_W124_11),
                               Mul(2, bNF[1], bNF[2], SS_W124_02), Mul(2, aNF[1], aNF[3], SS_W123_20),
                               Mul(Add(Mul(aNF[1], bNF[3]), Mul(aNF[3], bNF[1])), SS_W123_11),
                               Mul(2, bNF[1], bNF[3], SS_W123_02), Mul(Pow(aNF[2], 2), SS_W123_20),
                               Mul(aNF[2], bNF[2], SS_W123_11), Mul(Pow(bNF[2], 2), SS_W123_02),
                               Mul(Pow(aNF[1], 3), SS_W124_30),
                               Mul(Pow(aNF[1], 2), bNF[1], SS_W124_21),
                               Mul(aNF[1], Pow(bNF[1], 2), SS_W124_12),
                               Mul(Pow(bNF[1], 3), SS_W124_03),
                               Mul(3, Pow(aNF[1], 2), aNF[2], SS_W123_30),
                               Mul(2, aNF[1], aNF[2], bNF[1], SS_W123_21),
                               Mul(aNF[2], Pow(bNF[1], 2), SS_W123_12),
                               Mul(Pow(aNF[1], 2), bNF[2], SS_W123_21),
                               Mul(2, aNF[1], bNF[1], bNF[2], SS_W123_12),
                               Mul(3, Pow(bNF[1], 2), bNF[2], SS_W123_03),
                               Mul(Pow(aNF[1], 4), SS_W123_40),
                               Mul(Pow(aNF[1], 3), bNF[1], SS_W123_31),
                               Mul(Pow(aNF[1], 2), Pow(bNF[1], 2), SS_W123_22),
                               Mul(aNF[1], Pow(bNF[1], 3), SS_W123_13),
                               Mul(Pow(bNF[1], 4), SS_W123_04),
                               Mul(4, DS_phiW06_00), Mul(2, DS_phiW26_00),
                               Mul(4, DS_W12W05_00), Mul(2, DS_W12W25_00),
                               Mul(4, DS_W13W04_00), Mul(2, DS_W13W24_00),
                               Mul(4, DS_W14W03_00), Mul(2, DS_W14W23_00),
                               Mul(4, DS_W155W02_00), Mul(2, DS_W155W22_00),
                               Mul(4, aNF[1], DS_phiW05_10), Mul(2, aNF[1], DS_phiW25_10),
                               Mul(4, bNF[1], DS_phiW05_01), Mul(2, bNF[1], DS_phiW25_01),
                               Mul(4, aNF[1], DS_W12W04_10), Mul(2, aNF[1], DS_W12W24_10),
                               Mul(4, bNF[1], DS_W12W04_01), Mul(2, bNF[1], DS_W12W24_01),
                               Mul(4, aNF[1], DS_W13W03_10), Mul(2, aNF[1], DS_W13W23_10),
                               Mul(4, bNF[1], DS_W13W03_01), Mul(2, bNF[1], DS_W13W23_01),
                               Mul(4, aNF[1], DS_W14W02_10), Mul(2, aNF[1], DS_W14W22_10),
                               Mul(4, bNF[1], DS_W14W02_01), Mul(2, bNF[1], DS_W14W22_01),
                               Mul(4, aNF[2], DS_phiW04_10), Mul(2, aNF[2], DS_phiW24_10),
                               Mul(4, bNF[2], DS_phiW04_01), Mul(2, bNF[2], DS_phiW24_01),
                               Mul(4, aNF[2], DS_W12W03_10), Mul(2, aNF[2], DS_W12W23_10),
                               Mul(4, bNF[2], DS_W12W03_01), Mul(2, bNF[2], DS_W12W23_01),
                               Mul(4, aNF[2], DS_W13W02_10), Mul(2, aNF[2], DS_W13W22_10),
                               Mul(4, bNF[2], DS_W13W02_01), Mul(2, bNF[2], DS_W13W22_01),
                               Mul(4, aNF[3], DS_phiW03_10), Mul(2, aNF[3], DS_phiW23_10),
                               Mul(4, bNF[3], DS_phiW03_01), Mul(2, bNF[3], DS_phiW23_01),
                               Mul(4, aNF[3], DS_W12W02_10), Mul(2, aNF[3], DS_W12W22_10),
                               Mul(4, bNF[3], DS_W12W02_01), Mul(2, bNF[3], DS_W12W22_01),
                               Mul(4, aNF[4], DS_phiW02_10), Mul(2, aNF[4], DS_phiW22_10),
                               Mul(4, bNF[4], DS_phiW02_01), Mul(2, bNF[4], DS_phiW22_01),
                               Mul(4, Pow(aNF[1], 2), DS_phiW04_20),
                               Mul(2, Pow(aNF[1], 2), DS_phiW24_20),
                               Mul(4, aNF[1], bNF[1], DS_phiW04_11),
                               Mul(2, aNF[1], bNF[1], DS_phiW24_11),
                               Mul(4, Pow(bNF[1], 2), DS_phiW04_02),
                               Mul(2, Pow(bNF[1], 2), DS_phiW24_02),
                               Mul(4, Pow(aNF[1], 2), DS_W12W03_20),
                               Mul(2, Pow(aNF[1], 2), DS_W12W23_20),
                               Mul(4, aNF[1], bNF[1], DS_W12W03_11),
                               Mul(2, aNF[1], bNF[1], DS_W12W23_11),
                               Mul(4, Pow(bNF[1], 2), DS_W12W03_02),
                               Mul(2, Pow(bNF[1], 2), DS_W12W23_02),
                               Mul(4, Pow(aNF[1], 2), DS_W13W02_20),
                               Mul(2, Pow(aNF[1], 2), DS_W13W22_20),
                               Mul(4, aNF[1], bNF[1], DS_W13W02_11),
                               Mul(2, aNF[1], bNF[1], DS_W13W22_11),
                               Mul(4, Pow(bNF[1], 2), DS_W13W02_02),
                               Mul(2, Pow(bNF[1], 2), DS_W13W22_02),
                               Mul(8, aNF[1], aNF[2], DS_phiW03_20),
                               Mul(4, aNF[1], aNF[2], DS_phiW23_20),
                               Mul(4, Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), DS_phiW03_11),
                               Mul(2, Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), DS_phiW23_11),
                               Mul(8, bNF[1], bNF[2], DS_phiW03_02),
                               Mul(4, bNF[1], bNF[2], DS_phiW23_02),
                               Mul(8, aNF[1], aNF[2], DS_W12W02_20),
                               Mul(4, aNF[1], aNF[2], DS_W12W22_20),
                               Mul(4, Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), DS_W12W02_11),
                               Mul(2, Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), DS_W12W22_11),
                               Mul(8, bNF[1], bNF[2], DS_W12W02_02),
                               Mul(4, bNF[1], bNF[2], DS_W12W22_02),
                               Mul(8, aNF[1], aNF[3], DS_phiW02_20),
                               Mul(4, aNF[1], aNF[3], DS_phiW22_20),
                               Mul(4, Add(Mul(aNF[1], bNF[3]), Mul(aNF[3], bNF[1])), DS_phiW02_11),
                               Mul(2, Add(Mul(aNF[1], bNF[3]), Mul(aNF[3], bNF[1])), DS_phiW22_11),
                               Mul(8, bNF[1], bNF[3], DS_phiW02_02),
                               Mul(4, bNF[1], bNF[3], DS_phiW22_02),
                               Mul(4, Pow(aNF[2], 2), DS_phiW02_20),
                               Mul(2, Pow(aNF[2], 2), DS_phiW22_20),
                               Mul(4, aNF[2], bNF[2], DS_phiW02_11),
                               Mul(2, aNF[2], bNF[2], DS_phiW22_11),
                               Mul(4, Pow(bNF[2], 2), DS_phiW02_02),
                               Mul(2, Pow(bNF[2], 2), DS_phiW22_02),
                               Mul(4, Pow(aNF[1], 3), DS_phiW03_30),
                               Mul(2, Pow(aNF[1], 3), DS_phiW23_30),
                               Mul(4, Pow(aNF[1], 2), bNF[1], DS_phiW03_21),
                               Mul(2, Pow(aNF[1], 2), bNF[1], DS_phiW23_21),
                               Mul(4, aNF[1], Pow(bNF[1], 2), DS_phiW03_12),
                               Mul(2, aNF[1], Pow(bNF[1], 2), DS_phiW23_12),
                               Mul(4, Pow(bNF[1], 3), DS_phiW03_03),
                               Mul(2, Pow(bNF[1], 3), DS_phiW23_03),
                               Mul(4, Pow(aNF[1], 3), DS_W12W02_30),
                               Mul(2, Pow(aNF[1], 3), DS_W12W22_30),
                               Mul(4, Pow(aNF[1], 2), bNF[1], DS_W12W02_21),
                               Mul(2, Pow(aNF[1], 2), bNF[1], DS_W12W22_21),
                               Mul(4, aNF[1], Pow(bNF[1], 2), DS_W12W02_12),
                               Mul(2, aNF[1], Pow(bNF[1], 2), DS_W12W22_12),
                               Mul(4, Pow(bNF[1], 3), DS_W12W02_03),
                               Mul(2, Pow(bNF[1], 3), DS_W12W22_03),
                               Mul(12, Pow(aNF[1], 2), aNF[2], DS_phiW02_30),
                               Mul(6, Pow(aNF[1], 2), aNF[2], DS_phiW22_30),
                               Mul(8, aNF[1], aNF[2], bNF[1], DS_phiW02_21),
                               Mul(4, aNF[1], aNF[2], bNF[1], DS_phiW22_21),
                               Mul(4, aNF[2], Pow(bNF[1], 2), DS_phiW02_12),
                               Mul(2, aNF[2], Pow(bNF[1], 2), DS_phiW22_12),
                               Mul(4, Pow(aNF[1], 2), bNF[2], DS_phiW02_21),
                               Mul(2, Pow(aNF[1], 2), bNF[2], DS_phiW22_21),
                               Mul(8, aNF[1], bNF[1], bNF[2], DS_phiW02_12),
                               Mul(4, aNF[1], bNF[1], bNF[2], DS_phiW22_12),
                               Mul(12, Pow(bNF[1], 2), bNF[2], DS_phiW02_03),
                               Mul(6, Pow(bNF[1], 2), bNF[2], DS_phiW22_03),
                               Mul(4, Pow(aNF[1], 4), DS_phiW02_40),
                               Mul(2, Pow(aNF[1], 4), DS_phiW22_40),
                               Mul(4, Pow(aNF[1], 3), bNF[1], DS_phiW02_31),
                               Mul(2, Pow(aNF[1], 3), bNF[1], DS_phiW22_31),
                               Mul(4, Pow(aNF[1], 2), Pow(bNF[1], 2), DS_phiW02_22),
                               Mul(2, Pow(aNF[1], 2), Pow(bNF[1], 2), DS_phiW22_22),
                               Mul(4, aNF[1], Pow(bNF[1], 3), DS_phiW02_13),
                               Mul(2, aNF[1], Pow(bNF[1], 3), DS_phiW22_13),
                               Mul(4, Pow(bNF[1], 4), DS_phiW02_04),
                               Mul(2, Pow(bNF[1], 4), DS_phiW22_04),
                               Mul(9, TS_phiphiW155_00), Mul(18, TS_phiW12W14_00),
                               Mul(9, TS_phiW13W13_00), Mul(9, TS_W12W12W13_00),
                               Mul(9, aNF[1], TS_phiphiW14_10), Mul(9, bNF[1], TS_phiphiW14_01),
                               Mul(18, aNF[1], TS_phiW12W13_10), Mul(18, bNF[1], TS_phiW12W13_01),
                               Mul(3, aNF[1], TS_W12W12W12_10), Mul(3, bNF[1], TS_W12W12W12_01),
                               Mul(9, aNF[2], TS_phiphiW13_10), Mul(9, bNF[2], TS_phiphiW13_01),
                               Mul(9, aNF[2], TS_phiW12W12_10), Mul(9, bNF[2], TS_phiW12W12_01),
                               Mul(9, aNF[3], TS_phiphiW12_10), Mul(9, bNF[3], TS_phiphiW12_01),
                               Mul(3, aNF[4], TS_phiphiphi_10), Mul(3, bNF[4], TS_phiphiphi_01),
                               Mul(9, Pow(aNF[1], 2), TS_phiphiW13_20),
                               Mul(9, aNF[1], bNF[1], TS_phiphiW13_11),
                               Mul(9, Pow(bNF[1], 2), TS_phiphiW13_02),
                               Mul(9, Pow(aNF[1], 2), TS_phiW12W12_20),
                               Mul(9, aNF[1], bNF[1], TS_phiW12W12_11),
                               Mul(9, Pow(bNF[1], 2), TS_phiW12W12_02),
                               Mul(18, aNF[1], aNF[2], TS_phiphiW12_20),
                               Mul(9, Add(Mul(aNF[1], bNF[2]), Mul(aNF[2], bNF[1])), TS_phiphiW12_11),
                               Mul(18, bNF[1], bNF[2], TS_phiphiW12_02),
                               Mul(6, aNF[1], aNF[3], TS_phiphiphi_20),
                               Mul(3, Add(Mul(aNF[1], bNF[3]), Mul(aNF[3], bNF[1])), TS_phiphiphi_11),
                               Mul(6, bNF[1], bNF[3], TS_phiphiphi_02),
                               Mul(3, Pow(aNF[2], 2), TS_phiphiphi_20),
                               Mul(3, aNF[2], bNF[2], TS_phiphiphi_11),
                               Mul(3, Pow(bNF[2], 2), TS_phiphiphi_02),
                               Mul(9, Pow(aNF[1], 3), TS_phiphiW12_30),
                               Mul(9, Pow(aNF[1], 2), bNF[1], TS_phiphiW12_21),
                               Mul(9, aNF[1], Pow(bNF[1], 2), TS_phiphiW12_12),
                               Mul(9, Pow(bNF[1], 3), TS_phiphiW12_03),
                               Mul(9, Pow(aNF[1], 2), aNF[2], TS_phiphiphi_30),
                               Mul(6, aNF[1], aNF[2], bNF[1], TS_phiphiphi_21),
                               Mul(3, aNF[2], Pow(bNF[1], 2), TS_phiphiphi_12),
                               Mul(3, Pow(aNF[1], 2), bNF[2], TS_phiphiphi_21),
                               Mul(6, aNF[1], bNF[1], bNF[2], TS_phiphiphi_12),
                               Mul(9, Pow(bNF[1], 2), bNF[2], TS_phiphiphi_03),
                               Mul(3, Pow(aNF[1], 4), TS_phiphiphi_40),
                               Mul(3, Pow(aNF[1], 3), bNF[1], TS_phiphiphi_31),
                               Mul(3, Pow(aNF[1], 2), Pow(bNF[1], 2), TS_phiphiphi_22),
                               Mul(3, aNF[1], Pow(bNF[1], 3), TS_phiphiphi_13),
                               Mul(3, Pow(bNF[1], 4), TS_phiphiphi_04))).subs(extraparvals)

SS_W176_10 = first_order(1, 0, W176NF)
SS_W176_01 = first_order(0, 1, W176NF)
SS_W175_20 = first_order(2, 0, W175NF)
SS_W175_11 = first_order(1, 1, W175NF)
SS_W175_02 = first_order(0, 2, W175NF)
DS_phiW026_00 = second_order(0, 0, phiNF, W026NF)
DS_phiW226_00 = second_order(0, 0, phiNF, W226NF)
DS_W12W025_00 = second_order(0, 0, W12NF, W025NF)
DS_W12W225_00 = second_order(0, 0, W12NF, W225NF)
DS_W22W35_00 = second_order(0, 0, W22NF, W35NF)
DS_W03W124_00 = second_order(0, 0, W03NF, W124NF)
DS_W13W024_00 = second_order(0, 0, W13NF, W024NF)
DS_W13W224_00 = second_order(0, 0, W13NF, W224NF)
DS_W123W04_00 = second_order(0, 0, W123NF, W04NF)
DS_W123W24_00 = second_order(0, 0, W123NF, W24NF)
DS_W23W124_00 = second_order(0, 0, W23NF, W124NF)
DS_W23W34_00 = second_order(0, 0, W23NF, W34NF)
DS_W33W24_00 = second_order(0, 0, W33NF, W24NF)
DS_W165W02_00 = second_order(0, 0, W165NF, W02NF)
DS_W165W22_00 = second_order(0, 0, W165NF, W22NF)
DS_phiW025_10 = second_order(1, 0, phiNF, W025NF)
DS_phiW225_10 = second_order(1, 0, phiNF, W225NF)
DS_phiW025_01 = second_order(0, 1, phiNF, W025NF)
DS_phiW225_01 = second_order(0, 1, phiNF, W225NF)
DS_W12W024_10 = second_order(1, 0, W12NF, W024NF)
DS_W12W224_10 = second_order(1, 0, W12NF, W224NF)
DS_W12W024_01 = second_order(0, 1, W12NF, W024NF)
DS_W12W224_01 = second_order(0, 1, W12NF, W224NF)
DS_W123W03_10 = second_order(1, 0, W123NF, W03NF)
DS_W123W23_10 = second_order(1, 0, W123NF, W23NF)
DS_W123W03_01 = second_order(0, 1, W123NF, W03NF)
DS_W123W23_01 = second_order(0, 1, W123NF, W23NF)
DS_W124W02_10 = second_order(1, 0, W124NF, W02NF)
DS_W124W22_10 = second_order(1, 0, W124NF, W22NF)
DS_W124W02_01 = second_order(0, 1, W124NF, W02NF)
DS_W124W22_01 = second_order(0, 1, W124NF, W22NF)
DS_W22W34_10 = second_order(1, 0, W22NF, W34NF)
DS_W22W34_01 = second_order(0, 1, W22NF, W34NF)
DS_W23W33_10 = second_order(1, 0, W23NF, W33NF)
DS_W23W33_01 = second_order(0, 1, W23NF, W33NF)
DS_phiW024_20 = second_order(2, 0, phiNF, W024NF)
DS_phiW224_20 = second_order(2, 0, phiNF, W224NF)
DS_phiW024_11 = second_order(1, 1, phiNF, W024NF)
DS_phiW224_11 = second_order(1, 1, phiNF, W224NF)
DS_phiW024_02 = second_order(0, 2, phiNF, W024NF)
DS_phiW224_02 = second_order(0, 2, phiNF, W224NF)
DS_W22W33_20 = second_order(2, 0, W22NF, W33NF)
DS_W22W33_11 = second_order(1, 1, W22NF, W33NF)
DS_W22W33_02 = second_order(0, 2, W22NF, W33NF)
DS_W123W02_20 = second_order(2, 0, W123NF, W02NF)
DS_W123W22_20 = second_order(2, 0, W123NF, W22NF)
DS_W123W02_11 = second_order(1, 1, W123NF, W02NF)
DS_W123W22_11 = second_order(1, 1, W123NF, W22NF)
DS_W123W02_02 = second_order(0, 2, W123NF, W02NF)
DS_W123W22_02 = second_order(0, 2, W123NF, W22NF)
TS_phiphiW165_00 = third_order(0, 0, phiNF, phiNF, W165NF)
TS_phiphiW35_00 = third_order(0, 0, phiNF, phiNF, W35NF)
TS_phiW02W04_00 = third_order(0, 0, phiNF, W02NF, W04NF)
TS_phiW02W04_00 = third_order(0, 0, phiNF, W02NF, W04NF)
TS_phiW02W24_00 = third_order(0, 0, phiNF, W02NF, W24NF)
TS_phiW12W124_00 = third_order(0, 0, phiNF, W12NF, W124NF)
TS_phiW12W34_00 = third_order(0, 0, phiNF, W12NF, W34NF)
TS_phiW22W04_00 = third_order(0, 0, phiNF, W22NF, W04NF)
TS_phiW22W24_00 = third_order(0, 0, phiNF, W22NF, W24NF)
TS_phiW03W03_00 = third_order(0, 0, phiNF, W03NF, W03NF)
TS_phiW03W23_00 = third_order(0, 0, phiNF, W03NF, W23NF)
TS_phiW13W123_00 = third_order(0, 0, phiNF, W13NF, W123NF)
TS_phiW13W33_00 = third_order(0, 0, phiNF, W13NF, W33NF)
TS_phiW23W23_00 = third_order(0, 0, phiNF, W23NF, W23NF)
TS_W02W13W02_00 = third_order(0, 0, W02NF, W13NF, W02NF)
TS_W02W13W22_00 = third_order(0, 0, W02NF, W13NF, W22NF)
TS_W02W12W03_00 = third_order(0, 0, W02NF, W12NF, W03NF)
TS_W02W12W23_00 = third_order(0, 0, W02NF, W12NF, W23NF)
TS_W12W12W123_00 = third_order(0, 0, W12NF, W12NF, W123NF)
TS_W12W12W33_00 = third_order(0, 0, W12NF, W12NF, W33NF)
TS_W12W22W03_00 = third_order(0, 0, W12NF, W22NF, W03NF)
TS_W12W22W23_00 = third_order(0, 0, W12NF, W22NF, W23NF)
TS_W22W22W13_00 = third_order(0, 0, W22NF, W22NF, W13NF)
TS_phiphiW124_10 = third_order(1, 0, phiNF, phiNF, W124NF)
TS_phiphiW34_10 = third_order(1, 0, phiNF, phiNF, W34NF)
TS_phiphiW124_01 = third_order(0, 1, phiNF, phiNF, W124NF)
TS_phiphiW34_01 = third_order(0, 1, phiNF, phiNF, W34NF)
TS_phiW02W03_10 = third_order(1, 0, phiNF, W02NF, W03NF)
TS_phiW02W23_10 = third_order(1, 0, phiNF, W02NF, W23NF)
TS_phiW02W03_01 = third_order(0, 1, phiNF, W02NF, W03NF)
TS_phiW02W23_01 = third_order(0, 1, phiNF, W02NF, W23NF)
TS_phiW12W123_10 = third_order(1, 0, phiNF, W12NF, W123NF)
TS_phiW12W33_10 = third_order(1, 0, phiNF, W12NF, W33NF)
TS_phiW12W123_01 = third_order(0, 1, phiNF, W12NF, W123NF)
TS_phiW12W33_01 = third_order(0, 1, phiNF, W12NF, W33NF)
TS_phiW22W03_10 = third_order(1, 0, phiNF, W22NF, W03NF)
TS_phiW22W23_10 = third_order(1, 0, phiNF, W22NF, W23NF)
TS_phiW22W03_01 = third_order(0, 1, phiNF, W22NF, W03NF)
TS_phiW22W23_01 = third_order(0, 1, phiNF, W22NF, W23NF)
TS_W02W02W12_10 = third_order(1, 0, W02NF, W02NF, W12NF)
TS_W02W02W12_01 = third_order(0, 1, W02NF, W02NF, W12NF)
TS_W12W22W02_10 = third_order(1, 0, W12NF, W22NF, W02NF)
TS_W12W22W22_10 = third_order(1, 0, W12NF, W22NF, W22NF)
TS_W12W22W02_01 = third_order(0, 1, W12NF, W22NF, W02NF)
TS_W12W22W22_01 = third_order(0, 1, W12NF, W22NF, W22NF)
TS_phiphiW123_20 = third_order(2, 0, phiNF, phiNF, W123NF)
TS_phiphiW33_20 = third_order(2, 0, phiNF, phiNF, W33NF)
TS_phiphiW123_11 = third_order(1, 1, phiNF, phiNF, W123NF)
TS_phiphiW33_11 = third_order(1, 1, phiNF, phiNF, W33NF)
TS_phiphiW123_02 = third_order(0, 2, phiNF, phiNF, W123NF)
TS_phiphiW33_02 = third_order(0, 2, phiNF, phiNF, W33NF)
TS_phiW02W02_20 = third_order(2, 0, phiNF, W02NF, W02NF)
TS_phiW02W02_11 = third_order(1, 1, phiNF, W02NF, W02NF)
TS_phiW02W02_02 = third_order(0, 2, phiNF, W02NF, W02NF)
TS_phiW22W02_20 = third_order(2, 0, phiNF, W22NF, W02NF)
TS_phiW22W22_20 = third_order(2, 0, phiNF, W22NF, W22NF)
TS_phiW22W02_11 = third_order(1, 1, phiNF, W22NF, W02NF)
TS_phiW22W22_11 = third_order(1, 1, phiNF, W22NF, W22NF)
TS_phiW22W02_02 = third_order(0, 2, phiNF, W22NF, W02NF)
TS_phiW22W22_02 = third_order(0, 2, phiNF, W22NF, W22NF)
Q4S_phiphiphiW04_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W04NF)
Q4S_phiphiphiW24_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W24NF)
Q4S_phiphiW12W03_00 = fourth_order(0, 0, phiNF, phiNF, W12NF, W03NF)
Q4S_phiphiW12W23_00 = fourth_order(0, 0, phiNF, phiNF, W12NF, W23NF)
Q4S_phiphiW13W02_00 = fourth_order(0, 0, phiNF, phiNF, W13NF, W02NF)
Q4S_phiphiW13W22_00 = fourth_order(0, 0, phiNF, phiNF, W13NF, W22NF)
Q4S_phiW12W12W02_00 = fourth_order(0, 0, phiNF, W12NF, W12NF, W02NF)
Q4S_phiW12W12W22_00 = fourth_order(0, 0, phiNF, W12NF, W12NF, W22NF)
Q4S_phiphiphiW03_10 = fourth_order(1, 0, phiNF, phiNF, phiNF, W03NF)
Q4S_phiphiphiW23_10 = fourth_order(1, 0, phiNF, phiNF, phiNF, W23NF)
Q4S_phiphiphiW03_01 = fourth_order(0, 1, phiNF, phiNF, phiNF, W03NF)
Q4S_phiphiphiW23_01 = fourth_order(0, 1, phiNF, phiNF, phiNF, W23NF)
Q4S_phiphiW12W02_10 = fourth_order(1, 0, phiNF, phiNF, W12NF, W02NF)
Q4S_phiphiW12W22_10 = fourth_order(1, 0, phiNF, phiNF, W12NF, W22NF)
Q4S_phiphiW12W02_01 = fourth_order(0, 1, phiNF, phiNF, W12NF, W02NF)
Q4S_phiphiW12W22_01 = fourth_order(0, 1, phiNF, phiNF, W12NF, W22NF)
Q4S_phiphiphiW02_20 = fourth_order(2, 0, phiNF, phiNF, phiNF, W02NF)
Q4S_phiphiphiW22_20 = fourth_order(2, 0, phiNF, phiNF, phiNF, W22NF)
Q4S_phiphiphiW02_11 = fourth_order(1, 1, phiNF, phiNF, phiNF, W02NF)
Q4S_phiphiphiW22_11 = fourth_order(1, 1, phiNF, phiNF, phiNF, W22NF)
Q4S_phiphiphiW02_02 = fourth_order(0, 2, phiNF, phiNF, phiNF, W02NF)
Q4S_phiphiphiW22_02 = fourth_order(0, 2, phiNF, phiNF, phiNF, W22NF)
Q5S_phiphiphiphiW13_00 = fifth_order(0, 0, phiNF, phiNF, phiNF, phiNF, W13NF)
Q5S_phiphiphiW12W12_00 = fifth_order(0, 0, phiNF, phiNF, phiNF, W12NF, W12NF)
Q5S_phiphiphiphiW12_10 = fifth_order(1, 0, phiNF, phiNF, phiNF, phiNF, W12NF)
Q5S_phiphiphiphiW12_01 = fifth_order(0, 1, phiNF, phiNF, phiNF, phiNF, W12NF)
Q5S_phiphiphiphiphi_20 = fifth_order(2, 0, phiNF, phiNF, phiNF, phiNF, phiNF)
Q5S_phiphiphiphiphi_11 = fifth_order(1, 1, phiNF, phiNF, phiNF, phiNF, phiNF)
Q5S_phiphiphiphiphi_02 = fifth_order(0, 2, phiNF, phiNF, phiNF, phiNF, phiNF)

alpha143 = psiNF.dummy.dot(Add(Mul(aNF[1], SS_W176_10), Mul(bNF[1], SS_W176_01),
                               Mul(aNF[2], SS_W175_10), Mul(bNF[2], SS_W175_01),
                               Mul(Pow(aNF[1], 2), SS_W175_20), Mul(aNF[1], bNF[1], SS_W175_11),
                               Mul(Pow(bNF[1], 2), SS_W175_02), Mul(4, DS_phiW026_00),
                               Mul(2, DS_phiW226_00), Mul(4, DS_W12W025_00),
                               Mul(2, DS_W12W225_00), Mul(2, DS_W22W35_00),
                               Mul(4, DS_W03W124_00), Mul(4, DS_W13W024_00),
                               Mul(2, DS_W13W224_00), Mul(4, DS_W123W04_00),
                               Mul(2, DS_W123W24_00), Mul(2, DS_W23W124_00),
                               Mul(2, DS_W23W34_00), Mul(2, DS_W33W24_00),
                               Mul(4, DS_W165W02_00), Mul(2, DS_W165W22_00),
                               Mul(4, aNF[1], DS_phiW025_10), Mul(2, aNF[1], DS_phiW225_10),
                               Mul(4, bNF[1], DS_phiW025_01), Mul(2, bNF[1], DS_phiW225_01),
                               Mul(4, aNF[1], DS_W12W024_10), Mul(2, aNF[1], DS_W12W224_10),
                               Mul(4, bNF[1], DS_W12W024_01), Mul(2, bNF[1], DS_W12W224_01),
                               Mul(4, aNF[1], DS_W123W03_10), Mul(2, aNF[1], DS_W123W23_10),
                               Mul(4, bNF[1], DS_W123W03_01), Mul(2, bNF[1], DS_W123W23_01),
                               Mul(4, aNF[1], DS_W124W02_10), Mul(2, aNF[1], DS_W124W22_10),
                               Mul(4, bNF[1], DS_W124W02_01), Mul(2, bNF[1], DS_W124W22_01),
                               Mul(2, aNF[1], DS_W22W34_10), Mul(2, bNF[1], DS_W22W34_01),
                               Mul(2, aNF[1], DS_W23W33_10), Mul(2, bNF[1], DS_W23W33_01),
                               Mul(4, aNF[2], DS_phiW024_10), Mul(2, aNF[2], DS_phiW224_10),
                               Mul(4, bNF[2], DS_phiW024_01), Mul(2, bNF[2], DS_phiW224_01),
                               Mul(2, aNF[2], DS_W22W33_10), Mul(2, bNF[2], DS_W22W33_01),
                               Mul(4, aNF[2], DS_W123W02_10), Mul(2, aNF[2], DS_W123W22_10),
                               Mul(4, bNF[2], DS_W123W02_01), Mul(2, bNF[2], DS_W123W22_01),
                               Mul(4, Pow(aNF[1], 2), DS_phiW024_20),
                               Mul(2, Pow(aNF[1], 2), DS_phiW224_20),
                               Mul(4, aNF[1], bNF[1], DS_phiW024_11),
                               Mul(2, aNF[1], bNF[1], DS_phiW224_11),
                               Mul(4, Pow(bNF[1], 2), DS_phiW024_02),
                               Mul(2, Pow(bNF[1], 2), DS_phiW224_02),
                               Mul(2, Pow(aNF[1], 2), DS_W22W33_20),
                               Mul(2, aNF[1], bNF[1], DS_W22W33_11),
                               Mul(2, Pow(bNF[1], 2), DS_W22W33_02),
                               Mul(4, Pow(aNF[1], 2), DS_W123W02_20),
                               Mul(2, Pow(aNF[1], 2), DS_W123W22_20),
                               Mul(4, aNF[1], bNF[1], DS_W123W02_11),
                               Mul(2, aNF[1], bNF[1], DS_W123W22_11),
                               Mul(4, Pow(bNF[1], 2), DS_W123W02_02),
                               Mul(2, Pow(bNF[1], 2), DS_W123W22_02),
                               Mul(9, TS_phiphiW165_00), Mul(3, TS_phiphiW35_00),
                               Mul(24, TS_phiW02W04_00), Mul(12, TS_phiW02W24_00),
                               Mul(18, TS_phiW12W124_00), Mul(6, TS_phiW12W34_00),
                               Mul(12, TS_phiW22W04_00), Mul(12, TS_phiW22W24_00),
                               Mul(12, TS_phiW03W03_00), Mul(12, TS_phiW03W23_00),
                               Mul(18, TS_phiW13W123_00), Mul(6, TS_phiW13W33_00),
                               Mul(6, TS_phiW23W23_00), Mul(12, TS_W02W13W02_00),
                               Mul(12, TS_W02W13W22_00), Mul(24, TS_W02W12W03_00),
                               Mul(12, TS_W02W12W23_00), Mul(9, TS_W12W12W123_00),
                               Mul(3, TS_W12W12W33_00), Mul(12, TS_W12W22W03_00),
                               Mul(12, TS_W12W22W23_00), Mul(6, TS_W22W22W13_00),
                               Mul(9, aNF[1], TS_phiphiW124_10), Mul(3, aNF[1], TS_phiphiW34_10),
                               Mul(9, bNF[1], TS_phiphiW124_01), Mul(3, bNF[1], TS_phiphiW34_01),
                               Mul(24, aNF[1], TS_phiW02W03_10), Mul(12, aNF[1], TS_phiW02W23_10),
                               Mul(24, bNF[1], TS_phiW02W03_01), Mul(12, bNF[1], TS_phiW02W23_01),
                               Mul(18, aNF[1], TS_phiW12W123_10), Mul(6, aNF[1], TS_phiW12W33_10),
                               Mul(18, bNF[1], TS_phiW12W123_01), Mul(6, bNF[1], TS_phiW12W33_01),
                               Mul(12, aNF[1], TS_phiW22W03_10), Mul(12, aNF[1], TS_phiW22W23_10),
                               Mul(12, bNF[1], TS_phiW22W03_01), Mul(12, bNF[1], TS_phiW22W23_01),
                               Mul(12, aNF[1], TS_W02W02W12_10), Mul(12, bNF[1], TS_W02W02W12_01),
                               Mul(12, aNF[1], TS_W12W22W02_10), Mul(6, aNF[1], TS_W12W22W22_10),
                               Mul(12, bNF[1], TS_W12W22W02_01), Mul(6, bNF[1], TS_W12W22W22_01),
                               Mul(9, aNF[2], TS_phiphiW123_10), Mul(3, aNF[2], TS_phiphiW33_10),
                               Mul(9, bNF[2], TS_phiphiW123_01), Mul(3, bNF[2], TS_phiphiW33_01),
                               Mul(12, aNF[2], TS_phiW02W02_10), Mul(12, bNF[2], TS_phiW02W02_01),
                               Mul(12, aNF[2], TS_phiW22W02_10), Mul(6, aNF[2], TS_phiW22W22_10),
                               Mul(12, bNF[2], TS_phiW22W02_01), Mul(6, bNF[2], TS_phiW22W22_01),
                               Mul(9, Pow(aNF[1], 2), TS_phiphiW123_20),
                               Mul(3, Pow(aNF[1], 2), TS_phiphiW33_20),
                               Mul(9, aNF[1], bNF[1], TS_phiphiW123_11),
                               Mul(3, aNF[1], bNF[1], TS_phiphiW33_11),
                               Mul(9, Pow(bNF[1], 2), TS_phiphiW123_02),
                               Mul(3, Pow(bNF[1], 2), TS_phiphiW33_02),
                               Mul(12, Pow(aNF[1], 2), TS_phiW02W02_20),
                               Mul(12, aNF[1], bNF[1], TS_phiW02W02_11),
                               Mul(12, Pow(bNF[1], 2), TS_phiW02W02_02),
                               Mul(12, Pow(aNF[1], 2), TS_phiW22W02_20),
                               Mul(6, Pow(aNF[1], 2), TS_phiW22W22_20),
                               Mul(12, aNF[1], bNF[1], TS_phiW22W02_11),
                               Mul(6, aNF[1], bNF[1], TS_phiW22W22_11),
                               Mul(12, Pow(bNF[1], 2), TS_phiW22W02_02),
                               Mul(6, Pow(bNF[1], 2), TS_phiW22W22_02),
                               Mul(24, Q4S_phiphiphiW04_00), Mul(16, Q4S_phiphiphiW24_00),
                               Mul(72, Q4S_phiphiW12W03_00), Mul(48, Q4S_phiphiW12W23_00),
                               Mul(72, Q4S_phiphiW13W02_00), Mul(48, Q4S_phiphiW13W22_00),
                               Mul(72, Q4S_phiW12W12W02_00), Mul(48, Q4S_phiW12W12W22_00),
                               Mul(24, aNF[1], Q4S_phiphiphiW03_10),
                               Mul(16, aNF[1], Q4S_phiphiphiW23_10),
                               Mul(24, bNF[1], Q4S_phiphiphiW03_01),
                               Mul(16, bNF[1], Q4S_phiphiphiW23_01),
                               Mul(72, aNF[1], Q4S_phiphiW12W02_10),
                               Mul(48, aNF[1], Q4S_phiphiW12W22_10),
                               Mul(72, bNF[1], Q4S_phiphiW12W02_01),
                               Mul(48, bNF[1], Q4S_phiphiW12W22_01),
                               Mul(24, aNF[2], Q4S_phiphiphiW02_10),
                               Mul(16, aNF[2], Q4S_phiphiphiW22_10),
                               Mul(24, bNF[2], Q4S_phiphiphiW02_01),
                               Mul(16, bNF[2], Q4S_phiphiphiW22_01),
                               Mul(24, Pow(aNF[1], 2), Q4S_phiphiphiW02_20),
                               Mul(16, Pow(aNF[1], 2), Q4S_phiphiphiW22_20),
                               Mul(24, aNF[1], bNF[1], Q4S_phiphiphiW02_11),
                               Mul(16, aNF[1], bNF[1], Q4S_phiphiphiW22_11),
                               Mul(24, Pow(bNF[1], 2), Q4S_phiphiphiW02_02),
                               Mul(16, Pow(bNF[1], 2), Q4S_phiphiphiW22_02),
                               Mul(50, Q5S_phiphiphiphiW13_00), Mul(100, Q5S_phiphiphiW12W12_00),
                               Mul(50, aNF[1], Q5S_phiphiphiphiW12_10),
                               Mul(50, bNF[1], Q5S_phiphiphiphiW12_01),
                               Mul(10, aNF[2], Q5S_phiphiphiphiphi_10),
                               Mul(10, bNF[2], Q5S_phiphiphiphiphi_01),
                               Mul(10, Pow(aNF[1], 2), Q5S_phiphiphiphiphi_20),
                               Mul(10, aNF[1], bNF[1], Q5S_phiphiphiphiphi_11),
                               Mul(10, Pow(bNF[1], 2), Q5S_phiphiphiphiphi_02))).subs(extraparvals)

DS_phiW036_00 = second_order(0, 0, phiNF, W036NF)
DS_phiW236_00 = second_order(0, 0, phiNF, W236NF)
DS_W22W325_00 = second_order(0, 0, W22NF, W325NF)
DS_W123W024_00 = second_order(0, 0, W123NF, W024NF)
DS_W123W224_00 = second_order(0, 0, W123NF, W224NF)
DS_W33W224_00 = second_order(0, 0, W33NF, W224NF)
DS_W33W44_00 = second_order(0, 0, W33NF, W44NF)
DS_W175W02_00 = second_order(0, 0, W175NF, W02NF)
DS_W175W22_00 = second_order(0, 0, W175NF, W22NF)
TS_phiphiW175_00 = third_order(0, 0, phiNF, phiNF, W175NF)
TS_phiphiW325_00 = third_order(0, 0, phiNF, phiNF, W325NF)
TS_phiW02W024_00 = third_order(0, 0, phiNF, W02NF, W024NF)
TS_phiW02W224_00 = third_order(0, 0, phiNF, W02NF, W224NF)
TS_phiW22W024_00 = third_order(0, 0, phiNF, W22NF, W024NF)
TS_phiW22W224_00 = third_order(0, 0, phiNF, W22NF, W224NF)
TS_phiW22W44_00 = third_order(0, 0, phiNF, W22NF, W44NF)
TS_phiW123W123_00 = third_order(0, 0, phiNF, W123NF, W123NF)
TS_phiW123W33_00 = third_order(0, 0, phiNF, W123NF, W33NF)
TS_phiW33W33_00 = third_order(0, 0, phiNF, W33NF, W33NF)
TS_W02W02W123_00 = third_order(0, 0, W02NF, W02NF, W123NF)
TS_W02W22W123_00 = third_order(0, 0, W02NF, W22NF, W123NF)
TS_W02W22W33_00 = third_order(0, 0, W02NF, W22NF, W33NF)
TS_W22W22W123_00 = third_order(0, 0, W22NF, W22NF, W123NF)
TS_W22W22W33_00 = third_order(0, 0, W22NF, W22NF, W33NF)
Q4S_phiphiphiW024_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W024NF)
Q4S_phiphiphiW224_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W224NF)
Q4S_phiphiphiW44_00 = fourth_order(0, 0, phiNF, phiNF, phiNF, W44NF)
Q4S_phiphiW02W123_00 = fourth_order(0, 0, phiNF, phiNF, W02NF, W123NF)
Q4S_phiphiW02W33_00 = fourth_order(0, 0, phiNF, phiNF, W02NF, W33NF)
Q4S_phiphiW22W123_00 = fourth_order(0, 0, phiNF, phiNF, W22NF, W123NF)
Q4S_phiphiW22W33_00 = fourth_order(0, 0, phiNF, phiNF, W22NF, W33NF)
Q4S_phiW02W02W02_00 = fourth_order(0, 0, phiNF, W02NF, W02NF, W02NF)
Q4S_phiW02W02W22_00 = fourth_order(0, 0, phiNF, W02NF, W02NF, W22NF)
Q4S_phiW22W22W02_00 = fourth_order(0, 0, phiNF, W22NF, W22NF, W02NF)
Q4S_phiW22W22W22_00 = fourth_order(0, 0, phiNF, W22NF, W22NF, W22NF)
Q5S_phiphiphiphiW123_00 = fifth_order(0, 0, phiNF, phiNF, phiNF, phiNF, W123NF)
Q5S_phiphiphiphiW33_00 = fifth_order(0, 0, phiNF, phiNF, phiNF, phiNF, W33NF)
Q5S_phiphiphiW02W02_00 = fifth_order(0, 0, phiNF, phiNF, phiNF, W02NF, W02NF)
Q5S_phiphiphiW02W22_00 = fifth_order(0, 0, phiNF, phiNF, phiNF, W02NF, W22NF)
Q5S_phiphiphiW22W22_00 = fifth_order(0, 0, phiNF, phiNF, phiNF, W22NF, W22NF)
S6S_phiphiphiphiphiW02_00 = sixth_order(0, 0, phiNF, phiNF, phiNF, phiNF, phiNF, W02NF)
S6S_phiphiphiphiphiW22_00 = sixth_order(0, 0, phiNF, phiNF, phiNF, phiNF, phiNF, W22NF)
S7S_phiphiphiphiphiphiphi_00 = seventh_order(0, 0, phiNF, phiNF, phiNF, phiNF, phiNF, phiNF, phiNF)

alpha153 = psiNF.dummy.dot(Add(Mul(4, DS_phiW036_00), Mul(2, DS_phiW236_00),
                               Mul(2, DS_W22W325_00), Mul(4, DS_W123W024_00),
                               Mul(2, DS_W123W224_00), Mul(2, DS_W33W224_00),
                               Mul(2, DS_W33W44_00), Mul(4, DS_W175W02_00),
                               Mul(2, DS_W175W22_00), Mul(9, TS_phiphiW175_00),
                               Mul(3, TS_phiphiW325_00), Mul(24, TS_phiW02W024_00),
                               Mul(12, TS_phiW02W224_00), Mul(12, TS_phiW22W024_00),
                               Mul(12, TS_phiW22W224_00), Mul(6, TS_phiW22W44_00),
                               Mul(9, TS_phiW123W123_00), Mul(6, TS_phiW123W33_00),
                               Mul(6, TS_phiW33W33_00), Mul(12, TS_W02W02W123_00),
                               Mul(12, TS_W02W22W123_00), Mul(12, TS_W02W22W33_00),
                               Mul(6, TS_W22W22W123_00), Mul(3, TS_W22W22W33_00),
                               Mul(24, Q4S_phiphiphiW024_00), Mul(16, Q4S_phiphiphiW224_00),
                               Mul(4, Q4S_phiphiphiW44_00), Mul(72, Q4S_phiphiW02W123_00),
                               Mul(24, Q4S_phiphiW02W33_00), Mul(48, Q4S_phiphiW22W123_00),
                               Mul(36, Q4S_phiphiW22W33_00), Mul(32, Q4S_phiW02W02W02_00),
                               Mul(48, Q4S_phiW02W02W22_00), Mul(48, Q4S_phiW22W22W02_00),
                               Mul(12, Q4S_phiW22W22W22_00), Mul(50, Q5S_phiphiphiphiW123_00),
                               Mul(25, Q5S_phiphiphiphiW33_00), Mul(120, Q5S_phiphiphiW02W02_00),
                               Mul(160, Q5S_phiphiphiW02W22_00), Mul(70, Q5S_phiphiphiW22W22_00),
                               Mul(120, S6S_phiphiphiphiphiW02_00), Mul(90, S6S_phiphiphiphiphiW22_00),
                               Mul(35, S7S_phiphiphiphiphiphiphi_00))).subs(extraparvals)

alpha14 = Add(alpha53, - Mul(alpha5, alpha13, Pow(alpha1, - 1)),
              - Mul(alpha2, Pow(alpha1, - 1), Add(alpha23, Mul(alpha2, alpha13, Pow(alpha1, - 1)))))

alpha24 = Add(alpha63, - Mul(alpha2, alpha33, Pow(alpha1, - 1)), - Mul(2, alpha6, alpha13, Pow(alpha1, - 1)),
              - Mul(alpha3, Pow(alpha1, - 1), Add(alpha23, Mul(2, alpha2, alpha13, Pow(alpha1, - 1)))))

alpha34 = Add(alpha73, - Mul(3, alpha7, alpha13, Pow(alpha1, - 1)),
              - Mul(alpha3, Pow(alpha1, - 1), Add(alpha33, Mul(alpha3, alpha13, Pow(alpha1, - 1)))),
              Mul(alpha4, Pow(alpha1, - 1), Add(alpha43, Mul(alpha4, alpha13, Pow(alpha1, - 1)))))

alpha44 = Add(alpha83, Mul(alpha2, alpha43, Pow(alpha1, - 1)), - Mul(alpha4, alpha23, Pow(alpha1, - 1)),
              - Mul(alpha6, alpha13, Pow(alpha1, - 1)))

alpha54 = Add(alpha93, Mul(alpha3, alpha43, Pow(alpha1, - 1)),
              - Mul(alpha4, alpha33, Pow(alpha1, - 1)),
              - Mul(2, alpha7, alpha13, Pow(alpha1, - 1)))

alpha64 = Add(alpha103, Mul(alpha3, alpha13, Pow(alpha1, - 1)))

alpha74 = Add(alpha113, Mul(Add(alpha3, Mul(2, alpha4)), alpha13, Pow(alpha1, - 1)))

alpha84 = Add(alpha123, - Mul(alpha5, Pow(alpha1, - 1),
                              Add(alpha23, Mul(alpha2, alpha13, Pow(alpha1, - 1)))))

alpha94 = Add(alpha133,
              - Mul(alpha5, Pow(alpha1, - 1),
                    Add(alpha33, alpha43, Mul(Add(alpha3, alpha4), alpha13, Pow(alpha1, - 1)))),
              - Mul(alpha6, Pow(alpha1, - 1), Add(alpha23, Mul(alpha2, alpha13, Pow(alpha1, - 1)))))

alpha104 = Add(alpha143,
               - Mul(alpha6, Pow(alpha1, - 1),
                     Add(alpha33, alpha43, Mul(Add(alpha3, alpha4), alpha13, Pow(alpha1, - 1)))),
               - Mul(alpha7, Pow(alpha1, - 1), Add(alpha23, Mul(alpha2, alpha13, Pow(alpha1, - 1)))))

alpha114 = Add(alpha153,
              - Mul(alpha7, Pow(alpha1, - 1),
                    Add(alpha33, alpha43, Mul(Add(alpha3, alpha4), alpha13, Pow(alpha1, - 1)))))

if simp=='y':
    alpha1 = simplify(evaluation_alpha(alpha1))
    alpha2 = simplify(evaluation_alpha(alpha2))
    alpha3 = simplify(evaluation_alpha(alpha3))
    alpha4 = simplify(evaluation_alpha(alpha4))
    alpha5 = simplify(evaluation_alpha(alpha5))
    alpha6 = simplify(evaluation_alpha(alpha6))
    alpha7 = simplify(evaluation_alpha(alpha7))
else:
    alpha1 = evaluation_alpha(alpha1)
    alpha2 = evaluation_alpha(alpha2)
    alpha3 = evaluation_alpha(alpha3)
    alpha4 = evaluation_alpha(alpha4)
    alpha5 = evaluation_alpha(alpha5)
    alpha6 = evaluation_alpha(alpha6)
    alpha7 = evaluation_alpha(alpha7)

# if len(solve(alpha52, dict = True))==1:
#     extraparvals = extraparvals | solve(alpha52, dict = True)[0]
# else:
#     print(solve(alpha52, dict = True))
#     solnum = input('You have to pick a solution: ')
#     while True:
#         try:
#             extraparvals = extraparvals | solve(alpha52, dict = True)[int(solnum - 1)]
#         except:
#             solnum = input('The number you input was not valid. You have to pick a valid solution: ')

beta1 = - Mul(Add(Pow(alpha2, 2), Mul(4, alpha1, alpha5)), Pow(Mul(4, Pow(alpha1, 2)), - 1))

beta3 = - Mul(Add(Mul(alpha2, Add(alpha3, - alpha4)), Mul(2, alpha1, alpha6)),
              Pow(Mul(4, Pow(alpha1, 2)), - 1))

beta5 = - Mul(Add(Mul(Add(alpha3, alpha4), Add(Mul(3, alpha3), Mul(- 5, alpha4))), Mul(16, alpha1, alpha7)),
              Pow(Mul(48, Pow(alpha1, 2)), - 1))

if modelname=='Swift-Hohenberg':
    beta = Pow(sqrt(734), - 1)
    
    negativeRHS.actualcoord = DS_phiphi_01
    
    Ws4NF = Vector('Ws4NF')
    Ws24NF = Vector('Ws24NF')
    
    Ws4NF = linearsolver(Ws4NF, negativeRHS, coefmat0)
    
    Ws24NF = linearsolver(Ws24NF, negativeRHS, coefmat2)
    
    Ws4NF.actualcoord = Ws4NF.actualcoord.subs(phiNF_eval).subs(extraparvals)
    Ws24NF.actualcoord = Ws24NF.actualcoord.subs(phiNF_eval).subs(extraparvals)
    
    Ws4NF_eval = evaluation_dict(Ws4NF)
    Ws24NF_eval = evaluation_dict(Ws24NF)
    
    DS_phiWs4_00 = second_order(0, 0, phiNF, Ws4NF)
    DS_phiWs24_00 = second_order(0, 0, phiNF, Ws24NF)
    
    deltaE = symbols('deltaE', real = True)
    
    quantity = psiNF.actualcoord.dot(Mul(2, deltaE,
                                         Add(Mul(2, DS_phiWs4_00), DS_phiWs24_00,
                                             Mul(2, DS_phiW02_01), DS_phiW22_01))).subs(extraparvals)
    
    quantity = quantity.subs(Ws4NF_eval).subs(Ws24NF_eval).subs(W02NF_eval).subs(W22NF_eval).subs(phiNF_eval)
elif modelname=='SHDM':
    mu = sqrt(- bNF[4]).subs(extraparvals)
    
    deltar = symbols('deltar', real = True)
    
    quantity = psiNF.actualcoord.dot(Mul(deltar, SS_phi_01)).subs(extraparvals).subs(phiNF_eval)
    
elif modelname=='Brusselator':
    negativeRHS.actualcoord = DS_phiphi_01
    
    Ws4NF = Vector('Ws4NF')
    Ws24NF = Vector('Ws24NF')
    
    Ws4NF = linearsolver(Ws4NF, negativeRHS, coefmat0)
    
    Ws24NF = linearsolver(Ws24NF, negativeRHS, coefmat2)
    
    Ws4NF.actualcoord = Ws4NF.actualcoord.subs(phiNF_eval).subs(extraparvals)
    Ws24NF.actualcoord = Ws24NF.actualcoord.subs(phiNF_eval).subs(extraparvals)
    
    Ws4NF_eval = evaluation_dict(Ws4NF)
    Ws24NF_eval = evaluation_dict(Ws24NF)
    
    DS_phiWs4_00 = second_order(0, 0, phiNF, Ws4NF)
    DS_phiWs24_00 = second_order(0, 0, phiNF, Ws24NF)
    
    deltagamma = symbols('deltagamma', real = True)
    
    quantity1 = psiNF.actualcoord.dot(Mul(2, deltagamma,
                                          Add(Mul(2, DS_phiWs4_00), DS_phiWs24_00,
                                              Mul(2, DS_phiW02_01), DS_phiW22_01))).subs(extraparvals)
    
    quantity2 = psiNF.actualcoord.dot(Mul(deltagamma, SS_phi_01)).subs(extraparvals)
    
    quantity1 = simplify(quantity1.subs(Ws4NF_eval).subs(Ws24NF_eval).subs(W02NF_eval).subs(W22NF_eval).subs(phiNF_eval).subs(muNF, muval).subs(parameters))
    
    quantity2 = simplify(quantity2.subs(phiNF_eval).subs(muNF, muval).subs(parameters))
    
elif modelname=='Bru 4':
    negativeRHS.actualcoord = DS_phiphi_01
    
    Ws4NF = Vector('Ws4NF')
    Ws24NF = Vector('Ws24NF')
    
    Ws4NF = linearsolver(Ws4NF, negativeRHS, coefmat0)
    
    Ws24NF = linearsolver(Ws24NF, negativeRHS, coefmat2)
    
    Ws4NF.actualcoord = Ws4NF.actualcoord.subs(phiNF_eval).subs(extraparvals)
    Ws24NF.actualcoord = Ws24NF.actualcoord.subs(phiNF_eval).subs(extraparvals)
    
    Ws4NF_eval = evaluation_dict(Ws4NF)
    Ws24NF_eval = evaluation_dict(Ws24NF)
    
    DS_phiWs4_00 = second_order(0, 0, phiNF, Ws4NF)
    DS_phiWs24_00 = second_order(0, 0, phiNF, Ws24NF)
    
    deltab = symbols('deltab', real = True)
    
    quantity1 = psiNF.actualcoord.dot(Mul(2, deltab,
                                          Add(Mul(2, DS_phiWs4_00), DS_phiWs24_00,
                                              Mul(2, DS_phiW02_01), DS_phiW22_01))).subs(extraparvals)
    
    quantity2 = psiNF.actualcoord.dot(Mul(deltab, SS_phi_01)).subs(extraparvals)
    
    quantity1 = simplify(quantity1.subs(Ws4NF_eval).subs(Ws24NF_eval).subs(W02NF_eval).subs(W22NF_eval).subs(phiNF_eval).subs(muNF, muval).subs(parameters))
    
    quantity2 = simplify(quantity2.subs(phiNF_eval).subs(muNF, muval).subs(parameters))

Maxwell_equation = simplify(Add(Pow(beta3, 2), Mul(- 4, beta1, beta5)))

if Maxwell_equation!=0:
    if len(solve(Maxwell_equation, dict = True))==0:
        print('There are no parameters to find a Maxwell point.')
        exit()
    elif len(solve(Maxwell_equation, dict = True))==1:
        if simp=='y':
            tempvar = solve(Maxwell_equation, dict = True)[0]
            for par in tempvar.keys():
                tempvar[par] = simplify(tempvar[par])
            extraparvals = extraparvals | tempvar
        else:
            extraparvals = extraparvals | solve(Maxwell_equation, dict = True)[0]
    elif len(solve(Maxwell_equation, dict = True))>1:
        print("You have to choose a solution:")
        print(solve(Maxwell_equation, dict = True))
        solchoice = input("Which solution would you like to use? ")
        while True:
            try:
                int(solchoice)
                if simp=='y':
                    tempvar = solve(Maxwell_equation, dict = True)[int(solchoice) - 1]
                    for par in tempvar.keys():
                        tempvar[par] = simplify(tempvar[par])
                    extraparvals = extraparvals | tempvar
                else:
                    extraparvals = extraparvals | solve(Maxwell_equation, dict = True)[int(solchoice) - 1]
                break
            except:
                solchoice = input("You did not provide a valid integer. Which solution would you like to use? ")

if simp=='y':
    beta1 = simplify(beta1.subs(extraparvals))
    beta3 = simplify(beta3.subs(extraparvals))
    beta5 = simplify(beta5.subs(extraparvals))
else:
    beta1 = beta1.subs(extraparvals)
    beta3 = beta3.subs(extraparvals)
    beta5 = beta5.subs(extraparvals)

xiNF = Mul(alpha2, Pow(Mul(4, alpha1, sqrt(beta1)), - 1)).subs(extraparvals)
etaNF = Mul(Add(alpha3, alpha4), sqrt(beta1), Pow(Mul(4, alpha1, beta3), - 1)).subs(extraparvals)

if simp=='y':
    xiNF = simplify(xiNF)
    etaNF = simplify(etaNF)

Maxwell_equation = simplify(Maxwell_equation.subs(extraparvals))

try:
    if simplify(beta1) > 0 and simplify(beta3) < 0 and simplify(beta5) > 0:
        print('A Maxwell point has been found for the parameter values you provided.')
    else:
        print('Something went wrong. Review extra parameters.')
except:
    print('Some information is required to determine all the parameters of the expansion')

C11NF = Mul(I, exp(Mul(pi, xiNF)), Pow(Mul(2, sqrt(beta1)), Rational(1, 2)),
            Pow(Mul(- 2, beta3), - Rational(1, 2)))

C11NF_conj = Mul(I, exp(Mul(- pi, xiNF)), Pow(Mul(2, sqrt(beta1)), Rational(1, 2)),
                  Pow(Mul(- 2, beta3), - Rational(1, 2)))

if simp=='y':
    C11NF = simplify(C11NF)
    C11NF_conj = simplify(C11NF_conj)

# C02NF = - Mul(2, sqrt(beta1), Pow(Mul(- 2, beta3), - 1))

# C22NF = - Mul(exp(Mul(2, pi, xiNF)), Pow(Mul(2, sqrt(beta1)), 1), Pow(Mul(- 2, beta3), - 1))

# C22NF_conj = - Mul(exp(Mul(- 2, pi, xiNF)), Pow(Mul(2, sqrt(beta1)), 1), Pow(Mul(- 2, beta3), - 1))

# C123NF = Mul(Pow(I, 3), exp(Mul(pi, xiNF)), Pow(Mul(2, sqrt(beta1)), Rational(3, 2)),
#              Pow(Mul(- 2, beta3), - Rational(3, 2)))

# C123NF_conj = Mul(Pow(I, 3), exp(Mul(- pi, xiNF)), Pow(Mul(2, sqrt(beta1)), Rational(3, 2)),
#                   Pow(Mul(- 2, beta3), - Rational(3, 2)))

# C133NF = Mul(I, Rational(1, 2), exp(Mul(pi, xiNF)), Pow(Mul(2, sqrt(beta1)), Rational(1, 2)),
#              Pow(Mul(- 2, beta3), - Rational(1, 2)), Add(Mul(2, etaNF), I))

# C133NF_conj = Mul(I, Rational(1, 2), exp(Mul(- pi, xiNF)),
#                   Pow(Mul(2, sqrt(beta1)), Rational(1, 2)),
#                   Pow(Mul(- 2, beta3), - Rational(1, 2)), Add(Mul(2, etaNF), - I))

# C33NF = Mul(Pow(I, 3), exp(Mul(3, pi, xiNF)), Pow(Mul(2, sqrt(beta1)), Rational(3, 2)),
#             Pow(Mul(- 2, beta3), - Rational(3, 2)))

# C33NF_conj = Mul(Pow(I, 3), exp(Mul(- 3, pi, xiNF)), Pow(Mul(2, sqrt(beta1)), Rational(3, 2)),
#                  Pow(Mul(- 2, beta3), - Rational(3, 2)))

# C024NF = Mul(Pow(Mul(2, sqrt(beta1)), 2), Pow(Mul(- 2, beta3), - 2))

# C034NF = - Mul(I, Rational(1, 2), Pow(Mul(2, sqrt(beta1)), 1), Pow(Mul(- 2, beta3), - 1),
#                Add(1, Mul(- 2, etaNF, I)))

# C034NF_conj = Mul(I, Rational(1, 2), 2, sqrt(beta1), Pow(Mul(- 2, beta3), - 1),
#                   Add(1, Mul(2, etaNF, I)))

# C224NF = Mul(exp(Mul(2, pi, xiNF)), Pow(Mul(2, sqrt(beta1)), 2), Pow(Mul(- 2, beta3), - 2))

# C224NF_conj = Mul(exp(Mul(- 2, pi, xiNF)), Pow(Mul(2, sqrt(beta1)), 2),
#                   Pow(Mul(- 2, beta3), - 2))

# C234NF = - Mul(I, Rational(1, 2), exp(Mul(2, pi, xiNF)), Pow(Mul(2, sqrt(beta1)), 1),
#                Pow(Mul(- 2, beta3), - 1), Add(1, Mul(- 2, etaNF, I)))

# C234NF_conj = Mul(I, Rational(1, 2), exp(Mul(- 2, pi, xiNF)), Pow(Mul(2, sqrt(beta1)), 1),
#                   Pow(Mul(- 2, beta3), - 1), Add(1, Mul(2, etaNF, I)))

# To compare

# N(Mul(C11NF, phiNF.actualcoord[0]))

# N(Mul(2, C02NF, W02NF.actualcoord[0]).subs(phiNF_eval).subs(extraparvals))

# N(Mul(C22NF, W22NF.actualcoord[0]).subs(phiNF_eval).subs(muNF, muval).subs(extraparvals))

# N(Mul(C33NF, W33NF.actualcoord[0]).subs(W22NF_eval).subs(phiNF_eval).subs(muNF, muval).subs(extraparvals))

# N(Add(Mul(2, C024NF, W024NF.actualcoord[0]), Mul(C034NF, W034NF.actualcoord[0])).subs(W123NF_eval).subs(W133NF_eval).subs(W02NF_eval).subs(W22NF_eval).subs(phiNF_eval).subs(muNF, muval).subs(extraparvals))

# N(Add(Mul(C224NF, W224NF.actualcoord[0]), Mul(C234NF, W234NF.actualcoord[0])).subs(W123NF_eval).subs(W133NF_eval).subs(W33NF_eval).subs(W02NF_eval).subs(W22NF_eval).subs(phiNF_eval).subs(muNF, muval).subs(extraparvals))

# Higher order

if simp=='y':
    alpha12 = simplify(evaluation_alpha(alpha12))
    alpha22 = simplify(evaluation_alpha(alpha22))
    alpha32 = simplify(evaluation_alpha(alpha32))
    alpha42 = simplify(evaluation_alpha(alpha42))
    alpha52 = simplify(evaluation_alpha(alpha52))
    alpha62 = simplify(evaluation_alpha(alpha62))
    alpha72 = simplify(evaluation_alpha(alpha72))
    alpha13 = simplify(evaluation_alpha(alpha13))
    alpha23 = simplify(evaluation_alpha(alpha23))
    alpha33 = simplify(evaluation_alpha(alpha33))
    alpha43 = simplify(evaluation_alpha(alpha43))
    alpha53 = simplify(evaluation_alpha(alpha53))
    alpha63 = simplify(evaluation_alpha(alpha63))
    alpha73 = simplify(evaluation_alpha(alpha73))
    alpha83 = simplify(evaluation_alpha(alpha83))
    alpha93 = simplify(evaluation_alpha(alpha93))
    alpha103 = simplify(evaluation_alpha(alpha103))
    alpha113 = simplify(evaluation_alpha(alpha113))
    alpha123 = simplify(evaluation_alpha(alpha123))
    alpha133 = simplify(evaluation_alpha(alpha133))
    alpha143 = simplify(evaluation_alpha(alpha143))
    alpha153 = simplify(evaluation_alpha(alpha153))
    alpha14 = simplify(evaluation_alpha(alpha14))
    alpha24 = simplify(evaluation_alpha(alpha24))
    alpha34 = simplify(evaluation_alpha(alpha34))
    alpha44 = simplify(evaluation_alpha(alpha44))
    alpha54 = simplify(evaluation_alpha(alpha54))
    alpha64 = simplify(evaluation_alpha(alpha64))
    alpha74 = simplify(evaluation_alpha(alpha74))
    alpha84 = simplify(evaluation_alpha(alpha84))
    alpha94 = simplify(evaluation_alpha(alpha94))
    alpha104 = simplify(evaluation_alpha(alpha104))
    alpha114 = simplify(evaluation_alpha(alpha114))
else:
    alpha12 = evaluation_alpha(alpha12)
    alpha22 = evaluation_alpha(alpha22)
    alpha32 = evaluation_alpha(alpha32)
    alpha42 = evaluation_alpha(alpha42)
    alpha52 = evaluation_alpha(alpha52)
    alpha62 = evaluation_alpha(alpha62)
    alpha72 = evaluation_alpha(alpha72)
    alpha13 = evaluation_alpha(alpha13)
    alpha23 = evaluation_alpha(alpha23)
    alpha33 = evaluation_alpha(alpha33)
    alpha43 = evaluation_alpha(alpha43)
    alpha53 = evaluation_alpha(alpha53)
    alpha63 = evaluation_alpha(alpha63)
    alpha73 = evaluation_alpha(alpha73)
    alpha83 = evaluation_alpha(alpha83)
    alpha93 = evaluation_alpha(alpha93)
    alpha103 = evaluation_alpha(alpha103)
    alpha113 = evaluation_alpha(alpha113)
    alpha123 = evaluation_alpha(alpha123)
    alpha133 = evaluation_alpha(alpha133)
    alpha143 = evaluation_alpha(alpha143)
    alpha153 = evaluation_alpha(alpha153)
    alpha14 = evaluation_alpha(alpha14)
    alpha24 = evaluation_alpha(alpha24)
    alpha34 = evaluation_alpha(alpha34)
    alpha44 = evaluation_alpha(alpha44)
    alpha54 = evaluation_alpha(alpha54)
    alpha64 = evaluation_alpha(alpha64)
    alpha74 = evaluation_alpha(alpha74)
    alpha84 = evaluation_alpha(alpha84)
    alpha94 = evaluation_alpha(alpha94)
    alpha104 = evaluation_alpha(alpha104)
    alpha114 = evaluation_alpha(alpha114)

beta33 = Add(Mul(alpha2, Add(alpha24, - alpha44), Pow(Mul(2, alpha1), - 1)),
             Mul(Add(alpha3, - alpha4), alpha14, Pow(Mul(2, alpha1), - 1)),
             Mul(Pow(alpha2, 2), Add(alpha74, - alpha64), Pow(Mul(4, Pow(alpha1, 2)), - 1)), alpha94)

beta53 = Add(Mul(alpha2, Add(alpha34, - alpha54), Pow(Mul(2, alpha1), - 1)),
             Mul(alpha3, Add(Mul(3, alpha24), - alpha44), Pow(Mul(8, alpha1), - 1)),
             - Mul(alpha4, Add(alpha24, Mul(5, alpha44)), Pow(Mul(8, alpha1), - 1)),
             Mul(alpha2, alpha3, Add(Mul(2, alpha74), - Mul(3, alpha64)), Pow(Mul(8, Pow(alpha1, 2)), - 1)),
             Mul(alpha2, alpha4, Add(alpha64, Mul(2, alpha74)), Pow(Mul(8, Pow(alpha1, 2)), - 1)), alpha104)

beta73 = Add(Mul(alpha3, Add(Mul(2, alpha34), - alpha54), Pow(Mul(6, alpha1), - 1)),
             - Mul(alpha4, alpha54, Pow(Mul(2, alpha1), - 1)),
             Mul(Pow(alpha3, 2), Add(Mul(3, alpha74), - Mul(5, alpha64)), Pow(Mul(48, Pow(alpha1, 2)), - 1)),
             Mul(alpha3, alpha4, Add(Mul(3, alpha74), - alpha64), Pow(Mul(24, Pow(alpha1, 2)), - 1)),
             Mul(Pow(alpha4, 2), Add(alpha64, alpha74), Pow(Mul(16, Pow(alpha1, 2)), - 1)), alpha114)

if simple=='y':
    if len(solve(alpha12, dict = True))==1:
        if simp=='y':
            tempvar = solve(alpha12, dict = True)[0]
            for par in tempvar.keys():
                tempvar[par] = simplify(tempvar[par])
            extraparvals = extraparvals | tempvar
        else:
            extraparvals = extraparvals | solve(alpha12, dict = True)[0]
    if len(solve(alpha22, dict = True))==1:
        if simp=='y':
            tempvar = solve(alpha22, dict = True)[0]
            for par in tempvar.keys():
                tempvar[par] = simplify(tempvar[par])
            extraparvals = extraparvals | tempvar
        else:
            extraparvals = extraparvals | solve(alpha22, dict = True)[0]
    if len(solve(alpha32, dict = True))==1:
        if simp=='y':
            tempvar = solve(alpha32, dict = True)[0]
            for par in tempvar.keys():
                tempvar[par] = simplify(tempvar[par])
            extraparvals = extraparvals | tempvar
        else:
            extraparvals = extraparvals | solve(alpha32, dict = True)[0]
    if len(solve(alpha42, dict = True))==1:
        if simp=='y':
            tempvar = solve(alpha42, dict = True)[0]
            for par in tempvar.keys():
                tempvar[par] = simplify(tempvar[par])
            extraparvals = extraparvals | tempvar
        else:
            extraparvals = extraparvals | solve(alpha42, dict = True)[0]
    if len(solve(alpha52, dict = True))==1:
        if simp=='y':
            tempvar = solve(alpha52, dict = True)[0]
            for par in tempvar.keys():
                tempvar[par] = simplify(tempvar[par])
            extraparvals = extraparvals | tempvar
        else:
            extraparvals = extraparvals | solve(alpha52, dict = True)[0]
    if len(solve(alpha62, dict = True))==1:
        if simp=='y':
            tempvar = solve(alpha62, dict = True)[0]
            for par in tempvar.keys():
                tempvar[par] = simplify(tempvar[par])
            extraparvals = extraparvals | tempvar
        else:
            extraparvals = extraparvals | solve(alpha62, dict = True)[0]
    if len(solve(alpha72, dict = True))==1:
        if simp=='y':
            tempvar = solve(alpha72, dict = True)[0]
            for par in tempvar.keys():
                tempvar[par] = simplify(tempvar[par])
            extraparvals = extraparvals | tempvar
        else:
            extraparvals = extraparvals | solve(alpha72, dict = True)[0]

omega1NF = symbols('omega1NF', real = True)
omega2NF = symbols('omega2NF', real = True)
omega3NF = symbols('omega3NF', real = True)
omega4NF = symbols('omega4NF', real = True)
zetaprime = symbols('zetaprime', real = True)

extraparvals[omega2NF] = - Mul(Pow(beta1, 3), Add(alpha64, alpha74),
                               Pow(Mul(6, alpha1, Pow(beta3, 2)), - 1))

extraparvals[omega3NF] = 0

auxeq1 = Add(Mul(- 3, Pow(alpha2, 2), Pow(beta3, 3)),
             Mul(alpha1, Add(Mul(- 6, Pow(beta3, 3), alpha84),
                             Mul(6, beta1, Pow(beta3, 2), beta33),
                             Mul(Pow(beta1, 2), beta3, Add(Mul(beta3, Add(alpha64, alpha74)),
                                                           Mul(- 8, beta53))),
                             Mul(12, Pow(beta1, 3), beta73))),
             Mul(- 3, alpha2, Pow(beta3, 3), Add(alpha14, Mul(2, alpha1, zetaprime)))).subs(extraparvals)

auxeq2 = Add(Mul(3, Pow(alpha2, 2), Add(alpha3, alpha4), Pow(beta3, 3)),
             Mul(2, alpha1, Add(Mul(2, alpha1, Pow(beta3, 2),
                                    Add(Mul(3, Pow(beta3, 2), alpha14),
                                        Mul(- 3, beta1, beta3, Add(alpha24, alpha44)),
                                        Mul(4, Pow(beta1, 2), Add(alpha34, alpha54)))),
                                Mul(Add(alpha3, alpha4),
                                    Add(Mul(3, Pow(beta3, 3), alpha84),
                                        Mul(- Pow(beta1, 2), beta3, Add(Mul(beta3,
                                                                            Add(Mul(5, alpha64), alpha74)),
                                                                        Mul(4, beta53))),
                                        Mul(12, Pow(beta1, 3), beta73))),
                                Mul(12, Pow(alpha1, 2), Pow(beta3, 4), zetaprime))),
             Mul(3, alpha2, Pow(beta3, 3), Add(Mul(Add(alpha3, alpha4), alpha14),
                                               Mul(2, alpha1, Add(Mul(2, alpha1, beta3),
                                                                  Mul(2, beta1, alpha64),
                                                                  Mul(Add(alpha3, alpha4),
                                                                      zetaprime)))))).subs(extraparvals)

if abs(alpha2)<tol:
    extraparvals[zetaprime] = list(solveset(auxeq2, zetaprime))[0]
    
    Maxwell_correction = auxeq1
    
else:
    extraparvals[zetaprime] = list(solveset(auxeq1, zetaprime))[0]
    
    Maxwell_correction = auxeq2

Maxwell_correction = Maxwell_correction.subs(extraparvals)

if Maxwell_correction!=0:
    if len(solve(Maxwell_correction, dict = True))==1:
        if simp=='y':
            tempvar = solve(Maxwell_correction, dict = True)[0]
            for par in tempvar.keys():
                tempvar[par] = simplify(tempvar[par])
            extraparvals = extraparvals | tempvar
        else:
            extraparvals = extraparvals | solve(Maxwell_correction, dict = True)[0]
    else:
        print('You have to choose a solution between:')
        print(solve(Maxwell_correction, dict = True))
        
extraparvals[zetaprime] = extraparvals[zetaprime].subs(extraparvals)
    
den = - Mul(192, Pow(alpha1, 3), Pow(beta1, Rational(5, 4)),
            Pow(- beta3, Rational(9, 2))).subs(extraparvals)

num = Add(Mul(- 3, Add(Mul(3, I), Mul(2, pi)), Pow(alpha2, 2), Pow(beta3, 3),
              Add(Mul(Add(alpha3, alpha4), sqrt(beta1)),
                  Mul(2, I, alpha1, beta3))),
          Mul(alpha1, I, Add(Mul(Add(alpha3, alpha4), sqrt(beta1),
                                 Add(Mul(6, Add(- 3, Mul(2, I, pi)), Pow(beta3, 3),
                                         alpha84),
                                     Mul(6, beta1, Pow(beta3, 2), beta33),
                                     Mul(Pow(beta1, 2), beta3,
                                         Add(Mul(beta3, Add(Mul(15, alpha64), - alpha74)),
                                             Mul(8, beta53))),
                                     Mul(12, Add(- 5, Mul(- 2, I, pi)), Pow(beta1, 3), beta73))),
                             Mul(2, alpha1, beta3,
                                 Add(Mul(beta3, Add(Mul(- 16, Pow(beta1, Rational(5, 2)),
                                                        Add(alpha34, alpha54)),
                                                    Mul(- 6, Add(Mul(3, I), Mul(2, pi)),
                                                        Pow(beta3, 2), alpha84),
                                                    Mul(6, I, beta1, beta3, beta33),
                                                    Mul(I, Pow(beta1, 2),
                                                        Add(Mul(5, beta3, Add(alpha64, alpha74)),
                                                            Mul(8, beta53))))),
                                     Mul(12, Add(Mul(- 3, I), Mul(2, pi)),
                                         Pow(beta1, 3), beta73))))),
          Mul(- 3, Add(Mul(3, I), Mul(2, pi)), alpha2,
              Pow(beta3, 3), Add(Mul(Add(alpha3, alpha4), sqrt(beta1)),
                                 Mul(2, I, alpha1, beta3)),
              Add(alpha14, Mul(2, alpha1, zetaprime)))).subs(extraparvals)

A3_coef = Mul(- I, exp(Mul(pi, xiNF)), num, Pow(den, - 1))

print('The value of ANF3 is given by ' + latex(N(A3_coef.subs(omega1NF, 0))))

phiinf = Mul(2, sqrt(beta1), Add(etaNF, - xiNF))

Delta11 = Mul(sqrt(Mul(- 2, beta1, Pow(beta3, - 1)))).subs(extraparvals)

Delta12 = Mul(sqrt(- Mul(2, beta1, Pow(beta3, - 1))),
              Add(- 1, Mul(2, etaNF, I)), Pow(2, - 1)).subs(extraparvals)

den = Mul(Pow(alpha1, 2), sqrt(beta1), Pow(- beta3, Rational(7, 2)))

num = Mul(Pow(2, Rational(- 5, 2)), Add(Mul(alpha2, Pow(beta3, 3),
                                            Add(alpha2, alpha14, Mul(2, alpha1, zetaprime))),
                                        Mul(2, alpha1,
                                            Add(Mul(Pow(beta3, 3), alpha84),
                                                Mul(- 8, Pow(beta1, 3), beta73),
                                                Mul(4, Pow(beta1, 2), beta3, beta53),
                                                Mul(- 2, beta1, Pow(beta3, 2), beta33),
                                                Mul(4, I, alpha1, beta1, Pow(beta3, 3),
                                                    omega4NF))))).subs(extraparvals)

Delta31 = Mul(num, Pow(den, - 1))

den = Mul(Pow(alpha1, 2), Pow(- beta3, Rational(7, 2))).subs(extraparvals)

num = Mul(Pow(2, - Rational(3, 2)), Add(1, Mul(- 2, etaNF, I)),
          Add(Mul(alpha2, Pow(beta3, 3), Add(alpha2, alpha14, Mul(2, alpha1, zetaprime))),
              Mul(alpha1, Add(Mul(- Pow(beta1, 2), Pow(beta3, 2), Add(alpha64, alpha74)),
                              Mul(2, Pow(beta3, 3), alpha84),
                              Mul(- 4, Pow(beta1, 3), beta73))))).subs(extraparvals)

Delta32 = Mul(num, Pow(den, - 1))

den = Mul(384, Pow(alpha1, 6), sqrt(beta1), Pow(beta3, 7))

num = Add(Mul(3, Pow(alpha2, 2), Pow(beta3, 3),
              Add(Mul(Pow(Add(alpha3, Mul(- 3, alpha4)), 2), Pow(Add(alpha3, alpha4), 2),
                      Pow(beta1, 2)), Mul(- 28, Pow(alpha1, 2), Add(alpha3, Mul(- 3, alpha4)),
                                          Add(alpha3, alpha4), beta1, Pow(beta3, 2)),
                                          Mul(- 144, Pow(alpha1, 4), Pow(beta3, 4)))),
          Mul(- 3, alpha2, Pow(beta3, 3),
              Add(Mul(- Pow(Add(alpha3, Mul(- 3, alpha4)), 2), Pow(Add(alpha3, alpha4), 2),
                      Pow(beta1, 2), alpha14),
                  Mul(8, Pow(alpha1, 2), Add(alpha3, Mul(- 3, alpha4)),
                      Add(alpha3, alpha4), beta1, beta3, Add(Mul(Add(alpha3, - alpha4), beta1),
                                                             Mul(4, beta3, alpha14))),
                  Mul(16, Pow(alpha1, 4), Pow(beta3, 3), Add(Mul(Add(alpha3, Mul(- 15, alpha4)), beta1),
                                                             Mul(8, beta3, alpha14))),
                  Mul(288, Pow(alpha1, 5), Pow(beta3, 4), zetaprime),
                  Mul(- 2, alpha1, Pow(Add(alpha3, Mul(- 3, alpha4)), 2), Add(alpha3, alpha4),
                      Pow(beta1, 2), Add(Mul(- 2, beta1, alpha64), Mul(Add(alpha3, alpha4), zetaprime))),
                  Mul(8, Pow(alpha1, 3), Add(alpha3, Mul(- 3, alpha4)), beta1, Pow(beta3, 2),
                      Add(Mul(12, beta1, alpha64), Mul(7, Add(alpha3, alpha4), zetaprime))))),
          Mul(2, alpha1, Add(Mul(- 2, alpha1, Pow(Add(alpha3, Mul(- 3, alpha4)), 2),
                                 Add(alpha3, alpha4), Pow(beta1, 2), Pow(beta3, 2),
                                 Add(Mul(3, Pow(beta3, 2), alpha14),
                                     Mul(- 3, beta1, beta3, Add(alpha24, alpha44)),
                                     Mul(4, Pow(beta1, 2), Add(alpha34, alpha54)))),
                             Mul(- 48, Pow(alpha1, 3), Add(alpha3, Mul(- 3, alpha4)),
                                 beta1, Pow(beta3, 4), Add(Mul(3, Pow(beta3, 2), alpha14),
                                                           Mul(- 3, beta1, beta3, Add(alpha24, alpha44)),
                                                           Mul(4, Pow(beta1, 2), Add(alpha34, alpha54)))),
                             Mul(- Pow(Add(alpha3, Mul(- 3, alpha4)), 2), Pow(Add(alpha3, alpha4), 2),
                                 Pow(beta1, 2), Add(Mul(- 3, Pow(beta3, 3), alpha84),
                                                    Mul(6, beta1, Pow(beta3, 2), beta33),
                                                    Mul(- 4, Pow(beta1, 2), beta3, Add(Mul(beta3, alpha64),
                                                                                       Mul(3, beta53))),
                                                    Mul(24, Pow(beta1, 3), beta73))),
                             Mul(- 4, Pow(alpha1, 2), Add(alpha3, Mul(- 3, alpha4)), Add(alpha3, alpha4),
                                 beta1, Pow(beta3, 2), Add(Mul(24, Pow(beta3, 3), alpha84),
                                                           Mul(- Pow(beta1, 2), beta3,
                                                               Add(Mul(beta3, Add(Mul(31, alpha64),
                                                                                  Mul(7, alpha74))),
                                                                   Mul(16, beta53))),
                                                           Mul(60, Pow(beta1, 3), beta73),
                                                           Mul(beta1, Pow(beta3, 2),
                                                               Add(Mul(- 6, beta33),
                                                                   Mul(3, Add(alpha3, Mul(- 3, alpha4)),
                                                                       zetaprime))))),
                             Mul(96, Pow(alpha1, 4), Pow(beta3, 4),
                                 Add(Mul(- 4, Pow(beta3, 3), alpha84),
                                     Mul(- 8, Pow(beta1, 2), beta3, beta53),
                                     Mul(8, Pow(beta1, 3), beta73),
                                     Mul(beta1, Pow(beta3, 2), Add(Mul(6, beta33),
                                                                   Mul(Add(Mul(- 3, alpha3),
                                                                           Mul(9, alpha4)),
                                                                       zetaprime)))))))).subs(extraparvals)

gamma111 = Mul(num, Pow(den, - 1))

# Order 2

# C02NF = simplify(Mul(C11NF, C11NF_conj))

# C22NF = Pow(C11NF, 2)

# C22NF_conj = Pow(C11NF_conj, 2)

# # Order 3

# C123NF = simplify(Mul(Pow(C11NF, 2), C11NF_conj))

# C123NF_conj = simplify(Mul(C11NF, C11NF_conj, C11NF_conj))

# C133NF = Mul(C11NF, Rational(1, 2), Add(Mul(2, eta), I))

# C133NF_conj = Mul(C11NF_conj, Rational(1, 2), Add(Mul(2, eta), - I))

# C33NF = Pow(C11NF, 3)

# C33NF_conj = Pow(C11NF_conj, 3)

# # Order 4

# C024NF = simplify(Mul(Pow(C11NF, 2), Pow(C11NF_conj, 2)))

# C034NF = simplify(Mul(I, Rational(1, 2), C11NF, C11NF_conj, Add(1, Mul(- 2, eta, I))))

# C034NF_conj = simplify(Mul(- I, Rational(1, 2), C11NF, C11NF_conj, Add(1, Mul(2, eta, I))))

# C224NF = simplify(Mul(Pow(C11NF, 3), C11NF_conj))

# C224NF_conj = simplify(Mul(C11NF, Pow(C11NF_conj, 3)))

# C234NF = simplify(Mul(I, Rational(1, 2), Pow(C11NF, 2), Add(1, Mul(- 2, eta, I))))

# C234NF_conj = simplify(Mul(- I, Rational(1, 2), Pow(C11NF_conj, 2), Add(1, Mul(2, eta, I))))

# alphaNF = symbols('alphaNF')

print('Here, eta = ' + str(N(etaNF)))

# gammavals1 = [Add(- 1, Mul(- eta, I)), Mul(- eta, I), Add(2, Mul(- eta, I)), Add(3, Mul(- eta, I))]
# gammavals2 = [Add(- 1, Mul(eta, I)), Mul(eta, I), Add(2, Mul(eta, I)), Add(3, Mul(eta, I))]
    
kappa = symbols('kappa')

kappaval = sqrt(Mul(I, Pow(sqrt(muNF), - 1)))

c00NF = symbols('c00')
c20NF = symbols('c20')

gammaNF = symbols('gamma')

h1 = Add(Mul(Add(Mul(4, gammaNF, Add(gammaNF, - 2)), 3), alpha1),
         Mul(12, Pow(C11NF, 2), Pow(C11NF_conj, 2), alpha7),
         Mul(4, C11NF, C11NF_conj, Add(Mul(Add(etaNF, Mul(Add(gammaNF, - 1), I)), alpha3),
                                       Mul(- Add(Mul(2, etaNF), I), alpha4))))

k1 = Mul(2, Pow(C11NF_conj, 2), Add(Mul(Add(Mul(- 2, gammaNF, I), Mul(4, etaNF), I), alpha4),
                                    Mul(Add(Mul(- 2, etaNF), I), alpha3),
                                    Mul(- 4, C11NF, C11NF_conj, alpha7)))

ratio = Mul(h1, Pow(k1, - 1)).subs(gammaNF, Add(3, Mul(- etaNF, I)))

coefK2 = - Mul(Pow(I, - 1), Rational(1, 6), Add(3, Mul(2, etaNF, I)), exp(Mul(pi, xiNF)),
               Pow(Mul(2, sqrt(beta1)), - Rational(1, 2)),
               Pow(Mul(- 2, beta3), Rational(1, 2)))
coefK22 = - Mul(Pow(I, - 1), Rational(1, 6), Add(3, Mul(- 2, etaNF, I)), exp(Mul(- pi, xiNF)),
                Pow(Mul(2, sqrt(beta1)), - Rational(1, 2)),
                Pow(Mul(- 2, beta3), Rational(1, 2)))

invcoefK2 = Pow(coefK2, - 1)

if simp=='y':
    coefK2 = simplify(coefK2)
    coefK22 = simplify(coefK22)
    invcoefK2 = simplify(invcoefK2)