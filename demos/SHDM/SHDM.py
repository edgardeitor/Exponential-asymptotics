parameters = {}

unevaluatedparameters = []

par1 = 's'
par2 = 'r'

muval = '1'

translate = False

equilibrium = [0, 0]

var = ['u',
       'v']

diffmatrix = [[0, - 1],
              [- 1, 0]]

kinetics = ["r*u - u - 2*v + s*u**3 - u**5",
            "v"]

tol = 1e-7

parameter_functions = {}

phiunit = 'n'

simple = 'y'

def extrapars(extraparvals):
    extraparvals[aNF[0]] = 0
    extraparvals[aNF[2]] = 1
    extraparvals[bNF[0]] = 0
    extraparvals[bNF[6]] = 0
    return extraparvals