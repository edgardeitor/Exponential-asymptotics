parameters = {}

unevaluatedparameters = []

par1 = 'C'
par2 = 'E'

muval = '1'

translate = True

equilibrium = [0, 0]

var = ['u',
       'v']

diffmatrix = [[0, - 1],
              [- 1, 0]]

kinetics = ["- (1 + C)*u - 3*E*u**2 - u**3 - 2*v",
            "v"]

tol = 1e-7

parameter_functions = {}

phiunit = 'n'

simp = 'y'

simple = 'y'

def extrapars(extraparvals):
    extraparvals[aNF[0]] = 0
    extraparvals[bNF[0]] = sqrt(Rational(3, 38))
    extraparvals[aNF[4]] = 1
    extraparvals[aNF[6]] = 0
    return extraparvals