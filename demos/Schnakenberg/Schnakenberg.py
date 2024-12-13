parameters = {'d': Add(3, sqrt(8))}

unevaluatedparameters = []

par1 = 'sigma'
par2 = 'lambda0'

muval = Add(sqrt(2), - 1)

translate = True

equilibrium = ['lambda0', '1/lambda0']

var = ['u',
       'v']

diffmatrix=[['1', 0],
            [0, 'd']]

kinetics=['- u + u**2*v + sigma*(u - 1/v)**2',
          'lambda0 - u**2*v - sigma*(u - 1/v)**2']

tol = 1e-7

parameter_functions = {}

phiunit = 'n'

simple = 'y'

def extrapars(extraparvals):
    extraparvals[aNF[0]] = Add(Rational(- 31, 11), Mul(- 15, sqrt(2), Pow(22, - 1)),
                               Mul(21, sqrt(2), sqrt(Add(9407, Mul(- 6651, sqrt(2)))), Pow(22, - 1)),
                               Mul(15, sqrt(Add(9407, Mul(- 6651, sqrt(2)))), Pow(11, - 1)))
    extraparvals[aNF[2]] = - 1
    extraparvals[aNF[3]] = 0
    extraparvals[aNF[4]] = 0
    extraparvals[aNF[5]] = 0
    extraparvals[aNF[6]] = 0
    extraparvals[bNF[0]] = 1
    extraparvals[bNF[1]] = 0
    extraparvals[bNF[2]] = 0
    extraparvals[bNF[3]] = 0
    return extraparvals