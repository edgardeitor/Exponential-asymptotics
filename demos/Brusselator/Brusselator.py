# parameters = {'delta': Rational(95, 100)}
# parameters = {}
parameters = {'delta': Mul(Add(21, sqrt(313)), Pow(48, - 1))}

# unevaluatedparameters = ['delta']
unevaluatedparameters = []

# par1 = 'alpha'
par1 = 'a'
par2 = 'c'

muval = '(21 + sqrt(313))/(16 * delta**2)'

translate = True

# equilibrium = ['alpha/delta', '(c - 1)*delta/alpha']
equilibrium = ['a', '(c - 1)/a']

var = ['u',
       'v']

diffmatrix=[['delta^2', 0],
            [0, 1]]

# kinetics=['alpha/delta - c * u + u**2*v',
#           '(c - 1)*u - u**2*v']


kinetics=['a - c * u + u**2*v',
          '(c - 1)*u - u**2*v']

tol = 1e-7

parameter_functions = {}

phiunit = 'n'

simp = 'y'

simple = 'y'

def extrapars(extraparvals):
#     extraparvals[aNF[0]] = Mul(Add(21, sqrt(313)), Pow(16, - 1))
    extraparvals[aNF[0]] = 3
    extraparvals[aNF[2]] = 1
    extraparvals[aNF[3]] = 0
    extraparvals[aNF[4]] = 0
    extraparvals[aNF[5]] = 0
    extraparvals[aNF[6]] = 0
#     extraparvals[bNF[0]] = Add(Pow(Add(Mul(Add(21, sqrt(313)), Pow(16, - 1)), 1), 2), 1)
    extraparvals[bNF[0]] = Mul(Add(969, Mul(37, sqrt(313))), Pow(128, - 1))
    extraparvals[bNF[3]] = 0
    return extraparvals