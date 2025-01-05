parameters = {'alpha': 1,
              'beta': Rational(1, 2),
              'delta': Mul(Add(21, sqrt(313)), Pow(48, - 1))}

unevaluatedparameters = []

par1 = 'a'
par2 = 'b'

muval = Mul(Rational(9, 8), Add(21, - sqrt(313)))

equilibrium = ['a', 'b/a', 'a', 'b/a']

var = ['u',
       'v',
       'w',
       'z']

diffmatrix=[['delta**2', 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 'delta**2', 0],
            [0, 0, 0, 1]]

kinetics=['a - (b + 1)*u + u**2*v + alpha*(w - u)',
          'b*u - u**2*v + beta*(z - v)',
          'a - (b + 1)*w + w**2*z + alpha*(u - w)',
          'b*w - w**2*z + beta*(v - z)']

tol = 1e-7

parameter_functions = {}

phiunit = 'n'

simp = 'y'

simple = 'y'

def extrapars(extraparvals):
    extraparvals[aNF[0]] = 3
    extraparvals[aNF[1]] = 0
    extraparvals[aNF[2]] = 1
    extraparvals[aNF[3]] = 0
    extraparvals[aNF[4]] = 0
    extraparvals[aNF[5]] = 0
    extraparvals[aNF[6]] = 0
    extraparvals[bNF[0]] = Mul(Add(841, Mul(37, sqrt(313))), Pow(128, - 1))
    extraparvals[bNF[1]] = 0
    return extraparvals