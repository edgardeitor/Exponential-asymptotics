parameters = {'a': 2.4182377602262224,
              'b': 11.684349092782858, 
              'alpha': 1.0,
              'beta': 0.5,
              'delta': 0.5}

unevaluatedparameters = []

par1 = 'a'
par2 = 'b'

muval = 9.673

equilibrium = ['a/delta', 'b*delta/a', 'a/delta', 'b*delta/a']

var = ['u',
       'v',
       'w',
       'z']

diffmatrix=[['delta**2', 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 'delta**2', 0],
            [0, 0, 0, 1]]

kinetics=['a/delta - (b + 1)*u + u**2*v + alpha*(w - u)',
          'b*u - u**2*v + beta*(z - v)',
          'a/delta - (b + 1)*w + w**2*z + alpha*(u - w)',
          'b*w - w**2*z + beta*(v - z)']

tol = 6e-4

parameter_functions = {}

phiunit = 'n'

simp = 'y'

simple = 'y'

def extrapars(extraparvals):
    extraparvals[aNF[0]] = 2.4182377602262224
    extraparvals[aNF[1]] = 0
    extraparvals[aNF[2]] = 1.0
    extraparvals[aNF[3]] = 0
    extraparvals[aNF[4]] = 0
    extraparvals[aNF[5]] = 0
    extraparvals[aNF[6]] = 0
    extraparvals[bNF[0]] = 11.684349092782858
    extraparvals[bNF[1]] = 0
    return extraparvals