# solver.py
# Aritra Biswas
#
# Implementatations of the bisection, Newton-Rhaphson, and secant methods
# of root-finding for arbitrary functions.

import numpy as np

TOLERANCE = 1e-8    # desired precision for the root

def bisection(f, x1, x2):
    '''Given a function f and two initial guesses x1 and x2,
    find a root up to TOLERANCE using the bisection method.'''

    assert x1 < x2, "initial guesses must have x1 < x2"
    assert f(x1) * f(x2) < 0, "function must cross x-axis between guesses"

    x1, x2 = float(x1), float(x2)

    while (np.abs(x1 - x2) > TOLERANCE):
        m = (x1 + x2) / 2
        if (f(m) * f(x1) > 0): # sgn f(m) = sgn f(x1)
            x1 = m
        else:
            x2 = m

    return m

def newton(f, fp, x):
    '''Given a function f, its derivative fp, and an initial guess x,
    find a root up to TOLERANCE using the Newton-Raphson method.'''

    while True:
        shift = f(x) / fp(x)
        x -= shift
        if (abs(shift) < TOLERANCE): break

    return x

def secant(f, x, xn):
    '''Given a function f and two initial guesses x, xn,
    find a root up to TOLERANCE using the Newton-Raphson method
    while approximating the derivative with a secant.'''

    while True:
        fp = (f(xn) - f(x)) / (xn - x)
        shift = f(xn) / fp
        x = xn
        xn -= shift
        if (abs(shift) < TOLERANCE): break

    return xn

if __name__ == "__main__":

    c = 0.03
    f = lambda x: np.sin(x) - c
    fp = lambda x: np.cos(x)
    sol_anal = np.arcsin(c)
    sol_bisect = bisection(f, 0, 1)
    sol_newton = newton(f, fp, 1)
    sol_secant = secant(f, 0, 1)
    sol_safesc = safe_secant(f, 0, 1)

    print "f(x) = sin(x) - {} = 0".format(c)
    print "Analytical:     {}".format(sol_anal)
    print "Bisection:      {}, error: {}".format(sol_bisect, sol_bisect - sol_anal)
    print "Newton-Raphson: {}, error: {}".format(sol_newton, sol_newton - sol_anal)
    print "Secant:         {}, error: {}".format(sol_secant, sol_secant - sol_anal)
    print "Safe secant:    {}, error: {}".format(sol_safesc, sol_safesc - sol_anal)

