# solver_errors.py
# Aritra Biswas
#
# Tracking errors of the bisection, Newton-Rhaphson, and secant methods
# of root-finding for arbitrary functions.

import numpy as np
import matplotlib.pyplot as plt

TOLERANCE = 1e-8  # desired precision for the root

def bisection(f, x1, x2, sol_anal):
    '''Given a function f and two initial guesses x1 and x2,
    find a root up to TOLERANCE using the bisection method.'''

    assert x1 < x2, "initial guesses must have x1 < x2"
    assert f(x1) * f(x2) < 0, "function must cross x-axis between guesses"

    x1, x2 = float(x1), float(x2)
    
    errors = np.array([])
    while (np.abs(x1 - x2) > TOLERANCE):
        m = (x1 + x2) / 2
        if (f(m) * f(x1) > 0): # sgn f(m) = sgn f(x1)
            x1 = m
        else:
            x2 = m
        errors = np.append(errors, np.abs(m - sol_anal))

    return errors

def newton(f, fp, x, sol_anal):
    '''Given a function f, its derivative fp, and an initial guess x,
    find a root up to TOLERANCE using the Newton-Raphson method.'''

    errors = np.array([])
    while True:
        shift = f(x) / fp(x)
        x -= shift
        errors = np.append(errors, np.abs(x - sol_anal))
        if (abs(shift) < TOLERANCE): break

    return errors

def secant(f, x, xn, sol_anal):
    '''Given a function f and two initial guesses x, xn,
    find a root up to TOLERANCE using the Newton-Raphson method
    while approximating the derivative with a secant.'''

    errors = np.array([])
    while True:
        fp = (f(xn) - f(x)) / (xn - x)
        shift = f(xn) / fp
        x = xn
        xn -= shift
        errors = np.append(errors, np.abs(xn - sol_anal))
        if (abs(shift) < TOLERANCE): break

    return errors

if __name__ == "__main__":

    c = 4e-3
    f = lambda x: np.sin(x) - c
    fp = lambda x: np.cos(x)
    sol_anal = np.arcsin(c)
    err_bisect = bisection(f, -0.5, 1, sol_anal)
    err_newton = newton(f, fp, 1, sol_anal)
    err_secant = secant(f, -0.5, 1, sol_anal)

    plt.hold(True)
    plt.plot(np.arange(err_bisect.size), err_bisect, label="Bisection")
    plt.plot(np.arange(err_newton.size), err_newton, label="Newton-Raphson")
    plt.plot(np.arange(err_secant.size), err_secant, label="Secant")
    plt.legend()
    plt.show()
