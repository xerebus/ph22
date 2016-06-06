# integrators.py
# Aritra Biswas
#
# Several ODE integration methods.

from Vector import *
import matplotlib.pyplot as plt
import math

def step_explicit_Euler(xi, t, h, f):
    '''Computes a single step of explicit Euler integration for an ODE
    given by (d/dt) xi(t) = f(xi(t), t).
    input   xi:     Vector, values of xi(t)
    input   t:      float
    input   h:      float, step size
    input   f:      function, for which f(xi(t)) = (d/dt) xi(t)
    returns xi_n:   Vector, values of xi(t + h).'''

    # step using derivative at starting point
    xi_n = xi + h * f(xi)

    return xi_n


def step_midpoint(xi, t, h, f):
    '''Computes a single step of midpoint integration for an ODE
    given by (d/dt) xi(t) = f(xi(t), t).
    input   xi:     Vector, values of xi(t)
    input   t:      float
    input   h:      float, step size
    input   f:      function, for which f(xi(t)) = (d/dt) xi(t)
    returns xi_n:   Vector, values of xi(t + h).'''

    # use Euler step to find midpoint
    xi_mp = xi + (h / 2) * f(xi)

    # step using derivative at midpoint
    xi_n = xi + h * f(xi_mp)

    return xi_n


def step_rk4(xi, t, h, f):
    '''Computes a single step of 4th-order Runge-Kutta integration for
    an ODE given by (d/dt) xi(t) = f(xi(t), t).
    input   xi:     Vector, values of xi(t)
    input   t:      float
    input   h:      float, step size
    input   f:      function, for which f(xi(t)) = (d/dt) xi(t)
    returns xi_n:   Vector, values of xi(t + h).'''

    # calculate intermediate values
    k1 = h * f(xi)
    k2 = h * f(xi + (k1 / 2))
    k3 = h * f(xi + (k2 / 2))
    k4 = h * f(xi + k3)

    # take weighted average
    xi_n = xi + (k1 + 2*k2 + 2*k3 + k4) / 6.0

    return xi_n


def stepper(xi0, t0, tf, h, f, step):
    '''Uses a specific integrator step repeatedly to determine
    xi(t) between t = t0 and t = tf in intervals dt = h.
    input   xi0:    Vector, values of xi(t0)
    input   t0:     float
    input   tf:     float, end time
    input   h:      float, step size dt
    input   f:      function, for which f(xi(t), t) = (d/dt) xi(t)
    input   step:   function, which step routine to use
    returns xxi:    list of vectors, xxi[n] = values of xi(t + n*h).'''

    xxi = [xi0]
    t = t0

    while t < tf:
        xi = xxi[-1]
        xi_n = step(xi, t, h, f)
        xxi.append(xi_n)
        t += h

    return xxi


def convergence_plot():

    # define ODE for an exponential
    f = lambda xi: xi

    # initial conditions for y(t)
    t0 = 0
    y0 = 1.0

    # final time
    tf = 30.0

    # varying step sizes for convergence plot
    h0 = 2.0
    hh = [h0 / s for s in [2, 4, 8, 16, 32, 64, 128]]

    # prepare error arrays
    err_Euler = []
    err_mp = []

    for h in hh:

        # compute step array
        nn = range(0, int(tf / h))
        tt = [t0 + n*h for n in nn]

        # perform Euler integration
        xi0 = Vector([y0])
        xxi = stepper(xi0, t0, tf, h, f, step_explicit_Euler)
        yy_Euler = [xxi[n][0] for n in nn]
        
        # perform midpoint integration
        xi0 = Vector([y0])
        xxi = stepper(xi0, t0, tf, h, f, step_midpoint)
        yy_mp = [xxi[n][0] for n in nn]

        # find errors
        t_last = tt[-1]
        y_last_Euler = yy_Euler[-1]
        y_last_mp = yy_mp[-1]
        y_last_true = math.exp(t_last)
        err_Euler.append(abs(y_last_Euler - y_last_true))
        err_mp.append(abs(y_last_mp - y_last_true))

    # plot errors
    plt.plot(hh, err_Euler, label="explicit Euler")
    plt.plot(hh, err_mp, label="midpoint")
    plt.legend()
    plt.show()


def orbit(pds):
    '''Plot pds full periods of the orbit using RK4 integration.'''

    # define ODEs for orbit
    def f(xi):
        
        x = xi[0]
        y = xi[1]
        vx = xi[2]
        vy = xi[3]

        xdot = vx
        ydot = vy
        vxdot = -x / math.pow(math.pow(x, 2) + math.pow(y, 2), 1.5)
        vydot = -y / math.pow(math.pow(x, 2) + math.pow(y, 2), 1.5)

        return Vector([xdot, ydot, vxdot, vydot])

    # initial conditions
    t0 = 0
    R = 100.0
    x0 = 0
    y0 = R
    vx0 = math.sqrt(1 / R)
    vy0 = 0

    # final time and step size
    tf = 2 * math.pi * math.pow(R, 1.5) * pds
    h = 10
        
    # compute step array
    nn = range(0, int(tf / h))
    tt = [t0 + n*h for n in nn]

    # perform RK4 integration
    xi0 = Vector([x0, y0, vx0, vy0])
    xxi = stepper(xi0, t0, tf, h, f, step_rk4)
    xx = [xxi[n][0] for n in nn]
    yy = [xxi[n][1] for n in nn]

    # plot
    plt.plot(xx, yy)
    plt.show()

if __name__ == "__main__":

    #convergence_plot()
    orbit(200)
