# general_3body.py
# Aritra Biswas
#
# Simulate a general 3-body problem for comparable masses.

import math
from Vector import *
from integrators import *
import matplotlib.pyplot as plt
from matplotlib import animation

# NUMERICAL PARAMETERS IN NATURAL UNITS

G = 1
M1 = 1
M2 = 1
M3 = 1
T = 8 * math.pi

def general(x10, y10, x20, y20, x30, y30,
            v1x0, v1y0, v2x0, v2y0, v3x0, v3y0):
    '''Return (xx1, yy1, xx2, yy2, xx3, yy3), the motion coordinates
    of 3 bodies under the gravitational influence of each other.'''

    # DIFFERENTIAL EQUATION

    def f(xi):

        # get components
        x1 = xi[0]
        y1 = xi[1]
        x2 = xi[2]
        y2 = xi[3]
        x3 = xi[4]
        y3 = xi[5]
        v1x = xi[6]
        v1y = xi[7]
        v2x = xi[8]
        v2y = xi[9]
        v3x = xi[10]
        v3y = xi[11]

        # couple position derivatives
        x1dot = v1x
        y1dot = v1y
        x2dot = v2x
        y2dot = v2y
        x3dot = v3x
        y3dot = v3y

        # 2d vectors to ease calculation
        r1 = Vector([x1, y1])
        r2 = Vector([x2, y2])
        r3 = Vector([x3, y3])

        # individual forces
        grav1on2 = -G*M1 * (r2 - r1) / math.pow(abs(r2 - r1), 3)
        grav1on3 = -G*M1 * (r3 - r1) / math.pow(abs(r3 - r1), 3)
        grav2on1 = -G*M2 * (r1 - r2) / math.pow(abs(r1 - r2), 3)
        grav2on3 = -G*M2 * (r3 - r2) / math.pow(abs(r3 - r2), 3)
        grav3on1 = -G*M3 * (r1 - r3) / math.pow(abs(r1 - r3), 3)
        grav3on2 = -G*M3 * (r2 - r3) / math.pow(abs(r2 - r3), 3)

        # get acceleration in pieces
        v1dot = grav2on1 + grav3on1
        v2dot = grav1on2 + grav3on2
        v3dot = grav1on3 + grav2on3

        # debug   
        a3 = abs(v3dot)
        if a3 > 5:
            print a3

        # get acceleration components
        v1xdot = v1dot[0]
        v1ydot = v1dot[1]
        v2xdot = v2dot[0]
        v2ydot = v2dot[1]
        v3xdot = v3dot[0]
        v3ydot = v3dot[1]

        return Vector([x1dot, y1dot, x2dot, y2dot, x3dot, y3dot,
                       v1xdot, v1ydot, v2xdot, v2ydot, v3xdot, v3ydot])

    # TIME
    t0 = 0.0
    tf = 3 * T
    h = T / 500

    # INTEGRATE

    # step array
    nn = range(0, int(tf / h))
    tt = [t0 + n*h for n in nn]

    # do RK4 integration
    xi0 = Vector([x10, y10, x20, y20, x30, y30,
                  v1x0, v1y0, v2x0, v2y0, v3x0, v3y0])
    xxi = stepper(xi0, t0, tf, h, f, step_rk4)
    xx1 = [xxi[n][0] for n in nn]
    yy1 = [xxi[n][1] for n in nn]
    xx2 = [xxi[n][2] for n in nn]
    yy2 = [xxi[n][3] for n in nn]
    xx3 = [xxi[n][4] for n in nn]
    yy3 = [xxi[n][5] for n in nn]

    return (xx1, yy1, xx2, yy2, xx3, yy3)

def movie(x10, y10, x20, y20, x30, y30,
          v1x0, v1y0, v2x0, v2y0, v3x0, v3y0):
    '''Make a movie of a 3-body simulation with given initial conditions.'''

    fig = plt.figure()
    ax = plt.axes(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5))
    (line1, line2, line3) = ax.plot([], [], 'bo',
                                    [], [], 'ro',
                                    [], [], 'go')

    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        line3.set_data([], [])
        return (line1, line2, line3)

    (xx1, yy1, xx2, yy2, xx3, yy3) = general(x10, y10, x20, y20, x30, y30,
                                        v1x0, v1y0, v2x0, v2y0, v3x0, v3y0)
    def update(i):
        line1.set_data([xx1[i]], [yy1[i]])
        line2.set_data([xx2[i]], [yy2[i]])
        line3.set_data([xx3[i]], [yy3[i]])
        return (line1, line2, line3)

    anim = animation.FuncAnimation(fig, update, init_func=init,
            frames=500, interval=30, blit=True)

    anim.save("movie.mp4", writer="mencoder", fps=30)

    plt.show()


def trace(x10, y10, x20, y20, x30, y30,
          v1x0, v1y0, v2x0, v2y0, v3x0, v3y0):
    '''Trace out a 3-body simulation with given initial conditions.'''

    fig = plt.figure()
    ax = plt.axes(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5))
    (line1, line2, line3) = ax.plot([], [], 'b-',
                                    [], [], 'r-',
                                    [], [], 'g-')

    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        line3.set_data([], [])
        return (line1, line2, line3)

    (xx1, yy1, xx2, yy2, xx3, yy3) = general(x10, y10, x20, y20, x30, y30,
                                        v1x0, v1y0, v2x0, v2y0, v3x0, v3y0)
    def update(i):
        line1.set_data([xx1[0:i]], [yy1[0:i]])
        line2.set_data([xx2[0:i]], [yy2[0:i]])
        line3.set_data([xx3[0:i]], [yy3[0:i]])
        return (line1, line2, line3)

    anim = animation.FuncAnimation(fig, update, init_func=init,
            frames=1500, interval=30, blit=True)

    plt.show()

def eight():
    
    x30 = 0
    y30 = 0
    x10 = 0.97000436
    y10 = -0.24308753
    x20 = -x10
    y20 = -y10
    v3x0 = -0.93240737
    v3y0 = -0.86473146
    v1x0 = -v3x0 / 2
    v1y0 = -v3y0 / 2
    v2x0 = -v3x0 / 2
    v2y0 = -v3y0 / 2

    movie(x10, y10, x20, y20, x30, y30,
          v1x0, v1y0, v2x0, v2y0, v3x0, v3y0)

def lagrange():
    
    d = 1.5                     # equilateral triangle side
    R = d / math.sqrt(3)        # radius of centripetal motion
    v = math.sqrt(G * M1 / d)   # speed of centripetal motion

    th1 = 0
    th2 = 2 * math.pi / 3
    th3 = 4 * math.pi / 3
    
    x10 = R * math.sin(th1)
    y10 = R * math.cos(th1)
    v1x0 = v * math.cos(th1)
    v1y0 = -v * math.sin(th1)
    
    x20 = R * math.sin(th2)
    y20 = R * math.cos(th2)
    v2x0 = v * math.cos(th2)
    v2y0 = -v * math.sin(th2)
    
    x30 = R * math.sin(th3)
    y30 = R * math.cos(th3)
    v3x0 = v * math.cos(th3)
    v3y0 = -v * math.sin(th3)

    trace(x10, y10, x20, y20, x30, y30,
          v1x0, v1y0, v2x0, v2y0, v3x0, v3y0)

def yin_yang():
    
    v1x0 = 0.51394
    v1y0 = 0.30474
    x10 = -1.0
    x20 = 1.0
    x30 = 0.0
    y10 = 0.0
    y20 = 0.0
    y30 = 0.0
    v2x0 = v1x0
    v3x0 = -2 * v1x0
    v2y0 = v1y0
    v3y0 = -2 * v1y0

    movie(x10, y10, x20, y20, x30, y30,
          v1x0, v1y0, v2x0, v2y0, v3x0, v3y0)

if __name__ == "__main__":

    lagrange()
