# restricted_3body.py
# Aritra Biswas
#
# Plot the orbit of an asteroid in the Sun-Jupiter system starting near
# Lagrange points L4, L5.

import math
from Vector import *
from integrators import *
    
# NUMERICAL PARAMETERS IN SI

M1 = 1.899e27       # mass of Jupiter
M2 = 1.989e30       # mass of Sun
R = 778.3e9         # semimajor axis of Jupiter's orbit around Sun
Tj = 3.743e8        # period of Jupiter's orbit around Sun
G = 6.6742e-11      # Newton's G

# ENVIRONMENTAL CONDITIONS

# Jupiter's initial and fixed position
X1 = (M2 * R) / (M1 + M2)
Y1 = 0

# Sun's initial and fixed position
X2 = -(M1 * R) / (M1 + M2)
Y2 = 0

# angular velocity of rotating frame
Om = math.sqrt(G*(M1 + M2) / math.pow(R, 3))

# asteroid's period
T = 2 * math.pi / Om

def restricted(lag_point, dalpha):
    '''Return (xx, yy), the motion coordinates of an asteroid in the
    Sun-Jupiter system starting at rest near Lagrange point L4 or L5,
    offset by an angle dalpha.'''

    # LAGRANGE POINTS WITH PERTURBATION

    if lag_point == "L4":
        alpha = (math.pi / 3) + dalpha
    elif lag_point == "L5":
        alpha = -(math.pi / 3) - dalpha
    else:
        raise ValueError("lag_point must be 'L4' or 'L5'")

    Lx = R * ((M2 - M1) / (M1 + M2)) * math.cos(alpha)
    Ly = R * math.sin(alpha)

    # INITIAL CONDITIONS

    x0 = Lx
    y0 = Ly
    vx0 = 0
    vy0 = 0

    # TIME

    t0 = 0.0
    tf = 30 * T
    h = T / 50

    # DIFFERENTIAL EQUATION

    def f(xi):

        # get components of xi
        x = xi[0]
        y = xi[1]
        vx = xi[2]
        vy = xi[3]

        # position derivative is simply the velocity
        xdot = vx
        ydot = vy

        # 2d vectors to ease calculation
        r = Vector([x, y])
        R1 = Vector([X1, Y1])
        R2 = Vector([X2, Y2])

        # calculate acceleration in pieces
        grav1 = -G*M1 * (r - R1) / math.pow(abs(r - R1), 3)
        grav2 = -G*M2 * (r - R2) / math.pow(abs(r - R2), 3)
        coriolis = 2*Om * Vector([vy, -vx])
        centrifugal = math.pow(Om, 2) * r

        # grab acceleration components
        vdot = grav1 + grav2 + coriolis + centrifugal
        vxdot = vdot[0]
        vydot = vdot[1]

        return Vector([xdot, ydot, vxdot, vydot])

    # INTEGRATE

    # step array
    nn = range(0, int(tf / h))
    tt = [t0 + n*h for n in nn]

    # do RK4 integration
    xi0 = Vector([x0, y0, vx0, vy0])
    xxi = stepper(xi0, t0, tf, h, f, step_rk4)
    xx = [xxi[n][0] for n in nn]
    yy = [xxi[n][1] for n in nn]

    return (xx, yy)

if __name__ == "__main__":

    for dalpha in [0.0, 0.01, 0.05, 0.1, 0.5, 1.5]:

        (xx4, yy4) = restricted("L4", dalpha)
        (xx5, yy5) = restricted("L5", dalpha)

        plt.plot([X1], [Y1], 'ro')  # Jupiter
        plt.plot([X2], [Y2], 'y*')  # Sun
        plt.plot(xx4, yy4, 'b')     # asteroid near L4
        plt.plot(xx5, yy5, 'b')     # asteroid near L5

        plt.show()
