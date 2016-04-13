# pulsars.py
# Aritra Biswas

# Use numerical methods to solve the orbit of the binary
# pulsar 1913+16.

import numpy as np
import matplotlib.pyplot as plt
import solver

if __name__ == "__main__":

    # constants for 1913+16
    e = 0.617139        # eccentricity
    T = 27906.98161     # orbital period [sec]
    a = 7.0207e8        # projected semimajor axis [m]

    drt = 0.001         # time step in reduced time (time / T)

    rtt = np.arange(0, 1 + drt, drt)     # time array [sec / T]
    xx = np.array([])                    # coords [m]
    yy = np.array([])

    eta_last = 0.       # to store last calculate eta for range

    for rt in rtt:

        # solve parametric equation to find eta
        t = rt * T
        f = lambda eta: (T / (2*np.pi)) * (eta - e*np.sin(eta)) - t
        eta = solver.secant(f, eta_last, eta_last + np.pi)

        # update x, y, and store last eta for range
        xx = np.append(xx, a*(np.cos(eta) - e))
        yy = np.append(yy, a*np.sqrt(1 - e*e)*np.sin(eta))
        eta_last = eta

    plt.plot(xx, yy)
    plt.show()

    vvx = np.array([])  # velocity arrays
    vvy = np.array([])

    for i in range(rtt.size - 1):

        dx = xx[i + 1] - xx[i]
        dy = yy[i + 1] - yy[i]
        vx = dx / drt
        vy = dy / drt

        vvx = np.append(vvx, vx)
        vvy = np.append(vvy, vy)

    # drop last time from rtt
    rtt = np.delete(rtt, rtt.size - 1)

    phi = 1.5 * np.pi
    
    # project velocity along certain angle
    pvv = vvx * np.cos(phi) + vvy * np.sin(phi)

    plt.plot(rtt, pvv)
    plt.show()
