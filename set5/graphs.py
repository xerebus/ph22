from __future__ import print_function

# graphs.py
# Aritra Biswas
#
# Radial position distributions of various particles.

from multibody import *
from Vector import *
import cPickle as pickle
import matplotlib.pyplot as plt
import sys

def radial_plot(cluster):
    '''Given a cluster object, plot a histogram of the radial positions
    with automatic binning.'''

    rads = [abs(particle.pos) for particle in cluster.particles]
    plt.hist(rads)
    plt.xlabel("$r = |\mathbf r|$ of particle")
    plt.ylabel("number of particles")
    plt.title("Positions at $t/dt = {}$".format(cluster.timestep))
    plt.show()

def velocity_plot(cluster):
    '''Given a cluster object, plot a histogram of the speeds
    with automatic binning.'''

    vels = [abs(particle.vel) for particle in cluster.particles]
    plt.hist(vels, color = "r")
    plt.xlabel("$v = |\mathbf v|$ of particle")
    plt.ylabel("number of particles")
    plt.title("Velocities at $t/dt = {}$".format(cluster.timestep))
    plt.show()

def get_energy(cluster):
    '''Given a cluster object, find the total mechanical energy, assuming
    that cluster.force is a gravitational force.'''

    energy = 0
    
    for p in cluster.particles:

        # kinetic energy
        energy += 0.5 * p.mass * math.pow(abs(p.vel), 2)

        # gravitational potential energy
        # ASSUMPTIONS: force is gravitational, G = 1
        # This is NOT for generalized cluster.force
        G = 1
        for o in cluster.particles:
            if o != p:
                sep = o.pos - p.pos
                energy -= G * p.mass * o.mass / abs(sep)

    return energy

if __name__ == "__main__":

    assert len(sys.argv) >= 2
    files = sys.argv[1:]

    for filename in files:

        f_input = open(filename, "rb")

        clusters = []
        while True:

            # stop when there are no more clusters in the pickle file
            try:
                clusters.append(pickle.load(f_input))
            except:
                break

        f_input.close()

        timesteps = [cluster.timestep for cluster in clusters]
        energies = [get_energy(cluster) for cluster in clusters]
        print("theta = ", clusters[0].theta, " has avg energy: ",
              sum(energies) / len(energies))
        pltlabel = "$\\theta$ = {}".format(clusters[0].theta)
        plt.plot(timesteps, energies, label = pltlabel)

    plt.xlabel("timestep $t/dt$")
    plt.ylabel("energy $T + V$")
    plt.title("Energy evolution for various $\\theta$")
    plt.legend()
    plt.show()
