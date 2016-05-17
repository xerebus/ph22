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

if __name__ == "__main__":

    assert len(sys.argv) == 2
    f_input = open(sys.argv[1], "rb")

    while True:

        # stop when there are no more clusters in the pickle file
        try:
            cluster = pickle.load(f_input)
        except:
            break

        radial_plot(cluster)
        del cluster

    f_input.close()
