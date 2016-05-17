# io_sims.py
# Aritra Biswas
#
# Some multibody simulations to output.

from multibody import *

def dense_sphere():
    '''Small particles randomly distributed in an unit sphere with random
    velocity directions, all with unit speed.'''

    # set number of particles to create
    N = 1000

    # initialize cluster
    cluster = Cluster()
    cluster.force = soft_grav_force

    for i in range(N):

        # determine random particle position
        r = random.uniform(0.0, 1.0)
        theta = random.uniform(0.0, math.pi)
        phi = random.uniform(0.0, 2 * math.pi)
        pos_sph = Vector([r, theta, phi])
        
        # determine random particle velocity
        v = 0.1
        theta = random.uniform(0.0, math.pi)
        phi = random.uniform(0.0, 2 * math.pi)
        vel_sph = Vector([v, theta, phi])

        # create particle
        p = Particle()
        p.set_pos_sph(pos_sph)
        p.set_vel_sph(vel_sph)
        p.radius = 0.01

        # add to cluster
        cluster.add_particle(p)

    # run
    cluster.io_run(0.01, 10000, 1000, "dense_sphere.pkl")

def dense_sphere_at_rest():
    '''Small particles randomly distributed in an unit sphere starting
    at rest.'''

    # set number of particles to create
    N = 1000

    # initialize cluster
    cluster = Cluster()
    cluster.force = soft_grav_force

    for i in range(N):

        # determine random particle position
        r = random.uniform(0.0, 1.0)
        theta = random.uniform(0.0, math.pi)
        phi = random.uniform(0.0, 2 * math.pi)
        pos_sph = Vector([r, theta, phi])
        
        # create particle
        p = Particle()
        p.set_pos_sph(pos_sph)
        p.radius = 0.01

        # add to cluster
        cluster.add_particle(p)

    # run
    cluster.io_run(0.01, 10000, 1000, "dense_sphere_at_rest.pkl")

if __name__ == "__main__":

    print("Pick a simulation:")
    print("(1) Unit sphere of masses and random velocities")
    print("(2) Unit sphere of masses at rest")
    sel = input("Selection: ")

    if sel == 1:
        dense_sphere()
    elif sel == 2:
        dense_sphere_at_rest()
    else:
        print("Invalid selection.")
