# vs_sims.py
# Aritra Biswas
#
# Some multibody simulations to visualize in VPython.

from multibody import *

def line_test():
    '''Three equal masses on a line.'''

    # set up 3 objects of unit mass and radius at the origin
    a = Particle()
    b = Particle()
    c = Particle()

    # place objects on the x-axis separated by 10 units
    a.set_pos(Vector([-10.0, 0.0, 0.0]))
    b.set_pos(Vector([0.0, 0.0, 0.0]))
    c.set_pos(Vector([10.0, 0.0, 0.0]))

    # place objects in a cluster
    cluster = Cluster()
    cluster.add_particle(a)
    cluster.add_particle(b)
    cluster.add_particle(c)

    # set the cluster's force
    cluster.force = soft_grav_force

    # visualize
    cluster.vs_init()
    cluster.vs_run(0.1, 15)

def dense_sphere():
    '''Small particles randomly distributed in an unit sphere with random
    velocity directions, all with unit speed.'''

    # set number of particles to create
    N = 30

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

    # visualize
    cluster.vs_init()
    cluster.vs_run(0.001, 1)

def dense_sphere_at_rest():
    '''Small particles randomly distributed in an unit sphere starting
    at rest.'''

    # set number of particles to create
    N = 30

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

    # visualize
    cluster.vs_init()
    cluster.vs_run(0.001, 1)

def choreographic():
    '''Choreographic 3-body orbit from Ph22.2.'''

    # create 3 unit-mass particles
    a = Particle()
    b = Particle()
    c = Particle()

    # get values of initial conditions
    (x3, y3) = (0.0, 0.0)
    (v3x, v3y) = (-0.93240737, -0.86473146)
    (x1, y1) = (0.97000436, -0.24308753)
    (v1x, v1y) = (-v3x / 2, -v3y / 2)
    (x2, y2) = (-x1, -y1)
    (v2x, v2y) = (-v3x / 2, -v3y / 2)

    # set initial conditions
    a.set_pos(Vector([x1, y1, 0.0]))
    a.set_vel(Vector([v1x, v1y, 0.0]))
    b.set_pos(Vector([x2, y2, 0.0]))
    b.set_vel(Vector([v2x, v2y, 0.0]))
    c.set_pos(Vector([x3, y3, 0.0]))
    c.set_vel(Vector([v3x, v3y, 0.0]))

    # build cluster
    cluster = Cluster()
    cluster.add_particle(a)
    cluster.add_particle(b)
    cluster.add_particle(c)
    cluster.force = grav_force

    # visualize
    cluster.vs_init()
    cluster.vs_run(0.01, 5)

if __name__ == "__main__":

    print("Pick a simulation:")
    print("(1) Three equal masses on a line")
    print("(2) Unit sphere of masses and random velocities")
    print("(3) Unit sphere of masses at rest")
    print("(4) Three-body choreographic orbit")
    sel = input("Selection: ")

    if sel == 1:
        line_test()
    elif sel == 2:
        dense_sphere()
    elif sel == 3:
        dense_sphere_at_rest()
    elif sel == 4:
        choreographic()
    else:
        print("Invalid selection.")
