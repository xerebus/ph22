from multibody import *
from QuadTree import *

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
        theta = math.pi / 2
        phi = random.uniform(0.0, 2 * math.pi)
        pos_sph = Vector([r, theta, phi])
        
        # determine random particle velocity
        v = 0.1
        theta = math.pi / 2
        phi = random.uniform(0.0, 2 * math.pi)
        vel_sph = Vector([v, theta, phi])

        # create particle
        p = Particle()
        p.set_pos_sph(pos_sph)
        p.set_vel_sph(vel_sph)
        p.radius = 0.01

        # add to cluster
        cluster.add_particle(p)

    # create quadtree
    qt = QuadTree(cluster)

    # visualize
    #cluster.vs_init()
    #cluster.vs_run(0.001, 1)

dense_sphere()
