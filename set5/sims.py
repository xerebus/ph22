from __future__ import print_function
from multibody import *
import time

def dense_circle(opening_angle):
    '''Small particles randomly distributed in an unit circle with random
    velocity directions, all with unit speed.'''

    # set number of particles to create
    N = 30

    # initialize cluster
    cluster = Cluster()
    cluster.force = soft_grav_force
    cluster.theta = opening_angle

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
        assert p.pos[2] < 1e-10     # check no z components
        assert p.vel[2] < 1e-10
        p.radius = 0.01

        # add to cluster
        cluster.add_particle(p)

    # visualize
    cluster.vs_init()
    cluster.vs_run(0.001, 1)

def dense_circle_out(opening_angle):
    '''Small particles randomly distributed in an unit circle with random
    velocity directions, all with unit speed.'''

    # set number of particles to create
    N = 100

    # initialize cluster
    cluster = Cluster()
    cluster.force = soft_grav_force
    cluster.theta = opening_angle

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
        assert p.pos[2] < 1e-10     # check no z components
        assert p.vel[2] < 1e-10
        p.radius = 0.01

        # add to cluster
        cluster.add_particle(p)

    # run
    outname = "dense_circle_theta_{}.pkl".format(opening_angle)
    cluster.io_run(0.01, 1000, 100, outname)

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

    print("Pick something:")
    print("(1) Visualize 3-body choreographic orbit")
    print("(2) Visualize dense circle with various theta")
    print("(3) Run and save dense circle with various theta")
    sel = input("Selection: ")

    if sel == 1:
        choreographic()
    elif sel == 2:
        while True:
            opening_angle = input("Give opening angle: ")
            assert opening_angle > 0
            dense_circle(opening_angle)
    elif sel == 3:
        random.seed(1738)
        for theta in [0.1, 0.2, 0.3]:
            print()
            print("Simulating for theta = ", theta)
            start_time = time.time()
            dense_circle_out(theta)
            print("Wall time: ", time.time() - start_time)
