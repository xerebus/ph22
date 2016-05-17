from __future__ import print_function

# multibody.py
# Aritra Biswas
#
# Framework for evolving motion of multiple particles under mutual
# interactions.

from Vector import *
import math
import visual as vs
import random
import time
import cPickle as pickle
import copy

class Particle:
    '''A particle in 3D space, holding a position and velocity.'''

    def __init__(self,
            pos = Vector([0.0, 0.0, 0.0]),
            vel = Vector([0.0, 0.0, 0.0]),
            mass = 1.0,
            radius = 0.1,
        ):

        assert isinstance(pos, Vector) and len(pos) == 3
        assert isinstance(vel, Vector) and len(vel) == 3

        self.pos = pos
        self.vel = vel
        self.mass = mass
        self.radius = radius
    
    def set_pos(self, pos):
        '''Set the position in Cartesian coordinates.'''
        
        assert isinstance(pos, Vector) and len(pos) == 3
        self.pos = pos
    
    def set_vel(self, vel):
        '''Set the velocity in Cartesian coordinates.'''

        assert isinstance(vel, Vector) and len(vel) == 3
        self.vel = vel

    def set_pos_sph(self, pos_sph):
        '''Set the position in spherical coordinates, where the
        input pos_ph = [r, theta, phi]. The physics convention is
        used, where phi is the azimuthal angle.'''

        assert isinstance(pos_sph, Vector) and len(pos_sph) == 3
        (r, theta, phi) = list(pos_sph)

        x = r * math.sin(theta) * math.cos(phi)
        y = r * math.sin(theta) * math.sin(phi)
        z = r * math.cos(theta)

        pos = Vector([x, y, z])
        self.set_pos(pos)

    def set_vel_sph(self, vel_sph):
        '''Set the velocity in spherical coordinates, where the
        input vel_sph = [v, theta, phi]. The physics convention is
        used, where phi is the azimuthal angle.'''

        assert isinstance(vel_sph, Vector) and len(vel_sph) == 3
        (v, theta, phi) = list(vel_sph)

        vx = v * math.sin(theta) * math.cos(phi)
        vy = v * math.sin(theta) * math.sin(phi)
        vz = v * math.cos(theta)

        vel = Vector([vx, vy, vz])
        self.set_vel(vel)

def grav_force(test, source):
    '''Returns the force vector of the particle test due to
    the gravitational force from particle source.'''

    assert isinstance(test, Particle)
    assert isinstance(source, Particle)

    G = 1

    sep = source.pos - test.pos
    force = G * test.mass * source.mass * sep / math.pow(abs(sep), 3)

    return force

def soft_grav_force(test, source):
    '''Similar to grav_force, except the effective distance between two
    particles is always incremented by a small offset, such that the force
    does not blow up when particles are near each other.'''
    
    assert isinstance(test, Particle)
    assert isinstance(source, Particle)

    G = 1
    a = 0.1

    sep = source.pos - test.pos
    force = G * test.mass * source.mass * sep / (math.pow(abs(sep) + a, 3))

    return force

class Cluster:
    '''A cluster of particles interacting with each other under some
    force.'''

    # INITIALIZATION AND PARTICLE MEMBERSHIP 

    def __init__(self):

        # start cluster without any member particles
        self.particles = []

        # interaction force between pairs of particles
        # force(test, source) should return force vector on test
        self.force = grav_force

        # keep track of evolution timesteps
        self.timestep = 0

    def add_particle(self, particle):
        '''Add a particle to the cluster.'''

        assert isinstance(particle, Particle)
        self.particles.append(particle)

    # EVOLUTION

    def evolve(self, dt):
        '''Evolve by flipping between two types of symplectic Euler
        steps. Output the current step and time.'''

        if self.timestep % 2 == 0:
            self.evolve_symplectic1(dt)
        else:
            self.evolve_symplectic2(dt)
        
        out = "step = {} \t time = {}".format(
                self.timestep,
                self.timestep * dt)
        print(out, end = "\r")
        
        self.timestep += 1


    def evolve_symplectic1(self, dt):
        '''Using the symplectic Euler method, update the positions of
        all particles in the cluster for a timestep dt. Update
        velocities before positions.'''

        for test in self.particles:

            # calculate the net acceleration on test
            force = Vector([0.0, 0.0, 0.0])
            for source in self.particles:
                if source != test:
                    force = force + self.force(test, source)
            accel = force / test.mass

            # use current acceleration to update velocity
            test.set_vel(test.vel + accel * dt)

            # use updated velocity to update position
            test.set_pos(test.pos + test.vel * dt)
    
    def evolve_symplectic2(self, dt):
        '''Using the symplectic Euler method, update the positions of
        all particles in the cluster for a timestep dt. Update positions
        before velocities.'''

        for test in self.particles:
            
            # use current velocity to update position
            test.set_pos(test.pos + test.vel * dt)

            # calculate the net acceleration on test
            force = Vector([0.0, 0.0, 0.0])
            for source in self.particles:
                if source != test:
                    force = force + self.force(test, source)
            accel = force / test.mass

            # use updated acceleration to update velocity
            test.set_vel(test.vel + accel * dt)

    # VPYTHON VISUALIZATION

    def vs_init(self):
        '''Set up the visualization with the current particles.'''

        # define a list of acceptable colors to use for particles
        colors = [vs.color.red, vs.color.yellow, vs.color.green,
                  vs.color.orange, vs.color.white, vs.color.blue,
                  vs.color.cyan, vs.color.magenta]

        # initialize a dictionary to link particles with visual objects,
        # such that vs_map[particle] = vs_object
        self.vs_map = {}

        for particle in self.particles:

            # create visual object using particle parameters
            vs_obj = vs.sphere()
            vs_obj.radius = particle.radius
            vs_obj.pos = list(particle.pos)
            vs_obj.color = random.choice(colors)

            # link visual object to particle
            self.vs_map[particle] = vs_obj

    def vs_update(self):
        '''Update the visualization with the current positions of
        the particles.'''

        for particle in self.particles:
            vs_obj = self.vs_map[particle]
            vs_obj.pos = list(particle.pos)

    def vs_run(self, dt, tx):
        '''Continuously evolve and update the visual simulation
        with timestep dt and visual speed tx relative to real-time.'''

        while True:
            self.evolve(dt)
            vs.rate(tx / dt)
            self.vs_update()

    # OUTPUT SAVING

    def io_run(self, dt, steps, wrstep, name = "out.pkl"):
        '''Evolve the cluster for given steps with timestep dt and write
        out the cluster every wrstep timesteps.'''

        f_output = open(name, "wb")

        while self.timestep < steps:

            self.evolve(dt)

            if self.timestep % wrstep == 0:

                obj = copy.deepcopy(self)
                pickle.dump(obj, f_output)
                del obj
                
                print("wrte @")

        f_output.close()
